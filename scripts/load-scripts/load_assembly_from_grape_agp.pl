#!/usr/local/bin/perl



BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/gramene_ensembl/'; 
}


# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

#use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } qw ( bioperl-live modules ensembl/modules conf  ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);

use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;

#you may not want all of these:
use Gramene::Config;
#use DBI;
use DBI qw(:sql_types);
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Registry;
use Text::RecordParser::Tab;

=head1 SYNOPSIS

load_assembly_from_grape_agp.pl  [options] 

for loading assembly from agp files from

http://www.genoscope.cns.fr/externe/English/Projets/Projet_ML/data/assembly/

Need to load chr.fasta and scaffold fasta files before hand using

 $ perl /usr/local/gramene_ensembl/load-scripts/load_assembly_from_fasta.pl \
	--species=grape --assembly_version=8X \
	--coord_system=chromosome -not_seq_level \
	/usr/local/data/genome/grape_10102007/grape.fa

and 

 $ perl /usr/local/gramene_ensembl/load-scripts/load_assembly_from_fasta.pl \
	-species=grape --assembly_version=8X \
	--coord_system=scaffold \
 	/usr/local/data/genome/grape_10102007/scaffolds/Vvinifera_v4.0.fa


 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for ensembl db adaptors
    --species		species in registry file
    --random            process random chromosomes
    
=head1 OPTIONS

=over 3

=item B<--registry>

The registry file

=item B<--species>

The species name in the registry file

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

=item B<--random> 

if set, only process random chrs, if not set, only process real chromosomes
    
=back

=head1 ARGUMENTS

list of agp files downloaded from http://www.genoscope.cns.fr/externe/English/Projets/Projet_ML/data/assembly/
    
The file format is like

chr1    1       7473787 S       scaffold_5      1       7473787 +       7473787

chr1    7473788 7483787 N       GAP_UC  1       10000   .       10000

chr1    7483788 9922561 S       scaffold_46     1       2438774 +       2438774

chr1    9922562 9932561 N       GAP_UC  1       10000   .       10000

Field1 : chromosome name

Field2 : Start position on the chromosome

Filed3 : End position on the chromosome

Field4 : Type of sequence (S=scaffold, N=gap)

Field5 : Name of sequence

          - SCAFxxxx : scaffold xxxx
          - GAP_IS   : estimated gap between scaffolds in ultracontigs
          - GAP_UC   : gap between two mapped ultracontigs (10Kb)
          - GAP_UN   : gap between 2 random ultracontigs, or between 
                       ultracontigs and scaffold, not mapped (1kb)

Field6 : Start position on scaffold

Field7 : End position on scaffold

Field8 : scaffold's strand (on the chromosome)

Field9 : scaffold's length

=cut

my %clone2acc;
my $ensembl_species=$ENV{ENSEMBL_SPECIES};
my $random;
my $registry_file;
    
{  #Argument Processing
  
  my $help=0;
  my $man=0;
  GetOptions( "help|?"     => \$help,
	      "man"        => \$man,
	      "registry=s" => \$registry_file,
	      "species=s"  => \$ensembl_species,
	      "random"     => \$random,
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  
}

pod2usage(-verbose => 2) unless @ARGV ;

###
### Get Ensembl DB adaptor
###

$ENV{'ENSEMBL_SPECIES'} = $ensembl_species;
$registry_file       ||= "$ENV{GrameneDir}/conf/ensembl.registry";

my $reg = "Bio::EnsEMBL::Registry";	#Use this to get adaptors 
$reg->load_all( $registry_file ) or die "load_all failed";
     
my $slice_adaptor =$reg->get_adaptor($ensembl_species,'core','Slice') 
  or die "can't get Slice adaptor for $ensembl_species";

my $attribute_adaptor =$reg->get_adaptor($ensembl_species,'core','Attribute')
     or die "can't get Attribute adaptor for $ensembl_species";

my $coord_system_adaptor = $reg->get_adaptor($ensembl_species,'core'
                                     ,'CoordSystem');

my $chr_cs      = $coord_system_adaptor->fetch_by_name('chromosome');
my $scaffold_cs = $coord_system_adaptor->fetch_by_name('scaffold');
$coord_system_adaptor->store_mapping_path($chr_cs, $scaffold_cs );

my $dba=$slice_adaptor->db; #DBAdaptor
my $dbc=$dba->dbc;	#DBConnection
warn "user ".$dbc->username
  .", db ".$dbc->dbname."\n";


#my $ens_dbh=$dbc->db_handle; #you can use this as a DBI database handle

my $assembly_sth = $dbc->prepare("insert into assembly 
      (asm_seq_region_id ,cmp_seq_region_id ,asm_start ,asm_end ,cmp_start ,cmp_end ,ori)
      values  (?, ?,?,?,?,?,?)");

#my $decruft_sth  = $dbc->prepare('DELETE sr, a 
#			       FROM   seq_region sr, assembly a 
#			       WHERE  sr.seq_region_id=a.asm_seq_region_id 
#			       AND    a.cmp_seq_region_id = ?
#			       AND    sr.name like ?');
###
### Get TIGR ASSEMBLY data
###
#>Count  cM      asmbl_id        clone_name      gb_acc  pseudo_lend     pseudo_rend     BAC_start       BAC_end seq_group       BAC_Len gb_phase        gb_gBAC_ori


my $p = new Text::RecordParser::Tab;

my $j = 1; # as suffix of gap names GAP_UC_1
for my $file (@ARGV){
  
  $p->filename($file);
  $p->bind_fields( qw[ chr chr_start chr_end seq_type seq_name scaffold_start scaffold_end scaffold_strand scaffold_len ] );
  #my $records = $p->fetchall_hashref('>Count');
  my $records = $p->fetchall_arrayref( { Columns => {} } );
  
  my @records_array = @{$records};


  my %chrs_hash;
  for my $record( @records_array ){
    
    my $chr = $record->{'chr'};
    next if ( ($chr =~ /_random/i && !$random) ||
	      ($chr !~ /_random/i && $random)
	    );
    push @{$chrs_hash{ $record->{'chr'} }}, $record;
    
  }

  for my $chr( keys %chrs_hash ){
				 
     my @chr_array = sort { $a->{'chr_start'} <=> $b->{'chr_start'} } @{$chrs_hash{$chr}};
     my $last_member = $chr_array[$#chr_array];

    # health check the chr assembly, make sure no gap and length adds up

     my $pseudo_len     = $last_member->{'chr_end'};
  
     my $pseudo_len_cal = 0;
     my ($pre, $gap)    = (0,0);
  
     for my $rd( @chr_array ){
 
       $pseudo_len_cal += ($rd->{'chr_end'} -  $rd->{'chr_start'} + 1);
       
       if($pre && ($rd->{'chr_start'} - $pre) > 1 ){
	 $gap += $rd->{'chr_start'} - $pre - 1;
	 die "There is a gap between $rd->{'chr_start'} and $pre\n";
       }
       $pre = $rd->{'chr_end'};
       
     }


     unless( $pseudo_len == ($pseudo_len_cal + $gap)){
       die "The pseudomolecule $chr length not match, $pseudo_len != $pseudo_len_cal (gap size = $gap)";
     }
  

     # now start to load assembly 
     #
     my $chr_num = $chr;
     $chr_num    =~ s/chr(omosome)?_?0*//i;
     my $chr_slice = get_slice($coord_system_adaptor,
				$slice_adaptor,
				$chr_num,
				$pseudo_len,
			      'chromosome') 
				 or die "Can't get slice for chromosome $chr_num";
  
     for ( my $i=0; $i < scalar @chr_array; $i++ ){
    
       my $r = $chr_array[$i];
       if( $r->{'chr_start'} < 1 || $r->{'chr_start'} > $r->{'chr_end'}
	   || $r->{'chr_end'} > $chr_slice->end) {
	 
	 die "bad coordinates for $chr_num: " . (join ", ", values %{$r} );
	 #next;
       }
    
       # We don't need to save gap in seq_region_table
       # so skip
       # otherwise, comment the if block out
       if( $r->{'seq_type'} eq 'N' ){
	 next;
       }

       my $scaffold_name = $r->{'seq_name'};
       $scaffold_name .=  $j++ if ( $r->{'seq_type'} eq 'N' );


       my $scaffold_len = $r->{'scaffold_len'};
       
       my $scaffold_seq = undef;
       $scaffold_seq    = 'N' x $scaffold_len if ( $r->{'seq_type'} eq 'N' );

       my $scaffold_slice = get_slice($coord_system_adaptor,
				      $slice_adaptor,
				      $scaffold_name,
				      $scaffold_len,
				      'scaffold',
				      $scaffold_seq,
				     ) #$slice_adaptor->fetch_by_region('scaffold', $scaffold_name)
      or die "Can't get slice for scaffold $scaffold_name";
    
       
       my ( $cmp_start, $cmp_end, $ori )=( $r->{'scaffold_start'}, 
					   $r->{'scaffold_end'},
					   ($r->{'scaffold_strand'} eq '-' ?
					    -1 : 1),
					 );
    
       printf "chr_seq_region_id = %s, scaffold_seq_region_id = %s\n",
	 $chr_slice->get_seq_region_id, $scaffold_slice->get_seq_region_id;
       $assembly_sth->execute($chr_slice->get_seq_region_id, 
			      $scaffold_slice->get_seq_region_id,
			      $r->{'chr_start'}, $r->{'chr_end'},
			      $cmp_start ,$cmp_end ,$ori);
       #$decruft_sth->execute($clone_slice->get_seq_region_id,"R%_$accn");
     }
  
   }			

}

exit;
    
# end of main program

########################## subroutines ######################################
sub get_slice { #add if necessary
  my ($csa,$sa,$chr,$length, $cs, $seq)=@_;
  my $chr_slice=$sa->fetch_by_region($cs,$chr) ;
  unless( $chr_slice) {
    
    my $chrcs = $csa->fetch_by_name($cs);
    $chr_slice = Bio::EnsEMBL::Slice->new(-SEQ_REGION_NAME => $chr
					  ,-COORD_SYSTEM    => $chrcs
					  ,-START           => 1
					  ,-END             => $length
					  ,-SEQ_REGION_LENGTH => $length);
    if($seq){
      
      $slice_adaptor->store($chr_slice, \$seq );
    }else{
      $slice_adaptor->store($chr_slice );
    }
  }
  return $chr_slice;
}


sub get_scaffold_slice { #add if necessary
  my ($csa,$sa,$chr,$length)=@_;
  my $chr_slice=$sa->fetch_by_region('chromosome',$chr) ;
  unless( $chr_slice) {
    
    my $chrcs = $csa->fetch_by_name('chromosome');
    $chr_slice = Bio::EnsEMBL::Slice->new(-SEQ_REGION_NAME => $chr
					  ,-COORD_SYSTEM    => $chrcs
					  ,-START           => 1
					  ,-END             => $length
					  ,-SEQ_REGION_LENGTH => $length);
    $slice_adaptor->store($chr_slice );
  }
  return $chr_slice;
}


__END__

=head1 OUTPUT

=item B<Standard Output>

=head1 NOTES
    
Assumes TIGR xml files have been loaded in same server as
core ensembl db for species
and are accessible by the same user.

=head1 AUTHOR

   Steven Schmidt (schmidt@cshl.edu)
   Gramene Project (www.gramene.org)
   Stein Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

