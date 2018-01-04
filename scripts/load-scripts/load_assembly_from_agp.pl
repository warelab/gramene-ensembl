#!/usr/local/bin/perl



BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/'; 
}


# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

#use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } qw ( bioperl-live modules ensembl/modules conf  ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);

use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;

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

Need to load chr.fasta and scaffold fasta files beforehand using

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
    --species		species in EnsEMBL Web conf to use for db
    --random            process random chromosomes
    --registry_file     registry file for the core database
    --csa               assembly coordinate system
    --csc               component coordinate system

=head1 OPTIONS

=over 3

=item B<--species>

The species name in the conf file, species.ini

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

OB_chr3S        1       853402  1       W       OB_3S_superscaffold_1   1       853402  +
OB_chr3S        853403  853412  2       N       100     contig  no

OB_3S_superscaffold_1   1       115173  1       W       3S_AB_Left_scaffold00019        1       115173  +
OB_3S_superscaffold_1   115174  115233  2       N       60      fragment        yes

3S_EF_Right_scaffold00001       1       692     1       W       3S_EF_Right_contig00004 1       692     +
3S_EF_Right_scaffold00001       693     2776    2       N       2084    fragment        yes

Field1 : chromosome name

Field2 : Start position on the chromosome

Filed3 : End position on the chromosome

Field4 : order of the component sequence

Field5 : Gap/non gap, possible values with sequences: A/D/F/G/O/P/W, Gaps: N/U 

Field6 : component seq name or number of Ns for gap region

Field7 : component seq start or fragment/contig for gap region

Field8 : component seq end or yes/no for gap region

Field9 : component seq strand


=cut

my %clone2acc;
my $ensembl_species=$ENV{ENSEMBL_SPECIES};
my $random;
my $csc;
my $csa;
my $registry_file;
    
{  #Argument Processing
  
  my $help=0;
  my $man=0;
  GetOptions( "help|?"    => \$help,
	      "man"       => \$man,
	      "species=s" => \$ensembl_species,
	      "random"    => \$random,
	      "csa=s"     => \$csa,
	      "csc=s"     => \$csc,
	      "registry_file=s" => \$registry_file,
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
$registry_file = "$ENV{GrameneEnsemblDir}/conf/ensembl.registry" 
  unless $registry_file;

my $reg = "Bio::EnsEMBL::Registry";	#Use this to get adaptors 
$reg->load_all( $registry_file ) or die "load_all failed";
     
my $slice_adaptor =$reg->get_adaptor($ensembl_species,'core','Slice') 
  or die "can't get Slice adaptor for $ensembl_species";

my $attribute_adaptor =$reg->get_adaptor($ensembl_species,'core','Attribute')
     or die "can't get Attribute adaptor for $ensembl_species";

my $coord_system_adaptor = $reg->get_adaptor($ensembl_species,'core'
                                     ,'CoordSystem');

my $chr_cs      = $coord_system_adaptor->fetch_by_name($csa);
my $scaffold_cs = $coord_system_adaptor->fetch_by_name($csc);
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
$p->comment( qr/^#/ ); 

my $j = 1; # as suffix of gap names GAP_UC_1
for my $file (@ARGV){
  
  $p->filename($file);
  $p->bind_fields( qw[ chr chr_start chr_end seq_order IsGap seq_name scaffold_start scaffold_end scaffold_strand ] );
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
  
     print STDERR "pseudo_len is $pseudo_len\n";
     my $pseudo_len_cal = 0;
     my ($pre, $gap)    = (0,0);
  
     for my $rd( @chr_array ){
	 print STDERR "DEBUG | $rd->{'chr_end'} -  $rd->{'chr_start'}\n";

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
     $chr_num    =~ s/.*chr(omosome)?_?0*//i;
     my $chr_slice = get_slice($coord_system_adaptor,
				$slice_adaptor,
				$chr_num,
				$pseudo_len,
			        $csa) 
				 or die "Can't get slice for $csa $chr_num";
  
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
       if( $r->{'IsGap'} !~ /^\s*[ADFGOPW]\s*$/i ){
	 next;
       }

       
       my ($scaffold_len, $scaffold_name);
       if ( $r->{'IsGap'}  !~ /^\s*[ADFGOPW]\s*$/i ){
	 $scaffold_len = $r->{'seq_name'};
	 $scaffold_name .=  $j++;
       }else{
	 $scaffold_name = $r->{'seq_name'};
	 $scaffold_len = $r->{'scaffold_end'} - $r->{'scaffold_start'} + 1;
       }

       

             
       my $scaffold_seq = undef;
       $scaffold_seq    = 'N' x $scaffold_len if ( $r->{'IsGap'}  !~ /^\s*[ADFGOPW]\s*$/i);

       my $scaffold_slice = get_slice($coord_system_adaptor,
				      $slice_adaptor,
				      substr($scaffold_name, 0, 40),
					# table collumn only has varchar(40) to store name
				      $scaffold_len,
				      $csc,
				      $scaffold_seq,
				     ) #$slice_adaptor->fetch_by_region('scaffold', $scaffold_name)
      or die "Can't get slice for $csc $scaffold_name";
    
       # set the default strand to be + if the orientation is undefined
       my ( $cmp_start, $cmp_end)=( $r->{'scaffold_start'}, 
					   $r->{'scaffold_end'},
				    );

       my $ori=$r->{'scaffold_strand'};
       if( !$ori || $ori eq '?'){
	   $ori = 1;
       }elsif( $ori eq '+' || $ori =~ /^plus$/i || $ori eq '?' ){
	   $ori = 1;
       }elsif( $ori eq '-' || $ori =~ /^minus$/i ){
	   $ori = -1;
       }else {
	   die "does not recognized orientation, $ori\n"; 
       }
       
    
       #map {print "$_ => $r->{$_}\n"} keys %$r;
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

  #print "in get_slice => cs=$cs, chr=$chr, len=$length\n";
  #$chr_slice=$sa->fetch_by_region('seqlevel', $chr) unless $chr_slice;

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

      warn("ERROR: $cs:$chr not found in the database\n");
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

