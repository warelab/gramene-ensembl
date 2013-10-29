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
    --acs               assembly coordinate system
    --ccs               optional, major component coordinate system, sometimes there are two or more component coordinate system
    --forgiving          skip the component region if no component seq_region found

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
    
=item B<--acs> 

assembly coordinate system, 

=item B<--ccs>

Optional, major component system. In most cases, there is only one component coordinate system. In rare cases, there are two or more component coordinate system (For example: chr assembled from both scaffolds and contigs), if the component seq_region_name cannot be found in the major cs, the cs will be decided by the highest ranking cs that the component seq_region_name can be found. For example, if a component seq_region_name cannot be found on the ccs, which is superscaffold, but it can be found on both scaffold and contig, scaffold will be used as the true cs for that component seq_region. Most of the time, the seq_region_names are unique globally across all coordinate systems
 


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

Field5 : Gap/non gap, possible values W, N, U 

Field6 : component seq name or number of Ns for gap region

Field7 : component seq start or fragment/contig for gap region

Field8 : component seq end or yes/no for gap region

Field9 : component seq strand


=cut

my %clone2acc;
my $ensembl_species=$ENV{ENSEMBL_SPECIES};
my $random;
my $ccs;
my $acs;
my $registry_file;
my $forgiving;
    
{  #Argument Processing
  
  my $help=0;
  my $man=0;
  GetOptions( "help|?"    => \$help,
	      "man"       => \$man,
	      "species=s" => \$ensembl_species,
	      "random"    => \$random,
	      "acs=s"     => \$acs,
	      "ccs=s"     => \$ccs,
	      "registry_file=s" => \$registry_file,
	      "forgiving" => \$forgiving,
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $ccs = uc $ccs;
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

my $coord_system_adaptor = $reg->get_adaptor($ensembl_species,'core'
                                     ,'CoordSystem');

# global vairable %component_css
our %component_css;
$component_css{$ccs} = 1;



my $dba=$slice_adaptor->db; #DBAdaptor
my $dbc=$dba->dbc;	#DBConnection
warn "user ".$dbc->username.", db ".$dbc->dbname."\n";


my $assembly_sth = $dbc->prepare("insert into assembly 
      (asm_seq_region_id ,cmp_seq_region_id ,asm_start ,asm_end ,cmp_start ,cmp_end ,ori)
      values  (?, ?,?,?,?,?,?)");

my $p = new Text::RecordParser::Tab;

my $j = 1; # as suffix of gap names GAP_UC_1
for my $file (@ARGV){
  
  $p->filename($file);
  $p->bind_fields( qw[ chr chr_start chr_end seq_order IsGap seq_name scaffold_start scaffold_end scaffold_strand ] );

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
     my $chr_slice = get_slice(
				$slice_adaptor,
				$chr_num,
			        $acs) 
				 or die "Can't get slice for $acs $chr_num";
  
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
       if( $r->{'IsGap'} ne 'W' ){
	 next;
       }

       
       my ($scaffold_len, $scaffold_name);
       if ( $r->{'IsGap'} ne 'W' ){
	 $scaffold_len = $r->{'seq_name'};
	 $scaffold_name .=  $j++;
       }else{
	 $scaffold_name = $r->{'seq_name'};
	 $scaffold_len = $r->{'scaffold_end'} - $r->{'scaffold_start'} + 1;
       }

       

       #Get component slice, no need to create a false slice for gaps      
       #my $scaffold_seq = undef;
       #$scaffold_seq    = 'N' x $scaffold_len if ( $r->{'IsGap'} ne 'W' );

       my $scaffold_slice;

	$scaffold_slice = get_slice(
				      $slice_adaptor,
				      $scaffold_name,
				      $ccs
				      ); 

	if( ! $scaffold_slice ){ 
	    die "Can't get slice for $ccs $scaffold_name" unless $forgiving;
	    next;
	}
    
       # set the default strand to be + if the orientation is undefined
       my ( $cmp_start, $cmp_end)=( $r->{'scaffold_start'}, 
					   $r->{'scaffold_end'},
				    );

       my $ori=$r->{'scaffold_strand'};
       #according to AGP specification, 
       #http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
       #unknown orientation usually is treated as if they had + orientation
      
       if( $ori eq '+' || $ori =~ /^\s*plus\s*$/i ){
	   $ori = 1;
       }elsif( $ori eq '-' || $ori =~ /^\s*minus\s*$/i ){
	   $ori = -1;
       }else {
	   $ori = 1;
	   warn "does not recognized orientation, $ori, treated as +\n"; 
       }
       
    
       #map {print "$_ => $r->{$_}\n"} keys %$r;
       printf "chr_seq_region_id = %s, scaffold_seq_region_id = %s\n",
	 $chr_slice->get_seq_region_id, $scaffold_slice->get_seq_region_id;
       $assembly_sth->execute($chr_slice->get_seq_region_id, 
			      $scaffold_slice->get_seq_region_id,
			      $r->{'chr_start'}, $r->{'chr_end'},
			      $cmp_start ,$cmp_end ,$ori);

     }
  
   }			

  
}


#store assembly.mapping 
#
my $chr_cs      = $coord_system_adaptor->fetch_by_name($acs);

for my $cs (keys %component_css){
    
    my $comp_cs = $coord_system_adaptor->fetch_by_name($cs);
    $coord_system_adaptor->store_mapping_path($chr_cs, $comp_cs );

}

exit;
    
# end of main program

########################## subroutines ######################################
sub get_slice { #add if necessary

  my ($sa,$chr,$cs)=@_;
  
  my $chr_slice;

  eval{
     $chr_slice=$sa->fetch_by_region($cs,$chr) ;
   };

  if( $@ ){ print "sequence $chr cannot be found on coordinate system $cs"; }

  #if the component seq_region_name cannot be found in the specified component cs
  #try to guess the most likely component cs
  unless( $chr_slice) {
      
      $chr_slice = $sa->fetch_by_region( undef, $chr) or return undef;
      my $cs_name = $chr_slice->coord_system_name;
      $component_css{ uc $cs_name } = 1; #add to the glabal registry for component coord_system
  }
  return $chr_slice;
}


# obsolete
sub get_scaffold_slice { #add if necessary
  my ($acs,$sa,$chr,$length)=@_;
  my $chr_slice=$sa->fetch_by_region('chromosome',$chr) ;
  unless( $chr_slice) {
    
    my $chrcs = $acs->fetch_by_name('chromosome');
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

