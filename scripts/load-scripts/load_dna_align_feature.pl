#!/lab/bin/perl -w

BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/'; 
}

# The first shall be last...
#use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

#use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use lib map { "$ENV{'GrameneEnsemblDir'}/ensembl-live/$_" } 
        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
             ensembl-draw/modules ensembl-compara/modules );

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;

#use Gramene::Config;
use DBI;
use DBI qw(:sql_types);

use Data::Dumper qw(Dumper); # For debug

use FindBin qw( $Bin );
use File::Basename qw( dirname );


# Import EnsEMBL modules
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

=head1 SYNOPSIS

load_dna_align_feature.pl  [options] 
 
 Options:
    --species		species in EnsEMBL registry to use for db [required]
    --v 		makes more verbose
    --help		help message
    --man		full documentation
    --coord_system	coord system to use
    --registry_file	Default is $GrameneEnsemblDir/conf/ensembl.registry
    --file_prefix	Remove from filename before interpreting as logic name (default=empty)
    --logic_name	Allowed logic name
    --replace		Remove all hits for each track being loaded


=head1 OPTIONS

=over 4

=item B<--species> 

A species in the EnsEMBL registry

=item B<--coord_system>

Coordinate system name (e.g. 'clone', 'chromosome' )

Defaults to the sequence level coordinate system 

Can be over-ridden by data file: if region name is chr_NUMBER
perhaps followed by _NUMBER then chromosome coordinates are used
and the 1st NUMBER is the chromosome name.

=item B<--logic_name>

Allowed logic name (can be repeated).

If none given then all are allowed.

=item B<--replace>

Before loading a track, remove all existing hits

If you manage to load several files into the same track,
the remove happens only once.


=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

files to process (produced by Gramene Pipeline)

filename is prefix+logic_name+'.whatever'

=cut

use vars qw($ENS_DBA $species @logic_names $file_prefix $coord_system
		$replace $verbose);
$verbose=0;
$file_prefix='';
{  #Argument Processing
  my $help=0;
  my $man=0;
  my $registry_file;
  GetOptions
      ( 
        "help|?"          => \$help,
        "man"             => \$man,
        "species=s"       => \$species,
	"coord_system=s"  => \$coord_system,
        "registry_file=s" => \$registry_file,
        "file_prefix=s"   => \$file_prefix,
        "logic_name=s"    => \@logic_names,
        "replace"         => \$replace,
        "v+"		  => \$verbose,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Validate the input files
  $registry_file    ||= $ENV{GrameneEnsemblDir}.'/conf/ensembl.registry';
  @ARGV || ( warn( "Need some input files to process\n" ) && pod2usage(1) );
  map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
  } $registry_file, @ARGV;
  warn( "Found ",scalar @ARGV," files to process...\n" );

  # Load the ensembl file
  $species || ( warn( "Need a --species\n" ) && pod2usage(1) );
  Bio::EnsEMBL::Registry->load_all( $registry_file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( warn( "No core DB for $species set in $registry_file\n" ) &&
                pod2usage(1) );

}
my %logic_names = map{$_=>1} @logic_names;

my $slice_adaptor    = $ENS_DBA->get_adaptor('Slice');
my $analysis_adaptor = $ENS_DBA->get_adaptor('Analysis');
my $feature_adaptor  = $ENS_DBA->get_adaptor('DnaAlignFeature');

unless ($coord_system) {
    my $coord_system_adaptor=  $ENS_DBA->get_adaptor('CoordSystem');
    $coord_system=$coord_system_adaptor->fetch_sequence_level->name;
    warn("$coord_system coordinates\n");
}


my %replaced;

foreach my $file( @ARGV ){
  warn( "Processing $file\n" );
  my( $logic_name ) = $file =~ m|([^/]+?)(_\d*)?\..*$|;
  warn ($logic_name) if $verbose;
  if( $file_prefix ){ 
      warn "removing $file_prefix" if $verbose;
      $logic_name =~ s/$file_prefix// 
  }
  if( %logic_names and ! $logic_names{$logic_name} ){
    warn( "  Skipping analysis logic_name: $logic_name\n" );
  }

  my $analysis = &fetch_analysis($logic_name); 
  if($replace && !$replaced{$logic_name}++ && $analysis->dbID) {
      $feature_adaptor->db->dbc->do(
       "delete from dna_align_feature where analysis_id=".$analysis->dbID)
       or die("delete $logic_name = ".$analysis->dbID.": $DBI::errstr");
  } 
  
  open( FEATURES, $file ) || die ("Cannot open $file: $!");
  while( my $line = <FEATURES> ){
    next if ($line =~ /^HIT_ID/);
    chomp $line;
    my( $hit_num, 
        $hit_name, 
        $seq_region_name_with_offset, 
        $seq_region_strand, 
        $hit_start, 
        $hit_end, 
	$seq_region_name,  
        $seq_region_start, 
        $seq_region_end ,
	$score,
	$perc_id,
	$cigar_string,
	$evalue );

    my $this_coord_system=$coord_system;

    if( $file_prefix eq 'final_filtered_offset_' ) {
	( $hit_num, 
	  $hit_name, 
	  $seq_region_name_with_offset, 
	  $seq_region_strand, 
	  $hit_start, 
	  $hit_end, 
	  $seq_region_name,  #not in 'final_filtered_' format
	  $seq_region_start, 
	  $seq_region_end ,
	  $score,
	  $perc_id,
	  $cigar_string,
	  $evalue )
	    =  split( /\t/, $line );
    } else {	# final_filtered_
	( $hit_num, 
	  $hit_name, 
	  $seq_region_name_with_offset, #or just clone name
	  $seq_region_strand, 
	  $hit_start, 
	  $hit_end, 
	  $seq_region_start, 
	  $seq_region_end ,
	  $score,
	  $perc_id,
	  $cigar_string,
	  $evalue )   #evalue has bogus value
	    =  split( /\t/, $line );
	  $seq_region_name=$seq_region_name_with_offset;
    }
    #It's stupid hack time again down at the Bar-G Corral
    if($seq_region_name_with_offset =~ /^chr_(\d+)(.*)$/) {
	$seq_region_name=$1;
	$this_coord_system='chromosome';
    }
    warn "$this_coord_system:$seq_region_name\n" if $verbose;

    my $slice = $slice_adaptor->fetch_by_region($this_coord_system,$seq_region_name) ||
        ( warn( "  Cannot fetch $this_coord_system $seq_region_name - skipping\n" ) && next );

    #my $cigar = ( $seq_region_end-$seq_region_start+1 ) . "M";

    $hit_name =~ s/\|.*//;
    my $feature = Bio::EnsEMBL::DnaDnaAlignFeature->new
        ( 
          -start        => $seq_region_start,
          -end          => $seq_region_end,
          -strand       => $seq_region_strand,
          -hstart       => $hit_start,
          -hend         => $hit_end,
          -hstrand      => 1,
          -hseqname     => $hit_name,
          -analysis     => $analysis,
          -slice        => $slice,
          -cigar_string => $cigar_string,
  #        -p_value 	=> $evalue,
          -score        => $score,
          -percent_id        => $perc_id, 
          );

    $feature_adaptor->store($feature);
    warn( "    Created feature $hit_name at ",
          $feature->feature_Slice->name , "\n") if $verbose;

  }


}

exit;

    
# end of main program

########################## subroutines ######################################

#======================================================================
# Returns an Analysis object; either fetched from the DB if one exists,
# or created fresh, in which case it is stored to the DB.
sub fetch_analysis{
  my $logic_name = shift || die("Need a logic_name" );
  my $db_file    = shift || '';
  my $adaptor = $ENS_DBA->get_adaptor('Analysis'); 

  my %args = ( -logic_name=>$logic_name, 
               $db_file ? (-db_file=>$db_file) : () );

  # Hard-coded nastyness to make analyses correspond to Ensembl
  if( $logic_name eq 'MyAnalysis' ){
    $args{-logic_name} = 'MyLogic';
    $args{-db}         = 'MyDB';
    $args{-program}    = 'MyProgram';
  }

  my $analysis;
  if( $analysis = $adaptor->fetch_by_logic_name($args{-logic_name}) ){
    # Analysis found in database already; use this.
    return $analysis;
  }

  # No analysis - create one from scratch
  $analysis = Bio::EnsEMBL::Analysis->new(%args);
  warn("created analysis for $logic_name.\n");
  return $analysis;
}

#======================================================================
1;

__END__

=head1 OUTPUT

Explanatory only

=head1 NOTES
    
Creates entry in analysis table if doesn't exist


=head1 AUTHOR

   Will Spooner, Kiran Kumar, Steven Schmidt
   Gramene Project (www.gramene.org)
   Stein Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

