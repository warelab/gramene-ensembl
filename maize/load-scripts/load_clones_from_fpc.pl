#!/usr/local/bin/perl -w

=head1 NAME

load_clones_from_fpc.pl - Populates a core Ensembl DB with FPC clone data

=head1 SYNOPSIS

perl load_clones_from_fpc.pl [options] fpc_file

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -f --fasta_file
 -c --cb_to_bp
 -n --no_insert

=head1 OPTIONS

Reads the B<fpc_file>, and uses its clone data to load feature tables
(dna_align_feature at the moment). Should be run after the fpc assembly
has been loaded using the load_assembly_from_fpc.pl script.

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

B<-f --fasta_file>
  Fasta file (GenBank dump) containing sequences for all accessioned
  clones [REQUIRED]. TODO: Get sqeuences from GenBank at runtime using
  Entrez?

B<-c --cb_to_bp>
  Conversion factor to convert native (cb) fpc coordinated into basepairs.
  Calculated by hand by comparison with coordinates shown on AGI browser.
  Default: 4096.

B<-n --no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Populates a core Ensembl DB with FPC clone mappings.
  Could do with additional description!
  
  A good place to look for clone sequences is;
  /usr/local/data/fpc/<species>/<species>.fpc
  
Maintained by Will Spooner <whs@ebi.ac.uk>

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper); # For debug 

use DBI;
use FindBin qw( $Bin );
use File::Basename qw( dirname );

use vars qw( $BASEDIR );
BEGIN{
  # Set the perl libraries
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/bioperl-live';
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::MapIO;
use Bio::SeqIO;

use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::Map::MarkerFeature;

use vars qw( $I $ENS_DBA $FPC_MAP $SCALING_FACTOR $SEQ_IO );
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $file, $cb_to_bp, $fasta_file, $no_insert );
  GetOptions
      (
       "help|?"          => \$help,
       "man"             => \$man,
       "species=s"       => \$species,
       "registry_file=s" => \$file,
       "fasta_file=s"    => \$fasta_file,
       "cb_to_bp=s"      => \$cb_to_bp,
       "no_insert"       => \$no_insert,
       ) or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $I = $no_insert ? 0 : 1; # Put stuff in the database?

  # CB to BP scaling foctor. Calculated from AGI site. 
  # Can we get this programatically?
  $SCALING_FACTOR = $cb_to_bp || 4096;

  # Validate file paths
  $file    ||= $BASEDIR.'/conf/ensembl.registry';
  my $fpc_file = shift @ARGV;
  $fpc_file   || ( warn( "Need the path to an FPC file\n" ) && pod2usage(1));
  $fasta_file || ( warn( "Need a --fasta_file" ) && pod2usage(1));

  map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
  } $file, $fpc_file, $fasta_file;

  # Set the FASTA_FILE
  warn( "Found a FASTA file: $fasta_file. Loading...\n" );
  $SEQ_IO = Bio::SeqIO->new(-file=>$fasta_file, -format=>'fasta');

  # Load the ensembl file
  $species || ( warn( "Need a --species\n" ) && pod2usage(1) );
  Bio::EnsEMBL::Registry->load_all( $file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( warn( "No core DB for $species set in $file\n" ) &&
                pod2usage(1) );

  # Load the FPC file
  warn( "Found an FPC file: $fpc_file. Loading...\n" );
  my $mapio = new Bio::MapIO(-format  => "fpc",
                             -file    => "$fpc_file",
                             -readcor => 0,
                             -verbose => 0);
  $FPC_MAP = $mapio->next_map(); # Single map per FPC file
}

#warn Dumper( $fpcmap );
#&list_contigs_by_chromosome($fpcmap);
my $meta    = $ENS_DBA->get_MetaContainer();
my $species = $meta->get_Species || 
    die( "Cannot find the species in the meta table of the DB" );
my $common_name = $species->common_name ||
    die( "Cannot find the species common name in the meta table of the DB" );

$common_name = ucfirst( $common_name );

###########
# Prepare some Ensembl adaptors
#my $analysis = &fetch_analysis( $common_name."_BAC_FPC" );
my $marker_analysis = &fetch_analysis( $common_name."_marker" );
my $sl_adapt = $ENS_DBA->get_adaptor('Slice');
my $f_adapt  = $ENS_DBA->get_adaptor('SimpleFeature');
my $mf_adapt = $ENS_DBA->get_adaptor('MiscFeature');


###########
# Process the BAC sequence file from GenBank to map clone names to
# accessions
my %accessions_by_name;
my %accessioned_length;
my %all_accessions; # Keep tally for reporting
warn( "Processing fasta file containing accessioned clones" );
while ( my $seq = $SEQ_IO->next_seq() ) {
  my $name = $seq->display_name;
  my $description = $seq->description;
  my( $accession, $version, $clone_name );
  if( $name =~ m/^gi\|\d+\|gb\|(\w+)\.(\d+)\|.*$/ ){ 
    # Parse genbank header
    $accession = $1;
    $version   = $2;
    $all_accessions{$accession} = 0;
    $accessioned_length{$accession} = $seq->length;
    $accessions_by_name{$accession} = $accession;
    $accessions_by_name{$version}   = $accession;
  } else { 
    warn( "Cannot determine accession from $name" );
    next;
  }
  foreach my $clone ( $description =~ /clone ([A-Z]+(\w+))/ ){
    # Clone names can be prefixed with library. Try both with and without
    $accessions_by_name{$1} = $accession;
    $accessions_by_name{$2} = $accession;
  }
}

###########
# Establish the different classes of feature. Similar to analysis
my $misc_set_fpc = Bio::EnsEMBL::MiscSet->new
    (-code        => 'superctgs', # Must conform to Ensembl for display
     -name        => 'FPC Contig',
     -description => '',
     -longest_feature => 8000000 );

my $misc_set_bac = Bio::EnsEMBL::MiscSet->new
    (-code        => 'bac_map', # Must conform to Ensembl for display
     -name        => 'BAC map',
     -description => 'Full list of FPC BAC clones',
     -longest_feature => 250000 );

my $misc_set_accbac =  Bio::EnsEMBL::MiscSet->new
    (-code        => 'acc_bac_map', # Must conform to Ensembl for display
     -name        => 'Accessioned BAC map',
     -description => 'List of mapped and accessioned BAC clones',
     -longest_feature => 250000 );

###########
# Loop through each FPC contig
foreach my $ctgid( sort $FPC_MAP->each_contigid ){
  $ctgid || next; # Skip ctg 0
#  last;
  # Create a slice from the FPC contig
  my $ctg_name = "ctg$ctgid";
  my $ctg_slice = $sl_adapt->fetch_by_region(undef,$ctg_name);
  $ctg_slice || die( "FPC Contig $ctg_name was not found in the Ensembl DB" );
  
  warn( "Processing $ctg_name at ",$ctg_slice->name );

  # Load the FPC as a misc feature
  my $fpc_feature = Bio::EnsEMBL::MiscFeature->new
      (
       -start  => $ctg_slice->start,
       -end    => $ctg_slice->end,
       -strand => 1,
       -slice  => $ctg_slice,
       );
  my $fpcname_attrib = Bio::EnsEMBL::Attribute->new
      (
       -VALUE => $ctg_slice->seq_region_name,
       -CODE  => 'name',
       -NAME  => 'name'
       );
  
  $fpc_feature->add_MiscSet($misc_set_fpc);
  $fpc_feature->add_Attribute($fpcname_attrib);
  $fpc_feature = $fpc_feature->transform('chromosome');
  $mf_adapt->store( $fpc_feature ) if $I;

  # Process each clone on FPC
  my $ctg = $FPC_MAP->get_contigobj( $ctgid );
  my @clones = $ctg->each_cloneid;
  warn( "  Num clones on $ctg_name: ",scalar(@clones) );
  foreach my $cloneid( sort $ctg->each_cloneid ){
    my $clone = $FPC_MAP->get_cloneobj( $cloneid );
    #my $clone_tmpl = "    Clone %s at Ctg %s:%s-%s\n";
    #printf( $clone_tmpl, $clone->name, $clone->contigid, 
    #        $clone->range->start,  $clone->range->end );

    my $start = $SCALING_FACTOR * $clone->range->start;
    my $end   = $SCALING_FACTOR * $clone->range->end;
    my $clone_name = $clone->name;
    my $feature = Bio::EnsEMBL::MiscFeature->new
        (
         -start  => $start,
         -end    => $end,
         -strand => 1,
         -slice  => $ctg_slice,
         );

    # Populate the feature with a variety of attributes
    $feature->add_MiscSet($misc_set_bac);
    $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                             ( -VALUE => $clone_name,
                               -CODE  => 'name',
                               -NAME  => 'name' ) );
    $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                             ( -VALUE => $ctg_name,
                               -CODE  => 'superctg',
                               -NAME  => 'FPC contig name' ) );
    $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                             ( -VALUE => $start,
                               -CODE  => 'inner_start',
                               -NAME  => 'Max start value' ) );
    $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                             ( -VALUE => $end,
                               -CODE  => 'inner_end',
                               -NAME  => 'Min end value' ) );
    $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                             ( -VALUE => $end-$start+1,
                               -CODE  => 'fp_size',
                               -NAME  => 'FP size' ) );

    if( my $acc = $accessions_by_name{$clone_name} ){
      # This is an accessioned clone
      $all_accessions{$acc} ++;     
      $feature->add_MiscSet($misc_set_accbac);
      $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                               ( -CODE  => 'state',
                                 -NAME  => 'Current state of clone',
                                 -VALUE => '12:Accessioned' ) );
      $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                               ( -CODE  => 'embl_acc',
                                 -NAME  => 'Accession number',
                                 -VALUE => "$acc" ) );
      $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                               ( -CODE  => 'seq_len',
                                 -NAME  => 'Accession length',
                                 -VALUE => $accessioned_length{$acc} ) );
    } else {
      # This is a free-state clone
      $feature->add_Attribute( Bio::EnsEMBL::Attribute->new
                               ( -CODE  => 'state',
                                 -NAME  => 'Current state of clone',
                                 -VALUE => '06:Free' ) );
    }

    # Store
    $feature = $feature->transform('chromosome');
    $mf_adapt->store( $feature ) if $I;
  }
}

###########
# Loop through each marker. 
my $marker_feature_adaptor = $ENS_DBA->get_adaptor('MarkerFeature');
my $marker_adaptor         = $ENS_DBA->get_adaptor('Marker');
my @markers = $FPC_MAP->each_markerid();
warn( "  Num markers: ",scalar(@markers) );
foreach my $marker( @markers ){
#  last;
  my $markerobj = $FPC_MAP->get_markerobj($marker);
  my $type = $markerobj->type; # STS, eMRK, Probe, OVERGO
  
  my $ensmarker = Bio::EnsEMBL::Map::Marker->new();
  $ensmarker->display_MarkerSynonym
      ( Bio::EnsEMBL::Map::MarkerSynonym->new(undef,$type,$marker) );
  $ensmarker->priority(100); #Ensure marker displays

  # get all the contigs where this marker hit
  my @contigs = $markerobj->each_contigid();
  foreach my $ctgid( @contigs ) {
    $ctgid || next; # Skip empty contig

    # Create a slice from the FPC contig
    my $ctg_name = "ctg$ctgid";
    my $ctg_slice = $sl_adapt->fetch_by_region(undef,$ctg_name);
    $ctg_slice || die( "FPC $ctg_name was not found in the Ensembl DB" );
    
    # Create the marker feature
    my $maploc = $markerobj->position($ctgid);
    my $marker_feature = Bio::EnsEMBL::Map::MarkerFeature->new
        (
         undef, #DBid
         undef, #MarkerFeatureAdaptor
         $maploc * $SCALING_FACTOR,   # Start
         ($maploc+1) * $SCALING_FACTOR, # End - assume length=SCALING_FACTOR
         $ctg_slice,         # Slice
         $marker_analysis,   # Analysis
         undef,              # Marker ID?
         scalar( @contigs ), # map weight
         $ensmarker,         # Ensembl Marker
         );
    $marker_feature = $marker_feature->transform('chromosome') ||
        ( warn("No chr mapping for $marker!") && next );
    $marker_feature_adaptor->store( $marker_feature ) if $I;
  }
}

warn( "Found ", scalar( grep{$_} values %all_accessions ), 
      " accessioned clones in FPC\n" );
warn( "Lost ", scalar( grep{! $_} values %all_accessions ),
      " accessioned clones from FPC\n" );

exit;

#======================================================================
sub list_contigs_by_chromosome{
  my $FPC_MAP = shift;

  my %chr_data;
  my %lengths;
  my $scaler = $SCALING_FACTOR || 4096;
  foreach my $ctgname( $FPC_MAP->each_contigid ){
    $ctgname || next; # Skip empty
    my $ctgobj   = $FPC_MAP->get_contigobj($ctgname); 
    my $ctgpos   = $ctgobj->position+0   || 0;
    my $ctgstart = $ctgobj->range->start || 0;
    my $ctgend   = $ctgobj->range->end   || 0;
    my $ctgbasepair = $ctgend * $scaler;

    my $chr = $ctgobj->group || '';
    $chr_data{$chr} ||= [];
    push @{$chr_data{$chr}}, [ $ctgname, 
                               sprintf( "%.2f", $ctgpos ), 
                               #$ctgstart, 
                               $ctgend ,
                               $ctgbasepair,
                               ];
    my $l_n = length($ctgname);
    my $l_p = length($ctgpos);
    my $l_s = length($ctgstart);
    my $l_e = length($ctgend);
    my $l_b = length($ctgbasepair);

    if( $l_n > ($lengths{name} ||0) ){ $lengths{name} = $l_n }
    if( $l_p > ($lengths{pos}  ||0) ){ $lengths{pos}  = $l_p }
    if( $l_s > ($lengths{start}||0) ){ $lengths{start}= $l_s }
    if( $l_e > ($lengths{end}  ||0) ){ $lengths{end}  = $l_e }
    if( $l_b > ($lengths{bp}   ||0) ){ $lengths{bp}   = $l_b }
  }
  
  my $num_per_line = 3;
  $lengths{name} += 3; # Allow for ctg prefix
  my $entry_t = join( " ", 
                      "|", 
                      "%-${lengths{name}}s",
                      "%${lengths{pos}}s",
                      #"%${lengths{start}}s",
                      "%${lengths{end}}s",
                      "%${lengths{bp}}s",
                      "|" );

  my $tmpl = "  ". ( $entry_t x $num_per_line ) . "\n"; 

  # Loop through each chromosome and print FPC data
  foreach my $chr( sort{$a<=>$b} keys %chr_data ){
    my @ctgs = sort{$a->[0]<=>$b->[0]} @{$chr_data{$chr}};
    map{ $_->[0] = "Ctg".$_->[0] } @ctgs;

    my $total_length = 0;
    map{ $total_length += $_->[3] } @ctgs;
    my $num_ctgs = scalar( @ctgs );
    print( "Chr: $chr ($num_ctgs FPCs, $total_length bp)\n" );

    while( my @entries = splice( @ctgs, 0, $num_per_line ) ){
      my @data = ( map{defined($_) ? $_ : ''} 
                   map{@{$entries[$_-1]||[]}[0..3]} 
                   (1..$num_per_line) );
      #warn Dumper( @data );
      printf( $tmpl, @data );
    }
  }
}

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
  # None needed for now
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
  $adaptor->store($analysis) if $I;
  return $analysis;
}

#======================================================================
1;
