#!/usr/local/bin/perl -w 

=pod

=head1 NAME

convert_coordinates.pl - convert the coordinates to a different coordinate system. 

=head1 SYNOPSIS

  convert_coordinates.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -from_cs|--from_coord_system       the coord system the features are on in the original input file
  -to_cs|--to_coord_system        the coord system the features are to be converted to

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s|--species>
  Use this species entry from the registry file [REQUIRED].

B<-l|--from_cs>
	the coord system the features are on in the original input file, such as contig or clone

B<-n|--to_cs>
  the coord system the features are to be converted to, such as chromosome or toplevel
  


=head1 DESCRIPTION


  example of the gff file:
scaffold_1      Ptrichocarpav2_0        gene    12632   13612   .       +       .       ID=POPTR_0001s00200;Name=POPTR_0001s00200
scaffold_1      Ptrichocarpav2_0        mRNA    12632   13612   .       +       .       ID=POPTR_0001s00200.1;Name=POPTR_0001s00200.1;PACid=17327967;Parent=POPTR_0001s00200
scaffold_1      Ptrichocarpav2_0        exon    12632   12650   .       +       .       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        5'-UTR  12632   12638   .       +       .       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        CDS     12639   12650   .       +       0       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        exon    12768   12891   .       +       .       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        CDS     12768   12891   .       +       0       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        exon    13117   13226   .       +       .       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        CDS     13117   13226   .       +       1       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        exon    13310   13612   .       +       .       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        3'-UTR  13385   13612   .       +       .       Parent=POPTR_0001s00200.1;PACid=17327967
scaffold_1      Ptrichocarpav2_0        CDS     13310   13384   .       +       0       Parent=POPTR_0001s00200.1;PACid=17327967


B<The Ensembl Registry>

  The database connection details for both Ensembl and interproscan
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'maize', 
    '-group'   => 'core', 
    '-dbname'  => 'zea_mays_core_30_bac20', );
  ---


Maintained by Sharon Wei <weix@cshl.edu>

=cut


use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Readonly;
use List::Util qw( first );
use List::MoreUtils;


use vars qw( $BASEDIR );

BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;
use Date::Calc;
use Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor;

use vars qw( $ENS_DBA $SLICE_ADAPTOR $AS_MAPPER_ADAPTOR  $GFF_HANDLE $FROM_CS $TO_CS $CS_ADAPTOR );

my $date = time(); 

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $from_cs, $to_cs );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "from_cs=s"       => \$from_cs,
        "to_cs=s"      => \$to_cs,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args

  my $gff_file = $ARGV[0] || pod2usage("\nNeed the path to a gff file\n");
  $species    || pod2usage("\nNeed a --species\n");
  $from_cs    || pod2usage("\nNeed a --from coordinate system name such as contig\n");
  $to_cs    || pod2usage("\nNeed a --to coordinate system name such as chromosome\n");
  $reg        ||= $BASEDIR.'/conf/ensembl.registry';

  map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
  } $reg, $gff_file;

  $FROM_CS = $from_cs;
  $TO_CS  = $to_cs;


  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) &&
                pod2usage() );

  $AS_MAPPER_ADAPTOR = Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "assemblymapper" );

  $CS_ADAPTOR = Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "coordsystem" );
 
  $SLICE_ADAPTOR = Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "slice" );
  
  # Create a GFF stream
  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
  }
}

my $from_cs_obj = $CS_ADAPTOR->fetch_by_name($FROM_CS);
my $to_cs_obj = $CS_ADAPTOR->fetch_by_name($TO_CS);
my $asm_mapper = $AS_MAPPER_ADAPTOR->fetch_by_CoordSystems( $from_cs_obj, $to_cs_obj );

while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  print && next if ( $line =~ /\#/ ); 
  print && next if ( $line =~ /^\s+/ );
  #chomp $line;

  #print "$line";  #####
  #warn(Dumper($TRPT2GENE));

  # Split gff line,\
  # start always <= end even for - strand genes
  my( $seqname, $source, $feature, 
      $start, $end, $score, $strand, $frame, 
      $attribute ) = split( /\s+/, $line, 9 );
  $feature = uc($feature);
  $seqname =~ s/chr0*//i;
  $strand = $strand eq '-' ? -1 : 1;

  my $from_slice = $SLICE_ADAPTOR->fetch_by_region( $FROM_CS, $seqname );

  if ( !$from_slice->is_toplevel ){
       my @to_coords = $asm_mapper->map( $seqname, $start, $end, $strand,
                                  $from_cs_obj );
	die "should not map to more than one segments" if @to_coords > 1; #we only deal with contig to chr now

	my $to_coord = pop @to_coords;
  	my $to_seqregion_id   = $to_coord->id;
        my $to_seqregion_slice = $SLICE_ADAPTOR->fetch_by_seq_region_id ( $to_seqregion_id );

        my $to_sr_name = $to_seqregion_slice->seq_region_name();
	my $to_sr_start = $to_coord->start;
        my $to_sr_end  = $to_coord->end;
        my $to_sr_strand = $to_coord->strand == -1 ? '-' : '+';

#print "converted $FROM_CS:$seqname:$start:$end:$strand to $TO_CS:$to_sr_name:$to_sr_start:$to_sr_end:$to_sr_strand\n";  
	print join "\t", ($to_sr_name, $source, $feature, $to_sr_start, $to_sr_end, $score, $to_sr_strand, $frame, $attribute); 
 }else{
	print "$line"; 
 }
  
}
