#!/usr/local/bin/perl -w 

=pod

=head1 NAME

load_genes_from_grape_gff.pl - load grape genes into EnsemblDB

=head1 SYNOPSIS

  load_genes_from_grape_gff.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -l|--logic_name       the logic name for the gene track
  -n|--no_insert        Do not make changes to the database. For debug.

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

B<-l|--logic_name>
  the logic name for these set of genes, the track name shown on the contigview, default = ensembl


B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.


=head1 DESCRIPTION

load repeat features from the gff files



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

B<Restoring the database>

  If the script bombs out half-way through, your database will be in a
  partial loaded state (i.e. in a bit of a mess). Here is some SQL to
  help put it right;

  delete rf.*, rc.* from repeat_feature rf join repeat_consensus rc using (repeat_consensus_id); 
  
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

use vars qw( $BASEDIR );
BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
  #unshift @INC, $BASEDIR.'/bioperl-live';
  unshift @INC, '/usr/local/ensembl-live/ensembl/modules';
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::RepeatConsensus;

use Date::Calc;


use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS );

my $date = time(); ;

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species,  $reg, $logic_name, $no_insert );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        
        "ensembl_registry=s" => \$reg,
        "logic_name=s"       => \$logic_name,
        "no_insert"          => \$no_insert,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args

  my $gff_file = $ARGV[0] || pod2usage("\nNeed the path to a gff file\n");
  $species    || pod2usage("\nNeed a --species\n");
  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  $logic_name ||= 'ensembl';

  map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
  } $reg, $gff_file;

  # Put stuff in the database?
  $I= $no_insert ? 0 : 1; 

  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) &&
                pod2usage() );

  # Create the analysis
  $ANALYSIS = Bio::EnsEMBL::Analysis->new
      ( -logic_name => $logic_name ); #'GeneModel_RiceIndica_BGI' 
  my $analAdaptor = $ENS_DBA->get_adaptor('Analysis');
  $analAdaptor->store ($ANALYSIS);

  # Create a GFF stream
  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


  print( "Loading repeatFeature for $species\n" );
  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
    $pre_text = "  Old DB: ";
  }
}

my $sa  = $ENS_DBA->get_adaptor('Slice');
my $rfa = $ENS_DBA->get_adaptor('RepeatFeature');
my $rca = $ENS_DBA->get_adaptor('RepeatConsensus');
my $count=0;

while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /\#/ ); 
  next if ( $line =~ /^\s+/ );
  chomp $line;

  #print "$line\n";  #####
  #warn(Dumper($TRPT2GENE));

  # Split gff line,\
  # start always <= end even for - strand genes
  my( $seqname, $source, $feature, 
      $start, $end, $score, $strand, $frame, 
      $attribute ) = split( /\s+/, $line, 9 );
  $feature = uc($feature);
  $seqname =~ s/chr(omosome)?_?//i;

  my %attribs;
  foreach my $id( split( /;/, $attribute ) ){

    if( $id  =~ m/([^=]*)=(.*)/ ){
      #ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
      #ID=11562.t00100;Name="ORF249";Note="pub_locus=LOC_Osp1g01070"

      my $key   = uc $1;
      my $value = $2;
      $value =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
      $value =~ s/^"//;
      $value =~ s/"$//;
      $attribs{$key} = $value;
      
    }
  }


  
  if ( $feature eq 'MASK_REGION'){    
  
    my $rc = Bio::EnsEMBL::RepeatConsensus->new
      (
       -NAME => $source,
       -REPEAT_CLASS => 'JGI-Sbi1.4',
       -ADAPTOR => $rca,
      );
    
    my $slice = $sa->fetch_by_region( undef, $seqname ) or 
      die "cannot fetch_by_region for $seqname\n";
    
    my $repeat_feature = Bio::EnsEMBL::RepeatFeature->new 
      (
       -start   => $start,
       -end     => $end,
       -strand  => 1,
       -slice   => $slice,
       -analysis => $ANALYSIS,
       -repeat_consensus => $rc,
       -hstart  => 1,
       -hend    => $end-$start+1,
       -hstrand => 1,
       #-score => 83.2
       
      );

    $rfa->store( $repeat_feature ) if $I;

    $count++;
    print "stored $count repeat features: $start, $end, $attribs{ID}, $attribs{MASK_LEN}\n";
  }
  else{
    die( "Unrecognized feature $feature" );
  }

}


