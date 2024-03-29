#! /usr/local/bin/perl -w 

=pod

=head1 NAME

load_repeats_from_gff_agi.pl - load AGI repeats from gff into EnsemblDB

=head1 SYNOPSIS

  load_repeats_from_gff_agi.pl [options] gff_file

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

load repeat features from the gff files similar to the following, it is generated by AGI for Oryza glaberrima whole genome

##gff-version 3
Oglab01_0232    agi_genomes_db  LTR_retrotransposon     80957   92887   .       +       .       ID=1;rclass=I;subclass=I;rorder=LTR;superfamily=Gypsy;family=RIRE2;is_complete=true;pbs_x=81402;pbs_y=81416;ppt_x=92365;ppt_y=92379;ltr_comparison_insertion_time=0.89;tsd_sequence=ATATA;presence_in_sativa=present
Oglab01_0232    agi_genomes_db  LTR_retrotransposon     112336  125583  .       +       .       ID=2;rclass=I;subclass=I;rorder=LTR;superfamily=Gypsy;family=ATLANTYS;pbs_x=113246;pbs_y=113264;ppt_x=124637;ppt_y=124651
Oglab01_0232    agi_genomes_db  LTR_retrotransposon     125593  134413  .       +       .       ID=3;rclass=I;subclass=I;rorder=LTR;superfamily=Gypsy;family=RIRE3;pbs_x=126220;pbs_y=126239;ppt_x=133703;ppt_y=133717
Oglab01_0232    agi_genomes_db  LTR_retrotransposon     126345  127017  .       -       .       ID=6;rclass=I;subclass=I;rorder=LTR;superfamily=Gypsy;family=RIRE3;tsd_sequence=unknown

repeat_name= "LTR_retrotransposon-1" (feature-ID)
repeat_class= "I/I/LTR/Gypsy/RIRE2" (rclas/subclass/rorder/superfamily/family)
repeat_type= "LTR_retrotransposon" (feature)
repeat_consensus= "ATATA" (tsd_sequence)


those attributes summarize extra info associated to ltr sequences. These info have been obtained by in silico analyses.
pbs_x and pbs_y are the start-end coordinates of the Primer Binding Site of the elements. Ltr_comparison_insertion_time is the estimated insertion time for that element as it was calculated using the molecular paleontology method devised by SanMiguel (Nat genetics, 1989). We used a mutation rate of 1.8 x 10-8. tsd_sequence is the 4-6 bp (usually) inverted repeat that could be found at the 5' and 3' ends of a complete not rearranged ltr-retroelement (and solo LTRs too). is_complete indicate that the element is complete i.e. (fort this annotation purposes) it has both LTRs.
Let me know if other info are needed. 


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
use Readonly;

Readonly my $repeat_name_key => "id";
Readonly my @repeat_class_keys => qw(name rclass subclass rorder superfamily family);
Readonly my $repeat_consensus_key => "tsd_sequence";



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


use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS $LOGIC_NAME );

my $date = time(); ;
my $SPECIES_NAME;

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my ( $species,  $reg, $logic_name, $no_insert);
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

  $LOGIC_NAME = $logic_name;

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

   $SPECIES_NAME = $species;
}

my $sa  = $ENS_DBA->get_adaptor('Slice');
my $rfa = $ENS_DBA->get_adaptor('RepeatFeature');
my $rca = $ENS_DBA->get_adaptor('RepeatConsensus');
my $count=0;
my $rc_name = $LOGIC_NAME || 'no-name';

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
  #$feature = uc($feature);
  #$seqname =~ s/.*chr(omosome)?_?0*//i;
  $seqname =~ s/(.*chr0?|\D*0?)//i unless $SPECIES_NAME =~ /longistaminata/i;
  $strand = $strand eq '-' ? -1: 1;
  
  my %attribs;
  foreach my $id( split( /;/, $attribute ) ){

    if( $id  =~ m/([^=]*)=(.*)/ ){
      #ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
      #ID=11562.t00100;Name="ORF249";Note="pub_locus=LOC_Osp1g01070"

      #ID=1;rclass=I;subclass=I;rorder=LTR;superfamily=Gypsy;family=RIRE2;is_complete=true;pbs_x=81402;pbs_y=81416;ppt_x=92365;ppt_y=92379;ltr_comparison_insertion_time=0.89;tsd_sequence=ATATA;presence_in_sativa=present

      my $key   = lc $1;
      my $value = $2;
      $value =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
      $value =~ s/^"//;
      $value =~ s/"$//;
      $attribs{$key} = $value;
      #print "$key => $value\n";
    }
  }


  #print "feature=$feature, $repeat_name_key => $attribs{$repeat_name_key}, $repeat_consensus_key=>$attribs{$repeat_consensus_key} \n";
  #my $repeat_name = join '-', ($feature, $attribs{$repeat_name_key});
  my $repeat_name = $attribs{$repeat_name_key};
  my $repeat_class = join '/', map{ $attribs{$_} } @repeat_class_keys;
  my $repeat_type = $feature;
  my $repeat_consensus = $attribs{$repeat_consensus_key} || 'unknown';
  
  #print "repeat_name=$repeat_name repeat_class=$repeat_class repeat_type=$repeat_type repeat_consensus=$repeat_consensus\n";


  my $rc = Bio::EnsEMBL::RepeatConsensus->new
      (
       -NAME => $repeat_name,
       -REPEAT_CLASS => $repeat_class,
       -REPEAT_TYPE => $repeat_type,
       -REPEAT_CONSENSUS => $repeat_consensus,
       -ADAPTOR => $rca,
       );
  $rca->store( $rc ) if $I;
 
  my $slice = $sa->fetch_by_region( undef, $seqname ) or 
      die "cannot fetch_by_region for $seqname\n";

  if( $start > $end ) {
      print "ERROR: start bigger than end! start=$start, end=$end, skip\n$line\n"; 
      next;
  }
  my $repeat_feature = Bio::EnsEMBL::RepeatFeature->new 
      (
       -start   => $start,
       -end     => $end,
       -strand  => $strand,
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
    print STDERR "stored $count repeat features: $start, $end, $attribs{id}\n";


}

print "Stored $count repeats from $ARGV[0]\n";
