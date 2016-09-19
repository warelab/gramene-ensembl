#! /usr/local/bin/perl -w 

=pod

=head1 NAME

load_repeats_from_gff.pl - load AGI repeats from gff into EnsemblDB

=head1 SYNOPSIS

  load_repeats_from_gff.pl [options] gff_file

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
B73V4_ctg2660   repeatmasker    match   72914   72957   12      +       .       ID=utg1472_utg114375_utg114381:hit:17618:1.3.0.0;Name=species:%28TGA%29n|genus:Simple_repeat;Target=species:%28TGA%29n|genus:Simple_repeat 1 45 +
B73V4_ctg2692   repeatmasker    match   99560   99994   369     +       .       ID=utg147687:hit:33065:1.3.0.0;Name=species:ZM_CACTA_31|genus:Unspecified;Target=species:ZM_CACTA_31|genus:Unspecified 6186 6644 +

repeat_name= "ZM_CACTA_31" (feature-ID)
repeat_class= "xxx" or the first two fields of species name (ZM_CACTA)_31
repeat_type= "Simple_repeat/Transposable Element if unspecified" (name:genus)
repeat_consensus= "(TGA)n" (Name:species:)

for repeatRunner
repeat_name=name
repeat_type='Transposable Element'
repeat_class='MobileElement'

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

our $BASEDIR ="/usr/local/";
use Readonly;

Readonly my $repeat_key_name => "name";
Readonly my $repeat_key_target => "target";
Readonly my $protein_match => "protein_match";
Readonly my $repeat_type_TE => "Transposable Element";
Readonly my $repeat_type_simple => "Simple repeats";
Readonly my $repeat_class_TE => "MobileElement";
Readonly my $repeat_class_simple => "Simple_repeat";

BEGIN{
  # Set the perl libraries 
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
my $rc_name = $LOGIC_NAME || 'no-name';

my $cnt = 0;

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
  next if $feature eq 'match_part';

  $strand = $strand eq '-' ? -1: 1;

  my ($repeat_name, $repeat_type, $repeat_class, $repeat_consensus) = ('', '', '', 'unknown');
  my ($repeat_target, $target_start, $target_end, $target_strand) = ('', 0, 0, 0); 
  my %attribs;
  ATTRIB: foreach my $id( split( /;/, $attribute ) ){

    if( $id  =~ m/([^=]*)=(.*)/ ){
	##gff-version 3
        #B73V4_ctg2660   repeatmasker    match   72914   72957   12      +       .       ID=utg1472_utg114375_utg114381:hit:17618:1.3.0.0;Name=species:%28TGA%29n|genus:Simple_repeat;Target=species:%28TGA%29n|genus:Simple_repeat 1 45 +
        #B73V4_ctg2692   repeatmasker    match   99560   99994   369     +       .       ID=utg147687:hit:33065:1.3.0.0;Name=species:ZM_CACTA_31|genus:Unspecified;Target=species:ZM_CACTA_31|genus:Unspecified 6186 6644 +
        #repeat_name= "logicName-1" (feature-ID)
        #repeat_class= "xxx" or the first two fields of species name (ZM_CACTA)_31
	#repeat_type= "Simple_repeat/Transposable Element if unspecified" (name:genus)
	#repeat_consensus= "(TGA)n" (Name:species:)
	#for repeatRunner
	#repeat_type='Transposable Element'
	#repeat_class='MobileElement'

      	my $key   = lc $1;
	my $flag_keep = 0;
      	map{  $flag_keep = 1 if ($key eq $_) } map ($repeat_key_name, $repeat_key_target);
      	next ATTRIB unless $flag_keep;	
	
	my $value = $2;
      	$value =~ s/\%([0-9]{2,4})/pack('C', hex($1))/seg;
      	$value =~ s/^"//;
      	$value =~ s/"$//;
      	$attribs{$key} = $value;
      	print "$key => $value\n";
    }
  }

  $repeat_target = $attribs{$repeat_key_target} ;
  $repeat_name = $attribs{$repeat_key_name} ;

  if ( $feature =~ /$protein_match/i){
	$repeat_type = $repeat_type_TE;
	$repeat_class = $repeat_class_TE;
  }else{
	
	my $repeat_key_value = $repeat_name;
	if( $repeat_key_value =~ /species:([()\w]+)\|genus:Simple_repeat/i){
		$repeat_type = $repeat_type_simple;
		$repeat_class = $repeat_class_simple;
		$repeat_consensus=$1;
		$repeat_name = $repeat_consensus;
	}elsif( $repeat_key_value =~ /species:(.+)\|genus:Unspecified/i ){
		$repeat_type = $repeat_type_TE;
		$repeat_name = $1;
		#if($s =~ /^([A-Z]{3})_/i)
		print "repeat_name=$repeat_name\n";
		if ( $repeat_name =~ /^([A-Z]{3})_/i ){
			$repeat_class = $1;
		}elsif( $repeat_name =~ /([A-Z]{2}_[A-Z0-9]+)/i ){
			$repeat_class = $1;
		}else{
			#RLG_cinful-zeon_AC208228_9805
			die "ERROR: wrong type of repeat_name $repeat_name \n";
		}
		
	}else{ 
  	      die "ERROR: wrong value for $repeat_key_name: $repeat_key_value  \n";
  	}
  }

  my ($tn, $tstart, $tend, $tstrand) = split ' ', $repeat_target;

  print "feature=$feature, name=$repeat_name, repeat_type=$repeat_type, repeat_class=$repeat_class, repeat_consensus=$repeat_consensus, target $tstart-$tend:$tstrand \n";
  
  $tstart = 1 unless $tstart >0;
  $tend = $end-$start+1 unless $tend > 0;
  $tstrand = $tstrand eq '-'? -1: 1;

  my $rc = $rca->fetch_by_name_class($repeat_name, $repeat_class);
	
  unless ($rc){	
    $rc = Bio::EnsEMBL::RepeatConsensus->new
      (
       -NAME => $repeat_name,
       -REPEAT_CLASS => $repeat_class,
       -REPEAT_TYPE => $repeat_type,
       -REPEAT_CONSENSUS => $repeat_consensus,
       -ADAPTOR => $rca,
       );
    $rca->store( $rc ) if $I;
  }
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
       -hstart  => $tstart,
       -hend    => $tend,
       -hstrand => $tstrand,
       -score => $score 
       
       );


    $rfa->store( $repeat_feature ) if $I;

    $cnt++;
    print STDERR "stored $cnt repeat features: $start, $end, $attribs{id}\n";


}

print "Stored $cnt repeats from $ARGV[0]\n";