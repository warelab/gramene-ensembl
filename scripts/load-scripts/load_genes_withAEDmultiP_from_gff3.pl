#!/usr/local/bin/perl -w 

=pod

=head1 NAME

load_genes_withMakerAED_gff3.pl - parse genes out from JGI gff3 file and load them into EnsemblDB
                              use Poplar v2.0 gff3 file as example

=head1 SYNOPSIS

  load_genes_withMakerAED_gff3.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -l|--logic_name       the logic name for the gene track
  -n|--no_insert        Do not make changes to the database. For debug.
  -c|--coding		nonCoding genes tRNA, miRNA ...
  -b|--biotype		nonCoding biotyoe such as miRNA, tRNA, lincRNA

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

  Parses the specified JGI gff3 file and uses the data therein to populate 
  the gene, exon, transcript and translation tables (inc stable_id tables) 
  in the ensembl DB indicated by the --species and --ensembl_registry.
  It will also load all the identifiers into xref, object_to_xref tables

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

B<Restoring the database>

  If the script bombs out half-way through, your database will be in a
  partial loaded state (i.e. in a bit of a mess). Here is some SQL to
  help put it right;

  DELETE FROM gene, transcript, translation, exon_transcript, exon,
              gene_stable_id, transcript_stable_id,
              translation_stable_id, exon_stable_id, analysis
  WHERE       gene.gene_id = transcript.gene_id
  AND         transcript.transcript_id=translation.transcript_id
  AND         transcript.transcript_id=exon_transcript.transcript_id
  AND         exon_transcript.exon_id=exon.exon_id
  AND         gene.gene_id = gene_stable_id.gene_id
  AND         transcript.transcript_id=transcript_stable_id.transcript_id
  AND         translation.translation_id=translation_stable_id.translation_id
  AND         exon.exon_id=exon_stable_id.exon_id
  AND         gene.analysis_id=analysis.analysis_id
  AND         analysis.logic_name="fgenesh_gene"

  The above could be improved with some left joins so that the
  constraints do not fail for partly-loaded gene models.

  Here is the brute-force approach that removes all gene models
  regardless of analysis;

  delete from gene; 
  delete from gene_stable_id; 
  delete from transcript; 
  delete from transcript_stable_id; 
  delete from translation; 
  delete from translation_stable_id;
  delete from exon_transcript; 
  delete from exon; 
  delete from exon_stable_id;


Maintained by Sharon Wei <weix@cshl.edu>

=cut

#delete from gene; delete from transcript; delete from translation; delete from exon; delete from exon_transcript; delete from gene_stable_id; delete from transcript_stable_id; delete from translation_stable_id; delete from exon_stable_id;

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

#Readonly my @NAME_FILEDS => qw(NAME ALIAS ID);
Readonly my @NAME_FILEDS => qw(ID NAME ALIAS);

Readonly my $UTR_REGEX   => qr{UTR}xms;

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


use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS $NONCODING $BIOTYPE );

my $date = time(); 

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $logic_name, $no_insert, $non_coding, $biotype );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "logic_name=s"       => \$logic_name,
        "no_insert"          => \$no_insert,
	"coding"             => \$non_coding,
	"biotype=s"	     => \$biotype,
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
  $NONCODING = $non_coding ? 1 : 0;
 
  $BIOTYPE = $biotype ? $biotype : undef;

warn ("noncoding flag is $NONCODING\n"); #exit;

  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) &&
                pod2usage() );

  # Create the analysis
  $ANALYSIS = Bio::EnsEMBL::Analysis->new
      ( -logic_name => $logic_name ); #'GeneModel_RiceIndica_BGI' 

  # Create a GFF stream
  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


  print( "Loading genes for $species\n" );
  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
  }
}


my $GENES = {};
my %TRPT2GENE = ();


while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /^\s*\#/ ); 
  next if ( $line =~ /^\s+/ );
  chomp $line;

  print "$line\n";  #####
  #warn(Dumper($TRPT2GENE));

  # Split gff line,\
  # start always <= end even for - strand genes
  my( $seqname, $source, $feature, 
      $start, $end, $score, $strand, $frame, 
      $attribute ) = split( /\s+/, $line, 9 );
  $feature = uc($feature);
  $seqname =~ s/^chr0*//i;

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


  
  #get gene name for this feature
  my ($gene_id, $gene_name, $gene_description); 
  my ($transcript_id, $transcript_name, $parent_id);
  my $organell = 0; #a flag for mitochondrion/chloroplast gff3,
                    #they were converted from xml using tigr tools

  if( $feature eq 'GENE' ){

    $gene_id   = $attribs{'ID'};
    $gene_name = get_first_func(\%attribs, @NAME_FILEDS); 
    
    $gene_description = $attribs{DESCRIPTION} || $gene_name;
    $attribs{DESCRIPTION} ||= $gene_description;
    
    print "Gene $gene_id $gene_name $gene_description\n";#####
    unless( $GENES->{$gene_id} ){
      $GENES->{$gene_id}->{GENE_NAME} = $gene_name;
      $GENES->{$gene_id}->{DESCRIPTION} = $gene_description;
      $GENES->{$gene_id}->{SEQ_NAME}  = $seqname;
      $GENES->{$gene_id}->{START}     = $start; #always smaller than end
      $GENES->{$gene_id}->{END}       = $end;
      $GENES->{$gene_id}->{STRAND}    = $strand eq '+' ? 1 : -1;
      $GENES->{$gene_id}->{ATTRIBS}   = \%attribs;
      
    }
  }
  elsif( $feature eq 'TRANSCRIPT' || $feature =~  /RNA/){
   
    $transcript_id   = $attribs{ID};
    $transcript_name = get_first_func->(\%attribs, @NAME_FILEDS);
    $gene_id         = $attribs{PARENT};

    next if ( $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id} );

    $attribs{DESCRIPTION} ||= $transcript_name;

    print "TRANSCRIPT_ID = $transcript_id, transcript_name = $transcript_name\n";
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{TRPT_NAME} 
      = $transcript_name;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{SEQ_NAME} 
      = $seqname;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{START} 
      = $start; 
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{END}  
      = $end;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{EXON_COUNT} 
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_COUNT} 
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_COUNT} 
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{STRAND} 
      =  $strand eq '+' ? 1 : -1;

    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{DESCRIPTION}
      = $attribs{DESCRIPTION};

    $TRPT2GENE{$transcript_id} = $gene_id;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{ATTRIBS}   
      = \%attribs;
    
  }
  elsif( $feature eq 'CDNA_MATCH'  ){
	$gene_id   = $attribs{'ID'};
    $gene_name = get_first_func(\%attribs, @NAME_FILEDS);

    $gene_description = $attribs{DESCRIPTION} || $gene_name;
    $attribs{DESCRIPTION} ||= $gene_description;

    print "Gene $gene_id $gene_name $gene_description\n";#####
    unless( $GENES->{$gene_id} ){
      $GENES->{$gene_id}->{GENE_NAME} = $gene_name;
      $GENES->{$gene_id}->{DESCRIPTION} = $gene_description;
      $GENES->{$gene_id}->{SEQ_NAME}  = $seqname;
      $GENES->{$gene_id}->{START}     = $start; #always smaller than end
      $GENES->{$gene_id}->{END}       = $end;
      $GENES->{$gene_id}->{STRAND}    = $strand eq '+' ? 1 : -1;
      $GENES->{$gene_id}->{ATTRIBS}   = \%attribs;

    	$GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{TRPT_NAME}
      = $gene_name.'.1';
    	$GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{SEQ_NAME}
      = $seqname;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{START}
      = $start;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{END}
      = $end;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{EXON_COUNT}
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{CDS_COUNT}
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{UTR_COUNT}
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{STRAND}
      =  $strand eq '+' ? 1 : -1;

    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{DESCRIPTION}
      = $gene_description;

    $TRPT2GENE{$gene_id} = $gene_id;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{ATTRIBS}
      = \%attribs;
    } 

    $GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{CDS_COUNT} ++;
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{CDS_EXON_START}}, 
      $start;
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$gene_id}->{CDS_EXON_END}}, 
      $end;

  }

  elsif( $feature =~ $UTR_REGEX  ){
      my @transcript_ids = split ',', $attribs{'PARENT'};
      
      for my $atranscript_id (@transcript_ids){
	  $gene_id = $TRPT2GENE{$atranscript_id};
	  print"UTR $atranscript_id, $gene_id\n";#####
	  unless ( $gene_id ){
	      die("cannot find gene for transcript $atranscript_id at ".
		  $GFF_HANDLE->input_line_number );
	  }
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{UTR_COUNT} ++;
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{UTR_EXON_START} ||= [];
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{UTR_EXON_END}   ||= [];
	  push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{UTR_EXON_START}}, 
	  $start;
	  push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{UTR_EXON_END}}, 
	  $end;
      }
  }
  elsif( $feature eq 'EXON'){
    #redundant
    #we can figure out the exon boundry from UTR and CDs
      my @transcript_ids = split ',', $attribs{'PARENT'};
      for my $atranscript_id (@transcript_ids){
	  $gene_id = $TRPT2GENE{$atranscript_id};
	  print"EXON $atranscript_id, $gene_id\n";#####
	  unless ( $gene_id ){
	      die("cannot find gene for transcript $atranscript_id at ".
		  $GFF_HANDLE->input_line_number );
	  }
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{EXON_COUNT} ++;
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{EXON_START} ||= [];
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{EXON_END}   ||= [];
	  push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{EXON_START}}, 
	  $start;
	  push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{EXON_END}}, 
	  $end;
      }
  }
  elsif( $feature eq 'CDS'){
    #redundant
      #we can figure out the exon boundry from UTR and CDs
      my @transcript_ids = split ',', $attribs{'PARENT'};

      for my $atranscript_id (@transcript_ids){
	  $gene_id = $TRPT2GENE{$atranscript_id};
	  print"CDS $atranscript_id, $gene_id\n";#####
	  unless ( $gene_id ){
	      die("cannot find gene for transcript $atranscript_id at ".
		  $GFF_HANDLE->input_line_number );
	  }
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{CDS_COUNT} ++;
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{CDS_EXON_START} ||= [];
	  $GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{CDS_EXON_END}   ||= [];
	  push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{CDS_EXON_START}},
	  $start;
	  push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$atranscript_id}->{CDS_EXON_END}},
	  $end;
      }
  }
  else{
    warn( "Unrecognized feature $feature" );
  }
}

# compute more features and store in the hash
#
foreach my $g( keys %$GENES ){
  my $genedata = $GENES->{$g};
  $genedata->{TRANSCRIPT_COUNT} = keys %{$genedata->{TRANSCRIPTS}};

  foreach my $t( keys %{$genedata->{TRANSCRIPTS}} ){

    my $trptdata = $genedata->{TRANSCRIPTS}->{$t};
    #reconstruct exons and save in the hash
    ## now we need to compute exons using CDS and UTRS

    my @trpt_ex_starts = sort{ $a<=>$b } (defined $trptdata->{EXON_START} ? @{$trptdata->{EXON_START}} : (),
					  defined $trptdata->{CDS_EXON_START} ? @{$trptdata->{CDS_EXON_START}} : (),
					  defined $trptdata->{UTR_EXON_START} ? @{$trptdata->{UTR_EXON_START}} : () 
					 );
    my @trpt_ex_ends   = sort{ $a<=>$b } (defined $trptdata->{EXON_END} ? @{$trptdata->{EXON_END}}:(),
					defined $trptdata->{CDS_EXON_END} ? @{$trptdata->{CDS_EXON_END}}:(),
				       defined $trptdata->{UTR_EXON_END} ? @{$trptdata->{UTR_EXON_END}} : ()
					 );
	
    my ($trpt_exon_starts, $trpt_exon_ends) = compute_exons(\@trpt_ex_starts,
							    \@trpt_ex_ends);
    
    die "cannot compute exons for $g, $t" 
      unless( $trpt_exon_starts &&
	      $trpt_exon_ends &&
	      scalar @{$trpt_exon_starts} == scalar @{$trpt_exon_ends});

    $trptdata->{EXON_COUNT} = scalar @{$trpt_exon_starts};
    $trptdata->{EXON_START} = $trpt_exon_starts;
    $trptdata->{EXON_END}   = $trpt_exon_ends;
        
    # Get exons into gene-order
    my $f = $trptdata->{STRAND} > 0 ? 1 : 0;
    foreach my $key qw( EXON_START EXON_END CDS_EXON_START CDS_EXON_END ){
      $trptdata->{$key} = [ sort{ $f ?$a<=>$b :$b<=>$a} @{$trptdata->{$key}} ] if defined $trptdata->{$key};
    }

    # Find the translation start/stop
    my( $start_codon ) = sort{ $a<=>$b } @{$trptdata->{CDS_EXON_START} || []};
    my( $stop_codon )  = sort{ $b<=>$a } @{$trptdata->{CDS_EXON_END} || []};
    if( $trptdata->{STRAND} < 0 ){ #gene oriented
      ($start_codon,$stop_codon) = ($stop_codon,$start_codon);
    }

    $trptdata->{START_CODON} = $start_codon;
    $trptdata->{STOP_CODON}  = $stop_codon;
warn ("start codon coord=$start_codon, stop codon coord =$stop_codon\n");

    # Flag exons where coding start/stops
    for( my $i=0; $i<$trptdata->{EXON_COUNT}; $i++ ){
      my $start = $trptdata->{EXON_START}->[$i];
      my $end   = $trptdata->{EXON_END}->[$i];
      if( $start_codon >= $start and $start_codon <= $end ){
	$trptdata->{CDS_START_EXON} = $i;
      }
      if( $stop_codon >= $start and $stop_codon <= $end ){
	$trptdata->{CDS_END_EXON} = $i;
      }
    }
    unless( $NONCODING ||  (defined $trptdata->{CDS_END_EXON} &&
          defined $trptdata->{CDS_START_EXON}) ){
      #warn Dumper( $trptdata );
      die( "Gene-Transcript ${g}-${t} has no CDS_END_EXON/CDS_START_EXON" );
    }
    
  }
}
 
warn ("Done parsing\n");

# Load the data
if( $NONCODING ){
        $BIOTYPE ||= 'nonCoding' ;
}

my $sa = $ENS_DBA->get_adaptor('Slice');
my $n = 1;
my $number = scalar( keys %$GENES );
foreach my $gene_id( keys %$GENES ){

  my $agenedata = $GENES->{$gene_id};
  warn("Load gene stable_id ", $agenedata->{GENE_NAME}, "\n");

  my $eGene = Bio::EnsEMBL::Gene->new();
  $eGene->analysis($ANALYSIS);
  $eGene->stable_id( $agenedata->{GENE_NAME} );
  $eGene->description( $agenedata->{DESCRIPTION} );
  $eGene->version(1);  
  $eGene->created_date( $date );
  $eGene->modified_date( $date );
  $eGene->biotype( $BIOTYPE ) if $BIOTYPE;

  # Add XREFs to gene
  print Dumper ($agenedata->{ATTRIBS});
  
  foreach my $trpt_id (keys %{$agenedata->{TRANSCRIPTS}} ){
    my $atrptdata = $agenedata->{TRANSCRIPTS}->{$trpt_id};
    my $trpt_name =  $atrptdata->{TRPT_NAME};

    warn("\tLoad transcript id $trpt_name\n");
    my $eTranscript = Bio::EnsEMBL::Transcript->new();
    $eTranscript->stable_id( $trpt_name );
    $eTranscript->analysis( $ANALYSIS );
    $eTranscript->version(1);
    $eTranscript->description( $atrptdata->{DESCRIPTION} );
    $eTranscript->created_date( $date );
    $eTranscript->modified_date( $date );
    $eTranscript->biotype( $BIOTYPE ) if $BIOTYPE;

	if(my $aed =$atrptdata->{ATTRIBS}->{'_AED'} ){
		
		warn ("aed=$aed\n");
		my $aed_attrib = Bio::EnsEMBL::Attribute->new
       		(-CODE => 'AED',
        	-VALUE => "$aed");
	$eTranscript->add_Attributes($aed_attrib);

	}
    my $transcript_xref;

    # Transcript xrefs...
    #print "For transcript\n";
    #print Dumper ($atrptdata->{ATTRIBS});
  
    my $seq_region_name = $atrptdata->{SEQ_NAME};
     #                         print "seq_region_name=$seq_region_name\n";
    $seq_region_name =~ s/^chr//i;
                              print "seq_region_name=$seq_region_name\n";
    my $slice = $sa->fetch_by_region( undef, $seq_region_name );
	#warn ("DEBUG sliceid = ", $slice->seq_region_name, "\n");

    my $start_exon;
    my $end_exon;
    my $trans_start_exon;
    my $trans_end_exon;
    my $translation_offset_left;
    my $translation_offset_right;
    my $exon_start_phase = -1;
    my $exon_end_phase   = -1;
    my $phase_diff       = 0;
    my $last_exon_phase  = -1;
    my $noncoding_exon_phase =  -1;

    my $transcript_start = $atrptdata->{START_CODON};
    my $transcript_stop  = $atrptdata->{STOP_CODON};
  
    if( $NONCODING ){
        $BIOTYPE ||= 'nonCoding' ;
        $eGene->biotype( $BIOTYPE );
    }else{ 
    	$transcript_start || die( "Trpt ${gene_id}-${trpt_id} has no start" );
    	$transcript_stop  || die( "Trpt ${gene_id}-${trpt_id} has no stop" );
    }
    #########
    # EXONS #
    #########
    
    for ( my $exon = 0; $exon < $atrptdata->{EXON_COUNT}; $exon++ ){
      
      my $exon_start  = $atrptdata->{EXON_START}[$exon];
      my $exon_end    = $atrptdata->{EXON_END}[$exon];

      warn("\t\tLoad exon $exon: $trpt_name.exon", $exon+1 , "\n");
      my $eExon = new Bio::EnsEMBL::Exon;
      $eExon->start($exon_start);
      $eExon->end($exon_end);
      $eExon->strand($atrptdata->{STRAND});	   
      $eExon->slice($slice);
      $eExon->stable_id( $trpt_name.'.exon'.($exon+1) );
      $eExon->version( 1 );
      
      # Phase calculations
      # start_phase of next exon is the end_phase of this exon by Ensembl API
      # which means the end phase of an exon is not the translation phase of the last base of the exon
      if($NONCODING){
	    $eExon->phase( $noncoding_exon_phase );
        $eExon->end_phase( $noncoding_exon_phase );
	    $eTranscript->add_Exon($eExon);
      }
      else{
        if ( $atrptdata->{STRAND} > 0 ) {
	        if ( ( $transcript_start > $exon_start ) &&
	            ( $transcript_start > $exon_end ) )  { # 5' exons, was different in original
	            $exon_start_phase = -1;
	            $exon_end_phase   = -1;
	        }
	        elsif(  $transcript_start == $exon_start  ) { # exon starts with a start
	            $phase_diff =  ( ( $exon_end - $exon_start + 1 ) % 3 );
	            $exon_start_phase = 0;
	            $exon_end_phase   = $transcript_stop < $exon_end ?  -1 : $phase_diff;
	        }
	        elsif ( ( $transcript_start > $exon_start  ) &&
		        ( $transcript_start <= $exon_end ) ) { # exon contains a start
	            $phase_diff = ( ( $exon_end - $transcript_start + 1 ) % 3 );
	            $exon_start_phase = - 1;
	            $exon_end_phase   = $transcript_stop < $exon_end ?  -1 : $phase_diff;
	        }
	        elsif ( ( $transcript_stop >= $exon_start ) &&
		            ( $transcript_stop < $exon_end ) ) { # exon contains a stop
	            $exon_start_phase = $last_exon_phase;
	            $exon_end_phase   = -1;
	        }
	        elsif ( ( $transcript_stop < $exon_start ) &&
		        ( $transcript_stop < $exon_end ) ) { # 3' exons
	            $exon_start_phase = -1;
	            $exon_end_phase   = -1;
	        }
        	else { # internal exon
	            $phase_diff = ( $exon_end - $exon_start + 1 ) % 3;
	            $exon_start_phase = $last_exon_phase;
	            $exon_end_phase   = ( $exon_start_phase + $phase_diff ) % 3;
	        }
	# set exon phase
	        $eExon->phase($exon_start_phase);
	        $eExon->end_phase($exon_end_phase);

	        $last_exon_phase = $exon_end_phase;
      }
      else{ # -ve strand
	# 5' exons
	        if    ( ( $transcript_start < $exon_start ) &&
		        ( $transcript_start < $exon_end ) )  {
	            $exon_start_phase = -1;
	            $exon_end_phase   = -1;
	        }
	# exon stops with a start 
	        elsif ( ( $transcript_start == $exon_end ) ) {
       	        $phase_diff = ( ( $exon_end - $exon_start + 1 ) % 3 );
	            $exon_start_phase = 0;
	            $exon_end_phase   = $transcript_stop <= $exon_start ? $phase_diff : -1;
	        }
	# exon contains a start codon
	        elsif ( ( $transcript_start >= $exon_start ) &&
		            ( $transcript_start < $exon_end ) ) {
        	    $phase_diff = ( ( $transcript_start - $exon_start + 1 ) % 3 );
	            $exon_start_phase = - 1;
	            $exon_end_phase   = $transcript_stop <= $exon_start ? $phase_diff : -1;
	        }
	# exon contains a stop codon
	        elsif ( ( $transcript_stop > $exon_start ) &&
		        ( $transcript_stop <= $exon_end ) ) {
	            $exon_start_phase = $last_exon_phase;
	            $exon_end_phase   = -1;
	        }
	# 3' exons
	        elsif ( ( $transcript_stop > $exon_start ) &&
		        ( $transcript_stop > $exon_end ) ) {
	            $exon_start_phase = -1;
	            $exon_end_phase   = -1;
	        }
	# internal exon
	        else {
	            $phase_diff = ( $exon_end - $exon_start + 1 ) % 3;
	            $exon_start_phase = $last_exon_phase;
	            $exon_end_phase   = ( $exon_start_phase + $phase_diff ) % 3;
	        }
	# set exon phase
	        $eExon->phase($exon_start_phase);
	        $eExon->end_phase($exon_end_phase);

    	    $last_exon_phase = $exon_end_phase;
      }
    #warn("Exon name2hashkey: ", $eExon->stable_id, "=", $eExon->hashkey, "\n");
      $eTranscript->add_Exon($eExon);
      
      # Start exon
      unless ( $start_exon ) {
	    if ( $exon == 0 ) {$start_exon = $eExon; }
      }
      # Final exon
      unless ( $end_exon ) {
	    if ( $exon == $atrptdata->{EXON_COUNT}-1 ) {$end_exon = $eExon;}
      }
      
      # Translation start exon
      unless ( $trans_start_exon ) {
    	if ( $exon == $atrptdata->{CDS_START_EXON} ) {$trans_start_exon = $eExon;}
      }
      
      # Translation stop exon
      unless ( $trans_end_exon ) {
	    if ( $exon == $atrptdata->{CDS_END_EXON} ) {$trans_end_exon = $eExon;}
      }
     }#end of noncoding
    } #_ END of exon loop _#
    
    unless ($NONCODING){
        $eTranscript->start_Exon($start_exon);
        $eTranscript->end_Exon($end_exon);
    
    ###############
    # TRANSLATION #
    ###############
    
    # Check for translation start being different to Exon start (i.e. UTR)
        if ( $atrptdata->{STRAND} > 0 ) {
            $translation_offset_left  =
                $transcript_start - $atrptdata->{EXON_START}[$atrptdata->{CDS_START_EXON}] + 1;
            $translation_offset_right =
                $transcript_stop  - $atrptdata->{EXON_START}[$atrptdata->{CDS_END_EXON}] + 1;
        }
        elsif ( $atrptdata->{STRAND} < 0 ) {
            $translation_offset_left  =
                $atrptdata->{EXON_END}[$atrptdata->{CDS_START_EXON}] - $transcript_start + 1;
            $translation_offset_right =
                $atrptdata->{EXON_END}[$atrptdata->{CDS_END_EXON}]   - $transcript_stop  + 1;
        }

        my $translation_name = $trpt_name;
        $translation_name =~ s/_T/_P/i;
        print "translation: $trans_start_exon, $trans_end_exon, $translation_offset_left, $translation_offset_right,   $translation_name\n";
        my $eTranslation = new  Bio::EnsEMBL::Translation
        (
            -START_EXON  => $trans_start_exon,
            -END_EXON    => $trans_end_exon,
            -SEQ_START   => $translation_offset_left,
            -SEQ_END     => $translation_offset_right,
            -STABLE_ID   => $translation_name,
            -VERSION     => 1,
        );
   # $eTranslation->stable_id( $trpt_name );
    
    #########################################
    # EnsEMBL add translation to transcript #
    #########################################
    
    $eTranscript->translation($eTranslation);

   }#end of else noncoding
    ##################################
    # EnsEMBL add transcript to gene #
    ##################################
    
    $eGene->add_Transcript($eTranscript);
  }
 
  if( $I ){
    my $dbid = $ENS_DBA->get_GeneAdaptor->store($eGene);
    print "// dbID = $dbid for gene = $gene_id ($n of $number)\n";
  }
  $n++;
}

#warn Dumper( $GENES );


exit();


#
# need to join the UTR and CDs that have continuous coordinates

#chr17   Gaze_filter     gene    5455880 5470620 178.7753        -       .       Gene GSVIVG00000008001 ; Complete 1
#chr17   Gaze_filter     mRNA    5455880 5470620 178.7753        -       .       mRNA GSVIVT00000008001 ; Gene GSVIVG00000008001 ; Complete 1
#chr17   Gaze_filter     UTR     5455880 5456137 -0.7006 -       .       UTR GSVIVT00000008001                
#chr17   Gaze_filter     UTR     5456161 5456247 -1.6200 -       .       UTR GSVIVT00000008001
#chr17   Gaze_filter     CDS     5456248 5457410 21.4375 -       1       CDS GSVIVT00000008001                chr17   Gaze_filter     CDS     5457519 5457702 3.6407  -       0       CDS GSVIVT00000008001
#chr17   Gaze_filter     CDS     5470468 5470620 8.4669  -       0       CDS GSVIVT00000008001    

#test at perl ~/scripts/test-merge-CDUTR.pl

sub compute_exons{ 

  my ($trpt_ex_starts, $trpt_ex_ends) = @_;
  
  my $start_cnt = scalar @{$trpt_ex_starts} ;
  my $end_cnt = scalar @{$trpt_ex_ends};
 
  if( $start_cnt != $end_cnt ){
    return ();
  }

  my @exon_starts = ( $trpt_ex_starts->[0] );
  my @exon_ends = ();
  for (my $i=1; $i<$start_cnt; $i++){
    if( $trpt_ex_starts->[$i] - $trpt_ex_ends->[$i-1] <= 1){
         
    }else{
      push @exon_ends, $trpt_ex_ends->[$i-1];
      push @exon_starts, $trpt_ex_starts->[$i];
    }
  }

  push @exon_ends, $trpt_ex_ends->[$start_cnt-1];

  return (\@exon_starts, \@exon_ends);
  
}

#----------------------------------------------------------------------

sub get_first_func  { 
  
  my $attribs = shift;

  first { defined } map { $attribs->{$_}  } @_;

}

1;
