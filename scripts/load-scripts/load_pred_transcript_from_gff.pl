#!/usr/local/bin/perl -w 

=pod

=head1 NAME

load_genes_from_gff3.pl - parse genes out from gff3 file and load them into EnsemblDB
                              use Poplar v2.0 and grape gff3 file as example

=head1 SYNOPSIS

  load_genes_from_gff3.pl [options] gff_file

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

Readonly my @NAME_FILEDS => qw(NAME ALIAS ID);
Readonly my $UTR_REGEX   => qr{UTR}xmsi;
Readonly my $GRAPE_REGEX => qr{ obsolete }xmsi;  #No need to deal with grape differently
Readonly my $GRAPE_FEATURE_NAME_REGEX => qr{ obsolete }xmsi;

Readonly my @FEATURES => qw( EXON_START EXON_END CDS_EXON_START CDS_EXON_END );
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

use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::PredictionExon;
use Bio::EnsEMBL::DBSQL::PredictionExonAdaptor;
use Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor;

use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS $SPECIES );

my $date = time(); 

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $logic_name, $no_insert, @xref_srcs );
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

  #my $PT_ADAPTOR = $ENS_DBA->get_PredictionTranscriptAdaptor();
  #my $PE_ADAPTOR = $ENS_DBA->get_PredictionExonAdaptor();

  # Create the analysis
  $ANALYSIS = Bio::EnsEMBL::Analysis->new
      ( -logic_name => $logic_name ); #'GeneModel_RiceIndica_BGI' 

  # Create a GFF stream
  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


  $SPECIES = $species;
  print( "Loading genes for $SPECIES\n" );
  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
  }
}


our $GENES = {};
our %TRPT2GENE = ();


while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /\#/ ); 
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
  $seqname =~ s/chr//i;

  my $attribs_ref = get_attrib( $attribute, $SPECIES) ;

  #print "The attribs is ", Dumper( $attribs_ref);
  #exit;

  
  #get gene name for this feature
  my ($gene_id, $gene_name, $gene_description); 
  my ($transcript_id, $transcript_name, $parent_id);
  my $organell = 0; #a flag for mitochondrion/chloroplast gff3,
                    #they were converted from xml using tigr tools

  if( $feature eq 'GENE' ){

    $gene_id   = $attribs_ref->{'ID'};
    $gene_name = get_first_func($attribs_ref, @NAME_FILEDS); 
    
    $gene_description = $attribs_ref->{DESCRIPTION} || $gene_name;
    $attribs_ref->{DESCRIPTION} ||= $gene_description;
    
    print "Gene $gene_id $gene_name $gene_description\n";#####
    unless( $GENES->{$gene_id} ){
      $GENES->{$gene_id}->{GENE_NAME} = $gene_name;
      $GENES->{$gene_id}->{DESCRIPTION} = $gene_description;
      $GENES->{$gene_id}->{SEQ_NAME}  = $seqname;
      $GENES->{$gene_id}->{START}     = $start; #always smaller than end
      $GENES->{$gene_id}->{END}       = $end;
      $GENES->{$gene_id}->{STRAND}    = $strand eq '+' ? 1 : -1;
      $GENES->{$gene_id}->{ATTRIBS}   = $attribs_ref;
      
    }
  }
  elsif( $feature eq 'MRNA' ){
   
    $transcript_id   = $attribs_ref->{ID};
    $transcript_name = get_first_func->($attribs_ref, @NAME_FILEDS);
    $gene_id         = get_gene_id( $attribs_ref, $SPECIES);

    next if ( $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id} );

    $attribs_ref->{DESCRIPTION} ||= $transcript_name;

    print "TRANSCRIPT_ID = $transcript_id, transcript_name = $transcript_name, gene=$gene_id\n";
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
      = $attribs_ref->{DESCRIPTION};

    $TRPT2GENE{$transcript_id} = $gene_id;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{ATTRIBS}   
      = $attribs_ref;
    
    #print " TRPT2GENE{$transcript_id} = $gene_id, $TRPT2GENE{$transcript_id} \n";
  }
  elsif( $feature eq 'CDS' ){ # use CDS and UTR as a proxy for exon

      #for (keys %TRPT2GENE){
#	  print "$_ => $TRPT2GENE{$_}\n";
#      }
    $transcript_id = $attribs_ref->{'PARENT'} || $attribs_ref->{ID};
    $gene_id       = $TRPT2GENE{$transcript_id};
    print"CDS for $transcript_id, $gene_id\n";#####
    unless ( $gene_id ){
      die("cannot find gene for transcript $transcript_id at " . 
	  $GFF_HANDLE->input_line_number );
      
    }
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_COUNT} ++;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_START} ||= [];
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_END}   ||= [];
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_START}}, 
      $start;
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_END}}, 
      $end;

  }

  elsif( $feature =~ $UTR_REGEX  ){
    $transcript_id = $attribs_ref->{'PARENT'} || $attribs_ref->{ID};
    $gene_id = $TRPT2GENE{$transcript_id};
    print"UTR $transcript_id, $gene_id\n";#####
    unless ( $gene_id ){
      die("cannot find gene for transcript $transcript_id at ".
	 $GFF_HANDLE->input_line_number );
      
    }
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_COUNT} ++;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_START} ||= [];
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_END}   ||= [];
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_START}}, 
      $start;
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_END}}, 
      $end;
    
  }
  elsif( $feature eq 'EXON'){
    #redundant
    #we can figure out the exon boundry from UTR and CDs
    $transcript_id = $attribs_ref->{'PARENT'} || $attribs_ref->{ID};
    $gene_id       = $TRPT2GENE{$transcript_id};
    print"EXON for $transcript_id, $gene_id\n";#####
    unless ( $gene_id ){
      die("cannot find gene for transcript $transcript_id at " .
	                $GFF_HANDLE->input_line_number );
    }
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{EXON_COUNT} ++;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{EXON_START} ||= [];
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{EXON_END}   ||= [];
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{EXON_START}},
	          $start;
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{EXON_END}},
	            $end;

  }
  else{
    die( "Unrecognized feature $feature" );
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

    my @trpt_ex_starts = ();
    my @trpt_ex_ends = ();

    if( defined $trptdata->{UTR_EXON_START} && defined $trptdata->{CDS_EXON_START}){
    	@trpt_ex_starts  = sort{ $a<=>$b } (@{$trptdata->{CDS_EXON_START}}, @{$trptdata->{UTR_EXON_START}});
    }else{
	    @trpt_ex_starts  = sort{ $a<=>$b }   defined $trptdata->{EXON_START} ? 
	    		@{$trptdata->{EXON_START}} : defined $trptdata->{CDS_EXON_START} ? 
		        $trptdata->{CDS_EXON_START} : ();	
    }

   if( defined $trptdata->{UTR_EXON_END} && defined $trptdata->{CDS_EXON_END}){ 
    	@trpt_ex_ends   = sort{ $a<=>$b } (@{$trptdata->{CDS_EXON_END}}, @{$trptdata->{UTR_EXON_END}} );
   }else{
	@trpt_ex_ends  = sort{ $a<=>$b }   defined $trptdata->{EXON_END} ? 
	                        @{$trptdata->{EXON_END}} : defined $trptdata->{CDS_EXON_END} ?
				$trptdata->{CDS_EXON_END} : ();

   }

    my ($trpt_exon_starts, $trpt_exon_ends) = compute_exons(\@trpt_ex_starts,
							    \@trpt_ex_ends);
    
    #print "gene=$g, trpt=$t, $trpt_exon_starts->[0], $trpt_exon_ends->[0]\n";
    
    die "cannot compute exons for $g, $t" 
      unless( $trpt_exon_starts &&
	      $trpt_exon_ends &&
	      scalar @{$trpt_exon_starts} == scalar @{$trpt_exon_ends});

    $trptdata->{EXON_COUNT} = scalar @{$trpt_exon_starts};
    $trptdata->{EXON_START} = $trpt_exon_starts;
    $trptdata->{EXON_END}   = $trpt_exon_ends;
        
    # Get exons into gene-order
    my $f = $trptdata->{STRAND} > 0 ? 1 : 0;
    foreach my $key ( @FEATURES ) {
      $trptdata->{$key} = [ sort{ $f ?$a<=>$b :$b<=>$a} @{$trptdata->{$key}} ];
    }
    


    # Find the translation start/stop
    my( $start_codon ) = sort{ $a<=>$b } @{$trptdata->{CDS_EXON_START}} || @{$trptdata->{EXON_START}} ;#
    my( $stop_codon )  = sort{ $b<=>$a } @{$trptdata->{CDS_EXON_END}} || @{$trptdata->{EXON_END}} ;
    if( $trptdata->{STRAND} < 0 ){ #gene oriented
      ($start_codon,$stop_codon) = ($stop_codon,$start_codon);
    }

    $trptdata->{START_CODON} = $start_codon;
    $trptdata->{STOP_CODON}  = $stop_codon;

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
    unless( defined $trptdata->{CDS_END_EXON} &&
          defined $trptdata->{CDS_START_EXON} ){
	  #warn Dumper( $trptdata );
      warn( "Gene-Transcript ${g}-${t} has no CDS_END_EXON/CDS_START_EXON" );
    }
   
  }
}


# Load the data
my $sa = $ENS_DBA->get_adaptor('Slice');
my $n = 1;
my $number = scalar( keys %$GENES );
foreach my $gene_id( keys %$GENES ){
  my $agenedata = $GENES->{$gene_id};
  #x#my $eGene = Bio::EnsEMBL::Gene->new();
  #x#$eGene->analysis($ANALYSIS);
  #x#$eGene->stable_id( $agenedata->{GENE_NAME} );
  #x#$eGene->description( $agenedata->{DESCRIPTION} );
  #x#$eGene->version(1);  
  #x#$eGene->created_date( $date );
  #x#$eGene->modified_date( $date );

  # Add XREFs to gene
  #print Dumper ($agenedata->{ATTRIBS});
  #add_xrefs( $eGene, $agenedata->{ATTRIBS}, 'GENE');
  
  foreach my $trpt_id (keys %{$agenedata->{TRANSCRIPTS}} ){
    my $atrptdata = $agenedata->{TRANSCRIPTS}->{$trpt_id};
    my $trpt_name =  $atrptdata->{TRPT_NAME};
    
    my $seq_region_name = $atrptdata->{SEQ_NAME};
    #print "seq_region_name=$seq_region_name\n";
    $seq_region_name =~ s/(chr)*0*//i;
    #                         print "seq_region_name=$seq_region_name\n";
    my $slice = $sa->fetch_by_region( undef, $seq_region_name );

    print "created slice for $seq_region_name\n";


    #print "constructing gene object, trpt_name=$trpt_name\n";
    #x#my $eTranscript = Bio::EnsEMBL::Transcript->new();
    my $eTranscript = Bio::EnsEMBL::PredictionTranscript->new
                  ( '-start'         =>  $atrptdata->{START},
                    '-end'           =>  $atrptdata->{END},
                    '-strand'        =>  $atrptdata->{STRAND},
                   # '-adaptor'       =>  $PT_ADAPTOR,
                    '-slice'         =>  $slice,
                    '-analysis'      =>  $ANALYSIS,
                    '-display_label' =>  $trpt_name);

    #X#$eTranscript->stable_id( $trpt_name );
    #X#$eTranscript->analysis( $ANALYSIS );
    #X#$eTranscript->version(1);
    #X#$eTranscript->description( $atrptdata->{DESCRIPTION} );
    #X#$eTranscript->created_date( $date );
    #X#$eTranscript->modified_date( $date );

   # my $transcript_xref;

    # Transcript xrefs...
    #print "For transcript\n";
    #print Dumper ($atrptdata->{ATTRIBS});
    #add_xrefs( $eTranscript, $atrptdata->{ATTRIBS}, 'TRANSCRIPT');
  

    my $start_exon;
    my $end_exon;
    my $trans_start_exon;
    my $trans_end_exon;
    my $translation_offset_left;
    my $translation_offset_right;
    my $exon_start_phase = 0;
    my $exon_end_phase   = 0;
    my $phase_diff       = 0;
    my $last_exon_phase  = 0;
    
    my $transcript_start = $atrptdata->{START_CODON};
    my $transcript_stop  = $atrptdata->{STOP_CODON};
    
    #$transcript_start || die( "Trpt ${gene_id}-${trpt_id} has no start" );
    #$transcript_stop  || die( "Trpt ${gene_id}-${trpt_id} has no stop" );
    
    #########
    # EXONS #
    #########
    
    for ( my $exon = 0; $exon < $atrptdata->{EXON_COUNT}; $exon++ ){
      
      my $exon_start  = $atrptdata->{EXON_START}[$exon];
      my $exon_end    = $atrptdata->{EXON_END}[$exon];
      
      
      my $eExon = new Bio::EnsEMBL::PredictionExon(
					       -START   => $exon_start,
					       -END     => $exon_end,
					       -STRAND  => $atrptdata->{STRAND},
					       -SLICE   => $slice,
					       );
      #X#my $eExon = new Bio::EnsEMBL::Exon;
      #X#$eExon->start($exon_start);
      #X#$eExon->end($exon_end);
      #X#$eExon->strand($atrptdata->{STRAND});	   
      #X#$eExon->slice($slice);
      #X#$eExon->stable_id( $trpt_name.'.exon'.($exon+1) );
      #X#$eExon->version( 1 );
      
      #print "constructed exon for  $trpt_name $exon+1, $exon_start, $exon_end\n";
      # Phase calculations 
      
      if ( $atrptdata->{STRAND} > 0 ) {
	if ( ( $transcript_start > $exon_start ) && 
	     ( $transcript_start > $exon_end ) )  { # 5' exons, was different in original
	  $exon_start_phase = -1;
	  $exon_end_phase   = -1;
	}
	elsif( ( $transcript_start == $exon_start ) && 
	       ( $transcript_start < $exon_end ) ) { # exon starts with a start
	  $phase_diff =  ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = 0;
	  $exon_end_phase   = $phase_diff;
	}
	elsif ( ( $transcript_start > $exon_start  ) && 
		( $transcript_start < $exon_end ) ) { # exon contains a start
	  $phase_diff = ( ( $exon_end - $transcript_start + 1 ) % 3 );
	  $exon_start_phase = - 1;
	  $exon_end_phase   = $phase_diff;
	}
	elsif ( ( $transcript_stop > $exon_start ) && 
		( $transcript_stop < $exon_end ) ) { # exon contains a stop
	  $phase_diff = ( ( $exon_end - $transcript_stop + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = -1;
	}
	elsif ( ( $transcript_stop == $exon_end ) && 
		( $transcript_stop > $exon_start ) ) { # exon stops with a stop
	  $phase_diff = ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = 0;
	}
	elsif ( ( $transcript_stop < $exon_start ) && 
		( $transcript_stop < $exon_end ) ) { # 3' exons
	  $exon_start_phase = -1;
	  $exon_end_phase   = -1;
	}
	elsif ( ( $transcript_start == $exon_start ) && 
		( $transcript_stop == $exon_end ) ) { # single exon genes
	  $phase_diff = 0;
	  $exon_start_phase = 0;
	  $exon_end_phase   = 0;
	}
	else { # internal exon
	  $phase_diff = ( $exon_end - $exon_start + 1 ) % 3;
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = ( $last_exon_phase + $phase_diff ) % 3;
	}
	
	# set exon phase
	$eExon->phase($exon_start_phase);
	#$eExon->end_phase($exon_end_phase);
	
	#$span = $exon_end - $exon_start + 1;
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
	  #was   $transcript_start == ($exon_end - 2)
	  $phase_diff = ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = 0;
	  $exon_end_phase   = $phase_diff;
	}
	# exon contains a start codon
	elsif ( ( $transcript_start > $exon_start ) && 
		( $transcript_start < $exon_end ) ) {
	  #$phase_diff = ( ( $exon_end - $transcript_start + 1 ) % 3 );
	  $phase_diff = ( ( $transcript_start - $exon_start + 1 ) % 3 );
	  $exon_start_phase = - 1;
	  $exon_end_phase   = $phase_diff;
	}
	# exon contains a stop codon
	elsif ( ( $transcript_stop > $exon_start ) && 
		( $transcript_stop < $exon_end ) ) {
	  $phase_diff = ( ( $exon_end - $transcript_stop + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = -1;
	}
	# exon stops with a stop
	elsif ( ( $transcript_stop == $exon_start ) ) { 
	  $phase_diff = ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = 0;
	}
	# 3' exons
	elsif ( ( $transcript_stop > $exon_start ) && 
		( $transcript_stop > $exon_end ) ) {
	  $exon_start_phase = -1;
	  $exon_end_phase   = -1;
	}
	# single exon genes
	elsif ( ( $transcript_start == $exon_start ) && 
		( $transcript_stop == $exon_end ) ) {
	  $phase_diff = 0;
	  $exon_start_phase = 0;
	  $exon_end_phase   = 0;
	}
	# internal exon
	else {
	  $phase_diff = ( $exon_end - $exon_start + 1 ) % 3;
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = ( $last_exon_phase + $phase_diff ) % 3;
	}
	
	# set exon phase
	$eExon->phase($exon_start_phase);
	#$eExon->end_phase($exon_end_phase);
	
	#$span = $exon_end - $exon_start + 1;      
	$last_exon_phase = $exon_end_phase;
      }
      
      eval{      $eTranscript->add_Exon($eExon); };
      if($@){print STDERR "$@"; die "failed adding exon to transcript" ;}
      
  } #_ END of exon loop _#
    
    
    #X# $eTranscript->start_Exon($start_exon);
    
    #X# $eTranscript->end_Exon($end_exon);

    if( $I ){
	$ENS_DBA->get_PredictionTranscriptAdaptor->store($eTranscript);
	print "// store prediction_transcript = $trpt_id of gene $gene_id, number $n\n";
    }
    $n++;
    
}#end of foreach transcript
    
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
    if( $trpt_ex_starts->[$i] - $trpt_ex_ends->[$i-1] == 1){
         
    }else{
      push @exon_ends, $trpt_ex_ends->[$i-1];
      push @exon_starts, $trpt_ex_starts->[$i];
    }
  }

  push @exon_ends, $trpt_ex_ends->[$start_cnt-1];

  return (\@exon_starts, \@exon_ends);
  
}


sub get_first_func  { 
  
  my $attribs = shift;

  first { defined } map { $attribs->{$_}  } @_;

}


sub get_attrib {

  my ($attribute_string, $species) =@_;
  
  #print "attribute_string=$attribute_string, species=$species\n";
  my %attribs;

  foreach my $id( split( /;/, $attribute_string ) ){
    
      $id =~ s/^\s+//;
      $id =~ s/\s+$//;
    if( $id  =~ m/([^=]*)=(.*)/ ){
      #ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
      #ID=11562.t00100;Name="ORF249";Note="pub_locus=LOC_Osp1g01070"
      
      my $key   = uc $1;
      my $value = $2;
      $value =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
      $value =~ s/^"//;
      $value =~ s/"$//;
      $attribs{$key} = $value;
  
    }else{
      #print "non Poplar\n$GRAPE_REGEX, $GRAPE_FEATURE_NAME_REGEX\n";
      if( $species =~ $GRAPE_REGEX && $id =~ $GRAPE_FEATURE_NAME_REGEX){
	$attribs{ID} = uc $& unless $attribs{ID};
      }
      
    }

  }

  return \%attribs;
}

sub get_gene_id{

  my ($attrib_ref, $species, ) = @_;
  my $parent_id;
  
  if($species =~ $GRAPE_REGEX ){
    
   $parent_id = $attrib_ref->{ID};
   $parent_id =~ s/GSVIV[A-Z]/GSVIVG/;
 
  }else{ 
    
    $parent_id = $attrib_ref->{PARENT};
  }

  return $parent_id;
  
}
1;
