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
  -o|--old_species      Species key for previous database UNUSED.
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

B<-o|--old_species>
  UNUSED 
  Use if you want to reuse gene ID from a previous database. IDs will
  be reused if the new gene model is in the same location on the same 
  seq_region as the old one.

B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.


=head1 DESCRIPTION

  Parses the specified 'grape' gff file from International Grape Genome Project, 
  and uses the data therein to populate the gene, exon, transcript and translation
  tables (inc stable_id tables) in the ensembl DB indicated by the
  --species and --ensembl_registry.

  The file can be got from here;
  http://www.genoscope.cns.fr/externe/English/Projets/Projet_ML/data/annotation/Vitis_vinifera_annotation_v1.gff
  

  At the moment there are no alternative splices.



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

use vars qw( $BASEDIR );
BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
  #unshift @INC, $BASEDIR.'/bioperl-live';
  #unshift @INC, '/usr/local/bioperl-HEAD';
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

# make it a function for reuse,
# could parse in an anonymous sub from command line 
# "sub {}", then $func->() to invoke
# for grape gff
  #Gene GSVIVG00000001001
  #Trpt GSVIVT00000001001
  #prtn GSVIVP00000001001

sub parse_gene_name {

  my $prefix = "GSVIV";
  if(/${prefix}[A-Z](\d+)/i){
    return "${prefix}G$1";
  }else{
    return undef;
  }
}


use vars qw( $ENS_DBA $ENS_DBA_OLD $I $INSERT $GFF_HANDLE $ANALYSIS );

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $old_species, $reg, $logic_name, $no_insert );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "old_species=s"      => \$old_species,
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

  if( $old_species ){
    $ENS_DBA_OLD = Bio::EnsEMBL::Registry->get_DBAdaptor($old_species,'core');
    $ENS_DBA_OLD || ( print( "No core DB for $species set in $reg\n" ) &&
                      pod2usage() );
  }

  # Create the analysis
  $ANALYSIS = Bio::EnsEMBL::Analysis->new
      ( -logic_name => $logic_name ); #'GeneModel_RiceIndica_BGI' 

  # Create a GFF stream
  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


  print( "Loading genes for $species\n" );
  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA, ( $ENS_DBA_OLD || () ) ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
    $pre_text = "  Old DB: ";
  }
}


my $GENES = {};
my $TRPT2GENE = {};
#my %EXONS;

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


  my %attribs;
  foreach my $id( split( /(;\s*)/, $attribute ) ){
    if( $id  =~ m/(\S+)\s+\"*([^\" ]+)/ ){
      $attribs{uc($1)} = $2;
    }
  }

  #my( $type, $feature_name ) = split( /\s+/, $attribute );

  
  #get gene name for this feature
  my ($gene, $transcript);

  if( $feature eq 'GENE' ){
    $gene = $attribs{'GENE'};
    print "Gene $gene\n";#####
    unless( $GENES->{$gene} ){
      $GENES->{$gene}->{GENE_NAME} = $gene;
      $GENES->{$gene}->{SEQ_NAME}  = $seqname;
      $GENES->{$gene}->{START} = $start; #always smaller than end
      $GENES->{$gene}->{END}  = $end;
      $GENES->{$gene}->{STRAND}     = $strand eq '+' ? 1 : -1;
    }
  }
  elsif( $feature eq 'MRNA' ){
    $gene = $attribs{'GENE'};
    unless( $GENES->{$gene} ){
      $GENES->{$gene}->{GENE_NAME} = $gene;
      $GENES->{$gene}->{SEQ_NAME}  = $seqname;
      $GENES->{$gene}->{START} = $start; #always smaller than end
      $GENES->{$gene}->{END}  = $end;
      $GENES->{$gene}->{STRAND}     = $strand eq '+' ? 1 : -1;
    }

    $transcript = $attribs{'MRNA'};
    unless( $GENES->{$gene}->{TRANSCRIPTS}->{$transcript} ){
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{TRPT_NAME} = $transcript;
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{SEQ_NAME} = $seqname;
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{START} = $start; 
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{END}  = $end;
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{EXON_COUNT} = 0;
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{CDS_COUNT} = 0;
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{UTR_COUNT} = 0;
      $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{STRAND} = $strand eq '+' ? 1 : -1;
    }
    
    #establish the mapping from transcript to gene
    $TRPT2GENE->{$transcript} = $gene;
    print"mRNA $transcript => $gene, $TRPT2GENE->{$transcript}\n";#####
    
  }
  #else{
    #$gene = parse_gene_name($feature_name);
  #}

  #$gene or die "Cannot parse gene name $transcript";


  

  #if( $feature !~ /^gene$/i ){
  #  $GENES->{$gene}->{TRANSCRIPT_IDS}->{$feature_name} ++;
  #}
  
    

  elsif( $feature eq 'CDS' ){ # use CDS and UTR as a proxy for exon
    #$GENES->{$gene}->{EXON_COUNT} ++;
    #$GENES->{$gene}->{EXON_START} ||= [];
    #$GENES->{$gene}->{EXON_END}   ||= [];
    #push @{$GENES->{$gene}->{EXON_START}}, $start;
    #push @{$GENES->{$gene}->{EXON_END}}, $end;

 #warn(Dumper($TRPT2GENE));
    $transcript = $attribs{'CDS'};
    $gene = $TRPT2GENE->{$transcript};
    print"CDS $transcript, $gene\n";#####
    unless ( $gene ){
      die("cannot find gene for transcript $transcript at " . 
	  $GFF_HANDLE->input_line_number );
      
    }
    $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{CDS_COUNT} ++;
    $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{CDS_EXON_START} ||= [];
    $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{CDS_EXON_END}   ||= [];
    push @{$GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{CDS_EXON_START}}, $start;
    push @{$GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{CDS_EXON_END}}, $end;

  }

  elsif( $feature eq 'UTR' ){
    $transcript = $attribs{'UTR'};
    $gene = $TRPT2GENE->{$transcript};
    print"UTR $transcript, $gene\n";#####
    unless ( $gene ){
      die("cannot find gene for transcript $transcript at ".
	 $GFF_HANDLE->input_line_number );
      
    }
    $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{UTR_COUNT} ++;
    $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{UTR_EXON_START} ||= [];
    $GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{UTR_EXON_END}   ||= [];
    push @{$GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{UTR_EXON_START}}, $start;
    push @{$GENES->{$gene}->{TRANSCRIPTS}->{$transcript}->{UTR_EXON_END}}, $end;
    
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
  #delete( $genedata->{TRANSCRIPT_IDS} );
  #delete( $genedata->{PROTEIN_IDS} );  
  if( $genedata->{TRANSCRIPT_COUNT} != 1 ){
    warn( "transcript count for gene $g != 1. ".
          "try handle this case ??" );
    #delete( $GENES->{$_} );
    #next;
  }
  
  
  foreach my $t( keys %{$genedata->{TRANSCRIPTS}} ){

    my $trptdata = $genedata->{TRANSCRIPTS}->{$t};

    #reconstruct exons and save in the hash
    ## now we need to compute exons using CDS and UTRS

    my @trpt_ex_starts = sort{ $a<=>$b } (@{$trptdata->{CDS_EXON_START}}, 
					  $trptdata->{UTR_EXON_START} ?
					  @{$trptdata->{UTR_EXON_START}}: () 
					 );
    my @trpt_ex_ends   = sort{ $a<=>$b } (@{$trptdata->{CDS_EXON_END}},
				       @{$trptdata->{UTR_EXON_END} || []}
				      );
    #my @trpt_utr_starts = sort{ $a<=>$b } @{$trptdata->{UTR_EXON_START}};
    #my @trpt_utr_ends = sort{ $a<=>$b } @{$trptdata->{UTR_EXON_END}};
	
    my ($trpt_exon_starts, $trpt_exon_ends) = compute_exons(\@trpt_ex_starts,
							    \@trpt_ex_ends);
    
    die "cannot compute exons for $g, $t" 
      unless( $trpt_exon_starts &&
	      $trpt_exon_ends &&
	      scalar @{$trpt_exon_starts} == scalar @{$trpt_exon_ends});
    #$GENES->{$g}->{TRANSCRIPTS}->{$t}->{EXON_COUNT} = scalar @{$trpt_exon_start};
    #$GENES->{$g}->{TRANSCRIPTS}->{$t}->{EXON_START} = $trpt_exon_starts;
    #$GENES->{$g}->{TRANSCRIPTS}->{$t}->{EXON_END}   = $trpt_exon_ends;
    
    #get updated transcript
    #$trptdata = $genedata->{TRANSCRIPTS}->{$t};


    $trptdata->{EXON_COUNT} = scalar @{$trpt_exon_starts};
    $trptdata->{EXON_START} = $trpt_exon_starts;
    $trptdata->{EXON_END}   = $trpt_exon_ends;
        
    # Get exons into gene-order
    my $f = $trptdata->{STRAND} > 0 ? 1 : 0;
    foreach my $key qw( EXON_START EXON_END CDS_EXON_START CDS_EXON_END ){
      $trptdata->{$key} = [ sort{ $f ?$a<=>$b :$b<=>$a} @{$trptdata->{$key}} ];
    }
    


    # Find the translation start/stop
    my( $start_codon ) = sort{ $a<=>$b } @{$trptdata->{CDS_EXON_START} || []};
    my( $stop_codon )  = sort{ $b<=>$a } @{$trptdata->{CDS_EXON_END} || []};
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
      warn Dumper( $trptdata );
      die( "Gene-Transcript ${g}-${t} has no CDS_END_EXON/CDS_START_EXON" );
    }
    
  }
}


# Load the data
my $sa = $ENS_DBA->get_adaptor('Slice');
my $n = 1;
my $number = scalar( keys %$GENES );
foreach my $gene( keys %$GENES ){
  my $genedata = $GENES->{$gene};
  my $eGene = Bio::EnsEMBL::Gene->new
      ( 
        -ANALYSIS  => $ANALYSIS,
        -STABLE_ID => $gene,
        -VERSION   => '1',
        #-TYPE      => $genedata->{TYPE} || die( "No TYPE for $gene" ) 
        );

  foreach my $trpt (keys %{$genedata->{TRANSCRIPTS}} ){
    my $trptdata = $genedata->{TRANSCRIPTS}->{$trpt};
    my $eTranscript = Bio::EnsEMBL::Transcript->new();
    $eTranscript->stable_id( $trpt );
    $eTranscript->version( 1 );
    $eTranscript->analysis( $ANALYSIS );
    
                              #warn Dumper( $trptdata );
    my $seq_region_name = $trptdata->{SEQ_NAME};
                              #print "seq_region_name=$seq_region_name\n";
    $seq_region_name =~ s/chr//i;
                              #print "seq_region_name=$seq_region_name\n";
    my $slice = $sa->fetch_by_region( undef, $seq_region_name );

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
    
    my $transcript_start = $trptdata->{START_CODON};
    my $transcript_stop  = $trptdata->{STOP_CODON};
    
    $transcript_start || die( "Trpt ${gene}-${trpt} has no start" );
    $transcript_stop  || die( "Trpt ${gene}-${trpt} has no stop" );
    
    #########
    # EXONS #
    #########
    
    for ( my $exon = 0; $exon < $trptdata->{EXON_COUNT}; $exon++ ){
      
      my $exon_start  = $trptdata->{EXON_START}[$exon];
      my $exon_end    = $trptdata->{EXON_END}[$exon];
      
      my $eExon = new Bio::EnsEMBL::Exon;
      $eExon->start($exon_start);
      $eExon->end($exon_end);
      $eExon->strand($trptdata->{STRAND});	   
      $eExon->slice($slice);
      $eExon->stable_id( $trpt.'.exon'.($exon+1) );
      $eExon->version( 1 );
      
      # Phase calculations 
      
      if ( $trptdata->{STRAND} > 0 ) {
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
	$eExon->end_phase($exon_end_phase);
	
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
	$eExon->end_phase($exon_end_phase);
	
	#$span = $exon_end - $exon_start + 1;      
	$last_exon_phase = $exon_end_phase;
      }
      
      $eTranscript->add_Exon($eExon);
      
      # Start exon
      unless ( $start_exon ) {
	if ( $exon == 0 ) {
	  $start_exon = $eExon; 
	}
      }
      
      # Final exon
      unless ( $end_exon ) {
	if ( $exon == $trptdata->{EXON_COUNT}-1 ) {
	  $end_exon = $eExon;
	}
      }
      
      # Translation start exon
      unless ( $trans_start_exon ) {
	if ( $exon == $trptdata->{CDS_START_EXON} ) {
	  $trans_start_exon = $eExon;
	}
      }
      
      # Translation stop exon
      unless ( $trans_end_exon ) {
	if ( $exon == $trptdata->{CDS_END_EXON} ) {
	  $trans_end_exon = $eExon;
	}
      }
    } #_ END of exon loop _#
    
    $eTranscript->start_Exon($start_exon);
    $eTranscript->end_Exon($end_exon);
    
    ###############
    # TRANSLATION #
    ###############
    
    # Check for translation start being different to Exon start (i.e. UTR)
    if ( $trptdata->{STRAND} > 0 ) {
      $translation_offset_left  = 
        $transcript_start 
	  - $trptdata->{EXON_START}[$trptdata->{CDS_START_EXON}] + 1;
      $translation_offset_right = 
        $transcript_stop  
	  #- $trptdata->{EXON_START}[$trptdata->{CDS_END_EXON}] + 3; ##???why +3 
	  - $trptdata->{EXON_START}[$trptdata->{CDS_END_EXON}] + 1;
    }
    elsif ( $trptdata->{STRAND} < 0 ) {
      $translation_offset_left  = 
        $trptdata->{EXON_END}[$trptdata->{CDS_START_EXON}] 
	  - $transcript_start + 1;
      $translation_offset_right = 
        $trptdata->{EXON_END}[$trptdata->{CDS_END_EXON}]   
	 # - $transcript_stop  + 3; ##???why +3
	   - $transcript_stop  + 1;
    }
    
    #  print "// [$gene] TranscriptLength ", $eTranscript->length, " [$translation_offset_left - $translation_offset_right]\n" if 1;#( $verbose);
    
    
    #$eTranslation_ID = $genedata->{TRANSCRIPT_NAME};
    #$eTranslation_ID =~ s/-R/-P/;
    
    my $eTranslation = new  Bio::EnsEMBL::Translation
      (
       -START_EXON  => $trans_start_exon,
       -END_EXON    => $trans_end_exon,
       -SEQ_START   => $translation_offset_left,
       -SEQ_END     => $translation_offset_right,
       -STABLE_ID   => $gene,
       -VERSION     => '1',
      );
    
    #########################################
    # EnsEMBL add translation to transcript #
    #########################################
    
    $eTranscript->translation($eTranslation);

    ##################################
    # EnsEMBL add transcript to gene #
    ##################################
    
    $eGene->add_Transcript($eTranscript);
  }
 
  if( $I ){
    my $dbid = $ENS_DBA->get_GeneAdaptor->store($eGene);
    print "// dbID = $dbid for gene = $gene ($n of $number)\n";
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
    if( $trpt_ex_starts->[$i] - $trpt_ex_ends->[$i-1] == 1){
         
    }else{
      push @exon_ends, $trpt_ex_ends->[$i-1];
      push @exon_starts, $trpt_ex_starts->[$i];
    }
  }

  push @exon_ends, $trpt_ex_ends->[$start_cnt-1];

  return (\@exon_starts, \@exon_ends);
  
}



1;
