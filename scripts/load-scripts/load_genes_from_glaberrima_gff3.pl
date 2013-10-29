##!/usr/local/bin/perl -w 
 
=pod

=head1 NAME

load_genes_from_jgi_brachy_gff3.pl - load poplar genes into EnsemblDB

=head1 SYNOPSIS

  load_genes_from_jgi_brachy_gff3.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -l|--logic_name       analysis.logic_name to use for the genes.
  -n|--no_insert        Do not make changes to the database. For debug.
  -v|--verbose          Print verbose output
  -source|--source	annotation source such as JGI, MIPS

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
  analysis.logic_name to use for the genes. Default is ensembl.

B<--source>
  annotation source, for xref

B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

B<-v|--verbose>
  Print verbose output

=head1 DESCRIPTION

  Parses the specified JGI genome gff file from JGI, and uses the data
  therein to populate the gene, exon, transcript and translation
  tables (inc stable_id tables) in the ensembl DB indicated by the
  --species and --ensembl_registry.

  Note; the 'poplar' gff is neither 'true' gff or gtf format. It is a
  hybrid. At the moment there are no alternative splices.

  The file can be got from here;
  ftp://ftp.jgi-psf.org/pub/JGI_data/Poplar/annotation/v1.1/Poptr1_1.JamboreeModels.gff.gz

B<The Ensembl Registry>

  The database connection details for both Ensembl and interproscan
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'Populus_trichocarpa', 
    '-group'   => 'core', 
    '-dbname'  => 'populus_trichocarpa_core_42_jgi2004', );
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


Maintained by Will Spooner <whs@ebi.ac.uk>

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
  unshift @INC, '/usr/local/ensembl-live/ensembl/modules';

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
use Readonly;

Readonly my $GENE_REGEX => qr/(PREDICTED_GENE|GENE)/i;
Readonly my $CDS_REGEX => qr/(CDS_PREDICTED|CD)/i;


use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS $V 
             $GENES  %PARENT $ANNOT_SOURCE );


#  0    1    2     3     4    5     6     7     8
#my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
#my $date=sprintf("%04d-%02d-%02d", $year, $mon, $mday);

my $date = time(); 

#print "date is $date\n";

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $logic_name, $no_insert );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "logic_name=s"       => \$logic_name,
        "no_insert"          => \$no_insert,
        "verbose"            => \$V,
	"source=s"	     => \$ANNOT_SOURCE,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args
  
  $ANNOT_SOURCE = 'JGI' unless $ANNOT_SOURCE;

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
  print( "[INFO] Genes will use analysis.logic_name of $logic_name\n" );
  $ANALYSIS = Bio::EnsEMBL::Analysis->new
      ( -logic_name => $logic_name );

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

#----------------------------------------------------------------------
# Main

while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /\#/ ); 
  next if ( $line =~ /^\s+/ );
  chomp $line;

  print "\n$line\n" if $V;

  # Split gff line
  my( $seqname, $source, $feature, 
      $start, $end, $score, $strand, $frame, 
      $attribute ) = split( /\s+/, $line, 9 );
  $feature = uc($feature);
  $seqname =~ s/Bd//i;

  my %attribs;
  foreach my $id( split( /;\s*/, $attribute ) ){
    #print "split attribute $id\n";
    if( $id  =~ m/(\S+)[\s\=]+\"*([^\"]+)/ ){
      $attribs{uc ($1)} = uc($2);
      #print "attrib mapping $1 => $2\n";
    }
  }

  my $id = $attribs{ID} ;
  my $parent = $attribs{PARENT} if defined $attribs{PARENT} ;
  $PARENT{ $id } = $parent if defined $id && defined $parent;
  my $gene_id = $parent ? &get_gene_id( $parent ) : $id;

#print "\n\n[DEBUG] get_gene_id $gene_id, feature=$feature xxxx\n\n";
  if( $feature =~ $GENE_REGEX){
    print "[INFO] Processing gene $gene_id" if $V;
    $GENES->{$gene_id}->{GENE_NAME}  = $id;
    $GENES->{$gene_id}->{SEQ_NAME}   = $seqname;
    $GENES->{$gene_id}->{TRPT_COUNT} = 0;
    $GENES->{$gene_id}->{STRAND}     = $strand eq '+' ? 1 : -1;
  
  }

  if( $feature eq 'MRNA' ){
    
    #$GENES->{$gene_id}->{TRPT}->{TRPT_COUNT}++;
    $GENES->{$gene_id}->{TRPT}->{$id}->{TRPT_NAME}  = $id; 
    $GENES->{$gene_id}->{TRPT}->{$id}->{EXON_COUNT} = 0;
    $GENES->{$gene_id}->{TRPT}->{$id}->{CDS_COUNT}  = 0;
    #$GENES->{$gene_id}->{TRPT}->{}{STRAND}     = $strand eq '+' ? 1 : -1;
  }
 
  if( $feature eq 'EXON' ||
      $feature =~ /UTR/i ||
      $feature =~ $CDS_REGEX
      #$feature eq 'CDS' ||
      #$feature eq 'THREE_PRIME_UTR' ||
      #$feature eq 'FIVE_PRIME_UTR'
    ){
    
    $GENES->{$gene_id}->{TRPT}->{$parent}->{EXON_COUNT} ++;
    #$GENES->{$gene_id}->{TRPT}->{$parent}->{EXON_START} ||= [];
    #$GENES->{$gene_id}->{TRPT}->{$parent}->{EXON_END}   ||= [];
    #push @{$GENES->{$gene_id}->{TRPT}->{$parent}->{EXON_START}}, $start;
    #push @{$GENES->{$gene_id}->{TRPT}->{$parent}->{EXON_END}}, $end;

    $GENES->{$gene_id}->{TRPT}->{$parent}->{EXON}{$start} = $end;
  }

  if( $feature =~ $CDS_REGEX ){
    $GENES->{$gene_id}->{TRPT}->{$parent}->{CDS_COUNT} ++;
    #$GENES->{$gene_id}->{TRPT}->{$parent}->{CDS_EXON_START} ||= [];
    #$GENES->{$gene_id}->{TRPT}->{$parent}->{CDS_EXON_END}   ||= [];
    #push @{$GENES->{$gene_id}->{TRPT}->{$parent}->{CDS_EXON_START}}, $start;
    #push @{$GENES->{$gene_id}->{TRPT}->{$parent}->{CDS_EXON_END}}, $end;
    $GENES->{$gene_id}->{TRPT}->{$parent}->{CDS_EXON}{$start} = $end;
  }

  if( $feature eq 'START_CODON' ){
    $GENES->{$gene_id}->{TRPT}->{$parent}->{START_CODON} = 
      $strand eq "+" ? $start : $end;
  }
  
  if( $feature eq 'STOP_CODON' ){
    $GENES->{$gene_id}->{TRPT}->{$parent}->{STOP_CODON} = 
      $strand eq "+" ? $end : $start;
  }

}

warn ( "[INFO] found " . scalar( keys %$GENES) . " genes in gff file\n" );


foreach my $geneid ( sort keys %$GENES ){

  my $genedata = $GENES->{$geneid};
  $genedata->{TRANSCRIPT_COUNT} = keys %{$genedata->{TRPT}};

  warn ( "[INFO] found $genedata->{TRANSCRIPT_COUNT} transcripts in gene $geneid\n" );

  foreach my $trptid ( sort keys %{$genedata->{TRPT}} ){
    
    my $trptdata = $genedata->{TRPT}->{$trptid};

    # Get exons into gene-order
    my $f = $genedata->{STRAND} > 0 ? 1 : 0;

    unless ( defined $trptdata->{CDS_EXON} ){
	%{$trptdata->{CDS_EXON}} = %{$trptdata->{EXON}};
    }
    
    foreach my $key qw( EXON CDS_EXON ){
      $trptdata->{"${key}_START"} = [ sort{ $f ?$a<=>$b :$b<=>$a} 
				      keys %{$trptdata->{$key}} ];
      $trptdata->{"${key}_END"} = [ sort{ $f ?$a<=>$b :$b<=>$a} 
				    values %{$trptdata->{$key}} ];
    }

#EXON_START EXON_END CDS_EXON_START CDS_EXON_END

  # Find the translation start/stop
    my( $start_codon ) = sort{ $a<=>$b } @{$trptdata->{CDS_EXON_START}};
    my( $stop_codon )  = sort{ $b<=>$a } @{$trptdata->{CDS_EXON_END}};
    if( $genedata->{STRAND} < 0 ){ 
      ($start_codon,$stop_codon) = ($stop_codon,$start_codon);
    }
    unless( $trptdata->{START_CODON} ){
      warn( "[WARN] Gene $geneid has no start_codon\n" ) if $V;
    }
    unless( $trptdata->{STOP_CODON} ){
      warn( "[WARN] Gene $geneid has no stop_codon\n" ) if $V;
    }
    
    if( ( $trptdata->{START_CODON} and 
	  $trptdata->{START_CODON} ne $start_codon ) or
	( $trptdata->{STOP_CODON}  and
	  $trptdata->{STOP_CODON}  ne $stop_codon ) ){
      warn( "[WARN]Start/stop for gene $_ inconsistent with CDS_EXONS\n" ) if $V;
    }
    $trptdata->{START_CODON} = $start_codon;
    $trptdata->{STOP_CODON}  = $stop_codon;

    # Flag exons where coding start/stops
    #my $coding_length = 0;
    for( my $i=0; $i<@{$trptdata->{EXON_START}}; $i++ ){
      my $start = $trptdata->{EXON_START}->[$i];
      my $end   = $trptdata->{EXON_END}->[$i];
      if( $start_codon >= $start and $start_codon <= $end ){
	$trptdata->{CDS_START_EXON} = $i;
      }
      if( $stop_codon >= $start and $stop_codon <= $end ){
	$trptdata->{CDS_END_EXON} = $i;
      }
      #$coding_length += ( $end - $start + 1 ); # Does not allow for UTRs
      #?????? Doesn't make sense here --Sharon
      # Test for overlapping exons
      if( $i > 0 ){
	my $last_start = $trptdata->{EXON_START}->[$i-1];
	my $last_end   = $trptdata->{EXON_END}->[$i-1];
	if( $last_end >= $start and $last_start <= $end ){
	  warn( "[WARN] Gene $_ has overlapping exons. Dubious. Skip gene.\n" );
	  delete( $genedata->{TRPT}->{$trptid} );
	  next;
	}
      }
    }


    #my $residues = $coding_length / 3;

    # Check for peptide length
    #if( $residues < 50 ){
    #  warn( "[WARN] Gene $_ has $residues residues. ".
    #        "Value should ideally be >50.\n" ) if $V;
    #  if( $sorghum ){
    #    delete $GENES->{$_};
    #    next;
    #  }
    #}
    
    # Check for phase
    #if( $residues - int( $residues ) > 0 ){
    #  warn( "[WARN] Gene $_ has $residues residues. ".
    #        "Value should be an integer.\n" ) if $V;
    #}
    
    unless( defined $trptdata->{CDS_END_EXON} &&
	    defined $trptdata->{CDS_START_EXON} ){
      warn Dumper( $trptdata );
      die( "[*DIE] Gene $_ has no CDS_END_EXON/CDS_START_EXON" );
    }
  }
}


#exit;

# Load the data
my $sa = $ENS_DBA->get_adaptor('Slice');
my $n = 1;
my $number = scalar( keys %$GENES );
warn( "[INFO] Storing $number genes to DB\n" );
foreach my $geneid( sort keys %$GENES ){

  my $genedata = $GENES->{$geneid};

  my $slice = $sa->fetch_by_region( undef, $genedata->{SEQ_NAME} );

  unless( $slice ){
      print STDERR "ERROR: No seq_region found for ", $genedata->{SEQ_NAME},"\n";
      next;
  }
  print "seq name is ",$genedata->{SEQ_NAME} , "\n";
  print "slice name is ", $slice->name, "\n";
  my $eGene = Bio::EnsEMBL::Gene->new();
  $eGene->analysis($ANALYSIS);
  $eGene->stable_id( $genedata->{GENE_NAME} );
  #$eGene->description( $genedata->{GENE_DESCRIPTION} || '');
  $eGene->version(1);  
  $eGene->created_date( $date );
  $eGene->modified_date( $date );
  $eGene->slice($slice);

  # Add XREFs to gene
  # The xref primary id cannot be null
  # make gene TU_feat_name as primary id and pub_locus as display_label

  my $gene_src = "${ANNOT_SOURCE}_GENE";

  my $gene_xref = [{
		    db                   => $gene_src,
		    primary_acc          => $genedata->{GENE_NAME},
		    display_label        => $genedata->{GENE_NAME},
		   # description          => $genedata->{GENE_DESCRIPTION}   || '',
		   },
                  ];
  
  foreach my $dbentry( make_dbentries($gene_xref) ){
    $eGene->add_DBEntry( $dbentry );
    if( $dbentry->dbname eq $gene_src ){ # Use as display_xref
      $eGene->display_xref( $dbentry );
    }
  }
  
  

  foreach my $trptid( sort keys %{$genedata->{TRPT}} ){
    
    my $trptdata = $genedata->{TRPT}->{$trptid};
    my $trpt_name =  $trptdata->{TRPT_NAME};

    my $eTranscript = Bio::EnsEMBL::Transcript->new();
    $eTranscript->stable_id( $trpt_name );
    $eTranscript->version( 1 );
    $eTranscript->analysis( $ANALYSIS );
    $eTranscript->created_date( $date );
    $eTranscript->modified_date( $date );
    $eTranscript->slice( $slice );

    my $trpt_src = "${ANNOT_SOURCE}_TRANSCRIPT"; 
    my $transcript_xref = [{
			    db            => $trpt_src,
			    primary_acc   => $trpt_name,
			    display_label => $trpt_name || '',
			    #description   => $genedata->{GENE_DESCRIPTION}  || '',
			   },
			  ];

    # Transcript xrefs...
    foreach my $dbentry( make_dbentries($transcript_xref) ){
      $eTranscript->add_DBEntry( $dbentry );
      if( $dbentry->dbname eq $trpt_src ){ # Use as display_xref
        $eTranscript->display_xref( $dbentry );
      }
    }

    #my $slice = $sa->fetch_by_region( undef, $genedata->{SEQ_NAME} );
    
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
    
    my $translation_start = $trptdata->{START_CODON};
    my $translation_stop  = $trptdata->{STOP_CODON};
    
    $translation_start || die( "Gene $genedata->{GENE_NAME}-$trpt_name has no start" );
    $translation_stop  || die( "Gene $genedata->{GENE_NAME}-$trpt_name has no stop" );
    
    #########
    # EXONS #
    #########
    
    for ( my $exon = 0; $exon < $trptdata->{EXON_COUNT}; $exon++ ){
      
      my $exon_start  = $trptdata->{EXON_START}[$exon];
      my $exon_end    = $trptdata->{EXON_END}[$exon];
      
      my $eExon = new Bio::EnsEMBL::Exon;
      $eExon->start($exon_start);
      $eExon->end($exon_end);
      $eExon->strand($genedata->{STRAND});	   
      $eExon->slice($slice);
      $eExon->stable_id( $trptid.'.exon'.($exon+1) );
      $eExon->version( 1 );
      
      # Phase calculations 
      
      if ( $genedata->{STRAND} > 0 ) {
	if ( ( $translation_start > $exon_start ) && 
	     ( $translation_stop > $exon_end ) )  { # 5' exons
	  $exon_start_phase = -1;
	  $exon_end_phase   = -1;
	}
	elsif( ( $translation_start == $exon_start ) && 
	       ( $translation_start < $exon_end ) ) { # exon starts with a start
	  $phase_diff =  ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = 0;
	  $exon_end_phase   = $phase_diff;
	}
	elsif ( ( $translation_start > $exon_start  ) && 
		( $translation_start < $exon_end ) ) { # exon contains a start
	  $phase_diff = ( ( $exon_end - $translation_start + 1 ) % 3 );
	  $exon_start_phase = - 1;
	  $exon_end_phase   = $phase_diff;
	}
	elsif ( ( $translation_stop > $exon_start ) && 
		( $translation_stop < $exon_end ) ) { # exon contains a stop
	  $phase_diff = ( ( $exon_end - $translation_stop + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = -1;
	}
	elsif ( ( $translation_stop == $exon_end ) && 
		( $translation_stop > $exon_start ) ) { # exon stops with a stop
	  $phase_diff = ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = 0;
	}
	elsif ( ( $translation_stop < $exon_start ) && 
		( $translation_stop < $exon_end ) ) { # 3' exons
	  $exon_start_phase = -1;
	  $exon_end_phase   = -1;
	}
	elsif ( ( $translation_start == $exon_start ) && 
		( $translation_stop == $exon_end ) ) { # single exon genes
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
	if    ( ( $translation_start < $exon_start ) && 
		( $translation_start < $exon_end ) )  {
	  $exon_start_phase = -1;
	  $exon_end_phase   = -1;
	}
	# exon stops with a start 
	elsif ( ( $translation_start == ($exon_end - 2) ) ) {  
	  $phase_diff = ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = 0;
	  $exon_end_phase   = $phase_diff;
	}
	# exon contains a start codon
	elsif ( ( $translation_start > $exon_start ) && 
		( $translation_start < $exon_end ) ) {
	  $phase_diff = ( ( $translation_start - $exon_start + 1 ) % 3 ); #
	  #( ( $exon_end - $translation_start + 1 ) % 3 );
	  $exon_start_phase = - 1;
	  $exon_end_phase   = $phase_diff;
	}
	# exon contains a stop codon
	elsif ( ( $translation_stop > $exon_start ) && 
		( $translation_stop < $exon_end ) ) {
	  $phase_diff = ( ( $exon_end - $translation_stop + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = -1;
	}
	# exon stops with a stop
	elsif ( ( $translation_stop == $exon_start ) ) { 
	  $phase_diff = ( ( $exon_end - $exon_start + 1 ) % 3 );
	  $exon_start_phase = $last_exon_phase;
	  $exon_end_phase   = 0;
	}
	# 3' exons
	elsif ( ( $translation_stop > $exon_start ) && 
		( $translation_stop > $exon_end ) ) {
	  $exon_start_phase = -1;
	  $exon_end_phase   = -1;
	}
	# single exon genes
	elsif ( ( $translation_start == $exon_start ) && 
		( $translation_stop == $exon_end ) ) {
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
    if ( $genedata->{STRAND} > 0 ) {
      $translation_offset_left  = 
        $translation_start 
	  - $trptdata->{EXON_START}[$trptdata->{CDS_START_EXON}] + 1;
      $translation_offset_right = 
        $translation_stop  
	  - $trptdata->{EXON_START}[$trptdata->{CDS_END_EXON}] + 1;
      #- $trptdata->{EXON_START}[$trptdata->{CDS_END_EXON}] + 3; 
    }
    elsif ( $genedata->{STRAND} < 0 ) {
      $translation_offset_left  = 
        $trptdata->{EXON_END}[$trptdata->{CDS_START_EXON}] 
	  - $translation_start + 1;
      $translation_offset_right = 
        $trptdata->{EXON_END}[$trptdata->{CDS_END_EXON}]   
	  - $translation_stop  + 1; #was + 3;
    }
    
      #print "// [g] TranscriptLength ", $eTranscript->length, ' ,', $eTranscript->slice->name, " [$translation_offset_left - $translation_offset_right]\n" if 1;#( $verbose);
    
    
    #$eTranslation_ID = $trptdata->{TRANSCRIPT_NAME};
    #$eTranslation_ID =~ s/-R/-P/;
    
    my $eTranslation = new  Bio::EnsEMBL::Translation
      (
       -START_EXON  => $trans_start_exon,
       -END_EXON    => $trans_end_exon,
       -SEQ_START   => $translation_offset_left,
       -SEQ_END     => $translation_offset_right,
       -STABLE_ID   => $trpt_name,
       -VERSION     => '1',
      );
    
    #########################################
    # EnsEMBL add translation to transcript #
    #########################################
    #print "before translation\n";
    $eTranscript->translation($eTranslation);
    
    ##################################
    # EnsEMBL add transcript to gene #
    ##################################
    
    $eGene->add_Transcript($eTranscript);
    
  }
  
  my $dbid = '???';
  if( $I ){
    $dbid = $ENS_DBA->get_GeneAdaptor->store($eGene);
    print "// dbID = $dbid for gene = $genedata->{GENE_NAME} ($n of $number)\n";
  }

  $n++;

}

warn( "[INFO] Stored $number genes to DB\n" );
warn( "[INFO] DONE!\n" );

sub get_gene_id{

  my $id = shift;

  #Find gene source
  
  my $gene = $id;
  while( defined $PARENT{ $id } and $gene = $PARENT{ $id }){
    $id = $gene;
  }    

  return $gene;
}



#----------------------------------------------------------------------
# Ensembl DBEntry objects are used for gene/transcript etc xrefs 
#
# There are only cDNA evidence xref in the tigrv4 xml
#
# create entry in xref table
# Before making it, need to check whether the same entry already exist, 
# if does, use the existing one instead of making new one

sub make_dbentries{

  my $xref  = shift;
  my @data  = @{ $xref || [] };
  my @xrefs = ();

  foreach my $entry( @data ){
    push @xrefs, Bio::EnsEMBL::DBEntry->new
        (
         -dbname      => $entry->{db} || '',
         -primary_id  => $entry->{primary_acc} || 'noid',
         -display_id  => $entry->{display_label} || '', # Filed required; Ens v28
         -description => $entry->{description} || '',
         -version     => 1,
         -release     => 1,
         );
  }
  return @xrefs;

  
}


1;
