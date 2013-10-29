##!/usr/local/bin/perl -w

=pod

=head1 NAME

load_genes_from_jgi_gff.pl - load gff genes into EnsemblDB, tailor for JGI poplarv2.0 gff3

=head1 SYNOPSIS

  load_genes_from_jgi_gff.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -l|--logic_name       analysis.logic_name to use for the genes.
  -n|--no_insert        Do not make changes to the database. For debug.
  -v|--verbose          Print verbose output

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

  An exmaple of the gff3 format is:
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


Written by Will Spooner <wspooner@cshl.edu>
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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;


use vars qw( $BASEDIR );
use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS $V );

BEGIN {

    # Set the perl libraries
    $BASEDIR = dirname($Bin);
    unshift @INC, $BASEDIR . '/ensembl-live/ensembl/modules';
    unshift @INC, $BASEDIR . '/ensembl-live/ensembl-compara/modules';
    #unshift @INC, '/usr/local/ensembl-live/ensembl/modules';

    #Argument Processing
    my $help = 0;
    my $man  = 0;
    my ( $species, $reg, $logic_name, $no_insert );
    GetOptions(
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "logic_name=s"       => \$logic_name,
        "no_insert"          => \$no_insert,
        "verbose"            => \$V
    ) or pod2usage(2);
    pod2usage( -verbose => 2 ) if $man;
    pod2usage(1) if $help;

    # Process args

    my $gff_file = $ARGV[0] || pod2usage("\nNeed the path to a gff file\n");
    $species                || pod2usage("\nNeed a --species\n");
    $reg                    ||= $BASEDIR . '/conf/ensembl.registry';
    $logic_name             ||= 'ensembl';

    map {
        -e $_ || pod2usage("\nFile $_ does not exist\n");
        -r $_ || pod2usage("\nCannot read $_\n");
        -f $_ || pod2usage("\nFile $_ is not plain-text\n");
        -s $_ || pod2usage("\nFile $_ is empty\n");
    } $reg, $gff_file;

    # Put stuff in the database?
    $I = $no_insert ? 0 : 1;

    # Load the ensembl file
    Bio::EnsEMBL::Registry->load_all($reg);
    $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
    $ENS_DBA
      || ( print("No core DB for $species set in $reg\n")
        && pod2usage() );

    # Create the analysis
    print("[INFO] Genes will use analysis.logic_name of $logic_name\n");
    $ANALYSIS = Bio::EnsEMBL::Analysis->new( -logic_name => $logic_name );

    # Create a GFF stream
    $GFF_HANDLE = IO::File->new("< $gff_file")
      or die("Could not read $gff_file: $!");

    print("Loading genes for $species\n");
    my $pre_text = "  Target DB: ";
    foreach my $dba ($ENS_DBA) {
        print(  $pre_text
              . $dba->dbc->dbname . "@"
              . $dba->dbc->host . ":"
              . $dba->dbc->port
              . "\n" );
    }
}

#----------------------------------------------------------------------
# Main

my $GENES   = {};

while ( my $line = $GFF_HANDLE->getline ) {

    # Skip comment and empty lines
    next if ( $line =~ /\#/ );
    next if ( $line =~ /^\s+/ );
    chomp $line;

    # Split gff line
    my (
        $seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split( /\s+/, $line, 9 );
    $feature = uc($feature);

    my %attribs;
    foreach my $id ( split( /;\s*/, $attribute ) ) {
        if ( $id =~ m/(\S+)[\s\=]+\"*([^\"]+)/ ) {
            $attribs{$1} = $2;
        }
    }

    if ($sorghum) {
        $attribs{name} ||= $attribs{ID};        # JGI Sorghum GFF mRNA
        $attribs{name} ||= $attribs{Parent};    # JGI Sorghum GFF CDS
    }

    my $gene = $attribs{name}
      || die( "[*DIE] No name attrib here: $attribute\n"
          . "       Complete line: $line" );

    if ($sorghum) {
        $gene =~ s/^Gene\://;                   # Strip any 'Gene:' prefix
    }

    unless ( $GENES->{$gene} ) {
        warn "[INFO] Processing gene $gene" if $V;
        $GENES->{$gene}->{GENE_NAME}  = $gene;
        $GENES->{$gene}->{SEQ_NAME}   = $seqname;
        $GENES->{$gene}->{EXON_COUNT} = 0;
        $GENES->{$gene}->{CDS_COUNT}  = 0;
        $GENES->{$gene}->{STRAND}     = $strand eq '+' ? 1 : -1;
    }

    if ( my $tran = $attribs{transcriptId} ) {
        $GENES->{$gene}->{TRANSCRIPT_IDS}->{$tran}++;
    }
    if ( my $prot = $attribs{proteinId} ) {
        $GENES->{$gene}->{PROTEIN_IDS}->{$prot}++;
    }

    if ($sorghum) {
        $GENES->{$gene}->{TRANSCRIPT_IDS}->{$gene}++;
        $GENES->{$gene}->{PROTEIN_IDS}->{$gene}++;
    }

    if ( $feature eq 'EXON' or ( $sorghum and $feature eq 'CDS' ) ) {
        $GENES->{$gene}->{EXON_COUNT}++;
        $GENES->{$gene}->{EXON_START} ||= [];
        $GENES->{$gene}->{EXON_END}   ||= [];
        push @{ $GENES->{$gene}->{EXON_START} }, $start;
        push @{ $GENES->{$gene}->{EXON_END} },   $end;
    }
    if ( $feature eq 'CDS' ) {
        $GENES->{$gene}->{CDS_COUNT}++;
        $GENES->{$gene}->{CDS_EXON_START} ||= [];
        $GENES->{$gene}->{CDS_EXON_END}   ||= [];
        push @{ $GENES->{$gene}->{CDS_EXON_START} }, $start;
        push @{ $GENES->{$gene}->{CDS_EXON_END} },   $end;
        if ( !$attribs{exonNumber} and !$sorghum ) {
            warn("[WARN] Gene $gene has CDS with no exonNumber\n");
        }
        if ( !$attribs{proteinId} and !$sorghum ) {
            warn("[WARN] Gene $gene has CDS with no proteinId\n");
        }
    }
    elsif ( $feature eq 'START_CODON' ) {
        $GENES->{$gene}->{START_CODON} = $strand eq "+" ? $start : $end;
    }
    elsif ( $feature eq 'STOP_CODON' ) {
        $GENES->{$gene}->{STOP_CODON} = $strand eq "+" ? $end : $start;
    }
    elsif ( $sorghum and $feature eq 'mRNA' and $attribs{Note} ) {

        # Sorghum gene mRNA
        $GENES->{$gene}->{DESC} = $attribs{Note};
    }

}

warn( "[INFO] found " . scalar( keys %$GENES ) . " genes in gff file\n" );

foreach ( sort keys %$GENES ) {
    my $genedata = $GENES->{$_};
    $genedata->{TRANSCRIPT_COUNT} = keys %{ $genedata->{TRANSCRIPT_IDS} };
    $genedata->{PROTEIN_COUNT}    = keys %{ $genedata->{PROTEIN_IDS} };

    #delete( $genedata->{TRANSCRIPT_IDS} );
    #delete( $genedata->{PROTEIN_IDS} );
    if (   $genedata->{TRANSCRIPT_COUNT} != 1
        or $genedata->{PROTEIN_COUNT} != 1 )
    {
        warn(   "[WARN] Transcript/protein count for gene $_ != 1. "
              . "Cannot yet handle this case\n" );
        delete( $GENES->{$_} );
        next;
    }

    # Get exons into gene-order
    my $f = $genedata->{STRAND} > 0 ? 1 : 0;
    foreach my $key qw( EXON_START EXON_END CDS_EXON_START CDS_EXON_END ) {
        $genedata->{$key} =
          [ sort { $f ? $a <=> $b : $b <=> $a } @{ $genedata->{$key} } ];
    }

    # Find the translation start/stop
    my ($start_codon) = sort { $a <=> $b } @{ $genedata->{CDS_EXON_START} };
    my ($stop_codon)  = sort { $b <=> $a } @{ $genedata->{CDS_EXON_END} };
    if ( $genedata->{STRAND} < 0 ) {
        ( $start_codon, $stop_codon ) = ( $stop_codon, $start_codon );
    }
    unless ( $genedata->{START_CODON} ) {
        warn("[WARN] Gene $_ has no start_codon\n") if $V;
    }
    unless ( $genedata->{STOP_CODON} ) {
        warn("[WARN] Gene $_ has no stop_codon\n") if $V;
    }

    if (
        (
                $genedata->{START_CODON}
            and $genedata->{START_CODON} ne $start_codon
        )
        or (    $genedata->{STOP_CODON}
            and $genedata->{STOP_CODON} ne $stop_codon )
      )
    {
        warn("[WARN]Start/stop for gene $_ inconsistent with CDS_EXONS\n")
          if $V;
    }
    $genedata->{START_CODON} = $start_codon;
    $genedata->{STOP_CODON}  = $stop_codon;

    # Flag exons where coding start/stops
    #my $coding_length = 0;
    for ( my $i = 0 ; $i < @{ $genedata->{EXON_START} } ; $i++ ) {
        my $start = $genedata->{EXON_START}->[$i];
        my $end   = $genedata->{EXON_END}->[$i];
        if ( $start_codon >= $start and $start_codon <= $end ) {
            $genedata->{CDS_START_EXON} = $i;
        }
        if ( $stop_codon >= $start and $stop_codon <= $end ) {
            $genedata->{CDS_END_EXON} = $i;
        }

        #$coding_length += ( $end - $start + 1 ); # Does not allow for UTRs
        #?????? Doesn't make sense here --Sharon
        # Test for overlapping exons
        if ( $i > 0 ) {
            my $last_start = $genedata->{EXON_START}->[ $i - 1 ];
            my $last_end   = $genedata->{EXON_END}->[ $i - 1 ];
            if ( $last_end >= $start and $last_start <= $end ) {
                warn(
"[WARN] Gene $_ has overlapping exons. Dubious. Skip gene.\n"
                );
                delete( $GENES->{$_} );
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

    unless ( defined $genedata->{CDS_END_EXON}
        && defined $genedata->{CDS_START_EXON} )
    {
        warn Dumper($genedata);
        die("[*DIE] Gene $_ has no CDS_END_EXON/CDS_START_EXON");
    }
}

# Load the data
my $sa     = $ENS_DBA->get_adaptor('Slice');
my $n      = 1;
my $number = scalar( keys %$GENES );
warn("[INFO] Storing $number genes to DB\n");
foreach my $gene ( sort keys %$GENES ) {
    my $genedata = $GENES->{$gene};
    my $eGene    = Bio::EnsEMBL::Gene->new(
        -ANALYSIS  => $ANALYSIS,
        -STABLE_ID => $gene,
        -VERSION   => '1',

        #-TYPE      => $genedata->{TYPE} || die( "No TYPE for $gene" )
    );
    if ( my $d = $genedata->{DESC} ) { $eGene->description($d) }
    my $eTranscript = Bio::EnsEMBL::Transcript->new();
    $eTranscript->stable_id($gene);
    $eTranscript->version(1);
    $eTranscript->analysis($ANALYSIS);

    my $slice = $sa->fetch_by_region( undef, $genedata->{SEQ_NAME} );

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

    my $translation_start = $genedata->{START_CODON};
    my $translation_stop  = $genedata->{STOP_CODON};

    $translation_start || die("Gene $gene has no start");
    $translation_stop  || die("Gene $gene has no stop");

    #########
    # EXONS #
    #########

    for ( my $exon = 0 ; $exon < $genedata->{EXON_COUNT} ; $exon++ ) {

        my $exon_start = $genedata->{EXON_START}[$exon];
        my $exon_end   = $genedata->{EXON_END}[$exon];

        my $eExon = new Bio::EnsEMBL::Exon;
        $eExon->start($exon_start);
        $eExon->end($exon_end);
        $eExon->strand( $genedata->{STRAND} );
        $eExon->slice($slice);
        $eExon->stable_id( $gene . '.exon' . ( $exon + 1 ) );
        $eExon->version(1);

        # Phase calculations

        if ( $genedata->{STRAND} > 0 ) {
            if (   ( $translation_start > $exon_start )
                && ( $translation_stop > $exon_end ) )
            {    # 5' exons
                $exon_start_phase = -1;
                $exon_end_phase   = -1;
            }
            elsif (( $translation_start == $exon_start )
                && ( $translation_start < $exon_end ) )
            {    # exon starts with a start
                $phase_diff       = ( ( $exon_end - $exon_start + 1 ) % 3 );
                $exon_start_phase = 0;
                $exon_end_phase   = $phase_diff;
            }
            elsif (( $translation_start > $exon_start )
                && ( $translation_start < $exon_end ) )
            {    # exon contains a start
                $phase_diff = ( ( $exon_end - $translation_start + 1 ) % 3 );
                $exon_start_phase = -1;
                $exon_end_phase   = $phase_diff;
            }
            elsif (( $translation_stop > $exon_start )
                && ( $translation_stop < $exon_end ) )
            {    # exon contains a stop
                $phase_diff = ( ( $exon_end - $translation_stop + 1 ) % 3 );
                $exon_start_phase = $last_exon_phase;
                $exon_end_phase   = -1;
            }
            elsif (( $translation_stop == $exon_end )
                && ( $translation_stop > $exon_start ) )
            {    # exon stops with a stop
                $phase_diff       = ( ( $exon_end - $exon_start + 1 ) % 3 );
                $exon_start_phase = $last_exon_phase;
                $exon_end_phase   = 0;
            }
            elsif (( $translation_stop < $exon_start )
                && ( $translation_stop < $exon_end ) )
            {    # 3' exons
                $exon_start_phase = -1;
                $exon_end_phase   = -1;
            }
            elsif (( $translation_start == $exon_start )
                && ( $translation_stop == $exon_end ) )
            {    # single exon genes
                $phase_diff       = 0;
                $exon_start_phase = 0;
                $exon_end_phase   = 0;
            }
            else {    # internal exon
                $phase_diff       = ( $exon_end - $exon_start + 1 ) % 3;
                $exon_start_phase = $last_exon_phase;
                $exon_end_phase   = ( $last_exon_phase + $phase_diff ) % 3;
            }

            # set exon phase
            $eExon->phase($exon_start_phase);
            $eExon->end_phase($exon_end_phase);

            #$span = $exon_end - $exon_start + 1;
            $last_exon_phase = $exon_end_phase;
        }

        else {    # -ve strand
                  # 5' exons
            if (   ( $translation_start < $exon_start )
                && ( $translation_start < $exon_end ) )
            {
                $exon_start_phase = -1;
                $exon_end_phase   = -1;
            }

            # exon stops with a start
            elsif ( ( $translation_start == ( $exon_end - 2 ) ) ) {
                $phase_diff       = ( ( $exon_end - $exon_start + 1 ) % 3 );
                $exon_start_phase = 0;
                $exon_end_phase   = $phase_diff;
            }

            # exon contains a start codon
            elsif (( $translation_start > $exon_start )
                && ( $translation_start < $exon_end ) )
            {
                $phase_diff = ( ( $exon_end - $translation_start + 1 ) % 3 );
                $exon_start_phase = -1;
                $exon_end_phase   = $phase_diff;
            }

            # exon contains a stop codon
            elsif (( $translation_stop > $exon_start )
                && ( $translation_stop < $exon_end ) )
            {
                $phase_diff = ( ( $exon_end - $translation_stop + 1 ) % 3 );
                $exon_start_phase = $last_exon_phase;
                $exon_end_phase   = -1;
            }

            # exon stops with a stop
            elsif ( ( $translation_stop == $exon_start ) ) {
                $phase_diff       = ( ( $exon_end - $exon_start + 1 ) % 3 );
                $exon_start_phase = $last_exon_phase;
                $exon_end_phase   = 0;
            }

            # 3' exons
            elsif (( $translation_stop > $exon_start )
                && ( $translation_stop > $exon_end ) )
            {
                $exon_start_phase = -1;
                $exon_end_phase   = -1;
            }

            # single exon genes
            elsif (( $translation_start == $exon_start )
                && ( $translation_stop == $exon_end ) )
            {
                $phase_diff       = 0;
                $exon_start_phase = 0;
                $exon_end_phase   = 0;
            }

            # internal exon
            else {
                $phase_diff       = ( $exon_end - $exon_start + 1 ) % 3;
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
        unless ($start_exon) {
            if ( $exon == 0 ) {
                $start_exon = $eExon;
            }
        }

        # Final exon
        unless ($end_exon) {
            if ( $exon == $genedata->{EXON_COUNT} - 1 ) {
                $end_exon = $eExon;
            }
        }

        # Translation start exon
        unless ($trans_start_exon) {
            if ( $exon == $genedata->{CDS_START_EXON} ) {
                $trans_start_exon = $eExon;
            }
        }

        # Translation stop exon
        unless ($trans_end_exon) {
            if ( $exon == $genedata->{CDS_END_EXON} ) {
                $trans_end_exon = $eExon;
            }
        }
    }    #_ END of exon loop _#

    $eTranscript->start_Exon($start_exon);
    $eTranscript->end_Exon($end_exon);

    ###############
    # TRANSLATION #
    ###############

    # Check for translation start being different to Exon start (i.e. UTR)
    if ( $genedata->{STRAND} > 0 ) {
        $translation_offset_left =
          $translation_start -
          $genedata->{EXON_START}[ $genedata->{CDS_START_EXON} ] + 1;
        $translation_offset_right =
          $translation_stop -
          $genedata->{EXON_START}[ $genedata->{CDS_END_EXON} ] + 1;

        #- $genedata->{EXON_START}[$genedata->{CDS_END_EXON}] + 3;
    }
    elsif ( $genedata->{STRAND} < 0 ) {
        $translation_offset_left =
          $genedata->{EXON_END}[ $genedata->{CDS_START_EXON} ] -
          $translation_start + 1;
        $translation_offset_right =
          $genedata->{EXON_END}[ $genedata->{CDS_END_EXON} ] -
          $translation_stop + 1;    #was + 3;
    }

#  print "// [$gene] TranscriptLength ", $eTranscript->length, " [$translation_offset_left - $translation_offset_right]\n" if 1;#( $verbose);

    #$eTranslation_ID = $genedata->{TRANSCRIPT_NAME};
    #$eTranslation_ID =~ s/-R/-P/;

    my $eTranslation = new Bio::EnsEMBL::Translation(
        -START_EXON => $trans_start_exon,
        -END_EXON   => $trans_end_exon,
        -SEQ_START  => $translation_offset_left,
        -SEQ_END    => $translation_offset_right,
        -STABLE_ID  => $gene,
        -VERSION    => '1',
    );

    #########################################
    # EnsEMBL add translation to transcript #
    #########################################

    $eTranscript->translation($eTranslation);

    ##################################
    # EnsEMBL add transcript to gene #
    ##################################

    $eGene->add_Transcript($eTranscript);

    my $dbid = '???';
    if ($I) {
        $dbid = $ENS_DBA->get_GeneAdaptor->store($eGene);

=stub

    # Desperate attempt to fix sorghum;
    my( $trn ) = @{$eGene->get_all_Transcripts};
    my $phase_nudge = 1;
    my $pepseq='*XX';
    while( $pepseq =~ /\*.*/ and $phase_nudge <= 3 ){
      $ENS_DBA->dbc->do
          ( sprintf( "
UPDATE translation SET seq_start=%s WHERE translation_id=%s", 
                     $phase_nudge, $trn->translation->dbID ) );
      $phase_nudge ++;
      $trn->translation->{seq} = undef; # Clear cache
      $pepseq = $trn->translation->seq;
    }
    if( $trn->translation->seq =~ /\*.*/ ){
      die sprintf( "[WARN] translation %s (%s) does not translate", 
                   $trn->translation->display_id,
                   $trn->translation->dbID );
    }

=cut

    }
    print "[INFO] dbID = $dbid for gene = $gene ($n of $number)\n" if $V;
    $n++;
}
warn("[INFO] Stored $number genes to DB\n");
warn("[INFO] DONE!\n");

#warn Dumper( $GENES );

exit();

1;
