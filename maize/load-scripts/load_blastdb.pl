#!/usr/local/bin/perl -w

=head1 NAME

load_blastdb.pl - creates/updates the blast databases for an ensembl database

=head1 SYNOPSIS

perl load_blastdb.pl [options]

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -o --output
 -f --formatter

=head1 OPTIONS

Creates blast databases from clones, gene models, and protein predictions

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

=head1 DESCRIPTION

B<This program> 

  Populates a core Ensembl DB with copy number features given 
  a user-input rlevel log threshold  

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper); # For debug 

use DBI;
use FindBin qw($Bin) ;
use File::Basename qw(dirname);

use vars qw($BASEDIR);
BEGIN{
  # Set the perl libraries
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/bioperl-live';
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::SimpleFeature;
use Bio::Seq;
use Bio::SeqIO;

use vars qw($ENS_DBA);
my $help=0;
my $man=0;
my($species, $file, $outdir, $wu_formatdb);
GetOptions
    (
     "help|?"          => \$help,
     "man"             => \$man,
     "species=s"       => \$species,
     "registry_file=s" => \$file,
     "output=s"        => \$outdir,
     "formatter=s"     => \$wu_formatdb,
     ) or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

# Validate file paths
$file    ||= $BASEDIR.'/conf/ensembl.registry';

map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
} $file;

# Load the ensembl file
$species || ( warn( "Need a --species\n" ) && pod2usage(1) );
Bio::EnsEMBL::Registry->load_all( $file );
$ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( warn( "No core DB for $species set in $file\n" ) &&
	      pod2usage(1) );

# Prepare some Ensembl adaptors and some static vars
my $slice_adaptor    = $ENS_DBA->get_adaptor('Slice');

# get all the contig slices
my @slices = @{ $slice_adaptor->fetch_all('contig') };

# cycle through the slices and dump the contig fasta file
my $seqio_contigs = Bio::SeqIO->new('-format' => 'Fasta', '-file' => ">$outdir/contigs.fasta");
my $seqio_predictions = Bio::SeqIO->new('-format' => 'Fasta', '-file' => ">$outdir/predictions.fasta");
foreach my $slice (@slices){
    my $seq_region = $slice->seq_region_name();
    my $strand     = $slice->strand();

    print STDERR "Working on $seq_region\n";
    
    # force the forward direction
    if ($strand < 1){ 
	$slice = $slice->invert;
    }

    # get stuff for desc
    my $location   = $slice->name();
    my $type       = $slice->coord_system->name();

    # get the sequence
    my $sequence   = $slice->seq();

    # compose
    my $seq = new Bio::Seq (-seq => $sequence, -display_id => $seq_region);
    $seq->description($location);
    $seqio_contigs->write_seq($seq);
    
    # get all the genes in this slice
    my @genes = @{ $slice->get_all_Genes };
    
    # cycle through the genes and get all of its transcripts
    foreach my $gene (@genes) {
        my $gene_id = $gene->dbID();
	
        foreach my $trans (@{ $gene->get_all_Transcripts }) {
            next if (!$trans->translation);
	    
            my $identifier = $trans->stable_id;
	    
	    # check to see if there is a stop codon
            my $tseq = $trans->translate();
            if ($tseq->seq =~ /\*/) {
                print STDERR "Translation of $identifier has stop codons ",
                "- Skipping! (in ", $trans->slice->name(), ")\n";
                next;
            }
	    
	    # compose
            $tseq->display_id($identifier);
            $tseq->desc("Translation id $identifier gene $gene_id");
	    $seqio_predictions->write_seq($tseq);
	}
    }
    
}

# system call: create the blast database on the contig fasta
`$wu_formatdb -i $outdir/contigs.fasta -p F -o T`;

# system call: create the blast database on the contig fasta                          
`$wu_formatdb -i $outdir/predictions.fasta -p T -o T`;
