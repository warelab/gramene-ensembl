#!/usr/local/bin/perl -w

=head1 NAME

get_ensembl_data.pl -- this program provides methods to extract data from 
    an ensembl database for the maize ftp site.

=head1 SYNOPSIS

perl get_ensembl_data.pl [options]

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -o --output

=head1 OPTIONS

gets data from an ensembl database as specified

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

B<-s --output>                                                       
  a location to place output for the switches specified.

=head1 DESCRIPTION

B<This program> 

Extracts information from an ensembl database as specified by the user

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

use vars qw($I $ENS_DBA);
my $help=0;
my $man=0;
my($species, $file, $outdir);
GetOptions
    (
     "help|?"           => \$help,
     "man"              => \$man,
     "species=s"        => \$species,
     "registry_file=s"  => \$file,
     "output=s"         => \$outdir,
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

####MAIN####

# Load the ensembl registry file
$species || ( warn( "Need a --species\n" ) && pod2usage(1) );
Bio::EnsEMBL::Registry->load_all( $file );
$ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( warn( "No core DB for $species set in $file\n" ) &&
	      pod2usage(1) );
my $slice_adaptor = $ENS_DBA->get_adaptor('Slice');

my $clones = latest_clones ($slice_adaptor);
bacs ($slice_adaptor, $clones);
contigs ($slice_adaptor, $clones); # contains the genes and proteins methods

####SUBS####
sub latest_clones {
    my $slice_adaptor = shift;

    # get all the 'current' clones                                                                          
    my @slices = grep {
        scalar @{ $_->get_all_Attributes('current-version') } > 0
    } @{ $slice_adaptor->fetch_all('clone') };

    my %clones = ();
    for my $slice (@slices){
        my $seq_region    = $slice->seq_region_name();
	    $clones{$seq_region} = 1;	
    }
    return (\%clones);
}

sub bacs {
    my $slice_adaptor = shift;
    my $clones = shift;

    # get all the clones
    my @slices = @{ $slice_adaptor->fetch_all('clone') };
    
    open (BACS, ">$outdir/BACS.fasta");
    foreach my $slice (@slices){
        my $seq_region = $slice->seq_region_name();
	next unless (exists ($clones->{$seq_region}));

#	next if $seq_region eq "AC198200.3";

        my $strand     = $slice->strand();

        warn "Working on $seq_region\n";
	
        # force the forward direction                                                
        if ($strand < 1){
            $slice = $slice->invert;
        }

        # get stuff                                                                          
        my $location   = $slice->name();
        my $type       = $slice->coord_system->name();
        my $sequence   = $slice->seq();

        print BACS ">$seq_region\n$sequence\n";
    }
    close (BACS);
}

sub contigs {
    my $slice_adaptor = shift;
    my $clones = shift;

    # get all contig slices
    my @slices = @{ $slice_adaptor->fetch_all('contig') };
    open (CONTIGS, ">$outdir/BAC_contigs.fasta");
 
    # cycle through the slices and dump the files of interest
    foreach my $slice (@slices){
	my $seq_region = $slice->seq_region_name();
	my ($accession, $contig) = split (/\-/, $seq_region);

	next unless (exists ($clones->{$accession}));

#	next if ($accession eq "AC198200.3");

	my $strand     = $slice->strand();
	
	warn "Working on $seq_region\n";
	
	# force the forward direction
	if ($strand < 1){ 
	    $slice = $slice->invert;
	}
	
	# get stuff
	my $location   = $slice->name();
	my $type       = $slice->coord_system->name();
	my $sequence   = $slice->seq();
	
	# output the contigs if desired
	print CONTIGS ">$seq_region\n$sequence\n";
	
	# output the genes and translations
	open (ALL_GENES, ">>$outdir/ALL_GENES.fasta");
	open (ALL_PTS, ">>$outdir/ALL_TRANSLATIONS.fasta");
	open (TE_GENES, ">>$outdir/TE-LIKE_GENES.fasta");
	open (G_GENES, ">>$outdir/NON-TE-LIKE_GENES.fasta");
	open (TE_PTS, ">>$outdir/TE-LIKE_TRANSLATIONS.fasta");
	open (G_PTS, ">>$outdir/NON-TE-LIKE_TRANSLATIONS.fasta");
	my @genes = @{ $slice->get_all_Genes };

	foreach my $gene (@genes){
	    my $gene_id    = $gene->dbID();
	    my $stable_id  = $gene->stable_id();
	    my $seq_region = $gene->slice->seq_region_name();
	    my $start      = $gene->start();
	    my $end        = $gene->end();
	    my $strand     = $gene->strand();
	    my $biotype    = $gene->biotype();
	    my $sequence   = $gene->seq();

	    warn "$stable_id\n";
	    
	    print ALL_GENES ">$stable_id:$seq_region:$start-$end:$strand:$biotype\n$sequence\n";

	    if ($biotype eq "transposon_pseudogene"){
		print TE_GENES ">$stable_id:$seq_region:$start-$end:$strand:$biotype\n$sequence\n";
	    }
	    else {
		print G_GENES ">$stable_id:$seq_region:$start-$end:$strand:$biotype\n$sequence\n";
	    }

	    # get the translations
	    foreach my $trans (@{ $gene->get_all_Transcripts }) {
		next if (!$trans->translation);
		
		my $identifier = $trans->stable_id;
		
		# check to see if there is a stop codon
		my $tseq = $trans->translate();
		my $aa   = $tseq->seq;
		if ($aa =~ /\*/) {
		    print STDERR "Translation of $identifier has stop codons ",
		    "- Skipping! (in ", $trans->slice->name(), ")\n";
		    next;
		}
		
		else {
		    # compose
		    
		    print ALL_PTS ">$identifier\n$aa\n";
		    
		    if ($biotype eq "transposon_pseudogene"){
			print TE_PTS ">$identifier\n$aa\n";
		    }
		    else {
			print G_PTS ">$identifier\n$aa\n";
		    }
		}
	    }
	}
	close (TE_GENES);
	close (G_GENES);    
	close (TE_PTS);
	close (G_PTS);
	close (ALL_GENES);
	close (ALL_PTS);
    }
    close (CONTIGS);
}

