#!/usr/local/bin/perl 

=head1 NAME

extend_translations.pl

If translation starts at beginning of transcript, search upstream (up to 2001 bp) for the furthest possible start codon (no in frame stops).
If such a start codon is found, extend the exon to that point and add an additional --prepend nucleotides

If translation ends at end of transcript, search downstream (up to 2001 bp) for the furthest possible stop codon.
If found, extend last exon to that point and add an additional --append nucleotides.

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use Scalar::Util qw (reftype);

#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Upstream;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

=head1 SYNOPSIS

extend_translation.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species to dump
    --exclude		genes to exclude
    --bylogicname	only work on genes udner this analysis logic_name
    --guide tab delimited file of translation stable ids and distance to extend
    --prepend UTR length to prepend to extended transcript
    --append UTR length to append to extended transcript
    --debug             debug mode

=head1 OPTIONS

=over 4

=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name whose transcripts are to be dumped

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=item B<--guide> 

   Tab delimited file with translation_stable_id and position of upstream M to use (if possible)

=item B<--prepend> 

   number of bases to prepend upstream of the new start codon (default 0)

=item B<--append> 

   number of bases to append downstream of the new stop codon (default 0)

=back

=head1 ARGUMENTS

    Gene Stable ids - only dump transcripts of these genes
    None=All.

=cut

my ($species, $registry);
my ($debug, $nowrite, $guide_file, $prepend, $append);
my $margin=undef;
{  							#Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
        ,"guide=s"=>\$guide_file,
        ,"prepend=i"=>\$prepend,
        ,"append=i"=>\$append,
	      ,"nowrite"=>\$nowrite
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  							#pod2usage(2) if $margin<0;
  $prepend ||= 0;
  $append ||= 0;
}

# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 
my $transcript_adaptor = $ENS_DBA->get_TranscriptAdaptor;
my $transl_adaptor = $ENS_DBA->get_TranslationAdaptor;
#my $assembly_type=$ENS_DBA->get_MetaContainer()->get_default_assembly();

my $dbh = $ENS_DBA->dbc->db_handle;
my $update_translation_sql = qq{update translation set 
                                start_exon_id=?, 
                                seq_start=?
                                where translation_id=?
                              };
			#We cannot change the exon phases, may exons were shared between transcripts
			#The other OK transcripts may rely on the phase of the exons. 
my $update_exon_start_sql = qq{update exon set phase=? where exon_id=?};
my $update_exon_sql = qq{update exon set phase=-1, end_phase=-1 where exon_id=?};

my $update_translation_sth;
my $update_exon_start_sth; 
#my $update_exon_sth;
 
unless ($nowrite){
    $update_translation_sth = $dbh->prepare($update_translation_sql) ||
							die "cannot prepare $update_translation_sql\n";
	$update_exon_start_sth = $dbh->prepare($update_exon_start_sql) ||
									die "cannot prepare $update_exon_start_sql\n";
#	$update_exon_sth = $dbh->prepare($update_exon_sql) ||
#									die "cannot prepare $update_exon_sql\n";
}

#print "DBbase connected is ", $ENS_DBA->dbname, "\n" if $debug;

my %guide;
if ($guide_file and -e $guide_file) {
  open (my $guide_fh, "<", $guide_file);
  while (<$guide_fh>) {
    chomp;
    my ($stable_id, $idmx) = split /\t/, $_;
    $guide{$stable_id} = $idmx;
  }
  close $guide_fh;
}
my @translations = map { $transl_adaptor->fetch_by_stable_id($_)} keys %guide;

print "found ". scalar @translations . " translations\n";

foreach my $translation (@translations) {
  my $trans = $translation->transcript;

  my $newCDSStartPos = seekUpstream($trans);
  my $newCDSEndPos = seekDownstream($trans);

  if ($newCDSStartPos or $newCDSEndPos) {
    my $ccs = $trans->cdna_coding_start();
    my $cce = $trans->cdna_coding_end();
    print "BEFORE:",$translation->stable_id," cdna_coding_start:$ccs cdna_coding_end:$cce translation:",$trans->translate->seq,"\n";

    my @exons = @{$trans->get_all_Exons()};
    if ($newCDSStartPos) {
      $exons[0] = $exons[0]->adjust_start_end(0 - $newCDSStartPos - $prepend, 0);
      $trans->{_trans_exon_array} = \@exons;
      $translation->start_Exon($exons[0]);
      $ccs = $trans->cdna_coding_start($ccs + $prepend);
      $cce = $trans->cdna_coding_end($cce + $newCDSStartPos + $prepend);
    }
    if ($newCDSEndPos) {
      $exons[-1] = $exons[-1]->adjust_start_end(0, $newCDSEndPos + $append);
      $trans->{_trans_exon_array} = \@exons;
      $translation->end_Exon($exons[-1]);
      $cce = $trans->cdna_coding_end($cce + $newCDSEndPos);
    }
    print "AFTER:",$translation->stable_id," cdna_coding_start:$ccs cdna_coding_end:$cce translation:",$trans->translate->seq,"\n";
  }
}

$update_translation_sth->finish if $update_translation_sth;
$update_exon_start_sth->finish if $update_exon_start_sth;
#$update_exon_sth->finish if $update_exon_sth;
$dbh->disconnect;
  
########################## subroutines ######################################
sub seekUpstream {
  my $trans = shift;
  $trans->cdna_coding_start == 1 or return 0;
  my $len = 2001;
  my $seqToCheck;
  if ($trans->strand == 1) {
    my $crs = $trans->coding_region_start;
    $seqToCheck = $trans->slice->subseq($crs-$len, $crs-1, 1);
  }
  else {
    my $cre = $trans->coding_region_end;
    $seqToCheck = $trans->slice->subseq($cre+1, $cre+$len, -1);
  }
  return findEarliestStart($seqToCheck);
}

sub seekDownstream {
  my $trans = shift;
  $trans->three_prime_utr and return 0;
  print "NO 3'UTR in ",$trans->stable_id,"\n";
  my $len = 21;
  my $seqToCheck;
  if ($trans->strand == 1) {
    my $cre = $trans->coding_region_end;
    $seqToCheck = $trans->slice->subseq($cre-2, $cre+$len, 1);
  }
  else {
    my $crs = $trans->coding_region_start;
    $seqToCheck = $trans->slice->subseq($crs-$len, $crs+2, -1);
  }
  return findFirstStop($seqToCheck);
}

sub findFirstStop {
  my $seq = shift;
  my %stop = (
    TAA => 1,
    TAG => 1,
    TGA => 1
  );
  my $len = length $seq;
  for(my $i=0;$i<$len;$i+=3) {
    my $codon = substr($seq,$i,3);
    print "checking $codon $i\n";
    return $i if $stop{$codon};
  }
  return 0;
}

sub findEarliestStart {
  my $seq = shift;
  my %codons = (
    TAA => 'stop',
    TAG => 'stop',
    TGA => 'stop',
    ATG => 'start'
  );
  my $len = length $seq;
  my $stop=-1;
  my @starts;
  for(my $i=0;$i<$len;$i+=3) {
    my $codon = substr($seq,$i,3);
    next unless exists $codons{$codon};
    if ($codons{$codon} eq 'stop') {
      $stop=$i;
    }
    elsif ($codons{$codon} eq 'start') {
      push @starts, $i;
    }
  }
  while (@starts and $starts[0] < $stop) {
    my $skip = shift @starts;
  }
  if (@starts) {
    return $len - $starts[0];
  }
  return 0;
}

__END__


=head1 OUTPUT


=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

