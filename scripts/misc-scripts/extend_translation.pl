#!/usr/local/bin/perl 

=head1 NAME

extend_translations.pl

If translation starts at beginning of transcript, search upstream for the furthest possible start codon (no in frame stops).
If such a start codon is found, extend the exon to that point and add an additional --prepend nucleotides

If translation ends at end of transcript, search downstream for the furthest possible stop codon.
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
use POSIX;
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
    --extendM extend even if protein starts with M
    --extend max number of bases to extend coding region
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

=item B<--extend> 

   max number of bases to extend coding region (default 2001)

=item B<--extendM> 

   if set, will try to extend upstream even if current translation starts with M

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
my ($debug, %exclude_gene, $bylogicname, $nowrite, $extendM, $extend, $prepend, $append, $guideFile);
my $margin=undef;
{  							#Argument Processing
  my $help=0;
  my $man=0;
  my @exclude_gene=();

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
        ,"bylogicname=s"=>\$bylogicname
        ,"exclude=s"=>\@exclude_gene
	      ,"debug"=>\$debug
        ,"extendM"=>\$extendM
        ,"extend=i"=>\$extend
        ,"prepend=i"=>\$prepend
        ,"append=i"=>\$append
	      ,"nowrite"=>\$nowrite
        ,"guide=s"=>\$guideFile
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  							#pod2usage(2) if $margin<0;
  $prepend ||= 0;
  $append ||= 0;
  $extend ||= 2001;
  $extend % 3 == 0 or die "extend needs to be a multiple of 3\n";
  %exclude_gene = map {$_,1} map {split /,/} @exclude_gene;
}

my %guide;
if (-e $guideFile) {
  open (my $fh, "<", $guideFile);
  while (<$fh>) {
    chomp;
    my ($transStableId,$distance) = split /\t/, $_;
    $guide{$transStableId} = $distance;
  }
  close $fh;
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

# need to be able to check whether fixing the first/last exon could break other transcripts
print STDERR "reading exon_transcript\n" if $debug;
my $sql = "select exon_id,transcript_id,rank from exon_transcript order by transcript_id, rank DESC";
my $sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
my %exonCheck;
my %lastExonLUT;
while (my $row = $sth->fetchrow_arrayref) {
  my ($e,$t,$rank) = @$row;
  if (exists $lastExonLUT{$t}) {
    $exonCheck{isNotLast}{$e} = $t;
  }
  else {
    $lastExonLUT{$t}=$e;
  }
  if ($rank > 1) {
    $exonCheck{isNotFirst}{$e} = $t;
  }
}
$sth->finish;

my $update_translation_sql = qq{update translation set 
                                seq_start=?,
                                seq_end=?
                                where translation_id=?
                              };
my $update_gene_sql = qq{update gene set seq_region_start=?, seq_region_end=? where gene_id=?};
my $update_transcript_sql = qq{update transcript set seq_region_start=?, seq_region_end=? where transcript_id=?};
my $update_exon_sql = qq{update exon set seq_region_start=?, seq_region_end=? where exon_id=?};

my $update_translation_sth;
my $update_gene_sth;
my $update_transcript_sth;
my $update_exon_sth;
 
unless ($nowrite){
  $update_gene_sth = $dbh->prepare($update_gene_sql) or die "cannot prepare $update_gene_sql\n";
  $update_transcript_sth = $dbh->prepare($update_transcript_sql) or die "cannot prepare $update_transcript_sql\n";
  $update_translation_sth = $dbh->prepare($update_translation_sql) or die "cannot prepare $update_translation_sql\n";
  $update_exon_sth = $dbh->prepare($update_exon_sql) or die "cannot prepare $update_exon_sql\n";
}

#print "DBbase connected is ", $ENS_DBA->dbname, "\n" if $debug;
my @genes = map { $gene_adaptor->fetch_by_stable_id($_) } @ARGV;
@genes or
  @genes = !$bylogicname ? @{$gene_adaptor->fetch_all()}:
  @{$gene_adaptor->fetch_all_by_logic_name($bylogicname)};
  
@genes or die "what, no genes?\n";

print STDERR "will process ",scalar @genes - scalar keys %exclude_gene," genes\n" if $debug;
my %updates = (
  genes => 0,
  transcripts => 0,
  translations => 0,
  exons => 0
);
foreach my $gene (@genes) {
  next unless $gene->biotype eq 'protein_coding';
  my $geneID = $gene->dbID;
  next if %exclude_gene and $exclude_gene{$gene->stable_id};
  print STDERR "starting gene ",$gene->stable_id,"\n" if $debug;
  my $geneStart = $gene->start;
  my $geneEnd = $gene->end;
  my $updateGene=0;
  my @transcripts = @{$gene->get_all_Transcripts()};
  my (%adjustedStart, %adjustedEnd);
  for my $trans (@transcripts) {
    my $transID = $trans->dbID;
    my $transStableId = $trans->stable_id;
    next unless $trans->biotype eq 'protein_coding';
    my $translation = $trans->translation;
    $translation or die "huh? I expected to get a translation of transcript " . $trans->stable_id;
    my $translationID = $translation->dbID;
    my $firstAA = substr($translation->seq,0,1);
    my @exons = @{$trans->get_all_Exons()};
    my $startExonID = $trans->start_Exon->dbID;
    my $updateFirstExon=0;
    my $updateLastExon=0;
    if (not exists $adjustedStart{$startExonID}) {
      if ($exonCheck{isNotFirst}{$startExonID}) {
        print STDERR "WARNING - Do not adjust this start exon (not always first)\n";
      }
      elsif ($extendM or $firstAA ne 'M') {
        print STDERR "checking upstream region for translation ",$translation->stable_id," because firstAA is $firstAA\n" if $debug;
        my $newCDSStartPos = seekUpstream($trans,$guide{$translation->stable_id});
        if ($newCDSStartPos) {
          print STDERR "extending $transStableId CDS upstream $newCDSStartPos\n" if $debug;
          $adjustedStart{$startExonID} = $newCDSStartPos;
          # update the first exon start pos
          my $exonStart = $exons[0]->start;
          my $exonEnd = $exons[0]->end;
          if ($trans->strand == 1) {
            $exonStart = $exons[0]->start($exonStart - $newCDSStartPos - $prepend);
          }
          else {
            $exonEnd = $exons[0]->end($exonEnd + $newCDSStartPos + $prepend);
          }
          $updateFirstExon=1;
        }
      }
    }
    my $endExonID = $trans->end_Exon->dbID;
    if (not exists $adjustedEnd{$endExonID}) {
      if ($exonCheck{isNotLast}{$endExonID}) {
        print STDERR "WARNING - Do not adjust this end exon (not always last)\n";
      }
      else {
        my $newCDSEndPos = seekDownstream($trans);
        if ($newCDSEndPos) {
          print STDERR "extending $transStableId CDS downstream $newCDSEndPos\n" if $debug;
          $adjustedEnd{$endExonID} = $newCDSEndPos;
          # update the last exon end pos
          my $exonStart = $exons[-1]->start;
          my $exonEnd = $exons[-1]->end;
          my $utr_length = $trans->length - $trans->cdna_coding_end;
          if ($trans->strand == 1) {
            $exonEnd = $exons[-1]->end($exonEnd + $newCDSEndPos + $append - $utr_length);
          }
          else {
            $exonStart = $exons[-1]->start($exonStart - $newCDSEndPos - $append + $utr_length);
          }
          $updateLastExon=1;
        }
      }
    }
    if ($updateFirstExon) {
      $update_exon_sth->execute($exons[0]->start,$exons[0]->end,$startExonID) if $update_exon_sth;
      print STDERR "update exon $startExonID ",$exons[0]->start," ",$exons[0]->end,"\n" if $debug;
      $updates{exons}++;
    }
    if ($updateLastExon and (not $updateFirstExon or @exons > 1)) {
      $update_exon_sth->execute($exons[-1]->start,$exons[-1]->end,$endExonID) if $update_exon_sth;
      print STDERR "update exon $endExonID ",$exons[-1]->start," ",$exons[-1]->end,"\n" if $debug;
      $updates{exons}++;
    }
  }
  for my $trans (@transcripts) {
    my $transID = $trans->dbID;
    my $transcriptStart = $trans->start;
    my $transcriptEnd = $trans->end;
    next unless $trans->biotype eq 'protein_coding';
    my $translation = $trans->translation;
    $translation or die "huh? I expected to get a translation of transcript " . $trans->stable_id;
    my $translationID = $translation->dbID;
    my $translationStart = $translation->start;
    my $translationEnd = $translation->end;
    my $startExonID = $translation->start_Exon->dbID;
    my $endExonID = $translation->end_Exon->dbID;
    my ($updateTranscript, $updateTranslation);
    if ($adjustedStart{$startExonID}) {
      if ($translationStart > 3) { # some other translation was extended, but this one wasn't
        $translationStart += $adjustedStart{$startExonID};
        $updateTranslation = 1;
      }
      if ($prepend > 0) {
        $translationStart += $prepend;
        $updateTranslation = 1;
      }
      if ($updateTranslation == 1 and $startExonID == $endExonID) {
        $translationEnd += $prepend + $adjustedStart{$startExonID};
      }
      # definitely update start/end of transcript depending on strand
      $updateTranscript=1;
      if ($trans->strand == 1) {
        $transcriptStart = $transcriptStart - $adjustedStart{$startExonID} - $prepend;
      }
      else {
        $transcriptEnd = $transcriptEnd + $adjustedStart{$startExonID} + $prepend;
      }
    }
    if ($adjustedEnd{$endExonID}) {
      $translationEnd += $adjustedEnd{$endExonID};
      $updateTranslation = 1;
    
      # definitely update start/end of transcript depending on strand
      $updateTranscript=1;
      if ($trans->strand == 1) {
        $transcriptEnd = $trans->coding_region_end + $adjustedEnd{$endExonID} + $append;
      }
      else {
        $transcriptStart = $trans->coding_region_start - $adjustedEnd{$endExonID} - $append;
      }
    }
    if ($updateTranscript) {
      $update_transcript_sth->execute($transcriptStart,$transcriptEnd,$transID) if $update_transcript_sth;
      print STDERR "update transcript $transID $transcriptStart $transcriptEnd\n" if $debug;
      $updates{transcripts}++;
      if ($transcriptStart < $geneStart) {
        $geneStart = $transcriptStart;
        $updateGene=1;
      }
      if ($transcriptEnd > $geneEnd) {
        $geneEnd = $transcriptEnd;
        $updateGene=1;
      }
    }
    if ($updateTranslation) {
      $update_translation_sth->execute($translationStart,$translationEnd,$translationID) if $update_translation_sth;
      print STDERR "update translation $translationID $translationStart $translationEnd\n" if $debug;
      $updates{translations}++;
    }
  }
  if ($updateGene) {
    $update_gene_sth->execute($geneStart,$geneEnd,$geneID) if $update_gene_sth;
    print STDERR "update gene $geneID $geneStart $geneEnd\n" if $debug;
    $updates{genes}++;
  }
}

$update_translation_sth->finish if $update_translation_sth;
$update_gene_sth->finish if $update_gene_sth;
$update_transcript_sth->finish if $update_transcript_sth;
$update_exon_sth->finish if $update_exon_sth;
$dbh->disconnect;

print STDERR "updated $updates{genes} genes $updates{transcripts} transcripts $updates{translations} translations $updates{exons} exons\n";
  
########################## subroutines ######################################
sub seekUpstream {
  my $trans = shift;
  my $estimated_distance = shift;
  $trans->cdna_coding_start <= 3 or return 0;
  my $len = $estimated_distance ? 3*ceil(1.5 * $estimated_distance) : $extend;
  print STDERR "len to check is $len\n" if $debug;
  my $seqToCheck;
  if ($trans->strand == 1) {
    if ($trans->start < $len) {
      $len = $trans->start - ($trans->start % 3) - 1;
      print STDERR "WARNING: + transcript starts near start of seq_region. ",$trans->start," new len to scan is $len\n";
    }
    my $crs = $trans->coding_region_start;
    $seqToCheck = $trans->slice->subseq($crs-$len, $crs-1, 1);
  }
  else {
    my $cre = $trans->coding_region_end;
    if ($trans->end + $len > $trans->slice->end) {
      $len = $trans->slice->end - $cre;
      $len -= $len % 3;
      print STDERR "WARNING: - transcript starts near end of seq_region. ",$trans->slice->end - $trans->end," new len to scan is $len\n";
    }
    $seqToCheck = $trans->slice->subseq($cre+1, $cre+$len, -1);
  }
  return findEarliestStart($seqToCheck);
}

sub seekDownstream {
  my $trans = shift;
  my $utr_length = $trans->length - $trans->cdna_coding_end;
  $utr_length < 3 or return 0;
  my $len = $extend;
  my $seqToCheck;
  if ($trans->strand == 1) {
    my $cre = $trans->coding_region_end;
    if ($cre + $len > $trans->slice->end) {
      $len = $trans->slice->end - $cre;
      $len -= $len % 3;
      print STDERR "WARNING: + transcript ends near end of seq_region. ",$trans->slice->end - $cre," new len to scan is $len\n";
    }
    $seqToCheck = $trans->slice->subseq($cre-2, $cre+$len, 1);
  }
  else {
    my $crs = $trans->coding_region_start;
    if ($crs < $len) {
      $len = $crs - ($crs % 3) - 1;
      print STDERR "WARNING: - transcript ends near start of seq_region. ",$crs," new len to scan is $len\n";
    }
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

   Andrew Olson <olson@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

