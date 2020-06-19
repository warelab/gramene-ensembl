#!/usr/local/bin/perl 

=head1 NAME

recalculate_translation.pl - Some translations contain stop codons or don't start with M.
If there is a problem, try to recalculate a valid translation.
If that's not possible, mark transcript as non-coding

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

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

=head1 SYNOPSIS

recalculate_translation.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species to process
    --exclude		genes to exclude
    --bylogicname	only work on genes udner this analysis logic_name
    --complete          require complete CDS
    --minORF            minimum length for an ORF
    --nowrite           test only , no change to the database
    --debug             debug mode

=head1 OPTIONS

=over 4


=item B<--exclude>

    Genes to ignore.  

=item B<--bylogicname>

    Only work on Genes under this logicname

=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name whose translations are to be recalculated

=item B<--nowrite>

    test run without writing to database

=item B<--complete>

    require complete CDS

=item B<--minorf>

    minimum length for an ORF

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=back

=head1 ARGUMENTS

    Gene Stable ids - only recalculate translations for these genes
    None=All.

=cut

my ($species, $registry);
my (%exclude_gene, $bylogicname, $debug, $nowrite, $complete, $minORF);
{  #Argument Processing
  my $help=0;
  my $man=0;
  my @exclude_gene=();

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"exclude=s"=>\@exclude_gene
	      ,"bylogicname=s"=>\$bylogicname
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
	      ,"nowrite"=>\$nowrite
        ,"complete"=>\$complete
        ,"minorf=i"=>\$minORF
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  %exclude_gene= map { $_,1 }  map { split /,/ } @exclude_gene;
  $minORF ||= 0;
}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;

print STDERR "#INFO DBbase connected is ", $ENS_DBA->dbc->dbname, "\n" if $debug;

my @genes = map{ $gene_adaptor->fetch_by_stable_id($_) } @ARGV;

#my @genes = ($gene_adaptor->fetch_by_stable_id('Opunc01g00010'));
@genes or  
    @genes = $bylogicname
    ? @{$gene_adaptor->fetch_all_by_logic_name()}
    : @{$gene_adaptor->fetch_all()};
print STDERR "#INFO will process ",scalar @genes," genes\n";

my $dbh = $ENS_DBA->dbc->db_handle;

my %exon_transcript_phase; # {exon}{transcript} = phase
my $select_exon_transcript_sth = $dbh->prepare(qq{
  select et.exon_id, et.transcript_id, e.phase from exon_transcript et, exon e
  where et.exon_id = e.exon_id
});
$select_exon_transcript_sth->execute();
while (my $row = $select_exon_transcript_sth->fetchrow_arrayref) {
  my ($e,$t,$phase) = @$row;
  $exon_transcript_phase{$e}{$t} = $phase;
}
$select_exon_transcript_sth->finish;
print STDERR "read exon transcript phase\n" if $debug;
my @exon_fields = qw(exon_id seq_region_id seq_region_start seq_region_end seq_region_strand phase end_phase is_current is_constitutive stable_id version);
my $exon_fields_string = join(',',@exon_fields);
my $question_marks = join(',', map {'?'} @exon_fields);
my $select_exon_sth = $dbh->prepare("select $exon_fields_string from exon where exon_id=?");
my $insert_exon_sth = $dbh->prepare("insert into exon ($exon_fields_string) VALUES ($question_marks)");
my $update_exon_phase_sth = $dbh->prepare("update exon set phase=? where exon_id=?");
my $update_exon_transcript_sth = $dbh->prepare("update exon_transcript set exon_id=? where exon_id=? and transcript_id=?");
my $update_translation_start_sth = $dbh->prepare("update translation set start_exon_id=?, seq_start=? where translation_id=?");
my $update_translation_end_sth = $dbh->prepare("update translation set end_exon_id=?, seq_end=? where translation_id=?");

my %sql = (
  delete_translation => "delete from translation where transcript_id=?",
  update_transcript => "update transcript set biotype='non_coding' where transcript_id=?",
  update_gene => "update gene set biotype='non_coding' where gene_id=?"
);

my %sth; 
unless ($nowrite){
  for my $statement (keys %sql) {
    $sth{$statement} = $dbh->prepare($sql{$statement}) or die "cannot prepare $sql{$statement}\n";
  }
}

my %count;
my %transcript_with_internal_stop;

foreach my $gene(@genes) {
  print STDERR "#INFO geneid = ", $gene->stable_id, " biotype = ", $gene->biotype, "\n" if $debug;
  $count{total_genes}++;
  next unless ($gene->biotype eq 'protein_coding');
  
  next if %exclude_gene and $exclude_gene{$gene->stable_id};
  
  $count{qualified_genes}++;
  
  my @transcripts;
  @transcripts = @{$gene->get_all_Transcripts};
  my %updates;
  my $good_translations=0;
  for my $transcript (@transcripts) {
    my $biotype = $transcript->biotype;
    next unless $biotype eq 'protein_coding';
    my $id = $transcript->dbID;
    my $stableid = $transcript->stable_id;
    my $translation = $transcript->translation();

    my $aa_seq;
    eval{$aa_seq= $translation->seq()};
    $@ and die "#ERROR getting translation->seq() for $stableid, $@\n";

    my $bad_translation=0;
    if ($aa_seq) {
      $count{qualified_transcripts}++;
      if ($aa_seq !~ /^M/i and $complete) {
        warn("#WARN aa does not start with M\n");
        $bad_translation=1;
      }
      if ($aa_seq =~ /\*[A-Z]/i) {
        warn("#WARN aa contains internal stop\n");
        $bad_translation=1;
      }
    } else {
      warn("#WARN no aa translation for $stableid\n");
      $bad_translation=1;
    }
    
    if (not $bad_translation) {
      $good_translations++;
      next;
    }

    # recalculate translation
    my $mrna = $transcript->spliced_seq();
    
    my @orfs = find_orfs($mrna);
    
    my ($orfStart, $orfEnd, $orfLen, $phase) = @{$orfs[0]};
    if ($orfLen < $minORF) {
      $bad_translation=1;
      # remove translation from db
      print STDERR "delete from translation where transcript_id=",$transcript->dbID,";\n" if $debug;
      $sth{delete_translation}->execute($transcript->dbID) unless $nowrite;
      
      # change transcript from protein-coding to non-coding
      print STDERR "update transcript set biotype='non_coding' where transcript_id=",$transcript->dbID,";\n" if $debug;
      $sth{update_transcript}->execute($transcript->dbID) unless $nowrite;
      next;
    }
    $good_translations = 1;
    update_translation($transcript, $orfStart, $orfEnd-3, \%updates);
    
  }
  if (not $good_translations) {
    # change gene->biotype to non-coding
    print STDERR "update gene set biotype='non_coding' where gene_id=",$gene->dbID,";\n" if $debug;
    $sth{update_gene}->execute($gene->dbID) unless $nowrite;
  }
  next if $nowrite;
  for my $eid (keys %updates) {
    # check if updates are consistent w.r.t. exon phase
    # if not, create a new exon for the updated translation
    my %phases;
    my $orig_phase;
    for my $tid (keys %{$exon_transcript_phase{$eid}}) {
      $orig_phase = $exon_transcript_phase{$eid}{$tid};
      if (exists $updates{$eid}{$tid}) {
        $phases{$updates{$eid}{$tid}{phase}}{$tid} = 1;
      }
      else {
        $phases{$orig_phase}{$tid} = 1;
      }
    }
    if (keys %phases > 1) {
      # create a new exon for each new phase
      $select_exon_sth->execute($eid);
      my ($eid2,$sr,$st,$en,$str,$ph,$ep,$cur,$con,$sid,$vers,@dates) = @{$select_exon_sth->fetchrow_arrayref};
      for my $phase (keys %phases) {
        $phase = $phase + 0; # hash keys are strings, but phase is an integer
        next if $phase == $orig_phase;
        print STDERR "creating new exon with phase $phase\n" if $debug;
        $vers++;
        $insert_exon_sth->execute(undef,$sr,$st,$en,$str,$phase,$ep,$cur,$con,$sid,$vers);
        my $newExonId = $dbh->last_insert_id(undef, undef, undef, undef);
        print STDERR "got new exon id $newExonId\n" if $debug;
        # update exon_transcript table and startExon of translation (update happens after this if block)
        for my $tid (keys %{$phases{$phase}}) {
          print STDERR "updating exon_transcript $newExonId,$eid,$tid\n" if $debug;
          $update_exon_transcript_sth->execute($newExonId,$eid,$tid);
          $updates{$eid}{$tid}{exonID} = $newExonId;
        }
      }
    }
    else {
      my ($phase) = keys %phases;
      $phase = $phase + 0;
      if ($phase != $orig_phase) {
        print STDERR "updating exon phase $eid $phase\n" if $debug;
        $update_exon_phase_sth->execute($phase, $eid);
      }
    }
    # update the translations
    for my $tid (keys %{$updates{$eid}}) {
      my $u = $updates{$eid}{$tid};
      if ($u->{seq_start}) {
        print STDERR "updating translation start\n" if $debug;
        $update_translation_start_sth->execute($u->{exonID}, $u->{seq_start}, $u->{translationID});
      }
      if ($u->{seq_end}) {
        print STDERR "updating translation end\n" if $debug;
        $update_translation_end_sth->execute($u->{exonID}, $u->{seq_end}, $u->{translationID});
      }
    }
  }
}

unless ($nowrite) {
  for my $handle (values %sth) {
    $handle->finish
  }
}
$select_exon_sth->finish;
$insert_exon_sth->finish;
$update_exon_phase_sth->finish;
$update_exon_transcript_sth->finish;
$update_translation_start_sth->finish;
$update_translation_end_sth->finish;

$dbh->disconnect;

  
########################## subroutines ######################################
sub find_orfs {
    my ( $sequence ) = @_;
    $sequence    = uc $sequence;
    my $codon_table = Bio::Tools::CodonTable->new( -id => 1);
    my $is_start = sub { shift eq 'ATG' };
 
    # stores the begin index of the currently-running ORF in each
    # reading frame
    my @current_orf_start = (-1,-1,-1);
 
    #< stores coordinates of longest observed orf (so far) in each
    #  reading frame
    my @orfs;
 
    # go through each base of the sequence, and each reading frame for each base
    my $seqlen = length $sequence;
    my @start_frame_order;
    for( my $j = 0; $j <= $seqlen-3; $j++ ) {
        my $frame = $j % 3;
 
        my $this_codon = substr( $sequence, $j, 3 );
 
        # if in an orf and this is either a stop codon or the last in-frame codon in the string
        if ( $current_orf_start[$frame] >= 0 ) {
            if ( $codon_table->is_ter_codon( $this_codon ) ||( my $is_last_codon_in_frame = ($j >= $seqlen-5)) ) {
                # record ORF start, end (half-open), length, and frame
                my @this_orf = ( $current_orf_start[$frame], $j+3, undef, $frame );
                my $this_orf_length = $this_orf[2] = ( $this_orf[1] - $this_orf[0] );
                push @orfs, \@this_orf;
                $current_orf_start[$frame] = -1;
            }
        }
        # if this is a start codon
        elsif ( $is_start->($this_codon) ) {
            $current_orf_start[$frame] = $j;
            push @start_frame_order, $frame;
        }
    }
 
    return sort { $b->[2] <=> $a->[2] } @orfs;
}

sub update_translation {
  my ($transcript, $orfStart, $orfEnd, $updatesRef) = @_;
  my %updates = %$updatesRef;
  my $translation = $transcript->translation();
  my @exons = @{$transcript->get_all_Exons};
  my ($seq_start, $seq_end, $start_exon_id, $end_exon_id);
  my $posInTranscript = 0;
  for my $exon (@exons) {
    my ($phase, $end_phase) = (-1, -1);
    if ($orfStart >= $posInTranscript and $orfStart < $posInTranscript + $exon->length) {
      $seq_start = $orfStart - $posInTranscript + 1;
      $start_exon_id = $exon->dbID;
    }
    if ($orfEnd > $posInTranscript and $orfEnd <= $posInTranscript + $exon->length) {
      $seq_end = $orfEnd - $posInTranscript;
      $end_exon_id = $exon->dbID;
    }
    if ($orfStart < $posInTranscript and $posInTranscript < $orfEnd) {
      $phase = ($posInTranscript - $orfStart) % 3;
    }
    $posInTranscript += $exon->length;
    if ($orfStart < $posInTranscript && $posInTranscript <= $orfEnd) {
      $end_phase = ($posInTranscript - $orfStart) % 3;
    }
    # update exon
    # print STDERR "update exon set phase=$phase, end_phase=$end_phase where exon_id=",$exon->dbID,";\n" if $debug;
    # $sth{update_exon}->execute($phase, $end_phase, $exon->dbID) unless $nowrite;

    $updates{$exon->dbID}{$transcript->dbID} = {
      exonID => $exon->dbID, # this could change if a new exon needs to be created
      phase => $phase,
      end_phase => $end_phase
    };
  }
  # update translation
  # print STDERR "update translation set start_exon_id=$start_exon_id, seq_start=$seq_start, end_exon_id=$end_exon_id, seq_end=$seq_end where translation_id=",$translation->dbID,";\n" if $debug;
  # $sth{update_translation}->execute($start_exon_id, $seq_start, $end_exon_id, $seq_end, $translation->dbID) unless $nowrite;
  my $startUpdate = $updates{$start_exon_id}{$transcript->dbID};
  $startUpdate->{translationID} = $translation->dbID;
  $startUpdate->{seq_start} = $seq_start;
  my $endUpdate = $updates{$end_exon_id}{$transcript->dbID};
  $endUpdate->{translationID} = $translation->dbID;
  $endUpdate->{seq_end} = $seq_end;
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

