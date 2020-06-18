#!/usr/local/bin/perl 

=head1 NAME

fix_internal_start_codon.pl - Some translation start before Methonine, Find all the protein coding genes, check if they all start from M, if not, find the 1st M and set start codon there.
	

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

fix_internal_start_codon.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --exclude		genes to exclude
    --bylogicname	only work on genes udner this analysis logic_name
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

    supply the species name whose transcripts are to be dumped

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=item B<--guide> 

   Tab delimited file with translation_stable_id and position of internal M to use

=back

=head1 ARGUMENTS

    Gene Stable ids - only dump transcripts of these genes
    None=All.

=cut

my ($species, $registry);
my (%exclude_gene, $bylogicname, $debug, $nowrite, $guide_file, $max_exon_idx);
my $margin=undef;
{  							#Argument Processing
  my $help=0;
  my $man=0;
  my @exclude_gene=();

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"exclude=s"=>\@exclude_gene
	      ,"bylogicname=s"=>\$bylogicname
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
        ,"guide=s"=>\$guide_file
	      ,"nowrite"=>\$nowrite
        ,"maxexonidx=i"=>\$max_exon_idx
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  							#pod2usage(2) if $margin<0;
  %exclude_gene= map { $_,1 }  map { split /,/ } @exclude_gene;
  $max_exon_idx ||= 99;
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

# print "DBbase connected is ", $ENS_DBA->dbname, "\n" if $debug;

my @genes = map{ $gene_adaptor->fetch_by_stable_id($_) } @ARGV;

				#my @genes = ($gene_adaptor->fetch_by_stable_id('Opunc01g00010'));
@genes or  
    @genes = !$bylogicname ? @{$gene_adaptor->fetch_all()}:
    @{$gene_adaptor->fetch_all_by_logic_name($bylogicname)};

my %count;

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

print STDERR "processing ". scalar @genes . " genes\n";

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
my @exon_fields = qw(seq_region_id seq_region_start seq_region_end seq_region_strand phase end_phase is_current is_constitutive stable_id version);
my $exon_fields_string = join(',',@exon_fields);
my $question_marks = join(',', map {'?'} @exon_fields);
my $select_exon_sth = $dbh->prepare("select $exon_fields_string from exon where exon_id=?");
my $insert_exon_sth = $dbh->prepare("insert into exon ($exon_fields_string) VALUES ($question_marks)");
my $update_exon_phase_sth = $dbh->prepare("update exon set phase=? where exon_id=?");
my $update_exon_transcript_sth = $dbh->prepare("update exon_transcript set exon_id=?, transcript_id=? where exon_id=?");
my $update_translation_sth = $dbh->prepare("update translation set start_exon_id=?, seq_start=? where translation_id=?");
use Data::Dumper;
#warn ( Dumper(\@genes) );

foreach my $gene(@genes) {
  print STDERR "geneid = ", $gene->stable_id, "\n" if $debug;
  
  $count{total_genes}++;
  
  next if %exclude_gene and $exclude_gene{$gene->stable_id};
  
  $count{qualified_genes}++;
  
  my @transcripts;
  @transcripts = @{$gene->get_all_Transcripts};
  my %updates;
  foreach my $trans (@transcripts) {
    my $cdna_seq=$trans->spliced_seq;
    my $id = $trans->dbID;
    my $stableid = $trans->stable_id;
    my $slice_name = $trans->slice->name;
    my $biotype = $trans->biotype;
    my $logic_name = $trans->analysis->logic_name;
    my $strand = $trans->strand;
    my $comp_id = join "|", ($id, $stableid, $strand, $slice_name, $logic_name, $biotype);
    			#print "processing transcript $comp_id\n";
    unless ( $cdna_seq ){
      print STDERR "No cDNA seq for :$comp_id\n";
      next;
    }

    $count{qualified_transcripts}++;

    my $trmapper = $trans->get_TranscriptMapper;
    			#my $translation = $trans->translation;
    my $aa = $trans->translate->seq;
 
    # print STDERR ">$comp_id\n$aa\n" if $debug;
    my $translation = $trans->translation;
    my $translation_id= $translation->dbID;
    my $translation_stable_id= $translation->stable_id;
    my $translation_old_start= $translation->start;
    
    my $idxm = $guide{$translation_stable_id} || index( uc($aa), 'M', 1);
    $idxm += 1;
    next if ($guide_file and not $guide{$translation_stable_id});
    
    if($aa =~ /^M/i and not $guide{$translation_stable_id}){
      $count{qualified_transcripts_with_M}++;
      next;
    } else {
      if( $aa =~ /M/i){
        # print STDERR "matched\n" if $debug;
        $count{qualified_transcripts_withInternal_M}++;

        # print STDERR "$comp_id: 1 based index of 1st M is $idxm\n" if $debug;

        my @genomic_coords = $trmapper->pep2genomic( $idxm, $idxm );
        my $Met_start_genomic;
        my $start_exon_start_phase;

        # map{ print $_->start .", " . $_->end . "\n"} @genomic_coords if $debug;

        if( scalar @genomic_coords == 0 ){
          warn("No genomic coord found for M at $idxm, skip\n");
          next;
        }

        if ($strand > 0){
          my @starts = sort map {$_->start} @genomic_coords;
          $Met_start_genomic = shift @starts;
        } else {
          my @ends = reverse sort map {$_->end} @genomic_coords;
          $Met_start_genomic = shift @ends;
		    }

        # print "Met start genomic coord is $Met_start_genomic (the genomic start fo 1st M may looks off by one codon for minus strand gene, but believe it it will end up giving the correct translation in the end)\n" if $debug;
	
		    my @fiveUTRexonIDs2update;
		    my ($met_start_ExonID, $start_exon_start, $exon_start_phase);
		    my @ordered_Exons = $strand>0 ? @{$trans->get_all_Exons}:
											sort {$b->seq_region_start <=> $a->seq_region_start} @{$trans->get_all_Exons};
		    my $Exon;
        my $exon_idx=0;
		    while ( $Exon = shift @ordered_Exons){
          $exon_idx++;
          my $exon_gstart = $Exon->seq_region_start;
          my $exon_gend = $Exon->seq_region_end;
          $exon_start_phase = $Exon->phase;
          my $exon_seq = $Exon->seq->seq;

          # print "$Met_start_genomic ? [$exon_gstart, $exon_gend,  $exon_start_phase]\n$exon_seq\n" if $debug;
          if( $Met_start_genomic >= $exon_gstart && $Met_start_genomic <= $exon_gend) {
				    $met_start_ExonID = $Exon->dbID;
				    $start_exon_start = $strand>0 ? $Met_start_genomic-$exon_gstart+1:$exon_gend-$Met_start_genomic+1;
				    if($strand > 0){
              $start_exon_start_phase = $Met_start_genomic == $exon_gstart ? 0 : -1;
				    }else{
				      $start_exon_start_phase = $Met_start_genomic == $exon_gend ? 0 : -1;
				    }
				    last;
				  }
		
	        push @fiveUTRexonIDs2update, $Exon->dbID if $strand>0;
	    	}
	    	
	   		@fiveUTRexonIDs2update = map{ $_->dbID } @ordered_Exons if $strand < 0;
              # print STDERR "met_start_ExonID=$met_start_ExonID, start_exon_start_phase=$start_exon_start_phase\n" if $debug;
		    unless( $met_start_ExonID && defined $start_exon_start_phase){
				  warn("ERROR: no valid exon found and start found for genomic coord   $Met_start_genomic, skip $comp_id\n");
				  next;
	    	}

	    	unless($nowrite or $exon_idx > $max_exon_idx){
          $updates{$met_start_ExonID}{$trans->dbID} = {
            startExon => $met_start_ExonID, # this could change if a new exon needs to be created
            phase => $start_exon_start_phase,
            seqStart => $start_exon_start,
            translationID => $translation_id
          };
        # print "$translation_stable_id ($translation_id), old start $translation_old_start, startExonID $met_start_ExonID, Met start in startExon $start_exon_start, startPhase ($exon_start_phase -> $start_exon_start_phase) exon_idx $exon_idx\n";
	    		$count{qualified_transcripts_withInternal_M_fixed}++;
	    	}
      } else {
	    	$count{qualified_transcripts_without_M}++;
	    	next;
		  }
    }
   

  }					#trans
 				 	#last;
          # warn(Dumper(\%updates)) if $debug;
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
        $insert_exon_sth->execute($sr,$st,$en,$str,$phase,$ep,$cur,$con,$sid,$vers);
        my $newExonId = $dbh->last_insert_id(undef, undef, undef, undef);
        print STDERR "got new exon id $newExonId\n" if $debug;
        # update exon_transcript table and startExon of translation (update happens after this if block)
        for my $tid (keys %{$phases{$phase}}) {
          print STDERR "updating exon_transcript $newExonId,$tid,$eid\n" if $debug;
          $update_exon_transcript_sth->execute($newExonId,$tid,$eid);
          $updates{$eid}{$tid}{startExon} = $newExonId;
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
      print STDERR "updating translation ",$u->{startExon},$u->{seqStart},$u->{translationID},"\n" if $debug;
      $update_translation_sth->execute($u->{startExon}, $u->{seqStart}, $u->{translationID});
    }
  }
}					#Gene


for my $k (sort keys %count){
  print STDERR "$k = $count{$k}\n";
}

$select_exon_sth->finish;
$insert_exon_sth->finish;
$update_exon_phase_sth->finish;
$update_exon_transcript_sth->finish;
$update_translation_sth->finish;

$dbh->disconnect;
  
########################## subroutines ######################################

__END__


=head1 OUTPUT


=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

