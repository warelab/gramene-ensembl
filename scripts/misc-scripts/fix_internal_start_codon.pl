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

=back

=head1 ARGUMENTS

    Gene Stable ids - only dump transcripts of these genes
    None=All.

=cut

my ($species, $registry);
my (%exclude_gene, $bylogicname, $debug, $nowrite);
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
	      ,"nowrite"=>\$nowrite
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  							#pod2usage(2) if $margin<0;
  %exclude_gene= map { $_,1 }  map { split /,/ } @exclude_gene;
 
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

print "@ARGV\n";
my @genes = map{ $gene_adaptor->fetch_by_stable_id($_) } @ARGV;

				#my @genes = ($gene_adaptor->fetch_by_stable_id('Opunc01g00010'));
@genes or  
    @genes = !$bylogicname ? @{$gene_adaptor->fetch_all()}:
    @{$gene_adaptor->fetch_all_by_logic_name($bylogicname)};

my %count;

foreach my $gene(@genes) {
  #print "geneid = ", $gene->stable_id, "\n";
  
  $count{total_genes}++;
  
  next if %exclude_gene and $exclude_gene{$gene->stable_id};
  
  $count{qualified_genes}++;
  
  my @transcripts;
  @transcripts = @{$gene->get_all_Transcripts};

  foreach my $trans (@transcripts) {
  
    			#print join "\t", ($trans->stable_id, $trans->spliced_seq, "\n");
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

    # check the 1st aa, if it is M, incr count on intact_genes, skip;
    # if not M, find the 1st M, and reset the start codon to that M
    # if no M found, return error and skip;
 
    # To replace translation 
    # create new translation object
    # then do 
    # $transl_adaptor->remove($translation);
    # $transl_adaptor->store($translation_new);

    my $trmapper = $trans->get_TranscriptMapper;
    			#my $translation = $trans->translation;
    my $aa = $trans->translate->seq;
 
    			print ">$comp_id\n$aa\n" if $debug;
    			#exit;
    if($aa =~ /^M/i){
		$count{qualified_transcripts_with_M}++;
		next;
    }else{
		if( $aa =~ /M/i){
	    	print "matched\n" if $debug;
	    	$count{qualified_transcripts_withInternal_M}++;

	    	my $translation = $trans->translation;
	    	my $translation_id= $translation->dbID;
	    	my $translation_stable_id= $translation->stable_id;
	    	my $translation_old_start= $translation->start;

	    	my $idxm = index( uc($aa), 'M', 1);
	    	$idxm += 1;	
	    	#$idxm += $strand>0 ? 1 : 2;
	    	#my $idxm = $index_of_M+2;
	    	print "1 based index of 1st M is $idxm\n$aa\n" if $debug;
	    
	    	my @genomic_coords = $trmapper->pep2genomic( $idxm, $idxm );
	    	my $Met_start_genomic;
	    	my $start_exon_start_phase;
	    	
	    	map{ print $_->start .", " . $_->end . "\n"} @genomic_coords if $debug;

		    if( scalar @genomic_coords == 0 ){
				warn("No genomic coord found for M at $idxm, skip\n");
				next;
		    }
	    
		    if ($strand > 0){
				my @starts = sort map {$_->start} 
						#grep { reftype $_ eq 'Bio::EnsEMBL::Mapper::Coordinate' }
					@genomic_coords;
				$Met_start_genomic = shift @starts;	 
		    }else{
				my @ends = reverse sort map {$_->end}
						#grep { reftype $_ =~ /Bio::EnsEMBL::Mapper::Coordinate/ } 
				@genomic_coords;
						#print 'ends are ', @ends;
				$Met_start_genomic = shift @ends;
		    }

					    #$Met_start_genomic = $Met_start_genomic + $translation_old_start -1;
	    	print "Met start genomic coord is $Met_start_genomic (the genomic start fo 1st M may looks off by one codon for minus strand gene, but believe it it will end up giving the correct translation in the end)\n" if $debug;
	    				#exit;
	
		    my @fiveUTRexonIDs2update;
		    my ($met_start_ExonID, $start_exon_start, $exon_start_phase);
		    my @ordered_Exons = $strand>0 ? @{$trans->get_all_Exons}:
											sort {$b->seq_region_start <=> $a->seq_region_start} @{$trans->get_all_Exons};
		    my $Exon;
		    while ( $Exon = shift @ordered_Exons){
		    	
				my $exon_gstart = $Exon->seq_region_start;
				my $exon_gend = $Exon->seq_region_end;
				$exon_start_phase = $Exon->phase;
				my $exon_seq = $Exon->seq->seq;
       		
				print "$Met_start_genomic ? [$exon_gstart, $exon_gend,  $exon_start_phase]\n$exon_seq\n" if $debug;
				if( $Met_start_genomic >= $exon_gstart &&
				    $Met_start_genomic <= $exon_gend){		    

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
            		
		    unless( $met_start_ExonID && $start_exon_start_phase){
				warn("ERROR: no valid exon found and start found for genomic coord   $Met_start_genomic, skip $comp_id\n");
				next;
	    	}
	    	print "$update_exon_start_sql with [$start_exon_start_phase, $met_start_ExonID]\n" if $debug;	
	    	print "$update_translation_sql for $met_start_ExonID, $start_exon_start, $translation_id\n" if $debug;
	    
	    	unless($nowrite){
				#$update_exon_start_sth->execute($start_exon_start_phase, $met_start_ExonID);
#				map{ $update_exon_sth->execute($_) }@fiveUTRexonIDs2update;
				print "$translation_stable_id ($translation_id), old start $translation_old_start, startExonID $met_start_ExonID, Met start in startExon $start_exon_start, startPhase ($exon_start_phase -> $start_exon_start_phase)\n";
				$update_translation_sth->execute($met_start_ExonID, $start_exon_start, $translation_id) or die "cannot execute the sql for $met_start_ExonID, $start_exon_start, $translation_id";
	    	}
		
			#check the resulting translation, won't work since the update not officially complete.
			#my $newaa = $transcript_adaptor->fetch_by_dbID($id)->translate->seq;
			#my $newcds = $transcript_adaptor->fetch_by_dbID($id)->translateable_seq;
			#print "$stableid(old)=$aa\n$stableid(new)=$newaa\n\n";

		}else{
	    	$count{qualified_transcripts_without_M}++;
	    	next;
		}
    }
   

  }					#trans
 				 	#last;
}					#Gene


for my $k (sort keys %count){
  print "$k = $count{$k}\n";
}

$update_translation_sth->finish if $update_translation_sth;
$update_exon_start_sth->finish if $update_exon_start_sth;
#$update_exon_sth->finish if $update_exon_sth;
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

