#!/usr/local/bin/perl 

=head1 NAME

dump_ensembl_proteins.pl - dump protein fasta from ensembl database based on user entered options, for example longest 2 proteins.
	

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
#use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } 
#        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
#             ensembl-draw/modules ensembl-compara/modules );
use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;


use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Pod::Usage;


=head1 SYNOPSIS

dump_ensembl_proteins.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --logicname	        the analysis logic_name of the gene set
    --longestn          only dump the longest n transcript for each gene

=head1 OPTIONS

=over 4


=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped

=item B<--logicname>

 the analysis logic_name of the gene set, for example: araport11 for arabidopsis_thaliana

=item B<--longestn>

dump only the longest n sequences

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back


=cut

my ($species, $registry, $longestn, $logicname);
{  #Argument Processing
  my $help=0;
  my $man=0;
  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"logicname=s"=>\$logicname
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"longestn=i"=>\$longestn
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 

my @genes=@{$gene_adaptor->fetch_all_by_logic_name($logicname)};

my $seqio = new Bio::SeqIO(-format => 'fasta',
    		        -file => ">${species}.prot$longestn.fasta",
			);
my %count;

warn ("Got genes, ", scalar @genes);

foreach my $gene (@genes) {
  
  $count{total_genes}++;
  my $sid = $gene->stable_id;
  my $type = $gene->biotype;
  unless( $type =~ /^protein\s+coding$/i || $type =~ /^protein_coding$/i || $type =~ /^protein\-coding$/i   ){
	warn ("gene $sid does not seem to be protein coding gene, it is $type\n");	
  	next;
  }

  $count{qualified_genes}++;
  
  my $translations = &get_longestn_translation( $gene, $longestn );
  
  foreach my $trans (@{$translations}) {
  

      my $trpt_stable_id = $trans->stable_id;
      my $trpt_obj = $trans->transcript;
      my $trpt_seq_region = join '', ($trpt_obj->seq_region_name, ':', $trpt_obj->start, '-', $trpt_obj->end, ':', $trpt_obj->strand); 
      my $stable_id = join '|', ($trpt_stable_id, $sid, $trpt_seq_region);
	
      my $seq = $trans->seq; 
      my $seq_obj;

      $seq_obj = Bio::Seq->new(
	      -display_id => $stable_id,
	      -seq => $seq,
	      );

      $seqio->write_seq($seq_obj);
      $count{qualified_transcripts}++;
    
  }
  #last;
}

for my $k (sort keys %count){
    print "$k = $count{$k}\n";
}


sub get_longestn_translation{
  my $g = shift;
  my $n = shift;

  my @translations_srt =

	sort {

		$b->length <=> $a->length
	} 
	map{
		$_->translation
	}
	grep{
		$_->translation
        } 
	@{$g->get_all_Transcripts};

   my $size = scalar @translations_srt;  
   if( $n && $n < $size ){
   
	return [ splice (@translations_srt, 0, $n)];
   }else{
	return \@translations_srt;
   }
	
}

  
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

