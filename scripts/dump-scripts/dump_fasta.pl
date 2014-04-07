#!/usr/local/bin/perl 

=head1 NAME

dump_transcripts.pl - Make fasta files of transcripts, and 3' and 5'
    regions, for all Genes or for Genes given as arguments.
	

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


#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Gene;
use Getopt::Long;
use Pod::Usage;


=head1 SYNOPSIS

dump_fasta.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --exclude		genes to exclude
    --bylogicname	make a fasta file for each analysis logic_name
    --analysispgm	this(these) analysispgm(s) only
    --exclude-analysispgm	skip this analysispgm
    --exclude-clone	skip this clone (when doing all genes)
    --margin		margin for 3', 5'
    --type            coding | cDNA | protein 
    --longest           only dump the longest transcript for each gene

=head1 OPTIONS

=over 4


=item B<--exclude>

    Genes to ignore.  


=item B<--exclude-clone>
    
    Clones to ignore 
    This may be a comma-separated list.
    
=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped

=item B<--margin>

size in basepairs of 3' and 5' regions to dump 
default 500. 0=don't dump these.

=item B<--type>

types could be coding | cDNA | protein

=item B<--longest>

dump only the longest cDNA sequences

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

Gene Stable ids - only dump transcripts of these genes
None=All.

=cut

my ($species, $registry, @types, $longest);
my (%exclude_gene,%exclude_analysispgm,%analysispgm,%exclude_clone,$bylogicname);
my $margin=undef;
{  #Argument Processing
  my $help=0;
  my $man=0;
  my @exclude_gene=();
  my @exclude_analysispgm=();
  my @analysispgm=();
  my @exclude_clone=();
  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"exclude=s"=>\@exclude_gene
	      ,"bylogicname"=>\$bylogicname
	      ,"exclude-analysispgm=s"=>\@exclude_analysispgm
	      ,"analysispgm=s"=>\@analysispgm
	      ,"exclude-clone=s"=>\@exclude_clone
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"margin=i"=>\$margin
	      ,"type=s"=>\@types
	      ,"longest"=>\$longest
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  #pod2usage(2) if $margin<0;
  %exclude_gene= map { $_,1 }  map { split /,/ } @exclude_gene;
  %exclude_analysispgm= map { $_,1 } map { split /,/ } @exclude_analysispgm;
  %analysispgm= map { $_,1 }  map { split /,/ } @analysispgm;
  %exclude_clone= map { $_,1 } map { split /,/ } @exclude_clone;
}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 
#my $assembly_type=$ENS_DBA->get_MetaContainer()->get_default_assembly();

my @genes=@ARGV;
@genes or  @genes=@{$gene_adaptor->list_dbIDs()};

my %count;

my (%seqio,%seqio5,%seqio3);

unless($bylogicname) {

    for my $type (@types){
	$seqio{$type} = new Bio::SeqIO(-format => 'fasta',
				       -file => ">${species}.$type.fasta",
				       #-fh     => \*STDOUT
	    );
    }
}

if($margin) {

    for my $type (@types){
	$seqio5{$type} = Bio::SeqIO->new(-file =>">5prime.$type.fasta", '-format'=>'fasta');
	$seqio3{$type} = Bio::SeqIO->new(-file =>">3prime.$type.fasta", '-format'=>'fasta');
    }
}


foreach my $geneid (@genes) {
  print "! $geneid\n";
  
  $count{total_genes}++;
  my $gene;
  if(@ARGV) {    # Stable ids GRMGnnnnnnn
    eval { $gene=$gene_adaptor->fetch_by_stable_id($geneid) };
    print STDERR "$@\n" and next if $@;
  } else {	#internal ids
    eval { $gene= $gene_adaptor->fetch_by_dbID($geneid); };
    print STDERR "gene_id $geneid:\n$@\n" and next if $@;
    # fails e.g. if gene has no exons
  }
  next unless $gene;
  
  my ($working_set_attrib) = @{$gene->get_all_Attributes('working-set')};
  #print "gene attrib is ", $working_set_attrib->code, "\n";
  next if ($species =~ /zea|mays|maize/i && !$working_set_attrib);
  
  
  next if %exclude_gene and $exclude_gene{$gene->stable_id}
    or %analysispgm and !$analysispgm{$gene->analysis->program}
      or %exclude_analysispgm and $exclude_analysispgm{$gene->analysis->program}
	;
  
  if (%exclude_clone) {
    my $slice;
    eval {	#have some genes on phase1 clones -- no assembly, no slice
      $slice=$slice_adaptor->fetch_by_gene_stable_id($gene->stable_id,0);
    };
    print STDERR "$@\n" and next if $@;
    foreach my $tile (@{$slice->get_tiling_path}) {
      next GENE if $exclude_clone{$tile->component_Seq->clone->embl_id};
    }
    
  }
  
  $count{qualified_genes}++;
  
  my $ln;
  if ($bylogicname) {
    $ln=$gene->analysis->logic_name;
    for my $type (@types){
	$seqio{"$ln.$type"} ||= 
	    Bio::SeqIO->new(-file =>">$ln.$type.fasta", '-format'=>'fasta');
    }
  }
  
  my @transcripts;
  
  if($longest){
    my $lognest_trpt = &get_longest_transcript( $gene );
    push @transcripts, $lognest_trpt if $lognest_trpt;
  }else{
    @transcripts = @{$gene->get_all_Transcripts};
  }
  
  
  foreach my $trans (@transcripts) {
  
      my (%seq, %id, %comp_id);

      my $slice_name = $trans->slice->name;
      my $biotype = $trans->biotype;
      my $logic_name = $trans->analysis->logic_name;
      
      #print join "\t", ($trans->stable_id, $trans->spliced_seq, "\n");
      for my $type (@types){
	  $seq{$type}=$trans->spliced_seq if $type =~ /cdna/i;
	  $seq{$type}=$trans->translateable_seq if $type =~ /coding/i;
	  $seq{$type}=$trans->translation->seq if $type =~ /protein/i;
      }

      for my $type (@types){
	  $id{$type}=$trans->stable_id if $type =~ /cdna|coding/i;
	  $id{$type}=$trans->translation->stable_id if $type =~ /protein/i;
	  $comp_id{$type} =  join " ", ($id{$type}, join "|", ($slice_name, $logic_name, $biotype, $type));
	  if($type =~ /coding/i ){
	      my $cdna_coding_start = $trans->cdna_coding_start;
	      my $cdna_coding_end   = $trans->cdna_coding_end;
	      $comp_id{$type} .= 
		  "|coding region $cdna_coding_start-$cdna_coding_end";
	  }
      }


      for my $type (@types){
	  unless ( $seq{$type} ){
	      print STDERR "No $type seq for $geneid:$comp_id{$type}\n";
	      next;
	  }
      }
    
      for my $type (@types){
	  my $seq_obj;

	  $seq_obj = Bio::Seq->new(
	      -display_id => $comp_id{$type},
	      -seq => $seq{$type},
	      );

	  $seqio{$type}->write_seq($seq_obj);
      }

      $count{qualified_transcripts}++;
    
  }
  #last;
}

for my $k (sort keys %count){
    print "$k = $count{$k}\n";
}


sub get_longest_transcript{

  my $g = shift;

  my $max_len=0;
  my $max_len_trpt=undef;
  for my $trpt ( @{$g->get_all_Transcripts} ){
    my $length = $trpt->length;
    if ($length > $max_len){
      $max_len = $length;
      $max_len_trpt = $trpt;
    }
  }

  return $max_len_trpt;
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

