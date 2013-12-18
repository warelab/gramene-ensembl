#!/usr/local/bin/perl 

=head1 NAME

cleanup_overkap_fgenesh - when fgenesh genes overlap each other, keep only the longest one
	

=cut


#use

#fetch_all_by_logic_name('foobar');

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw { lib/perl };

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw { bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules };

use strict;
use warnings;
use List::MoreUtils qw{ distinct };


use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Pod::Usage;
use Readonly;

Readonly my $remove => 'remove';
Readonly my $keep => 'keep';

=head1 SYNOPSIS

dump_gene_coord.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --logic_name        logic name of the fgenesh analysis

=head1 OPTIONS

=over 4


=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped

=item B<--logic_name> 

 logic name of the fgenesh analysis

=item B<--d> 

 delete genes if set

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

Input files with Gene Stable ids 

=cut

my ($species, $registry, $logic_name, $delete);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,
	      "man"=>\$man,
	      "species=s"=>\$species,
	      "registry=s"=>\$registry,
	      "logic_name=s"=>\$logic_name,
	      "d"=>\$delete,
	    )
    or pod2usage(2);

  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

}

my %ofile = ($remove => "${species}.gene2$remove",
	     $keep   => "${species}.gene2$keep");

for my $k (keys %ofile){
    my $f = $ofile{$k};
    if( -e $f ){ 
	unlink $f or warn "Could not unlink $f: $!"; 
    }
}

# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $exon_adaptor=$ENS_DBA->get_ExonAdaptor;
my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;

my %coord_exon_list_hash;
foreach my $exon (@{$exon_adaptor->fetch_all_by_logic_name($logic_name)}) {
  
#print "! $geneid\n";

    my $seq_name = $exon->seqname;
    my $strand   = $exon->strand;
    my $start    = $exon->start;
    my $end      = $exon->end;
    #my $genelist = $exon->get_overlapping_Genes;

    #1. store in a hash, and cluster by overlapping exons
    #2. pick the longest gene from each cluster, identify the ones need to be removed
    #3. delete the genes identified

    #print "exon->$seq_name:$strand:$start-$end, exonID=", $exon->dbID, "\n";
    my $chr_strand = "$seq_name:$strand";
    push @{$coord_exon_list_hash{$chr_strand}}, [$start, $end, $exon];

 }

for my $k (keys %coord_exon_list_hash){

    print "region $k\n";
    #next unless $k =~ /barthii/i;
    my $exon_coord_listref = $coord_exon_list_hash{$k};
    my $exon_clusters = cluster_exons(@{$exon_coord_listref});
    my $genes_classified = find_overlap_genes2delete($exon_clusters);
    if ($delete){
	print "start to delete redundant genes\n";
	delete_genes($genes_classified->{$remove}) ;
    }
    #last;
}

sub cluster_exons{
    
    my @exon_list = @_;
    
    my @sorted_exon_list = sort { $a->[0] <=> $b->[0]; $a->[1] <=> $b->[1]} @exon_list;
    
    my %cluster;
    my $big_cluster=0;
    my $i=1;
   

    my $first_exon = shift @sorted_exon_list;
    my $pre_start = $first_exon->[0];
    my $pre_stop = $first_exon->[1];
    
    push @{$cluster{$i}}, $first_exon->[2]; 
    
    for my $exon( @sorted_exon_list ){
	
	if($exon->[0] > $pre_stop){
	
	    my $member_count = scalar @{$cluster{$i}};
	    $big_cluster++ if $member_count >1 ;
	    #print join ',', (map{ $_->dbID } @{$cluster{$i}});
	    #print "\n";
	    $i++;
	}

	$pre_start = $exon->[0];
	$pre_stop = $exon->[1];
	push @{$cluster{$i}}, $exon->[2];
    }

    print "There are $i exons clusters, $big_cluster cluster have more than one exon members\n";
    return \%cluster;

}

sub find_overlap_genes2delete{

    my $clusters =shift;
    
    my %gene2remove;
    my %gene2keep;
    my $gene_classify;

    for my $key (keys %$clusters){
	
	my $exon_list = $clusters->{$key};
	my @gene_list;
	foreach my $exon( @$exon_list ){
	    #print "$key => ", $exon->dbID, "\n";
	    my $overlapping_genes = $exon->get_overlapping_Genes;
	    push @gene_list, @$overlapping_genes;
	}

	my @uniq_gene_ids =  distinct (map{ $_->dbID } @gene_list);
	
	if(scalar @uniq_gene_ids > 1){
	    my @len_gene_id =  map {my $gene=$gene_adaptor->fetch_by_dbID($_);
				 my $cdnalen = length ($gene->canonical_transcript->spliced_seq);
				 [$cdnalen, $gene, $_]} @uniq_gene_ids;
	    my @sorted_gene_id_by_len = sort {$b->[0] <=> $a->[0]} @len_gene_id;

	    map{ $gene2keep{$_->[2]} = $_->[1]}  shift @sorted_gene_id_by_len;
	    map{ $gene2remove{$_->[2]} = $_->[1]} @sorted_gene_id_by_len;
	    
	}
	
    }

    print "number of gene to remove is ", scalar keys %gene2remove, 
          "\nnumber of genes to keep is ", scalar keys %gene2keep, "\n";

    
    $gene_classify->{$keep} = [values %gene2keep];
    $gene_classify->{$remove} = [values %gene2remove];

    for my $k (keys %$gene_classify){
	open my $fh,  ">>", $ofile{$k} or die "cannot open gene2$k to write";
	print_genes($gene_classify->{$k}, $fh);
	close $fh;
    }

    return $gene_classify;
}

sub print_genes{

    my $genelist=shift;
    my $fh=shift;

    for my $g (@$genelist){
	print $fh (join "\t", (
				$g->dbID,
				$g->stable_id,
				$g->seq_region_name,
				$g->seq_region_start,
				$g->seq_region_end,
				$g->seq_region_strand,
				));
	print $fh "\n";
    }
}



sub delete_genes{

    my $gene_list =shift;

    for my $g ( @$gene_list ){
	$gene_adaptor->remove($g);
    }
}


=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

