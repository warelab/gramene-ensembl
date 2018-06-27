#! /usr/bin/env perl

use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Bio::EnsEMBL::Registry;

#our $ENS_DBA;
#our $CMP_DBA;

my( $sp1,  $reg,);
GetOptions ( 
    "species=s"          => \$sp1,
    "ensembl_registry=s" => \$reg,
    #"other=s"            => \$sp2,
    );

# Load the ensembl file and adapters
Bio::EnsEMBL::Registry->load_all( $reg );

#adaptor for homologies
my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('compara', 'compara', 'GeneMember');
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('compara', 'compara', 'Homology');

# Get all genes from refernce species
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $sp1, 'core', 'Gene' );
my $gene_id_arrayref = $gene_adaptor->list_stable_ids;

#my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $sp1, 'core', 'Transcript' );

my %orthologs;
my %target_species;
#my %ref_t2g_id;

#my $counter = 0;
while( my $gid = shift @$gene_id_arrayref ){ #a gene stable id
    print "$gid\n";
    next unless $gid;
    $orthologs{$gid} = undef;

    #my $gobj = $gene_adaptor->fetch_by_stable_id($gid);
    #for ( @{$gobj->get_all_Transcripts}){
#	$ref_t2g_id{$_->stable_id} = $gid;
#    }
    my $gm1 = $gene_member_adaptor->fetch_by_stable_id($gid); #gene_member
    next unless $gm1;
    my $homologies = $homology_adaptor->fetch_all_by_Member($gm1); #this should filter for target but doesn't work
    next unless $homologies;

    for my $homology (@{$homologies}) {
	my $hom_type = $homology->description;
	next unless $hom_type =~ /ortholog/;
	my $tsp;
	my $high_confidence = $homology->is_high_confidence;
	my %one_ortholog_align;
	$one_ortholog_align{high_confidence} = $high_confidence || 1;
	my %sp2gid;
	for my $gm ( @{$homology->get_all_Members()} ){ #gene_member  @{$homology->gene_list()}
	    my $sp = $gm->genome_db()->name;
	    my $geneID = $gm->gene_member->stable_id; #$gm->stable_id;     
	    print STDERR "gene member id = $geneID\n";
	    if($sp eq $sp1){
		die "ERROR: gene name not match for reference $gid ne $geneID " unless ($geneID eq $gid);
		$one_ortholog_align{ref} = $gm->perc_id;
	    }else{
		$tsp = $sp;
		$sp =~ s/([a-z])[a-z]*_/uc($1)/e;
		$target_species{$tsp} = $sp;
		push @{$one_ortholog_align{target}}, [$geneID, $gm->perc_id];
	    }
	
	}
	$orthologs{$gid}{$tsp} = \%one_ortholog_align;
    }

    #last if $counter++ > 100;
}

print_orthologs(\%orthologs, \%target_species);


#need to delete this func
#sub get_gid_from_tid{
#	my $tid = shift;
#	print STDERR "find gene stable id for $tid\n";	
#	my $tobj = $transcript_adaptor->fetch_by_stable_id( $tid );
#	return $tobj->get_Gene()->stable_id;
#}

sub print_orthologs {

	my $ortholg_hashrf = shift;
	my $target_species_hashrf = shift;

	my @header = qw/refGid targetGid refPecID targetPecID HighConf/;
    
	my $headerline = join "\t", @header;

	for my $aspecies (keys %{$target_species_hashrf}){
		
		my $filename = $target_species_hashrf->{$aspecies} . '.rtm';
		print STDERR "process species $aspecies and write to $filename\n";

		open my $fh, '>', "./$filename" or die "Cannot open ./$filename for write";
	  	print $fh "$headerline\n";

		my $ort_line = '';
		for my $ageneID(keys %{$ortholg_hashrf} ){
			if( my $ort_hashrf = $ortholg_hashrf->{$ageneID}{$aspecies} ){
				for my $ort_target( @{$ort_hashrf->{target}} ){
					$ort_line .= join "\t", (
							$ageneID,
							$ort_target->[0],
							$ort_hashrf->{ref},
							$ort_target->[1],
							$ort_hashrf->{high_confidence}."\n"		
							);
				}								
			}else{
				$ort_line = "$ageneID\n";
			}
			
			print $fh $ort_line;
		}
		close $fh;
	 	
	}

}
