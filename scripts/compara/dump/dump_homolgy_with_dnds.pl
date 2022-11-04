#!/usr/bin/env perl

=pod

=head1 NAME

dump_complete_homology.pl - Generates matrix of cross-species homologs

=head1 SYNOPSIS

  perl dump_complete_homology.pl [options] 
	will output two files
	homology_report from homology table
	missing	the pairs not in homology table

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -t|--test		testing mode
  -q|--query_species	query species
  -s|--subject_species

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

=head1 DESCRIPTION

Script that dumps a pretty-print matrix of cross-species holologs
based on gene tree data. Takes a while to complete...

B<The Ensembl Registry>

  The database connection details for both Ensembl species and compara
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::Registry;
  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
  ( '-species' => 'compara',
    '-group'   => 'compara',
    '-dbname'  => 'ensembl_compara_48_28', );
  ---

TODO: Complete this section
                   
Maintained by Sharon Wei <weix@cshl.edu>

=cut


use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Algorithm::Combinatorics qw(combinations);

use vars qw( $BASEDIR );
BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname(dirname(dirname(dirname($Bin))));
  -d $BASEDIR.'/ensembl-live'
    || die( "\n[*DIE] Need $BASEDIR/ensembl-live symlinked to ensembl\n" );
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}

use Bio::EnsEMBL::Registry;

our $CMP_DBA;
our $CMP_DBH;
our $FH_HOMO;
our $FH_MISS;
our $TEST;
our ($qspecies, $sspecies);

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $reg, $slice_name );

  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "ensembl_registry=s" => \$reg,
	"test"	     	     => \$TEST,
	"query_species=s"    => \$qspecies,
	"subject_species=s"  => \$sspecies,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  map{
    -e $_ || pod2usage( "\n[*DIE]File $_ does not exist\n" );
    -r $_ || pod2usage( "\n[*DIE]Cannot read $_\n" );
    -f $_ || pod2usage( "\n[*DIE]File $_ is not plain-text\n" );
    -s $_ || pod2usage( "\n[*DIE]File $_ is empty\n" );
  } $reg;
  
  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );

  $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara','compara')
      || pod2usage("\n[*DIE] No compara DB set in $reg\n" );
  $CMP_DBH = $CMP_DBA->dbc;

  my $out_file_homo = "homology_report_${qspecies}_${sspecies}";
  open $FH_HOMO, '>', $out_file_homo or die "Cannot open $out_file_homo to write";

  my $out_file_miss = "missed_homolgy_report";
  open $FH_MISS, '>', $out_file_miss or die "Cannot open $out_file_miss to write";
}


our %sqls;
 
our %sths;

$sqls{gene_tree_members} = qq(select  gtr.root_id, s.seq_member_id
                from gene_tree_node grn  join gene_tree_root gtr using (root_id)
                join seq_member s using (seq_member_id)
                join genome_db g using (genome_db_id)
                where tree_type='tree' and
                        clusterset_id='default' 
		#and
                #g.name not in ('ciona_savignyi', 'cyanidioschyzon_merolae',
                #'drosophila_melanogaster', 'homo_sapiens', 'saccharomyces_cerevisiae')
        	#and gtr.root_id=141598
        ); #testing mode: and gtr.root_id=253467  237901

$sqls{homology_pairs} = qq(select hm1.seq_member_id, hm2.seq_member_id, h.homology_id,
                g1.name, g2.name, sm1.stable_id, d1.name, sm1.dnafrag_start, sm1.dnafrag_end, sm1.dnafrag_strand, 
		sm2.stable_id, d2.name, sm2.dnafrag_start, sm2.dnafrag_end, sm2.dnafrag_strand, h.description,
                h.wga_coverage, h.is_tree_compliant, h.is_high_confidence,
                hm1.perc_id, hm1.perc_cov, hm2.perc_id, hm2.perc_cov, h.dn, h.ds, h.n, h.s, h.lnl
                from homology_member hm1 join homology_member hm2
                on (hm1.homology_id=hm2.homology_id)
                join homology h on (hm1.homology_id=h.homology_id)
                join seq_member sm1 on( hm1.seq_member_id=sm1.seq_member_id) 
		join dnafrag d1 on (d1.dnafrag_id=sm1.dnafrag_id)  
                join seq_member sm2 on( hm2.seq_member_id=sm2.seq_member_id)
		join dnafrag d2 on (d2.dnafrag_id=sm2.dnafrag_id)
                join genome_db g1 on( g1.genome_db_id=sm1.genome_db_id)
                join genome_db g2 on( g2.genome_db_id=sm2.genome_db_id)
                where h.gene_tree_root_id=? and
                hm1.seq_member_id < hm2.seq_member_id 
		and 
                ((g1.name = ? and g2.name = ?) or (g2.name = ? and g1.name = ?))
		);  #253467 rid example

$sqls{gids} = qq(

                select g.genome_db_id, g.name, s.stable_id from
                genome_db g join seq_member s using(genome_db_id)
                where s.seq_member_id = ?
        );


our @header = qw(
species
gene_stable_id
gene_chr:start-end:strand
otherspecies
other_gene_stable_id
other_gene_chr:start-end:strand
homology_type
wga_coverage
is_tree_compliant
query_percent_identity
query_percent_coverage
other_percent_identify
other_percent_coverage
confidence
dn
ds
n
s
lnl
root_id
);

for my $k (keys %sqls) {
	
	
	$sths{$k} = $CMP_DBH->prepare( $sqls{$k} );


}

MAIN:{

  # Get a hash of gene_tree_root_id => seq_member_ids
  my %genetree_uniq_seqids = &get_genetree_members;
  	#warn Dumper( \%genetree_uniq_seqids );

  # Loop through each gene tree by root_id 
  foreach my $rid( keys %genetree_uniq_seqids ){
   
	warn("[INFO] Processing tree root_id $rid\n"); 
	my $seqid_pairs          = &pairwise_combination( $genetree_uniq_seqids{$rid} ); #genetree_uniq_seqid_pairs( $genetree_uniq_seqids{$rid}  );
	#$rid = 315780; #Only for testing
	my %homology_seqid_pairs = &get_homology_pairs_by_rid( $rid, $qspecies, $sspecies );
    	
	warn("[INFO] number of seqid pairs is ", scalar @{$seqid_pairs}, "\n");
	warn("[INFO] number of homology is ", scalar keys %homology_seqid_pairs, "\n");
	if ($TEST){
		map{ warn "inspect:$_\n"; warn "Not in homology:$_\n" if !$homology_seqid_pairs{$_} } @{$seqid_pairs};
		map{ warn "homolgy pair:$_\n" } keys %homology_seqid_pairs;
	}
	my %seqid_pairs;

	warn "[INFO] report homologies \n";
	print $FH_HOMO join "\t", map{uc} @header;
	print $FH_HOMO "\n"; 
	map{ print $FH_HOMO join "\t", @{$homology_seqid_pairs{$_}{report}}; print $FH_HOMO "\t$rid\n" } 
				keys %homology_seqid_pairs;

	#finish all the sths
	for my $k(keys %sths){
		$sths{$k}->finish;
	}
	#last;
  }

#this is the format we want to report in
=stub

query_species
query_gene_stable_id
query_gene_chr:query_gene_start-query_gene_end:query_gene_strand
otherspecies
other_gene_stable_id
other_gene_chr:other_gene_start-other_gene_end:other_gene_strand
homology_type
wga_coverage
is_tree_compliant
query_percent_identity
query_percent_coverage
other_percent_identify
other_percent_coverage
confidence
dn
ds
n
s
lnl

=cut



}


sub pairwise_combination{

	my $seqids = shift;
	my @seqid_array = @{$seqids};
	my @seqid_2ndarray = @seqid_array;
	my %pairs;

	for my $tseqid (@seqid_array){  
		warn("PROCESS $tseqid\n") if $TEST;
		warn ("FOUND outer $tseqid\n") if $TEST && ($tseqid == 67932 || $tseqid == 81584);
		for my $s( @seqid_2ndarray ){
			warn ("FOUND inner $s\n") if $TEST && ($s == 67932 || $s == 81584);
			next if ($tseqid == $s) ;#67932-81584
			my $k = $tseqid < $s ? "${tseqid}-$s" : "${s}-$tseqid";
			warn ("FOUND IT 67932-81584\n") if $TEST && ($k eq '67932-81584');
			$pairs{$k} = 1;
		}
	}

	return [keys %pairs];
}

sub get_genetree_members{


  my $rv = $sths{gene_tree_members}->execute || die( $sths{gene_tree_members}->errstr );
  my %genetree_seqmember;
  while( my $row = $sths{gene_tree_members}->fetchrow_arrayref ){
    #warn("DEBUG fetched ", $row->[0], ", ", $row->[1], "\n");
    push @{$genetree_seqmember{$row->[0]}}, $row->[1];
  }

  # get uinq seqmember list
  my %genetree_seqmember_sorted;
  for my $rid(keys  %genetree_seqmember){
	my %gmembers_uniq;
	map{ $gmembers_uniq{$_}=1 } @{$genetree_seqmember{$rid}};
	$genetree_seqmember_sorted{$rid} = [sort { $a <=> $b } keys %gmembers_uniq]; 
  }
  return %genetree_seqmember_sorted;
}


#======================================================================
# 
sub get_homology_pairs_by_rid{

  my $rid =shift or return;
  my $query =shift or return;
  my $subject =shift or return;

  my $rv = $sths{homology_pairs}->execute( $rid, $query, $subject, $query, $subject) || die( $sths{homology_pairs}->errstr );
  my %homology;
  while( my $row = $sths{homology_pairs}->fetchrow_arrayref ){
	my $seq_id_pair = join '-', ($row->[0], $row->[1]);
        $homology{$seq_id_pair}{homology_id} = $row->[2];

	#query_species, query_seq_stable_id, otherspecies, other_seq_stable_id, homology_type,
	#wga_coverage, is_tree_compliant, query_percent_identity, query_percent_coverage,
	#other_percent_identify,  other_percent_coverage,  confidence

	$homology{$seq_id_pair}{report} = [
		$row->[3] || 'NA', $row->[5] || 'NA', "$row->[6]:$row->[7]-$row->[8]:$row->[9]",
		$row->[4] || 'NA', 
		$row->[10] || 'NA', "$row->[11]:$row->[12]-$row->[13]:$row->[14]",
		$row->[15] || 'NA',  
		$row->[16] || 'NA', $row->[17] || 'NA', $row->[19] || 'NA', $row->[20] || 'NA', 
		$row->[21] || 'NA', $row->[22] || 'NA', $row->[18] || 'NA',
		$row->[23] || 'NA', $row->[24] || 'NA', $row->[25] || 'NA',
		$row->[26] || 'NA', $row->[27] || 'NA', 
		];
 
  }
  return %homology;
}

sub report_missing_pairs{

	my $miss_paris_seqids = shift;
	my $rid = shift;

	my $cnter = 0;	
	for my $seqid_pair( @{$miss_paris_seqids}){
		my $sid_meta = &get_gids($seqid_pair);

		my @gid_pair;
		my @sid_pair;
		for my $sid(keys %$sid_meta){
			push @sid_pair, $sid;
			push @gid_pair, $sid_meta->{$sid}->[0];
		}
	
		my $paf_data = &get_paf_for_seqid_pairs( $sid_meta );
		
		for my $seqkey (keys %{$paf_data}){
		
			print $FH_MISS (join "\t", map{ $_ = 'NA' unless $_ } @{$paf_data->{$seqkey}} );
			print $FH_MISS "\t$rid\n";
			$cnter++;
		}
		last if $cnter >= 1000;
		
	} 
}

sub get_gids{
	
	my $seqid_pair = shift;

	my $sql_gids = qq(

		select g.genome_db_id, g.name, s.stable_id from
		genome_db g join seq_member s using(genome_db_id)
		where s.seq_member_id = ?
	);


	my %sid2genome;
	for my $sid(split '-', $seqid_pair){
		my $rv = $sths{gids}->execute($sid) || die( $sths{gids}->errstr );	
		my @row = $sths{gids}->fetchrow_array;
		$sid2genome{ $sid } = [@row]; #gid, gname, seq_name
	}

	return \%sid2genome;
}

#----------------------------------------------------------------------
sub get_paf_for_seqid_pairs{

  	my $sid_meta = shift ;

  	my @gid_pair;
  	my @sid_pair;
        for my $sid(keys %$sid_meta){
       		push @sid_pair, $sid;
                push @gid_pair, $sid_meta->{$sid}->[0];
        }

	my %paf_data;
	for my $gid ( @gid_pair){
  		$sqls{"paf$gid"} ||= qq(
  			select g1.name, s1.stable_id, g2.name, s2.stable_id,
			paf.score, paf.evalue, paf.align_length,
			paf.identical_matches, paf.perc_ident, 
			paf.positive_matches, paf.perc_pos,
			paf.hit_rank 
			from peptide_align_feature_$gid paf 
			join seq_member s1 on (s1.seq_member_id=paf.qmember_id) 
			join seq_member s2 on (s2.seq_member_id=paf.hmember_id) 
			join genome_db g1 on (g1.genome_db_id=paf.qgenome_db_id) 
			join genome_db g2 on (g2.genome_db_id=paf.hgenome_db_id)  
			where paf.qmember_id = ? and paf.hmember_id = ? 
			or paf.qmember_id = ? and paf.hmember_id = ?
     			and paf.evalue < 0.01 and paf.score > 100
			#order by paf.evalue limit 1000
		);	
		$sths{"paf$gid"} ||= $CMP_DBH->prepare( $sqls{"paf$gid"} );
  		my $rv = $sths{"paf$gid"}->execute($sid_pair[0], $sid_pair[1], $sid_pair[1], $sid_pair[0]) || die( $sths{"paf$gid"}->errstr );

  		while( my @row = $sths{"paf$gid"}->fetchrow_array ){
    			$paf_data{$gid} = [@row];
  		}
	}	

	my %paf_merged;
	for my $k(keys %paf_data){
		my $genome = $paf_data{$k}->[0];

		my $seqkey;
		if($paf_data{$k}->[0] le $paf_data{$k}->[2]){
			$seqkey = join '-', ($paf_data{$k}->[1], $paf_data{$k}->[3]);
			if (exists $paf_merged{$seqkey} ){ 
				push @{ $paf_merged{$seqkey} }, (join '|', $paf_data{$k}->[4],
								$paf_data{$k}->[5],$paf_data{$k}->[6],
								$paf_data{$k}->[7],$paf_data{$k}->[8],
								$paf_data{$k}->[9],$paf_data{$k}->[10],$paf_data{$k}->[11]);
			}else{	
				$paf_merged{$seqkey} =	[$paf_data{$k}->[0], $paf_data{$k}->[1], $paf_data{$k}->[2], $paf_data{$k}->[3],
							(join '|', $paf_data{$k}->[4],
                                                                $paf_data{$k}->[5],$paf_data{$k}->[6],
                                                                $paf_data{$k}->[7],$paf_data{$k}->[8],
                                                                $paf_data{$k}->[9],$paf_data{$k}->[10],$paf_data{$k}->[11])
							];
			}
		}else{ 
			$seqkey = join '-', ($paf_data{$k}->[3], $paf_data{$k}->[1]);
			if( exists $paf_merged{$seqkey}){ 
				push @{ $paf_merged{$seqkey} }, (join '|', $paf_data{$k}->[4],
                                                                $paf_data{$k}->[5],$paf_data{$k}->[6],
                                                                $paf_data{$k}->[7],$paf_data{$k}->[8],
                                                                $paf_data{$k}->[9],$paf_data{$k}->[10],$paf_data{$k}->[11]); 
			}else{
				$paf_merged{$seqkey} = [$paf_data{$k}->[2], $paf_data{$k}->[3], $paf_data{$k}->[0], $paf_data{$k}->[1],
                                                        (join '|', $paf_data{$k}->[4],
                                                                $paf_data{$k}->[5],$paf_data{$k}->[6],
                                                                $paf_data{$k}->[7],$paf_data{$k}->[8],
                                                                $paf_data{$k}->[9],$paf_data{$k}->[10],$paf_data{$k}->[11]) 
						
							];
			}
		}
	}
  	return \%paf_merged;
}


#======================================================================
1;



__END__
