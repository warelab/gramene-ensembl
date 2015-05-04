#!/usr/bin/perl -w

#transfer/dump genes. The source genes/transcripts can be all 
#transcripts in query db or use IDs from an input file

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';
use vars qw[ $VERSION ];
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::CoordSystem;


$VERSION = sprintf "%d.%02d", q$Revision: 1.10 $ =~ /(\d+)\.(\d+)/;
pod2usage(0) unless $ARGV[0];

my ($_help, $_version, $sdbhost, $sdbport, $sdbuser, $sdbpass, $sdbname, 
    $sdnadb, $tdbhost, $tdbport, $tdbuser, $tdbpass, $tdbname, $idfile, $type, 
    $istrans, $coord_system, $addid, $idbase, $keepid, $addscore, $newtype,
    $glen, $plen, $write, $usechr, $baselen, $logic_name);
GetOptions(
	   'h|help'        => \$_help,
	   'v|version'     => \$_version,
	   'dbhost=s'      => \$sdbhost,
	   'dbport=s'      => \$sdbport,
	   'dbuser=s'      => \$sdbuser,
	   'dbpass=s'      => \$sdbpass,
	   'dbname=s'     => \$sdbname,
	   'dnadb=s'       => \$sdnadb,
	   'tdbhost=s'      => \$tdbhost,
	   'tdbport=s'      => \$tdbport,
	   'tdbuser=s'      => \$tdbuser,
	   'tdbpass=s'      => \$tdbpass,
	   'tdbname=s'     => \$tdbname,
	   'write'          => \$write,
	   'logic_name=s'   => \$logic_name,
	   ) or die;
pod2usage(2) if $_help;
if ( $_version ) {
    my $prog = basename( $0 );
    print "$prog v$VERSION\n";
    exit(0);
}

$sdbport ||= '3306';
$tdbhost ||= $sdbhost;
$tdbuser ||= $sdbuser;
$tdbport ||= $sdbport;
$tdbpass ||= $sdbpass;
$tdbpass ||= $sdbname;

my $contig_projection_table='map_seq_region_v0t3'; #in to database

print STDERR "Target: $tdbhost $tdbport $tdbname\n";

my $sdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					       -user   => $sdbuser,
					       -pass   => $sdbpass,
					       -dbname => $sdbname,
					       -host   => $sdbhost,
					       -port => $sdbport,
					       -driver => 'mysql',
					       );
my $tdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					       -user   => $tdbuser,
					       -pass   => $tdbpass,
					       -dbname => $tdbname,
					       -host   => $tdbhost,
					       -port => $tdbport,
					       -driver => 'mysql',
					       );

$coord_system = 'chromosome';
print "coord_system is $coord_system\n";

my $csa = $tdba->get_CoordSystemAdaptor;
my $cs = $csa->fetch_by_name($coord_system);
my $slice_adaptor = $tdba->get_SliceAdaptor;
my $genes;
my $transcripts;


$genes = get_all_genes($sdba, $type, $coord_system, $logic_name);
print "Starting with this number of genes to copy: ", scalar @$genes, "\n";


my $dbc = $tdba->dbc;

my $mapping_sql = "select v13chr, v13chr_start, v13chr_end, v10chr, v10chr_start, v10chr_end from map_seq_region_v0t3 where v10chr = ? and v10chr_start <= ? and v10chr_end >= ?";
my $stmt_handler = $dbc->prepare($mapping_sql) 
    or die "cannot prepare\n$mapping_sql\n";

my $v3_offset_sql = "select  v13chr_start, v10chr_start from map_seq_region_v0t3 where v10chr = ? and v10chr_start <= ? and v10chr_end >= ?";
my $v3coord_stmt_handler = $dbc->prepare($v3_offset_sql) 
    or die "cannot prepare\n$v3_offset_sql\n";


our $BioseqIO = Bio::SeqIO->new(-file => ">v1transferredPep.fa",
			       -format => 'Fasta',
    ) or die "Cannot open v1transferredPep.fa to write";

open LOGF, ">GenesNotCopied.log" or die "Cannot open GenesNotCopied.log to write";

my $genedba = $tdba->get_GeneAdaptor();

my $cnt_transferred = 0;
my $cnt_not_transferred = 0;
my $num=0;
foreach my $i (0..@$genes-1) {
    my $gene = $genes->[$i];
    $gene->load;
    my $g_stable_id = $gene->stable_id;

    my $gene_v13 = project_2v13chr($gene, $stmt_handler);
    
    if( $gene_v13 ){
	#next;
	if ($write) {
	    $genedba->store($gene_v13);	
	    $num++;
	    print "$num genes copied\n";
	    #last if $num > 100;
	}
    }else{

	my $gene_v13 = change_coordinates($gene, $stmt_handler);

	if (! $gene_v13 ){
	    ++$cnt_not_transferred;
	    print LOGF "gene $g_stable_id cannot be easily transferred to v1.3 chromosome\n";
	    next;
	}
	if ($write) {
	    eval{
		$genedba->store($gene_v13);	
	    };
	    
	    if( $@ ){
		print LOGF "gene $g_stable_id cannot be easily transferred to v1.3 chromosome, $@\n";

		next;
	    }
	    $num++;
	    print "$num genes copied\n";
	    #last if $num > 100;
	}
	
    }
    
    #last;
}

print "$cnt_transferred genes transferred\n$cnt_not_transferred genes not copied\n" ;
exit;



###########################################################

sub project_2v13chr {

    my $gene = shift;
    my $sth = shift;

    my $gv0chr = $gene->seq_region_name;
    my $gv0chr_start = $gene->seq_region_start;
    my $gv0chr_end   = $gene->seq_region_end;
    my $g_stable_id  = $gene->stable_id;

    $sth->execute($gv0chr, $gv0chr_start, $gv0chr_end);

    my $result = $sth->fetchall_arrayref;
    
    unless (@$result == 1){
	#++$cnt_not_transferred;
	#print LOGF "gene $g_stable_id cannot be easily transferred to v1.3 chromosome, there is ", @$result+0, " mapped v1.3 region for $gv0chr:$gv0chr_start-$gv0chr_end \n" ;
	return undef;
	
    }

    ++$cnt_transferred;
    my ($cv3chr, $cv3chr_start, $cv3chr_end, $cv0chr, $cv0chr_start, $cv0chr_end) = @{$result->[0]};

    print "gene $g_stable_id ($gv0chr:$gv0chr_start-$gv0chr_end) can be transferred to v1.3 region $cv3chr, $cv3chr_start, $cv3chr_end (v0 contig:$cv0chr, $cv0chr_start, $cv0chr_end)\n";

    
    map{ my $pep_seq = $_->translate; $BioseqIO->write_seq($pep_seq) } @{$gene->get_all_Transcripts}; 

    if($cv3chr eq $cv0chr && $cv3chr_start == $cv0chr_start && $cv3chr_end == $cv0chr_end){
	return $gene;

    }

    my $v0_slice = Bio::EnsEMBL::Slice->new(
	-coord_system => $cs,
	-start =>  $cv0chr_start,
	-end   => $cv0chr_end,
	-seq_region_name => $cv0chr,
	-adaptor => $slice_adaptor,
	);

    my $v3_slice = Bio::EnsEMBL::Slice->new(
	-coord_system => $cs,
	-start =>  $cv3chr_start,
	-end   => $cv3chr_end,
	-seq_region_name => $cv3chr,
	-adaptor => $slice_adaptor,
	);
    
    my $gv3 = $gene->transfer($v0_slice);
    
    #replace the slice to the v3 slice
    
    $gv3->slice( $v3_slice);
    
    if( exists $gv3->{'_transcript_array'} ) {

	for my $old_transcript ( @{$gv3->{'_transcript_array'}} ) {
	    $old_transcript->slice( $v3_slice );

	    if( exists $old_transcript->{'_trans_exon_array'} ) {
		for my $old_exon ( @{$old_transcript->{'_trans_exon_array'}} ) {
		    $old_exon->slice( $v3_slice );
		}
	    }

	    if( exists $old_transcript->{'_supporting_evidence'} ) {
		for my $old_feature ( @{$old_transcript->{'_supporting_evidence'}} ) {
		    $old_feature->slice( $v3_slice );
		}
	    }
  
	    if(exists $old_transcript->{_ise_array}) {
		foreach my $old_feature ( @{$old_transcript->{_ise_array}} ) {
		    $old_feature->slice( $v3_slice );
		}
	    }

	    # flush cached internal values that depend on the exon coords
	    $old_transcript->{'transcript_mapper'}   = undef;
	    $old_transcript->{'coding_region_start'} = undef;
	    $old_transcript->{'coding_region_end'}   = undef;
	    $old_transcript->{'cdna_coding_start'}   = undef;
	    $old_transcript->{'cdna_coding_end'}     = undef;
#
#====
	}
	
    }

    return $gv3;
}


sub project_start_end_coords{

    my $f = shift;

    my $chr = $f->seq_region_name;
    my $stable_id = $f->stable_id;
    
    #should I include coding_region_start coding_region_end?

    my $coord_pt = undef;
    foreach $coord_pt qw (start end){

	print "Now project $coord_pt->";
	if( my $v0c = $f->$coord_pt ){

	    print "project single coord $coord_pt $chr, $v0c \n";
	    my $v3c = project_single_coord( $v3coord_stmt_handler, $chr, $v0c);
	    unless (defined $v3c && $v3c >0){
		print "$stable_id $coord_pt coord $chr:$v0c cannot be easily transferred to v1.3 chromosome, not covered by mapped regions\n";
		return undef;
	    }
	    print "project single coord  $chr, $v0c to $v3c\n";
	    $f->$coord_pt( $v3c );
	    
	}
    }

    return 1;
    
}

sub project_single_coord{

    my $stmt = shift;
    my $chr = shift;
    my $old_coord = shift;

    eval{
	$stmt->execute($chr, $old_coord, $old_coord);
    };

    if( $@ ){
	print ("project_single_coord failed with $chr, $old_coord, $old_coord, $@\n");
	return undef;
    }
    my $result = $stmt->fetchall_arrayref;
    
    unless (@$result == 1){

	print "Error: Found ". scalar @$result . "mapped regions, should be 1 region\n";
	return undef;	
    }
    
    my ($v13chr_st, $v10chr_st) = @{$result->[0]};

    return $v13chr_st-$v10chr_st+$old_coord;
}


sub change_coordinates{

    my $gene = shift;
    my $sth = shift;


    #map{ my $pep_seq = $_->translate; $BioseqIO->write_seq($pep_seq) } @{$gene->get_all_Transcripts}; 

    my $gv3 = $gene->transform( $coord_system );

    unless (project_start_end_coords( $gv3 )){	 
	return undef;
    }
    
    if( exists $gv3->{'_transcript_array'} ) {

	for my $old_transcript ( @{$gv3->{'_transcript_array'}} ) {

	    #The translation start end should stay the same, since 
	    # they are the coord within the exon
	    #my $old_translation = $old_transcript->translation;
	    
	    unless (project_start_end_coords( $old_transcript )){
		return undef;
	    }
    
	    if( exists $old_transcript->{'_trans_exon_array'} ) {
		for my $old_exon ( @{$old_transcript->{'_trans_exon_array'}} ) {
		    unless (project_start_end_coords( $old_exon )){
			return undef;		
		    }		    
		}
	    }

	    if( exists $old_transcript->{'_supporting_evidence'} ) {
		for my $old_feature ( @{$old_transcript->{'_supporting_evidence'}} ) {
		    unless (project_start_end_coords( $old_feature )){
			return undef;
		    }
		}
	    }
  
	    if(exists $old_transcript->{_ise_array}) {
		foreach my $old_feature ( @{$old_transcript->{_ise_array}} ) {
		    unless (project_start_end_coords( $old_feature )){
			return undef;
		    }
		}
	    }

	    # flush cached internal values that depend on the exon coords
	    $old_transcript->{'transcript_mapper'}   = undef;
	    $old_transcript->{'coding_region_start'} = undef;
	    $old_transcript->{'coding_region_end'}   = undef;
	    $old_transcript->{'cdna_coding_start'}   = undef;
	    $old_transcript->{'cdna_coding_end'}     = undef;

	     
#
#====
	}
	
    }    

    my $g_stable_id = $gene->stable_id;
    ++$cnt_transferred;
    print "\ngene $g_stable_id can be transferred to v1.3 region\n";    
    map{ my $pep_seq = $_->translate; $BioseqIO->write_seq($pep_seq) } @{$gene->get_all_Transcripts}; 

    return $gv3;
}

# retrieve all transcripts on rice all chromosomes
sub get_all_genes {
    my ($adba, $atype, $acoord_system, $alogic_name) = @_;
    my $aslice_adaptor = $adba->get_SliceAdaptor;
    my @genes;
    my @slices = @{$aslice_adaptor->fetch_all($coord_system)};
    print STDERR scalar @slices, " $coord_system(s)\n";
    for my $slice (@slices) {
#	my $repeats = $slice->get_all_RepeatFeatures();
#	print scalar @$repeats, " repeat regions\n";
	my @g;
	if ($atype && $atype =~ /^predict/) { 
	    @g = @{$slice->get_all_PredictionTranscripts()};
	}elsif ($atype) { 
	    @g = @{$slice->get_all_Genes_by_type($atype)};
	}else { @g = @{$slice->get_all_Genes($alogic_name)};}
	push @genes, sort {$a->seq_region_start <=> $b->seq_region_start} @g;
    }
    return \@genes;
    #return [grep {$_->stable_id eq 'OMERI01G03050' } @genes];
}


# ----------------------------------------------------

=head1 NAME

gene-copy.pl

=head1 SYNOPSIS

  gene-copy.pl -sdbhost host -sdbport dbport -sdbuser user -sdbname name 
    -sdnadb dnadb -tdbhost host -tdbport port -tdbuser user -tdbname name 
    -idfile file -gene -type type -coord_system coord

Options:

  -h|--help        Show brief help and exit
  -v|--version     Show version and exit
  -dbhost=s         source db host
  -dbport=s        source db port
  -dbuser=s         db user 
  -dbname=s        db for query genes
  -dnadb=s          db for dna seq 
  -tdbhost=s        target  db host
  -sdbport=s        target db port
  -tdbuser=s         db user 
  -tdbname=s        db for target genes
  -logic_name=s     logic_name of that gene analysis
  -write

=head1 DESCRIPTION

To copy genes from one db to another 

=head1 SEE ALSO

perl.

=head1 AUTHOR

Sharon Wei E<lt>weix@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2005 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

