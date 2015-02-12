#!/bin/env perl

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
use File::Path qw(make_path remove_tree);

#blastp parameters

our $E=10 ; #E-value
our $MAXTARGET=10; #Number of hits to output
our $THREADS=5;
our $makeblastdb_prog="/opt/hpc/bin/makeblastdb";
#makeblastdb -dbtype prot -in proteins.fasta -input_type fasta -parse_seqids &
our $blastp_prog="/opt/hpc/bin/blastp";
our $seq_load=30;

$VERSION = sprintf "%d.%02d", q$Revision: 1.10 $ =~ /(\d+)\.(\d+)/;
pod2usage(0) unless $ARGV[0];

my ($_help, $_version, $sdbhost, $sdbport, $sdbuser, $sdbpass, $sdbname, 
    $sdnadb, $tdbhost, $tdbport, $tdbuser, $tdbpass, $tdbname, $idfile, $type, 
    $istrans, $coord_system, $addid, $idbase, $keepid, $addscore, $newtype,
    $glen, $plen, $rerun, $usechr, $baselen, $logic_name, $dump_dir);
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
	   'rerun'          => \$rerun,
	   'logic_name=s'   => \$logic_name,
           'dump_dir=s'     => \$dump_dir,
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
$dump_dir ||= '.';

my $database_subdir = 'reference';
my $database_dir = "$dump_dir/$database_subdir";
my $database_basename = 'pepdb.fa';
my $database_targzed_path = $database_dir.".tar.gz";
my $query_dir = "$dump_dir/query";
my $query_basename = 'query.fa';
my $blastp_output_dir = "$dump_dir/results";

foreach my $d ($database_dir, $query_dir, $blastp_output_dir){

    if( ! -d $d){
	make_path ($d);
    }
}

unless( $rerun){
  
#fetch protome of the two version annotations
#the old version protome build as blastp database
#the new version protome split into query files

    foreach my $d ($database_dir, $query_dir, $blastp_output_dir){

	if( -d $d){
	    remove_tree( $d );
	}
	make_path ($d);
    }

   
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

    $coord_system = 'toplevel';
    print "coord_system is $coord_system\n";

    my $slice_adaptor = $tdba->get_SliceAdaptor;
    my $transcripts;

    my $transcripts_src = get_all_transcripts(
	$sdba, $type, $coord_system, $logic_name);
    print STDERR "source transcripts: ", scalar @$transcripts_src, "\n";
    
    my $transcripts_target = get_all_transcripts(
	$tdba, $type, $coord_system, $logic_name);
    print STDERR "target transcripts: ", scalar @$transcripts_target, "\n";


    prepare_blastp_database(
	$transcripts_target, $database_dir, $database_basename);

    prepare_blastp_query_files($transcripts_src, $query_dir, $query_basename);

}

my $query_count_r = qx{ls $query_dir/${query_basename}.* | wc };
my $query_count = 0;
if( $query_count_r =~ /(\d+)/){
    $query_count = $1;
}

print "There are $query_count query files\n";

#create sge shell script
my $sge_job_script = "$dump_dir/sge_job.sh";

open my $fh, '>', $sge_job_script or die "cannot open $sge_job_script to write\n"; 


my $script_codes = <<SGEJOBCODES;
#!/bin/sh

## -- begin embedded SGE arguments --
#\$ -cwd
#\$ -S /bin/sh
#\$ -j y
#\$ -o ./sge-out/\$JOB_NAME.o\$JOB_ID
#\$ -e ./sge-err/\$JOB_NAME.e\$JOB_ID
## -- end embedded SGE arguments --


QUERY="$query_dir/${query_basename}.\$SGE_TASK_ID"
TMPQUERY="\$TMPDIR/${query_basename}.\$SGE_TASK_ID"
BOUTPUTDIR="\$CWD/$blastp_output_dir/"


cp \$QUERY \$TMPQUERY 
cp $database_targzed_path \$TMPDIR
cd \$TMPDIR
gunzip ${database_subdir}.tar.gz
tar -xf $database_subdir.tar
rm $database_subdir.tar

DB="$database_dir/$database_basename"

BOUTPUT="\$TMPDIR/blastpout.\$SGE_TASK_ID"
$blastp_prog -query \$TMPQUERY -out \$BOUTPUT -db \$DB -evalue $E -outfmt '6 qseqid sseqid qlen slen evalue score qstart qend sstart send gapopen gaps pident length' -max_target_seqs $MAXTARGET -num_threads=$THREADS 

gzip \$BOUTPUT
ZBOUTPUT="\$BOUTPUT.gz"
mv \$ZBOUTPUT \$BOUTPUTDIR 

SGEJOBCODES

print $fh $script_codes;

close $fh;

#submit job arrays
print "Submitting a blastp array job";
my $cmd = <<CMD;
export CWD=\$PWD
qsub -pe threads 6 -l m_mem_free=3.5G -v CWD -t 1-${query_count}:1 -N SwBlastp -S /bin/sh -cwd $sge_job_script
CMD

print "\n$cmd\n";


my $r = system( $cmd );

$r == 0 ? print "Succeed":"Failed";




#====================================
# Function to prepare blastp database
#

sub prepare_blastp_database{
    
    my $trpts = shift;
    my $dir = shift;
    my $database_name = shift;
#=stub
    my $bioseq_io = Bio::SeqIO->new(
	-format => 'Fasta',
	-file     => ">$dir/$database_name");
    
    
    foreach my $atrpt (@$trpts){

	my $pepseq_obj = $atrpt->translate();
	my $trpt_chr   = $atrpt->seq_region_name;
	my $trpt_start = $atrpt->seq_region_start;
	my $trpt_end   = $atrpt->seq_region_end;
	my $trpt_strand = $atrpt->seq_region_strand;
	my $trpt_sid   = $atrpt->stable_id;
	my $gene_sid   = $atrpt->get_Gene->stable_id;
	
	my $newID = join '|', 
	( $pepseq_obj->display_id, 
	  (join ':', ( $trpt_chr,  $trpt_start,  $trpt_end, $trpt_strand)),
	  $trpt_sid,
	  $gene_sid
	);

	$pepseq_obj->display_id( $newID);

	$bioseq_io->write_seq($pepseq_obj);
	
    }

    $bioseq_io->close;
   
#=cut
 
    my $cmd = qq{
     ${makeblastdb_prog} -dbtype prot -in $dir/$database_name  -input_type fasta
     tar -cvf $dir.tar $dir
     gzip $dir.tar
    };
    
    print "build blastp database\n$cmd\n ";

    system( $cmd ) == 0 or die "Faild build blastp database";

    return;

}

sub prepare_blastp_query_files{

    my $trpts      = shift;
    my $dir        = shift;
    my $query_name = shift;

    my $index = 1;

    my $bioseq_io = Bio::SeqIO->new(
	-format => 'Fasta',
	-file   => ">$dir/${query_name}.".$index++);
    
    my $number = 1;
    foreach my $atrpt (@$trpts) {
	
	my $pepseq_obj = $atrpt->translate();
	my $trpt_chr   = $atrpt->seq_region_name;
	my $trpt_start = $atrpt->seq_region_start;
	my $trpt_end   = $atrpt->seq_region_end;
	my $trpt_strand = $atrpt->seq_region_strand;
	my $trpt_sid   = $atrpt->stable_id;
	my $gene_sid   = $atrpt->get_Gene->stable_id;
	
	my $newID = join '|', 
	( $pepseq_obj->display_id, 
	  (join ':', ( $trpt_chr,  $trpt_start,  $trpt_end, $trpt_strand)),
	  $trpt_sid,
	  $gene_sid
	);

	$pepseq_obj->display_id( $newID);

	$bioseq_io->write_seq($pepseq_obj);
	
	if ( $number++ >= 30){
	    
	    $number = 1;
	    $bioseq_io->close;
	    $bioseq_io = Bio::SeqIO->new(
		-format => 'Fasta',
		-file   => ">$dir/${query_name}.".$index++);
	}
	
    }

    $bioseq_io->close;

}




###########################################################


# retrieve all transcripts on rice all chromosomes
sub get_all_transcripts {
    my ($dba, $type, $coord_system, $logic_name) = @_;
    my $slice_adaptor = $dba->get_SliceAdaptor;
    my @transcripts;
    my @slices = @{$slice_adaptor->fetch_all($coord_system)};
    print STDERR scalar @slices, " $coord_system(s)\n";
    for my $slice (@slices) {
#	my $repeats = $slice->get_all_RepeatFeatures();
#	print scalar @$repeats, " repeat regions\n";
	my @g;
	if ($type && $type =~ /^predict/) { 
	    @g = @{$slice->get_all_PredictionTranscripts()};
	}elsif ($type) { 
	    @g = @{$slice->get_all_Transcripts(1, undef, $type)};
	}else { 
	    @g = @{$slice->get_all_Transcripts(1, $logic_name, undef)};
	}
	
	push @transcripts, sort {$a->seq_region_start <=> $b->seq_region_start} @g;
    }
    return \@transcripts;
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
  -rerun            reuse the existing blastDB and query files
  -dump_dir         the working directory

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

