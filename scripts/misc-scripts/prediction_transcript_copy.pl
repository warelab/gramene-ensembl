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
my $genes;
my $transcripts;

$genes = get_all_genes($sdba, $coord_system, $logic_name);
print STDERR "genes to copy: ", scalar @$genes, "\n";

my $genedba = $tdba->get_PredictionTranscriptAdaptor();
$genedba->store(@{$genes});


#my $num=0;
#foreach my $i (0..@$genes-1) {
#    my $gene = $genes->[$i];
#    $gene->load;
#    
#    if ($write) {
#	$genedba->store($gene);	
#	$num++;
#    }
    
#    print "$num prediction transcript copied\n" if $write;
   # last;
#}


exit;



###########################################################


# retrieve all transcripts on rice all chromosomes
sub get_all_genes {
    my ($dba, $coord_system, $logic_name) = @_;
    my $slice_adaptor = $dba->get_SliceAdaptor;
    my @genes;
    my @slices = @{$slice_adaptor->fetch_all($coord_system)};
    print STDERR scalar @slices, " $coord_system(s)\n";
    for my $slice (@slices) {
#	my $repeats = $slice->get_all_RepeatFeatures();
#	print scalar @$repeats, " repeat regions\n";
	my @g;
	
	@g = @{$slice->get_all_PredictionTranscripts()};
	
	push @genes, sort {$a->seq_region_start <=> $b->seq_region_start} @g;
    }
    return \@genes;
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

