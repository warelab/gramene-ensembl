#!/software/bin/perl

use strict;
use DBI;
use IO::Uncompress::Gunzip qw (gunzip);
use IO::Compress::Gzip qw (gzip);
use Digest::MD5 qw (md5);

# Params are ensembl release (44), species short (hsapiens), link version (36) dataset ID (68)
my ($release,$short,$linkVersion,$dsID) = @ARGV;

# Open sequence_mart_X
my $user =     $ENV{'ENSMARTUSER'};
my $password =     $ENV{'ENSMARTPWD'};
my $host =     $ENV{'ENSMARTHOST'};
my $port =     $ENV{'ENSMARTPORT'};
my $driver =     $ENV{'ENSMARTDRIVER'};
my $dsn = "DBI:$driver:database=sequence_mart_$release;host=$host;port=$port";
my $dbh = DBI->connect($dsn,$user,$password);

# Read ds conf and template conf for hsapiens
my $sql1 = "SELECT xml FROM sequence_mart_$release.meta_conf__xml__dm WHERE dataset_id_key=?";
my $sth1 = $dbh->prepare($sql1) || die $dbh->errstr;
my $dsxml;
$sth1->execute($dsID);
my $row = $sth1->fetchrow_arrayref();
$dsxml = $$row[0];

my $sql2 = "SELECT compressed_xml FROM sequence_mart_$release.meta_template__xml__dm WHERE template='${short}_genomic_sequence'";
my $sth2 = $dbh->prepare($sql2) || die $dbh->errstr;
my $txml;
$sth2->execute();
my $row = $sth2->fetchrow_arrayref();
my $comptxml = $$row[0];
gunzip \$comptxml => \$txml;

# Update by grepping species names and dataset IDs
my $oldLinkVersion = $linkVersion - 1;
$txml =~ s/${short}_${oldLinkVersion}/${short}_${linkVersion}/g;
$dsxml =~ s/${short}_${oldLinkVersion}/${short}_${linkVersion}/g;

# Write them back with gzip
my $comptxml;
my $compdsxml;
gzip \$txml => \$comptxml;
gzip \$dsxml => \$compdsxml;
my $dsmd5 = md5($dsxml);
my $tmd5 = md5($txml);

my $ins1 = "update sequence_mart_$release.meta_conf__xml__dm set xml=?,compressed_xml=?,message_digest=? where dataset_id_key=$dsID";
$sth1 = $dbh->prepare($ins1);
$sth1->execute($dsxml, $compdsxml, $dsmd5);
my $ins2 = "update sequence_mart_$release.meta_template__xml__dm set compressed_xml=? where template='${short}_genomic_sequence'";
$sth2 = $dbh->prepare($ins2);
$sth2->execute($comptxml);
