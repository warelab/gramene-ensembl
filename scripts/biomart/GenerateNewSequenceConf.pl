#!/usr/local/bin/perl

use strict;
use DBI;
use IO::Uncompress::Gunzip qw (gunzip);
use IO::Compress::Gzip qw (gzip);
use Digest::MD5 qw (md5);

# Params are ensembl release (44), species short (hsapiens), link version (36)
# dataset ID (68)
my ( $release, $short, $linkVersion, $dsID ) = @ARGV;

# Open sequence_mart_X
my $user     = $ENV{'ENSMARTUSER'};
my $password = $ENV{'ENSMARTPWD'};
my $host     = $ENV{'ENSMARTHOST'};
my $port     = $ENV{'ENSMARTPORT'};
my $driver   = $ENV{'ENSMARTDRIVER'};
my $dsn = "DBI:$driver:database=sequence_mart_$release;host=$host;port=$port";
my $dbh = DBI->connect( $dsn, $user, $password );

# Read ds conf and template conf for hsapiens
my $sql1
    = "SELECT xml FROM sequence_mart_$release.meta_conf__xml__dm WHERE dataset_id_key=800";
my $sth1 = $dbh->prepare( $sql1 ) || die $dbh->errstr;
my $dsxml;
$sth1->execute();
my $row = $sth1->fetchrow_arrayref();
$dsxml = $$row[0];

my $sql2 = qq[
    SELECT compressed_xml 
    FROM   sequence_mart_$release.meta_template__xml__dm 
    WHERE  template='athaliana_genomic_sequence'
];
my $sth2 = $dbh->prepare( $sql2 ) || die $dbh->errstr;
my $txml;
$sth2->execute();
my $row      = $sth2->fetchrow_arrayref();
my $comptxml = $$row[0];
gunzip \$comptxml => \$txml;

# Update by grepping species names and dataset IDs
$txml  =~ s/athaliana/${short}/g;
$txml  =~ s/datasetID="800"/datasetID="${dsID}"/g;
$txml  =~ s/${short}_36/${short}_${linkVersion}/g;
$dsxml =~ s/athaliana/${short}/g;
$dsxml =~ s/datasetID="800"/datasetID="${dsID}"/g;
$dsxml =~ s/${short}_36/${short}_${linkVersion}/g;

# Write them back with gzip
my $comptxml;
my $compdsxml;
gzip \$txml  => \$comptxml;
gzip \$dsxml => \$compdsxml;
my $dsmd5 = md5( $dsxml );
my $tmd5  = md5( $txml );

my $ins1
    = "insert into sequence_mart_$release.meta_conf__xml__dm(dataset_id_key,xml,compressed_xml,message_digest) values ($dsID,?,?,?)";
$sth1 = $dbh->prepare( $ins1 );
$sth1->execute( $dsxml, $compdsxml, $dsmd5 );
my $ins2
    = "insert into sequence_mart_$release.meta_template__xml__dm(template,compressed_xml) values ('${short}_genomic_sequence',?)";
$sth2 = $dbh->prepare( $ins2 );
$sth2->execute( $comptxml );

#
# Write other important meta table stuff.
my $ins3
    = "insert into sequence_mart_$release.meta_conf__dataset__main(dataset_id_key,dataset,type,visible) values ($dsID,'${short}_genomic_sequence','GenomicSequence',0)";
my $sth3 = $dbh->prepare( $ins3 );
$sth3->execute();
my $ins4
    = "insert into sequence_mart_$release.meta_conf__user__dm(dataset_id_key,mart_user) values ($dsID,'default')";
my $sth4 = $dbh->prepare( $ins4 );
$sth4->execute();
$ins4
    = "insert into sequence_mart_$release.meta_conf__interface__dm(dataset_id_key,interface) values ($dsID,'default')";
$sth4 = $dbh->prepare( $ins4 );
$sth4->execute();
my $ins5
    = "insert into sequence_mart_$release.meta_template__template__main(dataset_id_key,template) values ($dsID,'${short}_genomic_sequence')";
my $sth5 = $dbh->prepare( $ins5 );
$sth5->execute();

