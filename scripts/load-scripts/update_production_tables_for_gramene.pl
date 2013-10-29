#!/bin/env perl

use strict;
use warnings;
use Net::FTP;
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;
use Getopt::Long;
use FileHandle;

my ($dbhost, $dbport, $dbuser, $dbpass);
my $db_name = "ensembl_production";
my $help = 0;

GetOptions( "host=s",        \$dbhost,
            "user=s",        \$dbuser,
            "pass=s",        \$dbpass,
            "port=i",        \$dbport,
	    "help|h",        \$help,
          );

if (!($dbhost && $dbport && $dbuser && $dbpass) || $help) {
    usage();
    exit 1;
}

my @tables = ("biotype", "master_external_db", "master_attrib_type", "master_misc_set", "master_unmapped_reason", "analysis_description", "web_data");

my $host = "ftp.ensemblgenomes.org";
my $user = "anonymous";
my $pass = "";

my $ftp_connection = Net::FTP->new($host, Debug => 0, Passsive => 0)
      or die "Cannot connect to $host: $@";
$ftp_connection->login($user,$pass)
      or die "Cannot login ", $ftp_connection->message;

$ftp_connection->cwd("/.transfer/production_database")
      or die "Cannot change working directory ", $ftp_connection->message;

my $remote_sum_filename = "production_tables.sum";

my $temp_sum_fh = File::Temp->new( UNLINK => 0, DIR => "/tmp", SUFFIX => ".sum");
my $temp_sum_filename = $temp_sum_fh->filename;
$temp_sum_fh->close();
$ftp_connection->get($remote_sum_filename, $temp_sum_filename)
      		or die "get failed ", $ftp_connection->message;
$temp_sum_fh->close;

my $sums_href = {};
my $sum_fh = FileHandle->new();
$sum_fh->open("<$temp_sum_filename") or die "can't open checksum file, $temp_sum_filename!\n";
while (<$sum_fh>) {
    my $line = $_;
    chomp($line);
    $line =~ /^([^\t]+)\t(.+)/;
    my $table = $1;
    my $sum = $2;
    
    $sums_href->{$table} = $sum;
}
$sum_fh->close();

foreach my $table (@tables) {

    # Get the checksum for the current table
    
    print STDERR "running the sum command over table, $table...\n";

    my $sum_command = "mysql -h $dbhost -P $dbport -u $dbuser -p$dbpass $db_name --column-names=false -e 'checksum table $table'";

    my $sum_result = qx/$sum_command/;
    chomp ($sum_result);

    if (! defined $sum_result || length($sum_result) == 0) {
	print STDERR "$sum_command\n";
	die "Couldn't get the checksum for table, $table\n";
    }
    
    $sum_result =~ /^([^\t]+)\t(.+)/;
    my $gramene_sum = $2;
    
    my $ebi_sum = $sums_href->{"$db_name.".$table};
    
    # print STDERR "gramene_sum, ebi_sum, $gramene_sum, $ebi_sum\n";
    
    # Update the production database only the checksums differ

    if ("$ebi_sum" != "$gramene_sum") {

	# Todo: Replace by ftp, and onto a temp file

	my $remote_table_filename = $table . ".sql";

	my $temp_table_fh = File::Temp->new( UNLINK => 0, DIR => "/tmp", SUFFIX => ".sql");
	my $temp_table_filename = $temp_table_fh->filename;
	$temp_table_fh->close();
	$ftp_connection->get($remote_table_filename, $temp_table_filename)
	    or die "get failed ", $ftp_connection->message;

        if (! -f $temp_table_filename) {
            die "can't find file, $temp_table_filename, found!\n";
        }
        elsif (-z $temp_table_filename) {
            die "file, $temp_table_filename, is empty!\n";
        }
	
        my $command = "cat $temp_table_filename | mysql -h $dbhost -P$dbport -u $dbuser -p$dbpass $db_name";

	print STDERR "synching table, $table...\n";
	print STDERR "$command\n";

        my $return = qx/$command/;
        if ($return) {
            print STDERR "comamnd, $command, has failed!\n";
        }
    }
}


sub usage {
    my $script = undef;
    ($script = $0) =~ s!.*/!!;
    print STDERR "Usage: perl $script -host <db_host> -port <db_port> -user <db_user> -pass <db_pass> [-h|-help]\n";
}
