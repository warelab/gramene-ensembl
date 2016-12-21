#!/usr/bin/env perl

# Copyright [2009-2014] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
use Net::FTP;
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;
use Getopt::Long;
use FileHandle;
use Date::Calc qw(:all);

my ($dbhost, $dbport, $dbuser, $dbpass);
my $db_name = "ensembl_compara_master";
my $help = 0;
my ($year,$month,$day) = Today([localtime()]); 
my $today = sprintf "%4d%02d%02d", $year,$month,$day;
my $remote_db_dump_filename = "${db_name}_$today.sql.gz";
my $local_dbname = "${db_name}_today";

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


my $remote_file_url = "ftp://ftp.ensemblgenomes.org/.transfer/compara_master/$remote_db_dump_filename";

my $download_cmd = "cd /tmp; rm -rf $remote_db_dump_filename; wget $remote_file_url";
my $create_empty_db = qq(mysql -h $dbhost -P $dbport -u $dbuser -p$dbpass -e "drop database if exists $local_dbname; create database $local_dbname;");
my $sqlimport_command = "gunzip -c /tmp/$remote_db_dump_filename | mysql -h $dbhost -P $dbport -u $dbuser -p$dbpass $local_dbname";  

my $result = qx{$download_cmd};
die "Failed to download $remote_db_dump_filename: $download_cmd" unless (defined $result);
print "$result\n";

$result = qx{$create_empty_db};
die "Failed to create empty database: $create_empty_db" unless (defined $result);
print "$result\n";

$result = qx{$sqlimport_command};
die "Failed to gunzip and import today's compara master db: $sqlimport_command" unless (defined $result);
print "$result\n";

system("rm -rf /tmp/$remote_db_dump_filename");



sub usage {
    my $script = undef;
    ($script = $0) =~ s!.*/!!;
    print STDERR "Usage: perl $script -host <db_host> -port <db_port> -user <db_user> -pass <db_pass> [-h|-help]\n";
}
