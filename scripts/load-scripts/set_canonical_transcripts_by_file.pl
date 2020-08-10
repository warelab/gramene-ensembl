#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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


# Script that selects the best candidate for canonical transcripts on each
# gene.
# For usage instructions, run ./select_canonical_transcripts.pl -help


use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::TranscriptSelector;
use IO::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

my ($host, $port, $dbname, $user,$pass);
my ($log_path,$help,$idfile);

Bio::EnsEMBL::Registry->no_cache_warnings(1);

my $write = 0;
my $verbose = 0;

GetOptions( 'dbhost:s'            => \$host,
            'dbport:n'            => \$port,
            'dbname:s'            => \$dbname,
            'dbuser:s'            => \$user,
            'dbpass:s'            => \$pass,
	    'idfile:s'		  => \$idfile,
            'write!'              => \$write,
            'verbose!'            => \$verbose, 
            'log:s'               => \$log_path, # log file used for analysing choices in bulk
            'help!'               => \$help,
            'h!'                  => \$help);

if ($help) { &usage; exit;}
unless ($write) {
  print "You have not used the -write option "
      . "so results will not be written into the database\n";
}

my $fh;

open $fh, '<', $idfile or die "Cannot open $idfile to read";
my @ids; 
chomp(@ids = <$fh>);

my $dba =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -port   => $port,
                                      -dbname => $dbname,
                                      -pass   => $pass,
                                      -species => 'default',
                                      -no_cache => 1,
                                      );
                                      
my $log_fh;
if ($log_path) {
    $log_fh = IO::File->new($log_path,"w") or throw ("Could not open logging file.");
}    


print "The 1st 10 transcript ids are\n";
foreach my $idx (1..10) {
	print "$ids[$idx]\n";
}


print "Set Canonical transcript for ". scalar @ids ." genes\n";
if ($log_fh) {$log_fh->close;}

#__END__


## Change database entries.

if ($write) {
    my $gene_update_sql = "UPDATE gene g join transcript t using (gene_id) SET g.canonical_transcript_id = t.transcript_id where t.stable_id = ?";
    my $gene_sth = $dba->dbc->prepare($gene_update_sql);
    
    print "Updating database with new canonical transcripts...\n";
    foreach my $trpt_stable_id (@ids) {
        print "Set $trpt_stable_id as canonical transcript for its gene\n" if $verbose;
        $gene_sth->execute( $trpt_stable_id );
    }
    print "Done\n";
}



sub usage {
print "
Example usage: perl set_canonical_transcripts_by_file.pl -dbhost host -dbuser user 
     -dbpass *** -dbname dbname -dbport 3306 -idfile transcript_stable_ids_file -write

Script options:

    -dbname       Database name

    -dbhost       Database host

    -dbport       Database port

    -dbuser       Database user

    -dbpass       Database password

    -idfile      file contains the transcript stable ids that's determined to be canonical or primary transcript of the gene 

    -write                Specify if results should be written to the database

    -verbose              Increase verbosity of output messages

    -log                  Dump decision matrices into a log file for analysis

";
    
}
