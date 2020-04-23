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
my ($dnahost, $dnaport, $dnadbname, $dnauser, $dnapass);
my ($ccds_host, $ccds_dbname, $ccds_user, $ccds_port, $ccds_pass);
my ($log_path,$help);

Bio::EnsEMBL::Registry->no_cache_warnings(1);

my $coord_system_name;
my $seq_region_name;
my $logic_name; # keep as undefined unless you only want to run on a specific analysis
my $write = 0;
my $include_non_ref = 1;
my $include_duplicates;
my $verbose = 0;

GetOptions( 'dbhost:s'            => \$host,
            'dbport:n'            => \$port,
            'dbname:s'            => \$dbname,
            'dbuser:s'            => \$user,
            'dbpass:s'            => \$pass,
            'dnadbhost:s'         => \$dnahost,
            'coord_system_name:s' => \$coord_system_name,
            'seq_region_name:s'   => \$seq_region_name,
            'logic_name:s'        => \$logic_name,
            'write!'              => \$write,
            'include_non_ref!'    => \$include_non_ref,
            'include_duplicates'  => \$include_duplicates,
            'verbose!'            => \$verbose, 
            'log:s'               => \$log_path, # log file used for analysing choices in bulk
            'help!'               => \$help,
            'h!'                  => \$help);

if ($help) { &usage; exit;}
unless ($write) {
  print "You have not used the -write option "
      . "so results will not be written into the database\n";
}

my $dba =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -port   => $port,
                                      -dbname => $dbname,
                                      -pass   => $pass,
                                      -species => 'default',
                                      -no_cache => 1,
                                      );
                                      
my $dna = check_if_DB_contains_DNA($dba);
  if(!$dna) {
    throw ("Your gene DB contains no DNA. You must provide a DNA_DB to connect to");
}

my $log_fh;
if ($log_path) {
    $log_fh = IO::File->new($log_path,"w") or throw ("Could not open logging file.");
}    

#my $transcript_selector = Bio::EnsEMBL::Utils::TranscriptSelector->new($ccds_dba);

my $slice_adaptor = $dba->get_SliceAdaptor;
my $slices;
if ($seq_region_name) {
    $slices = [$slice_adaptor->fetch_by_region($coord_system_name,
                              $seq_region_name,
                              $include_non_ref,
                              $include_duplicates) ];
} else {
    if (!$coord_system_name) {throw 'Requires a coordinate system name to function in this mode';}
    $slices = $slice_adaptor->fetch_all($coord_system_name,
                                        '',
                                        $include_non_ref,
                                        $include_duplicates);
}
my $canonical_changes = 0;
my $total_genes = 0;

my @change_list;

foreach my $slice (@$slices) {
    my $genes = $slice->get_all_Genes($logic_name, undef, 1);
    while (my $gene = shift @$genes) {
        $total_genes++;
        my $new_canonical = select_canonical_transcript_for_Gene($gene);
        
            # Original canonical transcript is now absent, or never set.
        if ($log_fh) {
              print $log_fh "//\n";
              printf $log_fh "canonical_transcript '%s'\n", ($new_canonical->stable_id) ;
         }
         push @change_list,[$gene->dbID,$new_canonical->dbID];
         $canonical_changes++;
        
    }

}


print "Canonical transcript alterations: ".$canonical_changes." from ".$total_genes." genes\n";
if ($log_fh) {$log_fh->close;}


## Change database entries.

if ($write) {
    my $gene_update_sql = "UPDATE gene SET canonical_transcript_id = ? where gene_id = ?";
    my $gene_sth = $dba->dbc->prepare($gene_update_sql);
    
    print "Updating database with new canonical transcripts...\n";
    foreach my $change (@change_list) {
        print "Changin' ".$change->[1]." on ".$change->[0]."\n" if $verbose;
        $gene_sth->execute( $change->[1], $change->[0]);
    }
    print "Done\n";
}


sub check_if_DB_contains_DNA {
  my ($dba)        = @_;
  my $sql_command = "select count(*) from dna";
  my $sth         = $dba->dbc->prepare($sql_command);
  $sth->execute();
  my @dna_array = $sth->fetchrow_array;
  if ( $dna_array[0] > 0 ) {
    print "Your DB "
      . $dba->dbc->dbname
      . " contains DNA sequences. No need to attach a "
      . "DNA_DB to it.\n"
      if ( $verbose );
    return 1;
  } else {
    print "Your DB " . $dba->dbc->dbname . " does not contain DNA sequences.\n"
      if ( $verbose );
    return 0;
  }
}

sub select_canonical_transcript_for_Gene {
    my $gene = shift;

    if (!$gene) {throw ('Cannot select canonical transcript without a gene.');}

    my $transcript_array = $gene->get_all_Transcripts;
    my @transcripts;
    if ($transcript_array) {
        @transcripts = @$transcript_array;
    } else {
        warning('No transcripts attached to gene '.$gene->stable_id);
        return;
    }
    my @encoded; # array of encoded transcripts

    foreach my $transcript (@transcripts) {
        my $encoded_transcript = encode_transcript($transcript);
        push(@encoded, $encoded_transcript );
        if ( $verbose ) {
            printf "%s encoded to: [%s,%s,%s,%s,%s,%s]\n",$transcript->stable_id, @$encoded_transcript;
        }
    }

    my $sorted_ids = sort_into_canonical_order(\@encoded);
    if ( $verbose ) {
        print "Sorted order: \n";
        foreach my $dbID (@$sorted_ids) {
            print $dbID."\n";
        }
    }
    my $canonical_dbID = $sorted_ids->[0];

    foreach my $transcript (@transcripts) {
        if ($transcript->dbID == $canonical_dbID) {
            if ($verbose ) {print 'Chosen transcript: '.$transcript->stable_id."\n";}
            return $transcript;
        }
    }
    throw ('Run out of transcripts without finding selected canonical dbID.')
}



sub encode_transcript {
    my $transcript = shift;

    my $type;

    my $translation = $transcript->translate;
    my $translation_length = 0;
    if ($translation) {
        $translation_length = $translation->length;
        # Translations containing premature stops are undesirable.
        if ($translation->seq =~ /\*/) {$translation_length = 0;}
     }
        #                 # Zero-length/non-existent translations are ok. We sort by coding and non-coding first
        my $translates = 0;
        if ($translation_length != 0) {$translates = 1;}

        my $transcript_length = $transcript->length;
	my ($aed) = @{$transcript->get_all_Attributes('AED')};
	
        return [$transcript->dbID,
                $translates,
                $translation_length,
                $transcript_length,
                $aed->value,
		$transcript->stable_id];
 }

sub sort_into_canonical_order {
    my $encoded_list_ref = shift;

    my @sorted_ids = map { $_->[0] }
        sort {
            # [0] contains ID
            $b->[1] <=> $a->[1] ||    # translates
            $b->[2] <=> $a->[2] ||    # translation length (largest is best, $a and $b reversed)
            $b->[3] <=> $a->[3] ||    # transcript length               "
	    $a->[4] <=> $b->[4] ||    # aed score, the smaller the better
            $a->[5] cmp $b->[5]       # stable id. All other things being equal, 'lowest' transcript ID wins
               } @{$encoded_list_ref};
            
     return \@sorted_ids;
 }
           
          

sub usage {
print "
Example usage: perl set_canonical_transcripts.pl -dbhost host -dbuser user 
     -dbpass *** -dbname dbname -dbport 3306 -coord_system toplevel -write

Script options:

    -dbname       Database name

    -dbhost       Database host

    -dbport       Database port

    -dbuser       Database user

    -dbpass       Database password

Other optional arguments:

    -coord_system_name    Coordinate system to use

    -include_non_ref      Specify if the non_reference regions should 
                          be _excluded_. (default: include) 

    -include_duplicates   Specify if the duplicate regions should be 
                          _included_. eg. Human PAR on Y (default: exclude) 

    -seq_region_name      Chromosome name if running a single seq_region

    -write                Specify if results should be written to the database

    -verbose              Increase verbosity of output messages

    -log                  Dump decision matrices into a log file for analysis

To check the script has run correctly you can run the
CoreForeignKeys healthcheck:

./run-healthcheck.sh -d dbname -output problem CoreForeignKeys

A warning about not using CCDS is perfectly acceptible when not running on
Human, Mouse and Zebrafish.
";
    
}
