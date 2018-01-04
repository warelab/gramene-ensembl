#!/bin/env perl
  

=head1 Name

    parse_interpro2ensembl.pl - parse interproscan raw xml ouptut file and load into ensembl core database. 

=head1 SYNOPSIS

    parse_interpro2ensembl.pl OPTIONS  FILE1 FILE2 FILE3 ...

=head2 OPTIONS

    -dbhost         database hose
    -dbname	    database name
    -dbuser	    database user
    -dbpass	    database password
    -dbport	    database port default 3306
    -debug            debug mode, print more test messages
    -test             no insertions to the database
    -help
    
=head1 Author

Sharon Wei (weix@cshl.edu)


=cut


use lib map {"$ENV{ENS_CODE_ROOT}/$_" } 
	qw ( bioperl-live  
	     ensembl/modules 
             ensembl-compara/modules
	     Sharon/modules 
	     ensembl/misc-scripts/xref_mapping);


use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use XrefParser::Database;
use XrefMapper::Interpro;
use XrefParser::InterproParser;

use Readonly;

Readonly our $SRC => 'Interpro';
Readonly our $SPECIES_ID_METa_KEY => 'species.taxonomy_id';

my($help, $dbname, $dbhost, $dbuser, $dbpass, $dbport, $test, $debug);

GetOptions(
	   "help"              => \$help,
	   "dbname=s"          => \$dbname,
	   "dbhost=s"          => \$dbhost,
	   "dbuser=s"	       => \$dbuser,
	   "dbpass=s"          => \$dbpass,
	   "dbport=i"	       => \$dbport,
	   "test"              => \$test,
	   "debug"             => \$debug,
	   ) || pod2usage();

if($help || !$dbname || !$dbhost || !$dbuser || !$dbpass){

    pod2usage();
}


$dbport ||= 3306;

warn "get optionts\nhelp=>$help, debug=>$debug, $dbuser:$dbpass@$dbhost:$dbport/$dbname\n" if $debug ;

my $db_obj = XrefParser::Database->new( {
			'dbname' => $dbname,
			'host'   => $dbhost,
			'user'	=> $dbuser,
			'pass'	=> $dbpass,
			'port'	=> $dbport,
			'verbose' => $debug,
		});
die "Cannot create XrefParser::Database object for $dbuser:$dbpass@$dbhost:$dbport/$dbname" unless $db_obj;

my $interpro_parser = XrefParser::InterproParser->new($db_obj, $debug);
die "Cannot create XrefParser::InterproParser for the database" unless $interpro_parser;

my $dbi = $db_obj->dbi;

my $source_id = &get_source_id( $dbi );
my $species_id = &get_species_id($dbi );

print "source_id for $SRC is $source_id\nspecies_id is $species_id\n" if $debug;

for my $xmlfile (@ARGV){

	print "Process $xmlfile";
	
	$interpro_parser->run(
		{
		'source_id'  => $source_id,
		'species_id' => $species_id,
		'files'	     => [$xmlfile],
		'rel_file'   => "Dummy",
		'verbose'    => $debug,
		});
	

}

sub get_source_id {

	my $dbi=shift;
	my $sql = "select external_db_id from external_db where db_name = ?";

	my $sth = $dbi->prepare($sql) or return;

	print "get sth for $sql\n" if $debug;
	$sth->execute('Interpro');

	my ($src_id) = $sth->fetchrow_array();

	print "get src_id=$src_id\n" if $debug;	
	return $src_id;
}

sub get_species_id {

        my $dbi=shift;
        my $sql = "select meta_value from meta where meta_key = ?";

        my $sth = $dbi->prepare($sql) or return;

        $sth->execute($SPECIES_ID_METa_KEY);

        my ($species_id) = $sth->fetchrow_array();
        
        return $species_id;
}


__END__


# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $dba = $multi->get_DBAdaptor("core");

#
# Test get_ProteinFeatureAdaptor works
#
my $pfa = $dba->get_ProteinFeatureAdaptor();

ok($pfa && ref($pfa) && $pfa->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor'));

my $pfs = $pfa->fetch_all_by_translation_id(21724);

print_features($pfs);

ok(@$pfs == 15);

sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
	if (defined($f)) {
	  debug($f->start . '-' . $f->end . ' -> ' . $f->hseqname . ':' . $f->hstart . '-' . $f->hend);
	}
  }
}

# test adding and retrieving a new feature
my $start  = 10;
my $end    = 100;
my $hstart = 1;
my $hend   = 90;
my $hstrand = 1;
my $hseqname = 'RF1231';
my $percent_id = 90.8;
my $p_value = '1.52';
my $score   = 50;
my $species = 'Homo_sapiens';
my $hspecies = 'Mus_musculus';
my $hdes = "Hit description";

my $idesc = 'interpro description';
my $interpro_ac = 'interpro accession';

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test');
$multi->save('core', 'protein_feature', 'meta_coord');


my $f = Bio::EnsEMBL::ProteinFeature->new
  (-START       => $start,
   -END         => $end,
   -ANALYSIS    => $analysis,
   -HSTART      => $hstart,
   -HEND        => $hend,
   -HSEQNAME    => $hseqname,
   -PERCENT_ID  => $percent_id,
   -P_VALUE     => $p_value,
   -SCORE       => $score,
   -SPECIES     => $species,
   -HSPECIES    => $hspecies,
   -HDESCRIPTION=> $hdes,
   -IDESC       => $idesc,
   -INTERPRO_AC => $interpro_ac);
   
$pfa->store($f,21724);

my $pfs = $pfa->fetch_all_by_translation_id(21724);

ok(@$pfs == 16);

my @pfs = grep{$_->hdescription() eq $hdes} @$pfs;

ok(scalar @pfs > 0);

$multi->restore('core', 'protein_feature');

done_testing();
