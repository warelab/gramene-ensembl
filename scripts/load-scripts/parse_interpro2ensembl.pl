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


use lib map {"$ENS_CODE_ROOT/$_" } 
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


