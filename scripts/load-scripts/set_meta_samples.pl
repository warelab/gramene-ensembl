#!/usr/local/bin/perl -w

# POD documentation - main docs before the code

=pod

=head1 NAME

  set_meta_samples.pl

=head1 SYNOPSIS

 set the samples in the meta table of the core database

=head1 DESCRIPTION

 The script find out of there are samples in the meta table, if not use the chr1 and 
 the first gene/transcript in the table as sample.
 To not update analyses in the database you need to pass the -noupdate option.

=head1 OPTIONS

     Database options

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname      For RDBs, what name to connect to (dbname= in locator)
    -dbuser      For RDBs, what username to connect as (dbuser= in locator)
    -dbpass      For RDBs, what password to use (dbpass= in locator)
                 analysis.descriptions in this directory can be used and is 
                 an example of the format. Multiple -file args can be specified
    -pattern     check databases matching this PATTERN
                 Note that this is a database pattern of the form %core_53_%
    -help print out documentation

=head1 EXAMPLES

 perl set_meta_samples.pl -dbhost my_host -dbuser user -dbpass ***** 
 -dbname my_db 

if you want to update all databases for a type

perl load_analysis_descriptions.pl -dbhost my_host -dbuser user -dbpass ***** 
 -pattern '%_55_%' 

syntax errors found as the definition file is parsed (usually the use of spaces rather
than tabs), cause the script to exit - fix each one and rerun.

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use DBI;

use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

$| = 1;

my ($dsn,$dbh);

my $dbhost = '';
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname = '';
my @files = ();
my $noupdate;
my $help = 0;
my $pattern;


&GetOptions (
  'host|dbhost=s'       => \$dbhost,
  'dbname=s'            => \$dbname,
  'user|dbuser=s'       => \$dbuser,
  'pass|dbpass=s'       => \$dbpass,
  'port|dbport=s'       => \$dbport,
  'file|descriptions=s' => \@files,
  'pattern=s'           => \$pattern,
  'h|help!'             => \$help
);

if (!$dbhost){
  print ("Need to pass a dbhost\n");
  $help =1;
}
if (!$dbname and !$pattern){
  $help =1;
  throw("Need to enter either a database name in -dbname or a pattern in -pattern\n");
}

if ($dbname and $pattern) {
  $help =1;
  throw("You should only either enter a database name using -dbname OR a pattern using -pattern but not both.\n");
  }


if($help){
  usage();
}

#connect to database
$dsn = "DBI:mysql:host=" . $dbhost . ";port=" . $dbport;

eval{
  $dbh = DBI->connect($dsn, $dbuser, $dbpass, 
		      {'RaiseError' => 1,
		       'PrintError' => 0});
};

if( $@ ){
	print "failed to create dbh, $@\n";
}

# get all database names that match pattern
my ($sth, $sql);
my $sql_pattern = $pattern || $dbname;
$sql = "SHOW DATABASES LIKE '". $sql_pattern ."'";
warn("$sql");
$sth = $dbh->prepare($sql);
$sth->execute;
while (my ($dbname) = $sth->fetchrow_array){
  if ($pattern) {   # only check dbname if a DB pattern has been specified.
    next unless $dbname =~ /core|cdna|otherfeatures/;
    next if $dbname =~ /coreexpression/;
  }
  print "\n\nLooking at ... $dbname\n";
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $dbhost,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $dbpass,
    -port   => $dbport
  );

  my $dbc = $db->dbc;

  # Pre-fetch samples in meta table

  my @essentail_sample_keys = qw( sample.location_param sample.location_text
				  sample.gene_param sample.gene_text
				  sample.transcript_param sample.transcript_text
				  sample.search_text );
  my %meta_smaples;
  my $sth = $dbc->prepare("SELECT meta_key, meta_value FROM meta where meta_key like '%sample%'");

  $sth->execute();

  while ( my $result = $sth->fetchrow_arrayref ){
	$meta_smaples{ lc $result->[0]} = $result->[1]; 
  }

  $sth->finish;

  my @missing_sample_keys;
  for my $k( @essentail_sample_keys ){
	next if $meta_smaples{$k};
	push @missing_sample_keys, $k; 
  }

  my $loc_sql = "select name from seq_region sr join seq_region_attrib sra using(seq_region_id) where sra.attrib_type_id=6 limit 1";	
  my $gene_sql = "select g.stable_id, t.stable_id from gene g join transcript t using(gene_id) where g.biotype = 'protein_coding' limit 1";
  my $sth_loc;
  my $sth_gene ;  

  for my $mk( @missing_sample_keys ){
  	if ($mk =~ /sample\.location/ ){
		next if $meta_smaples{'sample.location_param'};
                $sth_loc = $dbc->prepare($loc_sql);
                $sth_loc->execute();;
                my @names =  $sth_loc->fetchrow_array ;
                $meta_smaples{'sample.location_param'} = "$names[0]:8001-18000";
                $meta_smaples{'sample.location_text'} = "$names[0]:8001-18000";
	} 
	if ($mk =~ /sample\.gene/ ){
		next if $meta_smaples{'sample.gene_param'};
		$sth_gene = $dbc->prepare($gene_sql);
		$sth_gene->execute();;
		my @stable_ids =  $sth_gene->fetchrow_array ;		
                $meta_smaples{'sample.gene_param'} = $stable_ids[0];
                $meta_smaples{'sample.gene_text'} = $stable_ids[0];
		$meta_smaples{'sample.transcript_param'} = $stable_ids[1];
                $meta_smaples{'sample.transcript_text'} = $stable_ids[1];
        }
	if ($mk eq 'sample.search_text' ){
                $meta_smaples{'sample.search_text'} = 'Carboxypeptidase';
        }

  }

  $sth_loc->finish;
  $sth_gene->finish;

  my $sql_insert = "insert into meta(meta_key, meta_value) values(?,?)";
  my $sth_insert = $dbc->prepare($sql_insert);

  for my $mk( @missing_sample_keys ){
	$sth_insert->execute($mk, $meta_smaples{$mk});
  }

  $sth_insert->finish;
}

$dbh->disconnect;

sub usage{
  exec('perldoc', $0);
  exit;
}
