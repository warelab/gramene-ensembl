#!/usr/local/bin/perl -w 

=pod

=head1 NAME

get_vardb_transcript_affected.pl - 

select tv.consequence_types, s.name,  feature_stable_id  from transcript_variation tv join variation_feature vf using(variation_feature_id) join source s using(source_id); 

=head1 SYNOPSIS

  get_vardb_transcript_affected.pl -e ensembl.registry 

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -n|--no_insert        Do not make changes to the database. For debug.
  -v|--verbose		be verbose and print more messages

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry


B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

Maintained by Sharon Wei <weix@cshl.edu>

=cut


use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Readonly;
use List::Util qw( first );
use List::MoreUtils;

use vars qw( $BASEDIR );

BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;


use vars qw( $ENS_DBA $dbh $I $INSERT $F_HANDLE  $IGNORE $V $STH );
our $sql = qq(
        	select tv.consequence_types, s.name,  tv.feature_stable_id  
		from transcript_variation tv 
		join variation_feature vf using(variation_feature_id) 
		join source s using(source_id);
        );

  my $help=0;
  my $man=0;
  my( $species, $reg, $no_insert, $ignore, $verbose );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "ensembl_registry=s" => \$reg,
	"species=s"	     => \$species,
        "no_insert"          => \$no_insert,
	"verbose"	     => \$verbose,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args

  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  $species    || pod2usage("\nNeed species\n");
  
  map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
  } $reg;



  # Put stuff in the database?
  $I= $no_insert ? 0 : 1; 
 
  $V = $verbose ? 1 : 0;

  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'variation' );
  $ENS_DBA || ( print( "No var DB set in $reg\n" ) &&
                pod2usage() );



  print( "Query transcripts affected by variation for $species\n" );
  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
  }


  $dbh = $ENS_DBA->dbc;
  $STH = $dbh->prepare($sql) or die $dbh->errstr;



my $found_cnt=0;
my $missing_cnt=0;
my $updated_cnt=0;

$STH->execute();

my %trpt_src_conseq;
my %src_trpt_conseq;
while (my $rowref = $STH->fetchrow_arrayref ){
  
  my( $consequence_types, $src_name,  $tstable_id   ) = @{$rowref}; 
  
  print "transcript=$tstable_id, source=$src_name, consequence_types=$consequence_types \n" if $V;

  $trpt_src_conseq{$tstable_id}{$src_name}{$consequence_types}++;
  $src_trpt_conseq{$src_name}{$tstable_id}++;

  $found_cnt++;

}


print "Fetched $found_cnt transcript effects\n";

my %src2trpt;
my %src_type_trptcnt;
for my $trpt( keys %trpt_src_conseq){

	for my $src(keys %{$trpt_src_conseq{$trpt}}){
		$src2trpt{$src}++;
		for my $cs_type(keys %{$trpt_src_conseq{$trpt}{$src}}){
			$src_type_trptcnt{$src}{$cs_type}++;	
		}
	}
}

print "\nReporting total number of affected transcripts...\n\n";
for my $s(sort keys %src2trpt){
	print "$s\t$src2trpt{$s}\n";
}

print "\n\nReporting total number of affected transcripts by consequence types\n";
for my $s(sort keys %src_type_trptcnt){
	for my $t( sort keys %{$src_type_trptcnt{$s}}){
	        print "$s\t$t\t$src_type_trptcnt{$s}{$t}\n";
	}
}


$STH->finish; 
