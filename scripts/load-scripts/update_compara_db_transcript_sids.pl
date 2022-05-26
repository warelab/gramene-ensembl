#!/usr/local/bin/perl -w 

=pod

=head1 NAME

update_compara_db_transcript_sids.pl - Parses the transcript ID mapping file, update the stable_id in seq_member table with the new name  
 

=head1 SYNOPSIS

  load_genes_withMakerAED_gff3.pl [options] gff_file

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

  Parses the transcript ID mapping file, update the stable_id with the new name 

  example of the stable_id_mapping file:
  new_name		old_name
  ------------------------------------
  SORBI_3001G000100.1     EER90453.2
  SORBI_3001G000200.1     EER93047.1
  SORBI_3001G000300.2     KXG37057.2

B<The Ensembl Registry>

  The database connection details for both Ensembl and interproscan
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'compara', 
    '-group'   => 'core', 
    '-dbname'  => 'ensembl_compara_1_87_sorghum0422', );
  ---

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
our $sql_update_sid = qq(
                update seq_member set stable_id=? where stable_id=?
        );
my $date = time(); 

  my $help=0;
  my $man=0;
  my( $species, $reg, $xref_src, $no_insert, $non_coding, $biotype, $ignore, $verbose );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "ensembl_registry=s" => \$reg,
        "no_insert"          => \$no_insert,
	"verbose"	     => \$verbose,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args

  my $gff_file = $ARGV[0] || pod2usage("\nNeed the path to a gff file\n");
  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  
  map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
  } $reg, $gff_file;



  # Put stuff in the database?
  $I= $no_insert ? 0 : 1; 
 
  $V = $verbose ? 1 : 0;

  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara', 'compara' );
  $ENS_DBA || ( print( "No compara DB set in $reg\n" ) &&
                pod2usage() );


  # Create a GFF stream
  $F_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


  print( "Updating transcripts for $species\n" );
  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
  }


  $dbh = $ENS_DBA->dbc;
  $STH = $dbh->prepare($sql_update_sid) or die $dbh->errstr;



my $found_cnt=0;
my $missing_cnt=0;
my $updated_cnt=0;

while( my $line = $F_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /^\s*\#/ ); 
  next if ( $line =~ /^\s+$/ );
  chomp $line;

  my( $new_sid, $old_sid ) = split( /\s+/, $line);
  $old_sid =~ s/\..+//;
  
  print "new_sid=$new_sid, old_sid=$old_sid \n" if $V;

  if ( $I ){
	$STH->execute( $new_sid, $old_sid);
	
  	$updated_cnt++;
  }
}

$STH->finish;

print "Updated $updated_cnt\n";

 
