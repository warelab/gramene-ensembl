#!/usr/local/bin/perl -w 

=pod

=head1 NAME

load_genes_withMakerAED_gff3.pl - Parses the transcript ID mapping file, update the stable_id with the new name  
  while save the old name as xref names. 

=head1 SYNOPSIS

  load_genes_withMakerAED_gff3.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -x|--xref_src		the source of the old names
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

B<-s|--species>
  Use this species entry from the registry file [REQUIRED].

B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

  Parses the transcript ID mapping file, update the stable_id with the new name 
  while save the old name as xref names. 

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
  ( '-species' => 'maize', 
    '-group'   => 'core', 
    '-dbname'  => 'zea_mays_core_30_bac20', );
  ---

Maintained by Sharon Wei <weix@cshl.edu>

=cut

#delete from gene; delete from transcript; delete from translation; delete from exon; delete from exon_transcript; delete from gene_stable_id; delete from transcript_stable_id; delete from translation_stable_id; delete from exon_stable_id;

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


use vars qw( $ENS_DBA $dbh $I $INSERT $F_HANDLE $XREF_SRC $IGNORE $V );

my $date = time(); 

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $xref_src, $no_insert, $non_coding, $biotype, $ignore, $verbose );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "xref_src=s"       => \$xref_src,
        "no_insert"          => \$no_insert,
	"verbose"	     => \$verbose,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args

  my $gff_file = $ARGV[0] || pod2usage("\nNeed the path to a gff file\n");
  $species    || pod2usage("\nNeed a --species\n");
  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  $xref_src   ||= pod2usage("\nNeed a --xref source for old name to be associated with\n");;
  
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
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) &&
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
  my $rv = $dbh->do("select * from external_db where db_name='$xref_src'") || die $dbh->errstr;
  
  $rv == 1 ? $XREF_SRC = $xref_src : die "The external_db $xref_src does not exist"; 

}


our $ta = $ENS_DBA->get_adaptor('Transcript');
our $DBEA = $ENS_DBA->get_adaptor('DBEntry');
our $sql_update_sid = qq(
		update transcript set stable_id=? where transcript_id=?
	);
our $sth_update = $dbh->prepare($sql_update_sid) or die $dbh->errstr; 

my $found_cnt=0;
my $missing_cnt=0;
my $updated_cnt=0;

while( my $line = $F_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /^\s*\#/ ); 
  next if ( $line =~ /^\s+$/ );
  chomp $line;

  #warn(Dumper($TRPT2GENE));

  # Split gff line,\
  # start always <= end even for - strand genes
  my( $new_sid, $old_sid ) = split( /\s+/, $line);
  $old_sid =~ s/\..+//;
  
  print "new_sid=$new_sid, old_sid=$old_sid \n" if $V;
  my $eTranscript = $ta->fetch_by_stable_id( $old_sid );
  unless( $eTranscript  ){  
	warn "MISSING: No Transcript found for $old_sid ($new_sid)\n";
	$missing_cnt++;
	next;
  }
  $found_cnt++;
  $eTranscript->stable_id( $new_sid );
  add_xrefs( $eTranscript, $old_sid, $XREF_SRC);

  if ( $I ){
	warn("DEBUG insert: the new id is ", $eTranscript->stable_id, "\n");
	#$ta->update( $eTranscript ) ; #not working as expected
	$sth_update->execute( $new_sid, $eTranscript->dbID);
	
  	$updated_cnt++;
  }
}

$sth_update->finish;

print "Missing $missing_cnt\nFound $found_cnt\nUpdated $updated_cnt\n";

sub add_xrefs{ # 1.ensembl_object
               # 2.old_sid
               # 3.external_db (gene or transcript)
  # Add XREFs to gene/transcript
  # The xref primary id cannot be null
  # for Rice make gene TU_feat_name as primary id and pub_locus as display_label

  my $eObj    = shift;
  my $old_sid = shift;
  my $external_db    = shift;


  my $xref = {
               db                   => $external_db,
               primary_acc          => $old_sid,
               display_label        => $old_sid,
               description          => "old stable_id of the transcript",
              };
  


  my $dbentry = make_dbentries($xref);
  my $transcript_id = $eObj->dbID;

  $DBEA->store($dbentry, $transcript_id, 'Transcript') if $I;

  #$eObj->add_DBEntry( $dbentry );
    #print "dbentry->dbname=", $dbentry->dbname, "\n";


}



#----------------------------------------------------------------------
# Ensembl DBEntry objects are used for gene/transcript etc xrefs
#
# There are only cDNA evidence xref in the tigrv4 xml
#
# create entry in xref table
# Before making it, need to check whether the same entry already exist,
# if does, use the existing one instead of making new one

sub make_dbentries{

  my $entry  = shift;

  my $xref_obj =  Bio::EnsEMBL::DBEntry->new
        (
         -dbname      => $entry->{db} || '',
         -primary_id  => $entry->{primary_acc} || 'noid',
         -display_id  => $entry->{display_label} || '', # Filed required; Ens v28
         -description => $entry->{description} || '',
         -version     => 1,
         -release     => 1,
         );
  
  return $xref_obj;


}
 
