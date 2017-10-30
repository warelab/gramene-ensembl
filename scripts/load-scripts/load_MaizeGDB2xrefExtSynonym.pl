#!/usr/local/bin/perl -w 

=pod

=head1 NAME

load_MaizeGDB2xref.pl - parse gene synonyms and paper reference out from MaizeGDB web data and load them into EnsemblDB
 

=head1 SYNOPSIS

  load_MaizeGDB2xref.pl [options]

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  --species          Species key in Ensembl registry file.
  --synonym_file

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<--species>
  Use this species entry from the registry file [REQUIRED].

B<--synonym_file>
    JSON file with MaizeGDB synonyms for genes

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
use List::MoreUtils qw( uniq );

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;

use JSON;
use File::Slurp;

use lib '/usr/local/ensembl-90/ensembl/modules';
use lib '/usr/local/ensembl-90/ensembl-compara/modules';

Readonly my $SOURCE => 'MaizeGDB';
Readonly my $INFO_TYPE => 'NONE';
Readonly my $VERSION => '1';

my( $SYNOMYM_FILE, $PAPER_FILE, $ENS_DBA, $ANALYSIS );
my $help=0;
my $man=0;
my( $species, $reg);
GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
	"synonym_file=s"     => \$SYNOMYM_FILE,
        )
    or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

  # Process args

$species    || pod2usage("\nNeed a --species\n");
  
map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
} $reg, $SYNOMYM_FILE;


  # Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $reg );
$ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) &&
	      pod2usage() );


print( "Loading genes for $species\n" );
my $pre_text = "  Target DB: ";
foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
	   ":". $dba->dbc->port ."\n");
}


my $gene_adaptor = $ENS_DBA->get_adaptor('Gene');
my $dbentry_adaptor = $ENS_DBA->get_adaptor('DBEntry');

my $json = JSON->new;
#
#The hash keys are:

## B73V3_GM   => external_synonym for this xref
## B73V4_GM   => xref.dbprimary_acc
## DESCRIPTION => xref.description
## FULL_NAME   => xref.display_label
## MGDB_ID     
## NAME        => external_synonym
## REFERENCES  
## SYNONYMS    => external_synonym


my $geneinfo_string = read_file ($SYNOMYM_FILE);
my $gene_xref; 
my $geneinfo_json_obj = decode_json( $geneinfo_string );
my $total_genes = scalar @{$geneinfo_json_obj};
print "Total genes is $total_genes\n";

#exit;

foreach my $a_geneinfo( @$geneinfo_json_obj ){
#    print Dumper ($a_geneinfo->{SYNONYMS});
#    last;
   # map { print "$_\n" }keys %$_;

   
    my $gene_stable_id = $a_geneinfo->{B73V4_GM};
    my $display_label  = $a_geneinfo->{FULL_NAME} || $a_geneinfo->{NAME};
    my $description    = join ';', @{$a_geneinfo->{DESCRIPTION}};
    my $synonyms       = [uniq (@{$a_geneinfo->{SYNONYMS}}, $a_geneinfo->{B73V3_GM}, $a_geneinfo->{NAME})];
    
    unless( $gene_stable_id){ 
	print STDERR join "\t", ($SOURCE, 'missingGeneStableID', $display_label,$INFO_TYPE, @$synonyms, "$description\n" );
	next;
    }

    my $eGene          = $gene_adaptor->fetch_by_stable_id($gene_stable_id);    
    
    unless ($eGene){
	print STDERR join "\t", ($SOURCE, $gene_stable_id, $display_label,$INFO_TYPE, @$synonyms, $description, "NotFound\n" );
	next;
    }

    print join "\t", ($SOURCE, $gene_stable_id, $display_label,$description, $INFO_TYPE, @$synonyms ), "\n";

    my $DBxrefs = $dbentry_adaptor->fetch_all_by_Gene( $eGene);
   
    foreach my $a_xref( @{$DBxrefs} ){
	next unless ( uc $a_xref->dbname eq uc $SOURCE);
    
	$a_xref->adaptor($dbentry_adaptor);
	$a_xref->primary_id( $gene_stable_id );
	$a_xref->display_id( $display_label);
	$a_xref->description( $description);
	$a_xref->version($VERSION);
	$a_xref->info_type( $INFO_TYPE);

	$dbentry_adaptor->update( $a_xref );
	my $xref_id = $a_xref->dbID;
	my $r = add_synonyms( $dbentry_adaptor, $xref_id, $synonyms );
	print STDERR "Fail to add_synonyms for $gene_stable_id\n" unless $r;

	$eGene->display_xref( $a_xref );
	$gene_adaptor->update( $eGene );
	
	print STDERR "Succeeded updating $gene_stable_id\n";
	last;
    }

}

sub add_synonyms{

    my  ($dbentry_adaptor, $xref_id, $synonyms) = @_;

    my $mysql_select_stmt = "select count(*) from external_synonym where xref_id=$xref_id and synonym = ?";
    my $mysql_insert_stmt = "insert into external_synonym values ($xref_id, ?)";

    my $select_stm_handle = $dbentry_adaptor->prepare($mysql_select_stmt);
    my $insert_stm_handle = $dbentry_adaptor->prepare($mysql_insert_stmt);
    
    return undef unless $select_stm_handle && $insert_stm_handle;

    for my $a_syn( @{$synonyms}){
	next if $a_syn =~ /^\s*$/;
	$select_stm_handle->execute($a_syn);
	my $count = $select_stm_handle->fetchrow_array;
	print "found $count row for $xref_id, $a_syn\n";
	next if $count > 0;
	$insert_stm_handle->execute($a_syn);
    }
    $select_stm_handle->finish;
    $insert_stm_handle->finish;
    return 1;
}
