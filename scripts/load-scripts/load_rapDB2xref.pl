#!/usr/local/bin/perl -w 

=pod

=head1 NAME

load_rapDB2xref.pl - parse gene synonyms and paper reference out from MaizeGDB web data and load them into EnsemblDB
 

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
    GFF3 file with rapGDB synonyms for genes

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
use URI::Escape;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;

use JSON;
use File::Slurp;

use lib '/usr/local/ensembl-90/ensembl/modules';
use lib '/usr/local/ensembl-90/ensembl-compara/modules';

Readonly my @SOURCES => qw ( CGSNL Oryzabase RAP-DB )  ;
Readonly my $INFO_TYPE => 'DIRECT';
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

my $DBsrcs = join '|', sort { $b cmp $a} @SOURCES;
print "regex is $DBsrcs\n"; 

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
#ID=Os01g0100900;Name=Os01g0100900;Note=Sphingosine-1-phosphate lyase%2C Disease resistance response (Os01t0100900-01);Transcript variants=Os01t0100900-01;RAP-DB Gene Symbol Synonym(s)=SPL1%2C OsSPL1;RAP-DB Gene Name Synonym(s)=SPHINGOSINE-1-PHOSPHATE LYASE 1%2C Sphingosine-1-Phoshpate Lyase 1;Oryzabase Gene Symbol Synonym(s)=OsSPL%2C OsSPL1;Oryzabase Gene Name Synonym(s)=Sphingosine-1-phosphate lyase%2C sphingosine-1-phosphate lyase 1
#CGSNL Gene Name
#CGSNL Gene Symbol
#ID
#Name
#Note
#Oryzabase Gene Name Synonym(s)
#Oryzabase Gene Symbol Synonym(s)
#RAP-DB Gene Name Synonym(s)
#RAP-DB Gene Symbol Synonym(s)
#chr01   irgsp1_locus    gene    89763   91465   .       -       .       ID=Os01g0102000;Name=Os01g0102000;Note=Phosphoesterase family protein. (Os01t0102000-01);Transcript variants=Os01t0102000-01;CGSNL Gene Symbol=NPC5;CGSNL Gene Name=NON-SPECIFIC PHOSPHOLIPASE C5;Oryzabase Gene Symbol Synonym(s)=OsNPC6%2C OsNPC5%2C NPC6;Oryzabase Gene Name Synonym(s)=Non-specific phospholipase C5

open my $fh, $SYNOMYM_FILE or die "Cannot open $SYNOMYM_FILE to read";
my @geneinfo_objs = 
	map{ my $extra=$_; my %obj; map { my @t = split '='; $t[0] =~ s/Synonym.*//; $t[0] =~ s/^ +//; $t[0] =~ s/ +$//; $t[0] =~ s/\s+/_/g; $obj{ $t[0] } = uri_unescape( $t[1] ); } (split ';', $extra); \%obj; }
	map{chomp; my @parts = split /\t/; $parts[8];} <$fh>;
close $fh;

for my $agene_obj( @geneinfo_objs ){
	
#	map{	print "$_=", $agene_obj->{$_},"\n"} keys %$agene_obj; 
}


#exit;

my $gene_xref; 
my $total_genes = scalar @geneinfo_objs;
print "Total genes is $total_genes";

#exit;

foreach my $a_geneinfo( @geneinfo_objs ){
#    print Dumper ($a_geneinfo->{SYNONYMS});
#    last;
    #map { print $_,$a_geneinfo->{$_},"\n" } keys %$a_geneinfo;

   
    my $gene_stable_id = $a_geneinfo->{ID};
    my $display_label  = $a_geneinfo->{Name};
    my $description    = $a_geneinfo->{Note};
    my %synonyms ;
   
    map{ my $k=$_."_Gene_Symbol"; my $v=$_."_Gene_Name"; $synonyms{$_} = [$a_geneinfo->{$k}, $a_geneinfo->{$v}] } 
				grep{ my $k=$_."_Gene_Symbol";  $a_geneinfo->{$k} } @SOURCES;
    
    next unless (scalar keys %synonyms >= 1);

    unless( $gene_stable_id){ 
	print STDERR join "\t", ( 'missingGeneStableID', $display_label, $INFO_TYPE, keys %synonyms, values %synonyms, "$description\n" );
	next;
    }

    my $eGene          = $gene_adaptor->fetch_by_stable_id($gene_stable_id);    
    
    unless ($eGene){
	print STDERR join "\t", ($gene_stable_id, $display_label,$INFO_TYPE, keys %synonyms, values %synonyms, $description, "NotFound\n" );
	next;
    }

    print join "\t", ( $gene_stable_id, $display_label, $description, $INFO_TYPE, keys %synonyms, map{ @{$_}} values %synonyms ), "\n";

#last;

  # Add XREFs to gene
  # The xref primary id cannot be null

    my $gene_xref;
    my $display = 0;
   
    for my $dbn(keys %synonyms){
	my ($synId, $synName) = @{$synonyms{$dbn}}; 

	push @$gene_xref, {
		    db                   => $dbn,
        	    primary_acc          => $synId,
                    display_label        => $synName,
                    description          => $description   || '',
		    info_type            => $INFO_TYPE,
			   };


    }#end of for loop to build xref array for eGene

	
    my $gid = $eGene->dbID;
    foreach my $dbentry( make_dbentries($gene_xref) ){
	$dbentry_adaptor->store( $dbentry, $gid, 'Gene'  );
    	$eGene->add_DBEntry( $dbentry );
      	$eGene->display_xref( $dbentry ) unless $display++;	
    }

    $gene_adaptor->update( $eGene );
    print STDERR "Succeeded updating $gene_stable_id\n";

}

sub make_dbentries{

  my $xref  = shift;
  my @data  = @{ $xref || [] };
  my @xrefs = ();

  foreach my $entry( @data ){
    push @xrefs, Bio::EnsEMBL::DBEntry->new
        (
	 -adaptor     => $dbentry_adaptor,
         -dbname      => $entry->{db} || '',
         -primary_id  => $entry->{primary_acc} || 'noid',
         -display_id  => $entry->{display_label} || '', # Filed required; Ens v28
         -description => $entry->{description} || '',
         -info_type   => $entry->{info_type},
         -version     => 1,
         -release     => 1,
         );
  }
  return @xrefs;


}

