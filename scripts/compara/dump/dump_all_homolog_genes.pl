#!/usr/local/bin/perl -w

# Pod is used to construct a guide of the correct use of the program

=pod

=head1 NAME

perl dump_homolog_gene.pl - Dumps the homolog genes

=head1 SYNOPSIS

  perl dump_homolog_gene.pl [options] > resultfile.dat

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -s|--ortho_species      Homolog Species bionomial name.

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry				
																# ENSEMBLHOME refers to the dir in which the script is run.
B<-s|--species>
  Use this species entry from the registry file [REQUIRED].

=head1 DESCRIPTION

Script that iterates through each gene in the given species and dumps
the subroot ID of any protein tree alongside the gene genomic
coordinates

B<The Ensembl Registry>

  The database connection details for both Ensembl species and compara   ##
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'Oryza_sativa', 
    '-group'   => 'core', 
    '-dbname'  => 'Oryza_sativa_japonica_core_48_28', );
  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
  ( '-species' => 'compara',
    '-group'   => 'compara',
    '-dbname'  => 'ensembl_compara_48_28', );
  ---

TODO: Complete this section
                   
Maintained by Will Spooner <whs@ebi.ac.uk>

=cut

####### Perl code starts here

BEGIN {
  $ENV{'GrameneEnsemblDir'} ||= '/usr/local/ensembl-56';
}

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use Data::Dumper qw(Dumper); #

use File::Basename;		# basename takes a long UNIX path name and returns the file name at the end.

use FindBin qw( $Bin ); # Locates the full path to the script bin directory to allow the use of paths relative to the bin directory.
						# allow the use of modules in the lib directory without knowing where the software tree is installed.
						# Maybe try  FindBin::again(); to make sure that it finds the pathway.
						
use Pod::Usage;			# Module that prints parts of the manual described above.

use Getopt::Long;		# Module used to read the options. options are written as follow: --name_opt=value
						# this module takes the options first and then the arguments. To specify when the options end, we need to write
						# a double dash --. Example: --opt1=value --opt2=value -- --arg1=value
						
use IO::File;			#Module used to create a file

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor;  # added by Tzitziki

our $ENS_DBA;
our $CMP_DBA;
our $ortho_species_id;
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $ortho_species,$reg, $logic_name, $no_insert, $homolog_type );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species, # s stands for string, i for numeric data.
	"ortho_species=s"    => \$ortho_species,
        "ensembl_registry=s" => \$reg,
        "homolog_type=s" => \$homolog_type,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $reg        ||= $ENV{'GrameneEnsemblDir'}.'/conf/ensembl.registry';
  map{  					#Perl's map() function runs an expression on each element of an array, and returns a new array with the results.
    -e $_ || pod2usage( "\n[*DIE]File $_ does not exist\n" );
    -r $_ || pod2usage( "\n[*DIE]Cannot read $_\n" );
    -f $_ || pod2usage( "\n[*DIE]File $_ is not plain-text\n" );
    -s $_ || pod2usage( "\n[*DIE]File $_ is empty\n" );
  } $reg;
  
  unless( $species ){
    $species    || pod2usage("\n[*DIE] Need a --species\n");
  }

  # restrict to a homolog type
  $homolog_type ||= '';
  $homolog_type = qr($homolog_type);

  my %species_to_taxon=(
    'arabidopsis_lyrata'=>59689,
    'arabidopsis_thaliana'=>3702,
    'oryza_brachyantha'=>4533,
    'oryza_brachyantha_3s'=>4533,
    'oryza_minutabb'=>63629,
    'oryza_minutacc' =>63629,
    'oryza_sativa'=>39947,
    'oryza_sativa_japonica'=>39947,
    'oryza_sativa_indica'=>39946,
    'sorghum_bicolor'=>4558,
    'populus_trichocarpa'=>3694,
    'vitis_vinifera'=>29760,
    'zea_mays'=>4577,
  );



  unless( $ortho_species ){
    pod2usage("\n[*DIE] Need the homologue species\n");
  }
  $ortho_species_id = $species_to_taxon{lc($ortho_species)};
  unless( $ortho_species_id ){
    pod2usage("\n[*DIE] Homologue species taxonomy id not found\n");
  }


  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' )  # conects to Ensembl, and to the  'core' database of the selected specie.
      || pod2usage("\n[*DIE]No core DB for $species set in $reg\n" );
  $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara','compara')
      || pod2usage("\n[*DIE]No compara DB set in $reg\n" );
}

# There are two ways to connect to the Ensembl Compara database. The old way uses the Bio::EnsEMBL::Compara::DBSQL::DBAdaptor explicitly.
# The new one uses the Bio::EnsEMBL::Registry module, which can read either a global or a specific configuration file or auto-configure itself.

MAIN:{
	
  # print the header
  print join( "\t",
              'gene_stable_id',
              'sequence_name',
              'gene1_start',
              'gene1_end',
              'specie_id_1',
              'homolog2_member_id',
              'sequence_name2',
              'gene2_start',
              'gene2_end',
              'homology_relationship',
              'specie_id_2',
              'root_node_id',
              'root_taxon_name',
              'interpro_union',
              'interpro_intersection',
              'shared_interpro_id',
              'shared_interpro_name',
            ) . "\n";

  # Create the ADAPTORS that we will be using to fetch data from the database
  my $gene_adaptor = $ENS_DBA->get_adaptor('Gene')
      || die( "[*DIE] Cannot ENS_DBA->get_adaptor('Gene')" );
  my $member_adaptor = $CMP_DBA->get_adaptor('Member')				# A member can be a gene or a protein. Most of them are in the core databases of their
      || die( "[*DIE] Cannot CMP_DBA->get_adaptor('Member')" );		# correspondient specie.
      
  
  my $HomologyAdaptor=  $CMP_DBA->get_adaptor('Homology') 
      || die ("[*DIE] Cannot CMP_DBA->get_adaptor('Homology')"); # added by Tzitziki 
  my $TreeAdaptor = $CMP_DBA->get_adaptor('ProteinTree') 
      || die ("[*DIE] Cannot CMP_DBA->get_adaptor('ProteinTree')");

  # Debug
  my $db = $ENS_DBA->dbc->dbname;
  warn( "[INFO] Collecting all genes from $db \n" );

  # Get a list of all gene IDs, and loop through them
  my $gene_id_arrayref = $gene_adaptor->list_stable_ids;
  my $total_genes = scalar @$gene_id_arrayref;

  warn( "[INFO] Collected genes $total_genes; processing...\n" );
  

  # Some counters
  my( $num_genes, $num_coding, $num_chr_genes,$chr_specie_genes,$ignored_specie) = (0,0,0,0,0);
  my $ignored=0;
  
  # loop through the genes
  
#  if( my $id = 'GRMZM2G018837' ){
  while( my $gene_id = shift @$gene_id_arrayref ){
    $num_genes ++; #COUNT

    # Get the gene object, and skip unless protein coding
    my $gene = $gene_adaptor->fetch_by_stable_id( $gene_id )
        || die( "[*DIE] Cannot fetch_by_stable_id gene $gene_id" );
    next unless $gene->biotype eq 'protein_coding';
    $num_coding ++; # COUNT
    #warn $gene_id;
   
    ######### TEST CODE
    #liya : comment this line to get all genes instead of chromosome genes only
    #unless($gene->seq_region_name =~ /^\d/){$ignored++; next;} # MANUAL MODIFICATION -turn this two "ifs" off when analizing maize
    
    $num_chr_genes++; #COUNT

    my $member 
        = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$gene_id);
    unless( $member ){
      warn( "[WARN] No ENSEMBLGENE member for $gene_id" );
      next;
    }
   
    my $species_id = $member->taxon_id;
    $gene->project('toplevel');
 
    # List all InterPro domains on this gene
    my %interpros;
    foreach my $translation( map{$_->translation||()} 
                             @{$gene->get_all_Transcripts} ){
      foreach my $pf( @{$translation->get_all_ProteinFeatures} ){
        $interpros{$pf->interpro_ac} = $pf->idesc;
      }
    }

    #if( my $member 
    #      = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$gene_id) ){
  
    #Which protein tree is the member in? Get list of genes in other species.
    my $root_node_id = '';
    my $root_taxon = '';
    my @ospecies_genes = ();
    if( my $gene_tree = $TreeAdaptor->fetch_by_Member_root_id($member, 0) ){
      my $root_node   = $gene_tree->subroot;
      $root_node_id = $root_node->node_id;
      $root_taxon   = $root_node->get_tagvalue('taxon_name') || '';
      @ospecies_genes 
          = ( map{ $_->get_Gene }
              grep{ $_->genome_db->taxon_id() == $ortho_species_id } 
              @{$root_node->get_all_leaves} );
      $gene_tree->release_tree;
    }
    @ospecies_genes || next;

    #Obtain the homology descriptions of the gene	
    my %homology_descriptions;
    foreach my $homology( @{$HomologyAdaptor->fetch_all_by_Member($member)} ){
      my( $ogene )
          = ( grep{ $_->taxon_id == $ortho_species_id } 
              @{$homology->gene_list} ); 
      $ogene || next;
      $homology_descriptions{$ogene->stable_id} = $homology->description;
    }

    my $homologies = $HomologyAdaptor->fetch_all_by_Member($member);
    foreach my $ogene( @ospecies_genes ){

      #Check for InterPro domain concordance
      my %ointerpros;
      my %shared_interpros;
      foreach my $translation( map{$_->translation||()}
                               @{$ogene->get_all_Transcripts} ){
        foreach my $pf( @{$translation->get_all_ProteinFeatures} ){
          my $interpro_ac = $pf->interpro_ac;
          $ointerpros{$interpro_ac} = $pf->idesc;
          if( $interpros{$interpro_ac} ){$shared_interpros{$interpro_ac}++}
        }
      }
      my %all_interpros = map{$_=>1} keys %interpros, keys %ointerpros; 
      my ( $shared_interpro ) = keys( %shared_interpros );

      # Check relationship involves orthology
      my $ogene_id = $ogene->stable_id;
      my $homology_type = $homology_descriptions{$ogene_id}
        || 'between_species_paralog';

      #next unless $homology_type =~ $homology_type;
      $ogene->project('toplevel');  
      $chr_specie_genes++;
      
      print join( "\t",
                  $gene_id,
                  $gene->seq_region_name,
                  $gene->seq_region_start,
                  $gene->seq_region_end,
                  $species_id,
                  $ogene_id,
                  $ogene->seq_region_name,
                  $ogene->seq_region_start,
                  $ogene->seq_region_end,
                  $homology_type,
                  $ortho_species_id,
                  $root_node_id,
                  $root_taxon,
                  scalar( keys( %all_interpros ) )    || 0,
                  scalar( keys( %shared_interpros ) ) || 0,
                  $shared_interpro || '',
                  $interpros{$shared_interpro||''} || '',
                  ), "\n";
    }
    ##########
  }

  
  # Update progress each 1000 genes
  if( $num_genes/1000-int($num_genes/1000) == 0 ){
    warn ("[INFO] processed $num_genes of $total_genes\n");
  }

# Debug
warn sprintf( "[INFO] processed %d coding genes of %d total genes; chr genes %d, ignored genes %d, ignored species %d, chr_specie genes %d \n",
              $num_coding, $num_genes,$num_chr_genes, $ignored, $ignored_specie, $chr_specie_genes);
print "coding genes:$num_coding,num_genes:$num_genes, num_tree_genes:$num_chr_genes,ignored:$ignored\n";

}

exit;
