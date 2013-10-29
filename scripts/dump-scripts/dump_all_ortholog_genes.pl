#!/usr/local/bin/perl -w

BEGIN {
    $ENV{'EnsemblDir'} ||= '/usr/local/ensembl-live/';
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/';
}

use lib map { $ENV{'EnsemblDir'} . "/$_" }
    qw ( bioperl-live modules ensembl/modules ensembl-external/modules
    ensembl-variation/modules ensembl-draw/modules ensembl-compara/modules );
use lib '/usr/local/gramene/lib/perl';

use strict;
use warnings;
use Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor;
use Bio::EnsEMBL::Registry;
use Data::Dumper qw(Dumper);
use File::Basename qw( fileparse basename );
use File::Path qw( mkpath );
use File::Spec::Functions;
use FindBin qw( $Bin );
use Getopt::Long;
use Gramene::Utils qw( commify );
use Gramene::Marker::DB;
use List::Util qw( max );
use Text::TabularDisplay;
use IO::File;
use IO::Prompt qw( prompt );
use Pod::Usage;
use Readonly;

$| = 1;

my $help          = 0;
my $man           = 0;
my $out_file      = '';
my $markers_file  = '';
my $species       = '';
my $ortho_species = '';
my $reg           = catfile( $ENV{'GrameneDir'}, 'conf', 'ensembl.registry' );

GetOptions(
    'help|?'               => \$help,
    'man'                  => \$man,
    's|species=s'          => \$species,   
    'o|ortho_species=s'    => \$ortho_species,
    'e|ensembl_registry:s' => \$reg,
    'f|out-file:s'         => \$out_file,
    'm|markers-file:s'     => \$markers_file,
) or pod2usage( 2 );

pod2usage( -verbose => 2 ) if $man;

pod2usage( 1 ) if $help;

map { 
    -e $_ || pod2usage( "\n[*DIE]File $_ does not exist\n" );
    -r $_ || pod2usage( "\n[*DIE]Cannot read $_\n" );
    -f $_ || pod2usage( "\n[*DIE]File $_ is not plain-text\n" );
    -s $_ || pod2usage( "\n[*DIE]File $_ is empty\n" );
} $reg;

if ( !$species ) {
    pod2usage('Missing species');
}

if ( !$ortho_species ) {
    pod2usage('Missing orthologue species');
}

if ( $out_file ) {
    if ( -e $out_file && -s _ ) {
        my $overwrite = prompt "$out_file exists.  Overwrite? ", -yn;
        if ( !$overwrite ) {
            print "Not OK, exiting.\n";
            exit 0;
        }
    }

    my ( $file_name, $path, $suffix ) = fileparse( $out_file );

    if ( !-d $path ) {
        mkpath( $path );
    }
} 
else {
    pod2usage('No output file specified');
}

my ( $out_fh, $markers_fh );
if ( $out_file ) {
    open $out_fh, '>', $out_file or die "Can't write '$out_file': $!\n";
}

if ( $markers_file ) {
    open $markers_fh, '>', $markers_file 
        or die "Can't write '$markers_file': $!\n";
}

if ( !$out_file && !$markers_fh ) {
    pod2usage('Must specify out-file and/or markers-file');
}

for ( $species, $ortho_species ) {
    $_ =  lc $_;
    $_ =~ s/ /_/g;
}

( my $markers_species       = $species       ) =~ s/_/ /g;
( my $markers_ortho_species = $ortho_species ) =~ s/_/ /g;
my $markers_analysis    = basename( $0 );
my $markers_marker_type = 'Gene Prediction';

Bio::EnsEMBL::Registry->load_all( $reg );

# connects to Ensembl, and to the  'core' database of the selected species. 
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' ) 
    || pod2usage("\n[*DIE]No core DB for $species set in $reg\n");

my $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara', 'compara' )
    || pod2usage("\n[*DIE]No compara DB set in $reg\n");

my $ortho_genome_db_adaptor = $CMP_DBA->get_GenomeDBAdaptor();

my $ortho_genome_db
    = $ortho_genome_db_adaptor->fetch_by_registry_name( $ortho_species )
    or pod2usage( join("\n",
        "[*DIE] Orthologue species '$ortho_species' taxonomy id not found",
        "Use '--list' for a list of valid species" 
    ) );

my $ortho_species_id = $ortho_genome_db->taxon_id;

print "Processing $species (=> $ortho_species ($ortho_species_id)\n";

# There are two ways to connect to the Ensembl Compara database. The
# old way uses the Bio::EnsEMBL::Compara::DBSQL::DBAdaptor explicitly.
# The new one uses the Bio::EnsEMBL::Registry module, which can read
# either a global or a specific configuration file or auto-configure
# itself.

# Create the ADAPTORS that we will be using to fetch data from the database
my $gene_adaptor    = $ENS_DBA->get_adaptor( 'Gene' )
                   || die( "[*DIE] Cannot ENS_DBA->get_adaptor('Gene')" );

#
# A member can be a gene or a protein. Most of them are in the
# core databases of their correspondent specie.
#
my $member_adaptor  = $CMP_DBA->get_adaptor('Member') 
                   || die("[*DIE] Cannot CMP_DBA->get_adaptor('Member')");
my $HomologyAdaptor = $CMP_DBA->get_adaptor( 'Homology' )
                   || die("[*DIE] Cannot CMP_DBA->get_adaptor('Homology')");    
my $db = $ENS_DBA->dbc->dbname;

print "[INFO] Collecting all genes from $db \n";

my $mdb = Gramene::Marker::DB->new;

#
# Get a list of all gene IDs, and loop through them
#
my $gene_id_arrayref = $gene_adaptor->list_stable_ids;
my $total_genes      = scalar @$gene_id_arrayref;

printf "[INFO] Collected genes %s; processing...\n", 
    commify($total_genes);

#
# Header
#
if ( $markers_fh ) {
    print $markers_fh join( "\t", qw[
        from_marker_id
        to_marker_id
        analysis_name
        analytical_correspondence_type
    ]), "\n";
}

if ( $out_fh ) {
    print $out_fh join( "\t", qw[
        gene_stable_id        sequence_name gene1_start
        gene1_end             species_id_1  homolog2_member_id
        sequence_name2        gene2_start   gene2_end
        homology_relationship species_id_2 
    ]), "\n";
}

my ( $num_genes, $num_coding, $num_chr_genes,
    $chr_species_genes, $ignored, $ignored_species ) = ( 0, 0, 0, 0, 0, 0 );

while ( my $id = shift @$gene_id_arrayref ) {
    $num_genes++; 

    # Get the gene object, and skip unless protein coding
    my $gene = $gene_adaptor->fetch_by_stable_id( $id )
        || die( "[*DIE] Cannot fetch_by_stable_id gene $id" );

    next unless $gene->biotype eq 'protein_coding';

    $num_coding++;

    $num_chr_genes++;

    if ( 
        my $member 
        = $member_adaptor->fetch_by_source_stable_id( 'ENSEMBLGENE', $id )
    ) {
        # Get the homologues of the gene
        my $homologies
            = $HomologyAdaptor->fetch_all_by_Member( $member );

        for my $homology ( @$homologies ) {
            if ( $homology->description =~ /ortholog/ ) {    
                my $description     = $homology->description;
                my $homologue_genes = $homology->gene_list();

                my ( $homologue1, $homologue2 ) = @$homologue_genes;

                # To make sure that the $homologue2 gene is different from
                # the gene that been analyzed
                unless ( $homologue1->stable_id eq $id ) {
                    ( $homologue1, $homologue2 ) =
                    ( $homologue2, $homologue1 ) ;
                }

                if ( $homologue2->taxon_id == $ortho_species_id ) {
                    # Make sure the gene is on the longest assembled sequence
                    $gene->project('toplevel');

                    $chr_species_genes++;

                    if ( $markers_fh ) {
                        my @marker_ids;
                        for my $gene (
                            {
                                marker_name => $id,
                                marker_type => $markers_marker_type,
                                taxonomy    => $markers_species,
                            },
                            {
                                marker_name => $homologue2->stable_id,
                                marker_type => $markers_marker_type,
                                taxonomy    => $markers_ortho_species,
                            },
                        ) {
                            my @markers = $mdb->marker_search( %$gene );
                            if ( scalar @markers == 1 ) {
                                push @marker_ids, $markers[0]->{'marker_id'};
                            }
                            else {
                                printf STDERR 
                                    "WARNING: %s markers found for %s %s %s\n",
                                    scalar @markers, 
                                    $gene->{'taxonomy'},
                                    $gene->{'marker_type'},
                                    $gene->{'marker_name'},
                                ;
                            }
                        }

                        if ( scalar @marker_ids == 2 ) {
                            print $markers_fh join( "\t",
                                @marker_ids, 
                                $markers_analysis,
                                $description,
                            ), "\n";
                        }
                        else {
                            printf STDERR 
                                "Couldn't find markers for %s => %s\n",
                                $id, $homologue2->stable_id;
                        }
                    }

                    if ( $out_fh ) {
                        print $out_fh join( "\t",
                            $id,                     
                            $gene->seq_region_name,
                            $gene->seq_region_start, 
                            $gene->seq_region_end,
                            $homologue1->taxon_id,             
                            $homologue2->stable_id,
                            $homologue2->chr_name,   
                            $homologue2->chr_start,
                            $homologue2->chr_end,    
                            $description,
                            $homologue2->taxon_id, 
                        ), "\n";
                    }
                }
                else { 
                    $ignored_species++; 
                }
            }
        }
    }

    printf "[INFO] processed %s of %s\r",
        map { commify( $_ ) } $num_genes, $total_genes;
}

if ( $out_fh ) {
    print $out_fh "coding genes:$num_coding,num_genes:$num_genes, ",
          "num_tree_genes:$num_chr_genes,ignored:$ignored\n";
}

my @pairs;
for ( 
    [ 'Coding'         , $num_coding        ],
    [ 'Chr. Genes'     , $num_chr_genes     ],
    [ 'Ignored Genes'  , $ignored           ],
    [ 'Ignored Species', $ignored_species   ],
    [ 'Chr. Species'   , $chr_species_genes ],
    [ 'Total'          , $num_genes         ],
) {
    push @pairs, [ $_->[0], commify($_->[1]) ];
}

my $max_num = max( map { length($_) } ( map { $_->[1] }@pairs ), 'Number' );

my $tab = Text::TabularDisplay->new( qw[ Type Number ] );
for my $pair ( @pairs ) {
    $tab->add( $pair->[0], sprintf( "%${max_num}s", $pair->[1] ) );
}

print join( "\n", '', $tab->render, '', '' );

exit 0;

# -----------------------------------------------------------

=pod

=head1 NAME

perl dump_all_ortholog_genes.pl - Dumps the ortholog genes

=head1 SYNOPSIS

  ./dump_all_ortholog_genes.pl -s Oryza_sativa -o Sorghum_bicolor \
  -f os_sb_orthologs.tab 2> ./os_sb_orthologs.err

Required arguments:

  -f|--out-file          The file to print the output
  -m|--markers-file      The file to print output for markers import
  -s|--species           Species key in Ensembl registry file
  -o|--ortho_species     Ortholog Species bionomial name

Options:

  -e|--ensembl_registry  Path to Ensembl registry file; the default 
                         is "/usr/local/gramene/conf/ensembl.registry"

  -h|--help              Show brief help and exit.
  --man                  Show detailed help

=head1 DESCRIPTION

Script that iterates through each gene in the given species and dumps
the subroot ID of any protein tree alongside the gene genomic
coordinates.

B<The Ensembl Registry>

  The database connection details for both Ensembl species and compara   ##
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::Registry;

  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new( 
    '-species' => 'Oryza_sativa', 
    '-group'   => 'core', 
    '-dbname'  => 'Oryza_sativa_japonica_core_48_28', 
  );

  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( 
    '-species' => 'compara',
    '-group'   => 'compara',
    '-dbname'  => 'ensembl_compara_48_28', 
  );
  ---

=head1 AUTHORS

Will Spooner E<lt>whs@ebi.ac.ukE<gt>,
Liya Ren E<lt>ren@cshl.eduE<gt>,
Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=cut
