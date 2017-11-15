# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::SyntenyBlock
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::SyntenyBlock

=head1 SYNOPSIS

  my $synblock = Bio::EnsEMBL::Analysis::RunnableDB::SyntenyBlock->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -analysis => $analysis,
     );
  $synblock->fetch_input;
  $synblock->run;
  $synblock->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable SyntenyBlock which wraps the program synteny_block.pl

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the marker_feature and marker tables

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::SyntenyBlock;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::SyntenyBlock;

use vars qw(@ISA $BAND_TO_BASEPAIR);

@ISA              = qw(Bio::EnsEMBL::Analysis::RunnableDB);
$BAND_TO_BASEPAIR = 4096;

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Arg [3]   : Bio::EnsEMBL::Analysis
  Function  : create a Bio::EnsEMBL::Analysis::RunnableDB object
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB
  Exceptions: throws if not passed either a dbadaptor, input id or
  an analysis object
  Example   : $rdb = $perl_path->new( -analysis => $self->analysis,
                                      -input_id => $self->input_id,
                                      -db => $self->adaptor->db );

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($target_db, $compara_db, $ort)
        = rearrange([ 'TARGET_DB', 'COMPARA_DB', 'ORT' ], @args);
    if (!$target_db) {
        throw("[*DIE] You must supply target db");
    }
    if (!$compara_db) {
        throw("[*DIE] You must supply compara db");
    }

    $self->compara_db($compara_db);
    $self->target_db($target_db);
    $self->{ort} = $ort;

    my ($query_species, $target_species, $query_taxon_id, $target_taxon_id);
    $query_species = $self->db->get_MetaContainer->get_scientific_name;
    $query_species =~ s/\s/_/g;
    $target_species
        = $self->target_db->get_MetaContainer->get_scientific_name;
    $target_species =~ s/\s/_/g;
    $query_taxon_id  = $self->db->get_MetaContainer->get_taxonomy_id;
    $target_taxon_id = $self->target_db->get_MetaContainer->get_taxonomy_id;

    $self->{species_meta} = {
        query => {
            name         => $query_species,    #'Zea_mays',
            core_adaptor => $self->db,
            taxon_id     => $query_taxon_id,
            abbr =>
                join('', map { substr $_, 0, 1 } split /_/, $query_species),
        },
        target => {
            name         => $target_species,    #'Oryza_sativa',
            core_adaptor => $self->target_db,
            taxon_id     => $target_taxon_id,
            abbr =>
                join('', map { substr $_, 0, 1 } split /_/, $target_species),
        },
    };
    return $self;
}

=head2 compara_db

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Compara::DBAdaptor
  Function  : container for dbadaptor
  Returntype: Bio::EnsEMBL::Compara::DBSQL::DBAdaptor
  Exceptions: throws if not passed a Bio::EnsEMBL::Compara::DBSQL::DBAdaptor
  Example   : 

=cut

sub compara_db {
    my $self       = shift;
    my $compara_db = shift;
    if ($compara_db) {
        throw(
            "[*DIE] Must pass RunnableDB:compara_db a Bio::EnsEMBL::Compara::DBSQL::DBAdaptor "
                . "not a "
                . $compara_db)
            unless (
            $compara_db->isa('Bio::EnsEMBL::Compara::DBSQL::DBAdaptor'));
        $self->{'compara_db'} = $compara_db;
    }
    return $self->{'compara_db'};
}

=head2 target_db

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : container for dbadaptor
  Returntype: Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions: throws if not passed a Bio::EnsEMBL::DBSQL::DBAdaptor
  Example   : 

=cut

sub target_db {
    my $self      = shift;
    my $target_db = shift;
    if ($target_db) {
        throw(
            "[*DIE] Must pass RunnableDB:target_db a Bio::EnsEMBL::DBSQL::DBAdaptor "
                . "not a "
                . $target_db)
            unless ($target_db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
        $self->{'target_db'} = $target_db;
    }
    return $self->{'target_db'};
}

=head2 fetch_input
    
  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::SyntenyBlock
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : ort file format
AT1G01020       1       6788    9130    3702At  PINA_004507     Un      7707405 7710294 ortholog_many2many      29760Vp
AT1G01020       1       6788    9130    3702At  PINA_004513     Un      7746086 7750113 ortholog_many2many      29760Vp
AT1G01020       1       6788    9130    3702At  PINA_021731     15      16526104        16531329        ortholog_many2many      29760Vp
 

=cut

sub fetch_input {
    my ($self) = @_;

    my $db         = $self->db;
    my $compara_db = $self->compara_db;

    my $query_species  = $self->{species_meta}{query}{query_species};
    my $target_species = $self->{species_meta}{target}{target_species};

    my @ort;
    if ($self->{ort}) {
        @ort = @{ $self->{ort} };
    } else {
        push @ort,
            join("\t",
            'gene_stable_id', 'sequence_name',
            'gene1_start',    'gene1_end',
            'specie_id_1',    'homolog2_member_id',
            'sequence_name2', 'gene2_start',
            'gene2_end',      'homology_relationship',
            'specie_id_2',);

        ##########################################################################
        # The following fetch data logic is adapted from Liya's script
        # brie:/home/liya/ensembl_script/dump_all_ortholog_genes.pl
        ##########################################################################

# Create the ADAPTORS that we will be using to fetch data from the database
# A member can be a gene or a protein. Most of them are in the core databases of their correspondient specie.

        my $gene_adaptor = $db->get_adaptor('Gene')
            || throw("[*DIE] Cannot db->get_adaptor('Gene')");
        my $member_adaptor = $compara_db->get_adaptor('GeneMember')
            || throw("[*DIE] Cannot compara_db->get_adaptor('Member')");
        my $HomologyAdaptor = $compara_db->get_adaptor('Homology')
            || throw("[*DIE] Cannot compara_db->get_adaptor('Homology')");

        my $dbname = $db->dbc->dbname;

        warn("[iNFO] Collecting all genes from $dbname\n");
        my $gene_id_arrayref = $gene_adaptor->list_stable_ids;
        my $total_genes      = scalar @$gene_id_arrayref;
        warn("[INFO] Collected genes $total_genes; processing...\n");

        my ($num_genes, $num_coding, $num_chr_genes, $chr_specie_genes,
            $ignored_species)
            = (0, 0, 0, 0, 0);
        my $ignored = 0;

        while (my $gene_id = shift @$gene_id_arrayref) {
            $num_genes++;    #COUNT

            # Get the gene object, and skip unless protein coding
            my $gene = $gene_adaptor->fetch_by_stable_id($gene_id)
                || warn("[*DIE] Cannot fetch_by_stable_id gene $gene_id");
            next unless $gene->biotype eq 'protein_coding';

            $num_coding++;    # COUNT
            my ($homologue1, $homologue2);    #added

            ######### TEST CODE
#liya : comment this line to get all genes instead of chromosome genes only
#unless($gene->seq_region_name =~ /^\d/){$ignored++; next;} # MANUAL MODIFICATION -turn this two "ifs" off when analizing maize
            $num_chr_genes++;                 #COUNT

			#if (my $member = $member_adaptor->fetch_by_source_stable_id(
			if (my $member = $member_adaptor->fetch_by_stable_id(
					#'ENSEMBLGENE', $gene_id
					$gene_id
                )
                )
            {

                #Obtain the homologues of the gene
                my $homologies
                    = $HomologyAdaptor->fetch_all_by_Member($member);
                foreach my $homology (@{$homologies}) {

                    # if the relationship involves orthology
                    if ($homology->description =~ /ortholog/) {
                        ($homologue1, $homologue2)
                            = @{ $homology->gene_list() };

                    # To make sure that the $homologue2 gene is different from
                    # the gene that it is been analyzed
                        unless ($homologue1->stable_id eq $gene_id) {
                            ($homologue1, $homologue2)
                                = ($homologue2, $homologue1);
                        }

                        if ($homologue2->taxon_id
                            == $self->{species_meta}{target}{taxon_id})
                        {

                     # Make sure the gene is on the longest assembled sequence
                            $gene->project('toplevel');

                            $chr_specie_genes++;
							my $homologue_gene2=$homologue2->get_Gene();
                            push @ort,
                                join("\t",
                                $gene_id,
                                $gene->seq_region_name,
                                $gene->seq_region_start,
                                $gene->seq_region_end,
                                $homologue1->taxon_id,
                                $homologue2->stable_id,
                                $homologue_gene2->seq_region_name,
                                $homologue_gene2->seq_region_start,
                                $homologue_gene2->seq_region_end,
                                $homology->description,
                                $homologue2->taxon_id,
                            );
                        } else {
                            $ignored_species++;
                        }
                    }
                }
            }

            # Update progress each 1000 genes
            if ($num_genes % 1000 == 0) {
                warn("[INFO] processed $num_genes of $total_genes\n");
            }
        }

        ##########################################################################
    }

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::SyntenyBlock->new(
        -analysis => $self->analysis,

        #-program => "/data/ware/gramene-pipeline/bin/synteny_block.pl",
        -program        => $self->analysis->program,
        -ort            => \@ort,
        -query_species  => $self->{species_meta}{query}{name},
        -target_species => $self->{species_meta}{target}{name},
    );
    $self->runnable($runnable);
    return 1;
}

=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : set analysis and slice on each feature
  Returntype: 1
  Exceptions: none
  Example   : 

=cut

sub write_output {
    my ($self)     = @_;
    my $compara_db = $self->compara_db;
    my $gdba       = $compara_db->get_adaptor('GenomeDB');

    my $meta = $self->{species_meta};
    foreach my $key (keys %$meta) {
        my $species      = $meta->{$key}->{'name'};
        my $core_adaptor = $meta->{$key}->{'core_adaptor'};
        $meta->{$key}->{'slice_adaptor'}
            = $core_adaptor->get_adaptor('Slice');
        my $species_name  
            = $core_adaptor->get_MetaContainer->list_value_by_key(
            'species.production_name')->[0]; #species.ensembl_alias_name')->[0];
        my $species_assembly
            = $core_adaptor->get_CoordSystemAdaptor->fetch_all->[0]->version;
        $meta->{$key}->{'genome_db'}
            = $gdba->fetch_by_name_assembly($species_name, $species_assembly)
            || throw("[*DIE] Species $species not in compara DB");
    }
    my $switched = $meta->{query}{abbr} le $meta->{target}{abbr} ? 0 : 1;
    local $compara_db->dbc->db_handle->{'AutoCommit'};
    eval {
        while (my $line = shift @{ $self->output })
        {
            next if $line =~ /blockID/;
            my @parts = (split(/\s+/, $line))[ 1 .. 4, 8 .. 11 ];
            my @empty = grep { !length($_) } @parts[ 0 .. 7 ];

            if (@empty) {

                # Sanity check
                warn("[WARN] Incomplete line in synteny file: $line\n");
            }
            my $syn_meta = {
                query => {
                    seqname => $parts[0],
                    start   => $parts[1],
                    end     => $parts[2],
                    strand  => $parts[3],
                },
                target => {
                    seqname => $parts[4],
                    start   => $parts[5],
                    end     => $parts[6],
                    strand  => $parts[7],
                },
            };
            if ($switched) {
                ($syn_meta->{query}, $syn_meta->{target})
                    = ($syn_meta->{target}, $syn_meta->{query});
            }
            my $synteny_id = $self->create_synteny_region(
                $meta->{query}->{'genome_db'},
                $meta->{target}->{'genome_db'},
            );
            foreach my $key (keys %$meta) {
                my $sp_meta     = $meta->{$key};
                my $sp_syn_meta = $syn_meta->{$key};

                my $species = $sp_meta->{name};
                my $seqname = $sp_syn_meta->{seqname};
                my $start   = $sp_syn_meta->{start};
                my $end     = $sp_syn_meta->{end};
                my $strand  = $sp_syn_meta->{strand};

                unless ($strand =~ /^[01]$/) {
                    $strand = $strand eq '+' ? 1 : 0;
                }

                # Convert all units into basepair
                if ($sp_meta->{units} && $sp_meta->{units} eq 'band') {
                    $start = $start * $BAND_TO_BASEPAIR;
                    $end   = $end * $BAND_TO_BASEPAIR;
                }

                # Convert syntenic region into toplevel
                my $slice
                    = $sp_meta->{slice_adaptor}
                    ->fetch_by_region(undef, $seqname, $start, $end, $strand)
                    or throw(
                    "[*DIE] Sequence $seqname not in core DB for $species");
                $slice   = $slice->project('toplevel')->[0]->to_Slice();
                $seqname = $slice->seq_region_name;
                $start   = $slice->start;
                $end     = $slice->end;

                #$strand  = $slice->strand;
                # Get the DnaFrag corresponding to the top-level slice
                my $dna_frag
                    = $self->get_dna_frag($sp_meta->{genome_db}, $slice);

                # Load the syntenic region into compara
                $self->create_dnafrag_region($synteny_id, $dna_frag->dbID,
                    $start, $end, $strand);
            }
        }
        $compara_db->dbc->db_handle->commit;
    };
    if ($@) {
        throw("[*DIE} Transaction aborted because $@");
        eval { $compara_db->dbc->db_handle->rollback };
    }
    return 1;
}

#----------------------------------------------------------------------
#
sub get_dna_frag {
    my $self      = shift;
    my $genome_db = shift || throw("[*DIE] Need a valid GenomeDB object");
    my $slice     = shift || throw("[*DIE] Need a Slice object");

    my $ENS_DBA = $self->compara_db;
    my $seqname = $slice->seq_region_name;

    my $dnafa = $ENS_DBA->get_adaptor('DnaFrag');
    my $dnafrag = $dnafa->fetch_by_GenomeDB_and_name($genome_db, $seqname);

    unless ($dnafrag) {
        $dnafrag = Bio::EnsEMBL::Compara::DnaFrag->new(
            -name              => $seqname,
            -length            => $slice->length,
            -coord_system_name => $slice->coord_system->name,
            -genome_db         => $genome_db,
        );
        $dnafa->store($dnafrag);
    }
    return $dnafrag;
}

#----------------------------------------------------------------------
# Inserts a record into the synteny_region table of the ensembl DB and
# returns it's dbID
sub create_synteny_region {
    my $self       = shift;
    my @genome_dbs = @_;
    @genome_dbs >= 2 || throw("[*DIE] Need at least two GenomeDBs");

    my $ENS_DBA = $self->compara_db;
    my $method_name = "SYNTENY";

    my $mlssa     = $ENS_DBA->get_adaptor('MethodLinkSpeciesSet');

    my $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_name,
        [@genome_dbs], 1);

    unless ($mlss) {
		my $method=$ENS_DBA->get_MethodAdaptor()->fetch_by_type($method_name);
		my $species_set=$ENS_DBA->get_SpeciesSetAdaptor()->fetch_by_GenomeDBs([@genome_dbs]); 
        $mlss = Bio::EnsEMBL::Compara::MethodLinkSpeciesSet->new(
            -method => $method,
            -species_set_obj => $species_set,
        );
        $mlssa->store($mlss);
    }

    my $sql = qq(
INSERT INTO synteny_region ( method_link_species_set_id )
VALUES (?) );

    my $sth = $ENS_DBA->dbc->prepare($sql);
    $sth->execute($mlss->dbID) || throw("[*DIE] " . $sth->errstr);
    return $sth->{'mysql_insertid'};
}

#----------------------------------------------------------------------
#
sub create_dnafrag_region {
    my $self           = shift;
    my $synteny_id     = shift || throw("[*DIE] Need a synteny_region_id");
    my $dnafrag_id     = shift || throw("[*DIE] Need a dnafrag_id");
    my $dnafrag_start  = shift || throw("[*DIE] Need a dnafrag_start");
    my $dnafrag_end    = shift || throw("[*DIE] Need a dnafrag_end");
    my $dnafrag_strand = shift || 0;

    my $ENS_DBA = $self->compara_db;

    my $sql = qq(
INSERT INTO dnafrag_region 
       ( synteny_region_id, dnafrag_id, dnafrag_start, dnafrag_end, dnafrag_strand )
VALUES ( ?, ?, ?, ?, ? ) );

    my $sth = $ENS_DBA->dbc->prepare($sql);
    $sth->execute(
        $synteny_id,  $dnafrag_id, $dnafrag_start,
        $dnafrag_end, $dnafrag_strand
    ) || throw("[[*DIE] " . $sth->errstr);
    return $sth->{'mysql_insertid'};
}

1;
