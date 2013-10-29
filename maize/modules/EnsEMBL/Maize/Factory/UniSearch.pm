package EnsEMBL::Maize::Factory::UniSearch;

use strict;

use EnsEMBL::Web::Factory;
use EnsEMBL::Web::Proxy::Object;
our @ISA = qw(EnsEMBL::Web::Factory);

sub createObjects {
    my $self = shift;
    my $idx = $self->param('type') || $self->param('idx') || 'all';
    ## Action search...
    my $search_method = "search_" . uc($idx);
    if ($self->param('q')) {
        if ($self->can($search_method)) {
            $self->{to_return}     = 30;
            $self->{_result_count} = 0;
            $self->{_results}      = [];
            $self->$search_method();
            $self->DataObjects(
                new EnsEMBL::Web::Proxy::Object(
                    'UniSearch',
                    {   'idx'     => $idx,
                        'q'       => $self->param('q'),
                        'results' => $self->{results}
                    },
                    $self->__data
                )
            );
        } else {
            $self->problem(
                'fatal', 'Unknown search method', qq(
      <p>
        Sorry do not know how to search for features of type "$idx"
      </p>)
            );
        }
    } else {
        $self->DataObjects(
            new EnsEMBL::Web::Proxy::Object(
                'UniSearch', { 'idx' => $idx, 'q' => '', 'results' => {} },
                $self->__data
            )
        );
    }
}

sub terms {
    my $self = shift;
    my @list = ();
    foreach my $kw (split /\s+/, join ' ', $self->param('q')) {
        $kw =~ s/(\w{5,})$/$1\*/;
        my $seq = $kw =~ s/\*+$/%/ ? 'like' : '=';

        #------- BEGIN MAIZE-SPECIFIC KEYWORDS (pasternak) -------#
        # Remove ZmmBB prefix
        $kw =~ s/^zmmbb//i;

        # Preformat clone names (zero-padding), e.g b35L3 -> b0035L03
        if ($kw =~ m/^([bc])(\d+)([a-z])(\d+)$/i) {
            $kw = sprintf('%s%.4d%s%.2d', lc($1), $2, uc($3), $4);
        }

        #------- END MAIZE-SPECIFIC KEYWORDS (pasternak) ---------#
        push @list, [ $seq, $kw ];
    }
    return @list;
}

sub count {
    my ($self, $db, $sql, $comp, $raw_keyword) = @_;
    my $result = +{
        'count'   => 0,
        'species' => undef,
    };
    for my $species (qw(Zea_mays Zea_mays2 Zea_mays_external)) {
        $result->{'species'} = $species;
        my $dbh = $self->database($db, $species);
        next unless $dbh;
        my $keyword = $dbh->dbc->db_handle->quote($raw_keyword);
        (my $t = $sql) =~ s/'\[\[KEY\]\]'/$keyword/g;
        $t =~ s/\[\[COMP\]\]/$comp/g;

        #my( $res ) = $dbh->db_handle->selectrow_array( $t );
        my ($count) = $dbh->dbc->db_handle->selectrow_array($t);
        if ($count > 0) {
            $result->{'count'} = $count;

            # Use the result from the first occurrence.
            last;
        }
    }
    return $result;
}

sub _fetch {
    my ($self, $species, $db, $search_SQL, $comparator, $raw_keyword, $limit)
        = @_;
    my $dbh = $self->database($db, $species);
    next unless $dbh;
    my $keyword = $dbh->dbc->db_handle->quote($raw_keyword);
    (my $t = $search_SQL) =~ s/'\[\[KEY\]\]'/$keyword/g;
    $t =~ s/\[\[COMP\]\]/$comparator/g;

    #my $res = $dbh->db_handle->selectall_arrayref( "$t limit $limit" );
    my $res = $dbh->dbc->db_handle->selectall_arrayref("$t limit $limit");
    push @{ $self->{_results} }, +{ 'species' => $species, 'rows' => $res };
}

sub search_ALL {
    my ($self, $species) = @_;
    my $package_space = __PACKAGE__ . '::';
    no strict 'refs';
    my @methods = map { /(search_\w+)/ && $1 ne 'search_ALL' ? $1 : () }
        keys %$package_space;
    my @ALL = ();

    foreach my $method (@methods) {
        $self->{_result_count} = 0;
        $self->{_results}      = [];
        $self->{to_return}     = 25;
        if ($self->can($method)) {
            $self->$method;
        }
    }
    return @ALL;
}

sub _fetch_results {
    my $self  = shift;
    my @terms = $self->terms();
    foreach my $query (@_) {
        my ($db, $subtype, $count_SQL, $search_SQL) = @$query;
        foreach my $term (@terms) {
            my $count_ref
                = $self->count($db, $count_SQL, $term->[0], $term->[1]);
            my $count_new = $count_ref->{'count'};
            if ($count_new) {
                if ($self->{to_return} > 0) {
                    my $limit
                        = $self->{to_return} < $count_new
                        ? $self->{to_return}
                        : $count_new;
                    $self->_fetch($count_ref->{'species'},
                        $db, $search_SQL, $term->[0], $term->[1], $limit);
                    $self->{to_return} -= $count_new;
                }
                $self->{'_result_count'} += $count_new;
            }
        }
    }
}

sub search_SNP {
    my $self = shift;
    $self->_fetch_results(
        [   'variation', 'SNP',
            "select count(*) from variation as v where name = '[[KEY]]'",
            "select s.name as source, v.name
      from source as s, variation as v
     where s.source_id = v.source_id and v.name = '[[KEY]]'"
        ],
        [   'variation', 'SNP',
            "select count(*) from variation as v, variation_synonym as vs
       where v.variation_id = vs.variation_id and vs.name = '[[KEY]]'",
            "select s.name as source, v.name
      from source as s, variation as v, variation_synonym as vs
     where s.source_id = v.source_id and v.variation_id = vs.variation_id and vs.name = '[[KEY]]'"
        ]
    );

    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'idx'     => 'SNP',
                'subtype' => "$row->[0] SNP",
                'ID'      => $row->[1],
                'URL'  => "/$species/snpview?source=$row->[0];snp=$row->[1]",
                'desc' => '',
                'species' => $species
                };
        }
    }
    $self->{'results'}{'SNP'} = [ \@rows, $self->{_result_count} ];
}

sub search_GENOMICALIGNMENT {
    my $self = shift;
    $self->_fetch_results(
        [   'core',
            'DNA',
            "select count(distinct analysis_id, hit_name) from dna_align_feature where hit_name [[COMP]] '[[KEY]]'",
            "select a.logic_name, f.hit_name, 'Dna', 'core',count(*)  from dna_align_feature as f, analysis as a where a.analysis_id = f.analysis_id and f.hit_name [[COMP]] '[[KEY]]' group by a.logic_name, f.hit_name"
        ],
        [   'core',
            'Protein',
            "select count(distinct analysis_id, hit_name) from protein_align_feature where hit_name [[COMP]] '[[KEY]]'",
            "select a.logic_name, f.hit_name, 'Protein', 'core',count(*) from protein_align_feature as f, analysis as a where a.analysis_id = f.analysis_id and f.hit_name [[COMP]] '[[KEY]]' group by a.logic_name, f.hit_name"
        ],
        [   'vega',
            'DNA',
            "select count(distinct analysis_id, hit_name) from dna_align_feature where hit_name [[COMP]] '[[KEY]]'",
            "select a.logic_name, f.hit_name, 'Dna', 'vega', count(*) from dna_align_feature as f, analysis as a where a.analysis_id = f.analysis_id and f.hit_name [[COMP]] '[[KEY]]' group by a.logic_name, f.hit_name"
        ],
        [   'est',
            'DNA',
            "select count(distinct analysis_id, hit_name) from dna_align_feature where hit_name [[COMP]] '[[KEY]]'",
            "select a.logic_name, f.hit_name, 'Dna', 'est', count(*) from dna_align_feature as f, analysis as a where a.analysis_id = f.analysis_id and f.hit_name [[COMP]] '[[KEY]]' group by a.logic_name, f.hit_name"
        ]
    );

    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'idx'     => 'GenomicAlignment',
                'subtype' => "$row->[0] $row->[2] alignment feature",
                'ID'      => $row->[1],
                'URL' =>
                    "/$species/featureview?type=$row->[2]AlignFeature;db=$row->[3];id=$row->[1]",
                'desc' =>
                    "This $row->[2] alignment feature hits the genome in $row->[4] place(s).",
                'species' => $species
                };
        }
    }
    $self->{'results'}{'GenomicAlignments'}
        = [ \@rows, $self->{_result_count} ];
}

sub search_DOMAIN {
    my $self = shift;
    $self->_fetch_results(
        [   'core', 'Domain',
            "select count(*) from xref as x, external_db as e 
        where e.external_db_id = x.external_db_id and x.dbprimary_acc [[COMP]] '[[KEY]]'",
            "select x.dbprimary_acc, x.description
         FROM xref as x, external_db as e
        WHERE e.db_name = 'Interpro' and e.external_db_id = x.external_db_id and
              x.dbprimary_acc [[COMP]] '[[KEY]]'"
        ],
        [   'core', 'Domain',
            "select count(*) from xref as x, external_db as e 
        where e.external_db_id = x.external_db_id and x.dbprimary_acc [[COMP]] '[[KEY]]'",
            "SELECT x.dbprimary_acc, x.description
         FROM xref as x, external_db as e
        WHERE e.db_name = 'Interpro' and e.external_db_id = x.external_db_id and
              x.description [[COMP]] '[[KEY]]'"
        ],
    );

    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'URL'     => "/$species/domainview?domain=$row->[0]",
                'idx'     => 'Domain',
                'subtype' => 'Domain',
                'ID'      => $row->[0],
                'desc'    => $row->[1],
                'species' => $species
                };
        }
    }
    $self->{'results'}{'Domain'} = [ \@rows, $self->{_result_count} ];
}

sub search_FAMILY {
    my ($self, $species) = @_;
    $self->_fetch_results(
        [   'compara',
            'Family',
            "select count(*) from family where stable_id [[COMP]] '[[KEY]]'",
            "select stable_id, description FROM family WHERE stable_id  [[COMP]] '[[KEY]]'"
        ],
        [   'compara',
            'Family',
            "select count(*) from family where description [[COMP]] '[[KEY]]'",
            "select stable_id, description FROM family WHERE description [[COMP]] '[[KEY]]'"
        ]
    );
    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'URL'     => "/$species/familyview?family=$row->[0]",
                'idx'     => 'Family',
                'subtype' => 'Family',
                'ID'      => $row->[0],
                'desc'    => $row->[1],
                'species' => $species
                };
        }
    }
    $self->{'results'}{'Family'} = [ \@rows, $self->{_result_count} ];
}

sub search_SEQUENCE {
    my $self               = shift;
    my $misc_feature_count = <<__END_COUNT_SQL__;
select count(distinct misc_feature_id)
  from misc_attrib where value [[COMP]] '[[KEY]]'
   and attrib_type_id in (select attrib_type_id
                            from attrib_type
                           where code = 'name' or code = 'embl_acc')
__END_COUNT_SQL__
    my $misc_feature_select = <<__END_SELECT_SQL__;
select ma.value, ms.name, seq_region_end-seq_region_start, 'miscfeature'
  from misc_attrib as ma, misc_set as ms,
       misc_feature as mf, misc_feature_misc_set as mfms, attrib_type as at
 where ma.misc_feature_id = mf.misc_feature_id
   and mfms.misc_feature_id = mf.misc_feature_id
   and mfms.misc_set_id = ms.misc_set_id
   and at.attrib_type_id = ma.attrib_type_id
   and at.code in ('name', 'embl_acc')
   and ma.value [[COMP]] '[[KEY]]'
group by mf.misc_feature_id
__END_SELECT_SQL__

    $self->_fetch_results(
        [   'core',
            'Sequence',
            "select count(*) from seq_region where name [[COMP]] '[[KEY]]'",
            "select sr.name, cs.name, sr.length, 'region' from seq_region as sr, coord_system as cs where cs.coord_system_id = sr.coord_system_id and sr.name [[COMP]] '[[KEY]]'"
        ],
        [ 'core', 'Sequence', $misc_feature_count, $misc_feature_select ]
    );

=pod
        [   'core', 'Sequence',
            "select count(distinct misc_feature_id) from misc_attrib where value [[COMP]] '[[KEY]]'",
            "select ma1.value, ms.name, seq_region_end-seq_region_start, 'miscfeature'
        from misc_attrib as ma1, misc_set as ms, misc_attrib as ma2, misc_feature as mf, misc_feature_misc_set as mfms, attrib_type as at
       where ma1.misc_feature_id = mf.misc_feature_id and ma2.misc_feature_id and mfms.misc_feature_id = mf.misc_feature_id and mfms.misc_set_id = ms.misc_set_id and
             ma1.attrib_type_id = at.attrib_type_id and at.code in ('name','clone_name','embl_acc','synonym','sanger_project') and ma2.value [[COMP]] '[[KEY]]'
       group by mf.misc_feature_id"
        ]
=cut

    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};

        # TODO: (shiran) Hack Alert!!! This is maize-specific
        my $view = $species eq 'Zea_mays' ? 'cytoview' : 'contigview';

        for my $row (@{$result}) {
            my ($feature_name, $feature_map, $feature_length, $feature_type)
                = @$row;
            my $url
                = (lc($feature_map) eq 'chromosome'
                    && length($feature_name) < 10)
                ? "/$species/mapview?chr=$feature_name"
                : "/$species/$view?$feature_type=$feature_name";

            push @rows,
                +{
                'URL'     => $url,
                'idx'     => 'Sequence',
                'subtype' => ucfirst($feature_map),
                'ID'      => $feature_name,
                'desc'    => '',
                'species' => $species
                };
        }
    }
    $self->{'results'}{'Sequence'} = [ \@rows, $self->{_result_count} ];
}

sub search_OLIGOPROBE {
    my $self = shift;
    $self->_fetch_results(
        [   'core', 'OligoProbe',
            "select count(distinct probeset) from oligo_probe where probeset [[COMP]] '[[KEY]]'",
            "select ap.probeset, group_concat(distinct aa.name order by aa.name separator ' ') from oligo_probe ap, oligo_array as aa
        where ap.probeset [[COMP]] '[[KEY]]' and ap.oligo_array_id = aa.oligo_array_id group by ap.probeset"
        ],
    );
    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'URL' => "/$species/featureview?type=OligoProbe;id=$row->[0]",
                'idx' => 'OligoProbe',
                'subtype' => 'OligoProbe',
                'ID'      => $row->[0],
                'desc'    => "Is a member of the following arrays: $row->[1]",
                'species' => $species
                };
        }
    }
    $self->{'results'}{'OligoProbe'} = [ \@rows, $self->{_result_count} ];
}

sub search_QTL {
    my $self = shift;
    $self->_fetch_results(
        [   'core', 'QTL',
            "select count(*)
  from qtl_feature as qf, qtl as q
 where q.qtl_id = qf.qtl_id and q.trait [[COMP]] '[[KEY]]'",
            "select q.trait, concat( sr.name,':', qf.seq_region_start, '-', qf.seq_region_end ),
       qf.seq_region_end - qf.seq_region_start
  from seq_region as sr, qtl_feature as qf, qtl as q
 where q.qtl_id = qf.qtl_id and qf.seq_region_id = sr.seq_region_id and q.trait [[COMP]] '[[KEY]]'"
        ],
        [   'core', 'QTL',
            "select count(*)
  from qtl_feature as qf, qtl_synonym as qs ,qtl as q
 where qs.qtl_id = q.qtl_id and q.qtl_id = qf.qtl_id and qs.source_primary_id [[COMP]] '[[KEY]]'",
            "select q.trait, concat( sr.name,':', qf.seq_region_start, '-', qf.seq_region_end ),
       qf.seq_region_end - qf.seq_region_start
  from seq_region as sr, qtl_feature as qf, qtl_synonym as qs ,qtl as q
 where qs.qtl_id = q.qtl_id and q.qtl_id = qf.qtl_id and qf.seq_region_id = sr.seq_region_id and qs.source_primary_id [[COMP]] '[[KEY]]'"
        ]
    );

    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'URL'     => "/$species/cytoview?l=$row->[1]",
                'idx'     => 'QTL',
                'subtype' => 'QTL',
                'ID'      => $row->[0],
                'desc'    => '',
                'species' => $species
                };
        }
    }
    $self->{'results'}{'QTL'} = [ \@rows, $self->{_result_count} ];
}

sub search_MARKER {
    my $self = shift;
    $self->_fetch_results(
        [   'core',
            'Marker',
            "select count(distinct name) from marker_synonym where name [[COMP]] '[[KEY]]'",
            "select distinct name from marker_synonym where name [[COMP]] '[[KEY]]'"
        ]
    );

    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'URL'       => "/$species/markerview?marker=$row->[0]",
                'URL_extra' => [
                    'C',
                    'View marker in ContigView',
                    "/$species/cytoview?marker=$row->[0]"
                ],
                'idx'     => 'Marker',
                'subtype' => 'Marker',
                'ID'      => $row->[0],
                'desc'    => '',
                'species' => $species
                };
        }
    }
    $self->{'results'}{'Marker'} = [ \@rows, $self->{_result_count} ];
}

sub search_GENE {
    my $self = shift;

    my @databases = ('core');
    push @databases, 'vega'
        if $self->species_defs->databases->{'ENSEMBL_VEGA'};
    push @databases, 'est'
        if $self->species_defs->databases->{'ENSEMBL_OTHERFEATURES'};
    foreach my $db (@databases) {
        $self->_fetch_results(
            [   $db,
                'Gene',
                "select count(*) from gene_stable_id WHERE stable_id [[COMP]] '[[KEY]]'",
                "SELECT gsi.stable_id, g.description, '$db', 'geneview', 'gene' FROM gene_stable_id as gsi, gene as g WHERE gsi.gene_id = g.gene_id and gsi.stable_id [[COMP]] '[[KEY]]'"
            ],
            [   $db,
                'Gene',
                "select count(*) from transcript_stable_id WHERE stable_id [[COMP]] '[[KEY]]'",
                "SELECT gsi.stable_id, g.description, '$db', 'transview', 'transcript' FROM transcript_stable_id as gsi, transcript as g WHERE gsi.transcript_id = g.transcript_id and gsi.stable_id [[COMP]] '[[KEY]]'"
            ],
            [   $db,
                'Gene',
                "select count(*) from translation_stable_id WHERE stable_id [[COMP]] '[[KEY]]'",
                "SELECT gsi.stable_id, x.description, '$db', 'protview', 'peptide' FROM translation_stable_id as gsi, translation as g, transcript as x WHERE g.transcript_id = x.transcript_id and gsi.translation_id = g.translation_id and gsi.stable_id [[COMP]] '[[KEY]]'"
            ],

            [   $db, 'Gene',
                "select count( * ) from object_xref as ox, xref as x
        where ox.ensembl_object_type = 'Gene' and ox.xref_id = x.xref_id and x.dbprimary_acc [[COMP]] '[[KEY]]'",
                "SELECT gsi.stable_id, concat( display_label, ' - ', g.description ), '$db', 'geneview', 'gene' from gene_stable_id as gsi, gene as g, object_xref as ox, xref as x
        where gsi.gene_id = ox.ensembl_id and ox.ensembl_object_type = 'Gene' and gsi.gene_id = g.gene_id and
              ox.xref_id = x.xref_id and x.dbprimary_acc [[COMP]] '[[KEY]]'"
            ],
            [   $db, 'Gene',
                "select count( * ) from object_xref as ox, xref as x
        where ox.ensembl_object_type = 'Gene' and ox.xref_id = x.xref_id and
              x.display_label [[COMP]] '[[KEY]]' and not(x.dbprimary_acc [[COMP]] '[[KEY]]')",
                "SELECT gsi.stable_id, concat( display_label, ' - ', g.description ), '$db', 'geneview', 'gene' from gene_stable_id as gsi, gene as g, object_xref as ox, xref as x
        where gsi.gene_id = ox.ensembl_id and ox.ensembl_object_type = 'Gene' and gsi.gene_id = g.gene_id and
              ox.xref_id = x.xref_id and x.display_label [[COMP]] '[[KEY]]' and
              not(x.dbprimary_acc [[COMP]] '[[KEY]]')"
            ],
            [   $db, 'Gene',
                "select count( * ) from object_xref as ox, xref as x
        where ox.ensembl_object_type = 'Transcript' and ox.xref_id = x.xref_id and x.dbprimary_acc [[COMP]] '[[KEY]]'",
                "SELECT gsi.stable_id, concat( display_label, ' - ', g.description ), '$db', 'transview', 'transcript' from transcript_stable_id as gsi, transcript as g, object_xref as ox, xref as x
        where gsi.transcript_id = ox.ensembl_id and ox.ensembl_object_type = 'Transcript' and gsi.transcript_id = g.transcript_id and
              ox.xref_id = x.xref_id and x.dbprimary_acc [[COMP]] '[[KEY]]'"
            ],
            [   $db, 'Gene',
                "select count( * ) from object_xref as ox, xref as x
        where ox.ensembl_object_type = 'Transcript' and ox.xref_id = x.xref_id and
              x.display_label [[COMP]] '[[KEY]]' and not(x.dbprimary_acc [[COMP]] '[[KEY]]')",
                "SELECT gsi.stable_id, concat( display_label, ' - ', g.description ), '$db', 'transview', 'transcript' from transcript_stable_id as gsi, transcript as g, object_xref as ox, xref as x
        where gsi.transcript_id = ox.ensembl_id and ox.ensembl_object_type = 'Transcript' and gsi.transcript_id = g.transcript_id and
              ox.xref_id = x.xref_id and x.display_label [[COMP]] '[[KEY]]' and
              not(x.dbprimary_acc [[COMP]] '[[KEY]]')"
            ],
            [   $db, 'Gene',
                "select count( * ) from object_xref as ox, xref as x
        where ox.ensembl_object_type = 'Translation' and ox.xref_id = x.xref_id and x.dbprimary_acc [[COMP]] '[[KEY]]'",
                "SELECT gsi.stable_id, concat( display_label ), '$db', 'protview', 'peptide' from translation_stable_id as gsi, object_xref as ox, xref as x
        where gsi.translation_id = ox.ensembl_id and ox.ensembl_object_type = 'Translation' and 
              ox.xref_id = x.xref_id and x.dbprimary_acc [[COMP]] '[[KEY]]'"
            ],
            [   $db, 'Gene',
                "select count( * ) from object_xref as ox, xref as x
        where ox.ensembl_object_type = 'Translation' and ox.xref_id = x.xref_id and
              x.display_label [[COMP]] '[[KEY]]' and not(x.dbprimary_acc [[COMP]] '[[KEY]]')",
                "SELECT gsi.stable_id, concat( display_label ), '$db', 'protview', 'peptide' from translation_stable_id as gsi, object_xref as ox, xref as x
        where gsi.translation_id = ox.ensembl_id and ox.ensembl_object_type = 'Translation' and 
              ox.xref_id = x.xref_id and x.display_label [[COMP]] '[[KEY]]' and
              not(x.dbprimary_acc [[COMP]] '[[KEY]]')"
            ]
        );
    }
    my @rows = ();
    foreach my $result_hash (@{ $self->{_results} }) {
        my $species = $result_hash->{'species'};
        my $result  = $result_hash->{'rows'};
        for my $row (@{$result}) {
            push @rows,
                +{
                'URL' =>
                    "/$species/$row->[3]?db=$row->[2];$row->[4]=$row->[0]",
                'URL_extra' => [
                    'C',
                    'View marker in ContigView',
                    "/$species/contigview?db=$row->[2];$row->[4]=$row->[0]"
                ],
                'idx'     => 'Gene',
                'subtype' => ucfirst($row->[4]),
                'ID'      => $row->[0],
                'desc'    => $row->[1],
                'species' => $species
                };
        }
    }
    $self->{'results'}{'Gene'} = [ \@rows, $self->{_result_count} ];
}

## Result hash contains the following fields...
##
## { 'URL' => ?, 'type' => ?, 'ID' => ?, 'desc' => ?, 'idx' => ?, 'species' => ?, 'subtype' =>, 'URL_extra' => [] }
1;
