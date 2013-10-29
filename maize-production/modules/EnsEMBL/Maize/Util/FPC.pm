package EnsEMBL::Maize::Util::FPC;

use strict;

use Readonly;
use EnsEMBL::Maize::Util::FileCache;
use Memoize;

use constant FEATURE_ATTRIBUTE_ACCESSION => 'embl_acc';

use constant FPC_SPECIES          => 'Zea_mays';
use constant BAC_SPECIES          => 'Zea_mays2';
use constant BAC_SPECIES_EXTERNAL => 'Zea_mays_external';

use constant FPC_CHROMOSOME_CACHE_SIZE => 100;

memoize('extract_version_from');
memoize('latest_version_of');
memoize('latest_fpc_version_of');
memoize('archive_for_accession');
memoize('_fetch_fpc_feature_by_accession');

Readonly my $CACHE_NAME => '_fpc_cache';

sub new {
    my $class = shift;
    my $self = bless {}, $class;

    return $self;
}

sub fetch_analyzed_bacs {
    my $self = shift;
    my ($object) = @_;

    my $bac_cache = $self->_get_cache($object, 'analyzed_bacs');
    if (defined $bac_cache) {
        return $bac_cache;
    }
    my $bac_list = $self->fetch_all_bac_clones();

    if (defined($object) && $object->can('chr_name')) {
        my %clone_features = map { $_->{'accession'} }
            @{ $self->fetch_all_fpc_clones($object) };
        my @filtered_bacs
            = grep { exists $clone_features{ $_->seq_region_name } }
            @$bac_list;

        $bac_list = \@filtered_bacs;
    }

    my $analyzed_bacs = [
        map {
            my $name   = $_->seq_region_name;
            my $length = $_->seq_region_length;
            +{
                'name'  => "$name [$length]",
                'value' => $name,
            }
        } @$bac_list
    ];
    $self->_add_to_cache($object, 'analyzed_bacs', $analyzed_bacs);
    return $analyzed_bacs;
}

sub fetch_all_bac_clones {
    my $self   = shift;
    my $clones = [];
    for my $species (BAC_SPECIES, BAC_SPECIES_EXTERNAL) {
        my $clone_adaptor = $self->_slice_adaptor(BAC_SPECIES);
        push @$clones, @{ $clone_adaptor->fetch_all('clone') };
    }
    return $clones;
}

sub fetch_all_fpc_clones {
    my $self = shift;
    my ($chromosome) = @_;

    my $cached_clones = $self->_get_cache($chromosome, 'all_clones');
    if (defined $cached_clones) {
        return $cached_clones;
    }
    my $clones = [
        map {
            +{  'clone' => $_->display_id,
                'accession' =>
                    $_->get_scalar_attribute(FEATURE_ATTRIBUTE_ACCESSION),
                }
            } @{ $chromosome->Obj->get_all_MiscFeatures('bac_map') }
    ];
    $self->_add_to_cache($chromosome, 'all_clones', $clones);
    return $clones;
}

sub fetch_accessioned_bacs {
    my $self            = shift;
    my ($chromosome)    = @_;
    my $cached_acc_bacs = $self->_get_cache($chromosome, 'accessioned_bacs');
    if (defined $cached_acc_bacs) {
        return $cached_acc_bacs;
    }
    my $acc_bacs = [ grep { $_->{'accession'} ne '' }
            @{ $self->fetch_all_fpc_clones($chromosome) } ];
    $self->_add_to_cache($chromosome, 'accessioned_bacs', $acc_bacs);
    return $acc_bacs;
}

sub fetch_contigs {
    my $self           = shift;
    my ($object)       = @_;
    my $cached_contigs = $self->_get_cache($object, 'contigs');
    if (defined $cached_contigs) {
        return $cached_contigs;
    }
    my $slice_adaptor = $self->_slice_adaptor(FPC_SPECIES);
    my $contigs       = undef;
    if (defined $object && $object->can('chr_name')) {
        my $chromosome_slice = $slice_adaptor->fetch_by_region('chromosome',
            $object->chr_name);
        $contigs
            = [ map { $_->to_Slice } @{ $chromosome_slice->project('fpc') } ];
    } else {
        $contigs = $slice_adaptor->fetch_all('fpc');
    }

    my $contigs = [
        map      { $_->[1] }
            sort { $a->[0] <=> $b->[0] }
            map  { [ substr($_, 3), $_ ] }
            map  { $_->seq_region_name } @$contigs
    ];
    $self->_add_to_cache($object, 'contigs', $contigs);
    return $contigs;
}

sub fetch_clone_by_accession {
    my $self = shift;
    my ($accession) = @_;
    return $self->_fetch_fpc_feature_by_accession($accession, 'name');

}

sub fetch_contig_by_accession {
    my $self = shift;
    my ($accession) = @_;

    my $value
        = $self->_fetch_fpc_feature_by_accession($accession, 'superctg');
    return $value;

}

sub fetch_chromosome_by_accession {
    my $self = shift;
    my ($accession) = @_;
    my $chromosome_name = "UNKNOWN";
    eval {
        my $contig        = $self->fetch_contig_by_accession($accession);
        my $slice_adaptor = $self->_slice_adaptor(FPC_SPECIES);
        my $contig_slice
            = $slice_adaptor->fetch_by_region('toplevel', $contig);
        my $chromosome = $contig_slice->project('chromosome')->[0]->to_Slice;
        $chromosome_name = $chromosome->seq_region_name;
    };
    if ($@) {
        warn "Unable to fetch chromosome for Accession $accession: $@";
    }
    return $chromosome_name;
}

sub _fetch_fpc_feature_by_accession {
    my $self = shift;
    my ($accession, $attribute_type) = @_;

    my $attribute_value = undef;

    eval {
        my $latest = $self->latest_fpc_version_of($accession);
        my $mfadaptor = $self->_misc_feature_adaptor(FPC_SPECIES);
        my $features
            = $mfadaptor->fetch_all_by_attribute_type_value('embl_acc',
            $latest);

        # expecting only a single feature
        my $feature = $features->[0];
        return undef unless (defined $feature);

        # get its contig
        my $attribute = $feature->get_all_Attributes($attribute_type);

        $attribute_value = $attribute->[0]->value();
    };
    if ($@) {
        warn "Unable to fetch feature for Accession $accession: $@";
        $attribute_value = '';
    }
    return $attribute_value;
}

sub fetch_corebins {
    my $self = shift;
    my ($object) = @_;

    my $cached_bins = $self->_get_cache($object, 'virtualbins');
    if (defined $cached_bins) {
        return $cached_bins;
    }
    my $slicer = $self->_slice_adaptor(FPC_SPECIES);
    my @slices = @{ $slicer->fetch_all('chromosome') };

    my $mfadaptor = $self->_misc_feature_adaptor(FPC_SPECIES);

    my %virtual_bins = ();
    for my $slice (@slices) {
        my $virtual_bins = $mfadaptor->fetch_all_by_Slice_and_set_code($slice,
            'core_bins');
        for my $vb (@$virtual_bins) {
            my $id    = $vb->display_id();
            my $start = $vb->seq_region_start();
            my $end   = $vb->seq_region_end();

            my $range = $start . "\t" . $end;

            $virtual_bins{$id} = $range;
        }
    }
    $self->_add_to_cache($object, 'virtualbins', \%virtual_bins);

    return (\%virtual_bins);
}

sub fetch_corebinmarkers {
    my $self = shift;
    my ($object) = @_;

    my $cached_cbms = $self->_get_cache($object, 'core_bin_markers');
    if (defined($cached_cbms)) {
        return $cached_cbms;
    }
    my $virtualbins = $self->fetch_corebins($object);
    my %virtualbins = %$virtualbins;

    my $corebinmarker_slice_adaptor = $self->_slice_adaptor(FPC_SPECIES);

    my @slices = @{ $corebinmarker_slice_adaptor->fetch_all('chromosome') };

    my @corebinmarkers;

    #    my %corebinmarkers;
    for my $slice (@slices) {
        my $corebinmarkers
            = $slice->get_all_MarkerFeatures('core_bin_marker');

        for my $cbm (@$corebinmarkers) {
            push(@corebinmarkers, $cbm->display_id());
        }
    }

    # alpha order for now. order by chrom??
    my @sorted_corebinmarkers = sort (@corebinmarkers);

    my $corebinmarkers = \@sorted_corebinmarkers;

    $self->_add_to_cache($object, 'core_bin_markers', $corebinmarkers);
    return $corebinmarkers;
}

sub latest_version_of {
    my $self = shift;
    my ($accession) = @_;

    $accession =~ s/\.\d+$//;
    for my $species (BAC_SPECIES, BAC_SPECIES_EXTERNAL) {
        my $slice_adaptor = $self->_slice_adaptor($species);
        my $slice = $slice_adaptor->fetch_by_region('clone', $accession);
        if (defined $slice) {
            return $slice->seq_region_name;
        }
    }
}

=pod

=head2 latest_fpc_version_of
    Retrieves the latest version of the given BAC in the FPC database

=cut
sub latest_fpc_version_of {
    my $self = shift;
    my ($accession) = @_;
    $accession =~ s/\.\d+$//;
    my $adaptor = $self->_misc_feature_adaptor(FPC_SPECIES);
    my $statement = $adaptor->prepare(<<SQL);
select misc_attrib.value
  from misc_attrib 
  left join attrib_type using (attrib_type_id)
 where attrib_type.code = 'embl_acc'
   and misc_attrib.value like ?
SQL
    $statement->bind_param(1, "$accession\%");
    $statement->execute();
    my $latest
        = $statement->fetchrow_arrayref()->[0];
    return $latest;
}


=pod

=head2 archive_for_accession
    Returns a list of other versions of this accession

=cut

sub archive_for_accession {
    my $self = shift;
    my ($versioned_accession) = @_;

    my ($accession, $version) = $versioned_accession =~ m/(\w+)(\.\d+)?$/;
    my $slice_adaptor = $self->_slice_adaptor(BAC_SPECIES);
    my $statement     = $slice_adaptor->prepare(<<SQL);
select seq_region.name
  from seq_region
       left join coord_system
              on seq_region.coord_system_id = coord_system.coord_system_id
 where seq_region.name like ?
   and coord_system.name = 'clone'
SQL
    $statement->bind_param(1, "$accession.\%");
    $statement->execute();
    my $accessions = [
        sort {
            $self->extract_version_from($b)
                <=> $self->extract_version_from($a)
            }
            grep {
            $_ ne $versioned_accession
            }
            map {
            $_->[0]
            } @{ $statement->fetchall_arrayref() }
    ];

    return $accessions;
}

=pod

=head2 extract_version_from
    Returns the version for a given accession

=cut

sub extract_version_from {
    my $self = shift;
    my ($accession) = @_;
    my ($version) = $accession =~ m/^\w+\.(\d+)$/;
    return $version;
}

sub _slice_adaptor {
    my $self = shift;
    my ($species) = @_;
    return Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Slice');
}

sub _marker_feature_adaptor {
    my $self = shift;
    my ($species) = @_;
    return Bio::EnsEMBL::Registry->get_adaptor($species, 'core',
        'MarkerFeature');
}

sub _misc_feature_adaptor {
    my $self = shift;
    my ($species) = @_;
    return Bio::EnsEMBL::Registry->get_adaptor($species, 'core',
        'MiscFeature');
}

sub _get_cache {
    my $self = shift;
    my ($chromosome, $cache_key) = @_;
    my $slice_name = $self->_parse_cache_slice_name($chromosome);

    my $chromosome_cache = $self->_cache->get_value($slice_name);
    if (!defined $chromosome_cache) {
        return undef;
    }
    my $cached_value = $chromosome_cache->{$cache_key};

    return $cached_value;
}

sub _add_to_cache {
    my $self = shift;
    my ($chromosome, $cache_key, $value) = @_;
    my $slice_name  = $self->_parse_cache_slice_name($chromosome);
    my $cache       = $self->_cache;
    my $slice_cache = $cache->get_value($slice_name);
    if (!defined $slice_cache) {
        $slice_cache = +{};
        $cache->set_value($slice_name, $slice_cache);
    }
    $slice_cache->{$cache_key} = $value;
    $cache->save;
}

sub _parse_cache_slice_name {
    my $self = shift;
    my ($chromosome) = @_;
    if (defined($chromosome) && $chromosome->can('chr_name')) {
        return $chromosome->chr_name;
    } else {
        return 'genome';
    }
}

sub _print_cache {
    my $self = shift;
    my (@messages) = @_;
    $self->_cache->print(@messages);
}

sub _cache {
    my $self = shift;
    if (!defined $self->{$CACHE_NAME}) {
        $self->{$CACHE_NAME} = EnsEMBL::Maize::Util::FileCache->new(
            { 'filename' => 'genome_entry_points.txt', });
    }
    return $self->{$CACHE_NAME};
}

1;
