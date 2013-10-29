# BioPerl module for DAS
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# Copyright Tony Cox
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

DAS - DESCRIPTION of Object

=head1 SYNOPSIS

use Data::Dumper;
use Bio::EnsEMBL::ExternalData::DAS::DASAdaptor;
use Bio::EnsEMBL::ExternalData::DAS::DAS;

$das_adaptor = Bio::EnsEMBL::ExternalData::DAS::DASdaptor->new(
                                             -url   => 'some_server',
                                             -dsn   => 'twiddly-bits',
                                             -ensdb => $ensembl_dbh,
                                            );

my $ext_das = Bio::EnsEMBL::ExternalData::DAS::DAS->new($das_adaptor)

$dbobj->add_ExternalFeatureFactory($ext_das);

This class implements only contig based method:

$dbobj->get_Ensembl_SeqFeatures_contig('AL035659.00001');

Also
my @features = $ext_das->fetch_SeqFeature_by_contig_id("AL035659.00001");

Method get_Ensembl_SeqFeatures_clone returns an empty list.

=head1 DESCRIPTION

Interface to an external DAS data source - lovelingly mangled into the Ensembl database
adaptor scheme.

interface for creating L<Bio::EnsEMBL::ExternalData::DAS::DAS.pm>
objects from an external DAS database. 

The objects returned in a list are
L<Bio::EnsEMBL::ExternalData::DAS::DASSeqFeature> objects which might possibly contain
L<Bio::Annotation::DBLink> objects to give unique IDs in various
DAS databases.

=head1 CONTACT

 Tony Cox <avc@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::ExternalData::DAS::DAS;

use strict;
use vars qw(@ISA);
use Bio::Das;
use Bio::EnsEMBL::Root;
use Data::Dumper;

use Bio::EnsEMBL::ExternalData::DAS::DASSeqFeature;

# Object preamble
@ISA = qw(Bio::EnsEMBL::Root);

sub new {
    my ($class, $adaptor) = @_;
    my $self;
    $self = {};
    bless $self, $class;

    $self->adaptor($adaptor);

    return $self;    # success - we hope!
}

#----------------------------------------------------------------------

=head2 adaptor

  Arg [1]   : Bio::EnsEMBL::ExternalData::DAS::DASAdaptor (optional)
  Function  : getter/setter for adaptor attribute
  Returntype: Bio::EnsEMBL::ExternalData::DAS::DASAdaptor
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub adaptor {
    my $key  = '_adaptor';
    my $self = shift;
    if (@_) { $self->{$key} = shift }
    return $self->{$key};
}

=head2 fetch_dsn_info

  Arg [1]   : none
  Function  : Retrieves a list of DSN objects from registered URL
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub fetch_dsn_info {
    my $self = shift;

    my @sources  = ();
    my $callback = sub {
        my $obj = shift;
        $obj->isa('Bio::Das::DSN') || return;
        my $data = {};
        $data->{url}         = $obj->url;
        $data->{base}        = $obj->base;
        $data->{id}          = $obj->id;
        $data->{dsn}         = $obj->id;
        $data->{name}        = $obj->name || $obj->id;
        $data->{description} = $obj->description;
        $data->{master}      = $obj->master;
        push @sources, $data;
    };
    my $dsn      = $self->adaptor->url;
    my $das      = $self->adaptor->_db_handle;
    my $dsn_hash = $das->dsns();
    foreach my $key (%{$dsn_hash}) {
        foreach my $obj (@{ $dsn_hash->{$key} }) {
            my $data = {};

            if ($dsn =~ m!(.+/das)/([^/]+)!) {
                $data->{base} = $1;
                $data->{url}  = $dsn;
            } else {
                $data->{base} = $dsn;
                $data->{url}  = "$dsn/$obj->{source_id}";
            }
            $data->{id}          = $obj->{source_id};
            $data->{dsn}         = $obj->{source_id};
            $data->{name}        = $obj->{source} || $obj->{source_id};
            $data->{description} = $obj->{description};
            $data->{master}      = $obj->{mapmaster};
            push @sources, $data;
        }
    }
    return [@sources];
}

=head2  fetch_all_by_Slice

  Arg[1]  : Slice 
  Example : $features = $adaptor->fetch_all_by_Slice($slice);
  Description : fetches DAS features for this DAS Adaptor and maps to 
                slice coordinates
  ReturnType: arrayref of Bio::Ensembl::External::DAS::DASSeqFeature
  Exceptions: ?
  Caller  : mainly Slice.pm

=cut

sub fetch_all_by_Slice {
    my ($self, $slice) = @_;

    # Examine cache
    my $CACHE_KEY = $slice->name;
    if ($self->{$CACHE_KEY}) {
        return (
            $self->{$CACHE_KEY},
            $self->{"_stylesheet_$CACHE_KEY"},
            $self->{"_segments_$CACHE_KEY"}
        );
    }

    # Get all coord systems this Ensembl DB knows about
    my $csa  = $slice->coord_system->adaptor;
    my $csa2 = $slice->adaptor()->db()->get_CoordSystemAdaptor();

# The following bit has been put to investigate why we get error messages in the error log saying that fetch_all can not be called on an undefined value
    if (!defined($csa)) {
        my @ca = caller(2);
        warn(
            "WARNING: Could not get a coord system adaptor for slice [$slice]\n @ca"
        );

        if (!defined($csa2)) {
            warn("CSA2 is empty");
            return [];
        } else {
            $csa = $csa2;
        }
    }
    my %coord_systems = map { $_->name, $_ } @{ $csa->fetch_all || [] };

    # Get the slice representation for each coord system.
    my @segments_to_request;    # The DAS segments to query
    my %slice_by_segment;       # tally of which slice belongs to segment
    foreach my $system (keys %coord_systems) {
        foreach my $segment (@{ $slice->project($system) || [] }) {
            my $slice        = $segment->to_Slice;
            my $slice_name   = $slice->name;
            my $slice_start  = $slice->start;
            my $slice_end    = $slice->end;
            my $region_name  = $slice->seq_region_name;
            my $coord_system = $slice->coord_system;

            if ($slice_name =~ /^clone/)
            {                   # Clone-specific hack for embl versions
                my ($id, $version) = split(/\./, $region_name);
                if ($version) {
                    push(@segments_to_request, "$id:$slice_start,$slice_end");
                    $slice_by_segment{$id} = $slice;
                } else {
                    my ($version_attrib)
                        = @{ $slice->get_all_Attributes('acc-version') };
                    print STDERR "[DAS] Getting Version for $region_name\n";

                    if ($version_attrib) {
                        my $versioned_acc = "$id." . $version_attrib->value;
                        print STDERR "Versioned ACC: $version_attrib\n";

                        push(@segments_to_request,
                            "$versioned_acc:$slice_start,$slice_end");
                        $slice_by_segment{$versioned_acc} = $slice;
                    }
                }
            }
            push(@segments_to_request,
                "$region_name:$slice_start,$slice_end");
            $slice_by_segment{$region_name} = $slice;
        }
    }

    #warn("SEGMENTS : ".scalar(@segments_to_request));
    # Run the DAS query
    my ($features, $style)
        = $self->get_Ensembl_SeqFeatures_DAS([@segments_to_request]);

    # Map the DAS results into the coord system of the original slice
    my @result_list;
    foreach my $das_sf (@$features) {
        my $segment = $das_sf->das_segment
            || (warn("No das_segment for $das_sf") && next);
        my $das_slice = $slice_by_segment{ $segment->ref }
            || (warn("No Slice for ", $segment->ref) && next);
        $self->_map_DASSeqFeature_to_slice($das_sf, $das_slice, $slice)
            && push @result_list, $das_sf;
    }

    # Return the mapped features
    warn "RETURNING FEATURES ... STYLES ... and SEGMENTS";
    return (
        ($self->{ $slice->name }                  = \@result_list),
        ($self->{ "_stylesheet_" . $slice->name } = $style),
        ($self->{ "_segments_" . $slice->name }   = \@segments_to_request)
    );
}

sub fetch_all_Features {
    my ($self, $slice, $source_type) = @_;

    # Examine cache
    my $CACHE_KEY = $slice->name;
    if ($self->{$CACHE_KEY}) {
        return (
            $self->{$CACHE_KEY},
            $self->{"_stylesheet_$CACHE_KEY"},
            $self->{"_segments_$CACHE_KEY"}
        );
    }

    if ($source_type =~ /^ensembl_location(.+)?/) {
        my %coord_systems;
        if (defined(my $cs = $1)) {
            $cs =~ s/^_//;
            $coord_systems{$cs} = 1;
        } else {

            # Get all coord systems this Ensembl DB knows about
            my $csa = $slice->coord_system->adaptor;
            if (!defined($csa)) {
                my @ca = caller(2);
                warn(
                    "WARNING: Could not get a coord system adaptor for slice [$slice]\n @ca"
                );
                my $csa2 = $slice->adaptor()->db()->get_CoordSystemAdaptor();
                if (!defined($csa2)) {
                    warn("CSA2 is empty");
                    return [];
                } else {
                    $csa = $csa2;
                }
            }
            %coord_systems = map { $_->name, $_ } @{ $csa->fetch_all || [] };
        }

        #      warn("CS:".join('*', sort keys %coord_systems));

        # Get the slice representation for each coord system.
        my @segments_to_request;    # The DAS segments to query
        my %slice_by_segment;       # tally of which slice belongs to segment
        foreach my $system (keys %coord_systems) {
            foreach my $segment (@{ $slice->project($system) || [] }) {
                my $slice        = $segment->to_Slice;
                my $slice_name   = $slice->name;
                my $slice_start  = $slice->start;
                my $slice_end    = $slice->end;
                my $region_name  = $slice->seq_region_name;
                my $coord_system = $slice->coord_system;
                if ($slice_name =~ /^clone/)
                {                   # Clone-specific hack for embl versions
                    my ($id, $version) = split(/\./, $region_name);
                    if ($version) {
                        push(@segments_to_request,
                            "$id:$slice_start,$slice_end");
                        $slice_by_segment{$id} = $slice;
                    } else {
                        my ($version_attrib)
                            = @{ $slice->get_all_Attributes('acc-version') };

                        if ($version_attrib) {
                            my $versioned_acc
                                = "$id." . $version_attrib->value;
                            push(@segments_to_request,
                                "$versioned_acc:$slice_start,$slice_end");
                            $slice_by_segment{$versioned_acc} = $slice;
                        }
                    }
                }
                push(@segments_to_request,
                    "$region_name:$slice_start,$slice_end");
                $slice_by_segment{$region_name} = $slice;
            }
        }

        # Run the DAS query
        my ($features, $style)
            = $self->get_Ensembl_SeqFeatures_DAS([@segments_to_request]);
        if (@$features && $features->[0]->das_type eq '__ERROR__') {
            return ($self->{ $slice->name } = $features,
                [], \@segments_to_request);
        }

        # Map the DAS results into the coord system of the original slice
        my @result_list;
        foreach my $das_sf (@$features) {
            my $segment = $das_sf->das_segment
                || (warn("No das_segment for $das_sf") && next);
            my $das_slice = $slice_by_segment{ $segment->ref }
                || (warn("No Slice for ", $segment->ref) && next);
            $self->_map_DASSeqFeature_to_slice($das_sf, $das_slice, $slice)
                && push @result_list, $das_sf;
        }

        # Return the mapped features
        return (
            ($self->{ $slice->name }                  = \@result_list),
            ($self->{ "_stylesheet_" . $slice->name } = $style),
            ($self->{ "_segments_" . $slice->name }   = \@segments_to_request)
        );
    }
}

#----------------------------------------------------------------------

sub _map_DASSeqFeature_to_pep {
    my $self   = shift;
    my $dblink = shift || die("Need a DBLink object");
    my $dsf    = shift || die("Need a DASSeqFeature object");

    if (!ref($dblink)) { return 1 }    # Ensembl id_type - mapping not needed

    # Check for 'global' feature - mapping not needed
    if (   $dsf->das_feature_id eq $dsf->das_segment->ref
        or !$dsf->das_start
        or !$dsf->das_end)
    {
        $dsf->start(0);
        $dsf->end(0);
        return 1;
    }

    # Check that dblink is map-able
    if (!$dblink->can('get_mapper')) { return 0 }

    # Map
    my @coords = ();
    eval { @coords = $dblink->map_feature($dsf) };
    if ($@) { warn($@) }

    @coords = grep {
        $_->isa('Bio::EnsEMBL::Mapper::Coordinate')
            || $_->isa('Bio::EnsEMBL::Mapper::Gap')
    } @coords;
    @coords || return 0;
    $dsf->start($coords[0]->start);
    $dsf->end($coords[-1]->end);

    #  warn( "Ensembl:".$dsf->start."-".$dsf->end );
    return 1;
}

#----------------------------------------------------------------------

=head2 _map_DASSeqFeature_to_slice

  Arg [1]   : DASSeqFeature object
  Arg [2]   : Slice with CoordSystem and seq_region_name for DASSeqFeature
  Arg [3]   : Slice with offsets and CoordSystem to map DASSeqFeature to
  Function  : Maps DASSeqFeature in one CoordSystem to Slice coords in 
              another CoordSystem
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub _map_DASSeqFeature_to_slice {
    my $self      = shift;
    my $das_sf    = shift;
    my $das_slice = shift;
    my $usr_slice = shift;

    my $fr_csystem = $das_slice->coord_system;
    my $to_csystem = $usr_slice->coord_system;

    my $db = $usr_slice->adaptor->db;
    my $ma = $db->get_AssemblyMapperAdaptor;

    # Map
    my ($slice_start, $slice_end, $slice_strand);
    unless ($fr_csystem->equals($to_csystem)) {
        my $mapper = $ma->fetch_by_CoordSystems($fr_csystem, $to_csystem);
        my @coords = ();

        eval {
            @coords
                = $mapper->map($das_slice->seq_region_name,
                $das_sf->das_start, $das_sf->das_end,
                $das_sf->das_orientation, $fr_csystem);
        };
        if ($@) { warn($@) }
        @coords
            = grep { $_->isa('Bio::EnsEMBL::Mapper::Coordinate') } @coords;
        scalar(@coords) || return 0;
        $slice_start  = $coords[0]->start - $usr_slice->start + 1;
        $slice_end    = $coords[-1]->end - $usr_slice->start + 1;
        $slice_strand = $coords[0]->strand;
    } else {    # No mapping needed
        $slice_start  = $das_sf->das_start - $usr_slice->start + 1;
        $slice_end    = $das_sf->das_end - $usr_slice->start + 1;
        $slice_strand = $das_sf->das_orientation;
    }

    $das_sf->seqname($usr_slice->seq_region_name);
    $das_sf->start($slice_start);
    $das_sf->end($slice_end);
    $das_sf->strand($slice_strand);

    return 1;
}

=head2 get_Ensembl_SeqFeatures_clone

 Title   : get_Ensembl_SeqFeatures_clone (not used)
 Function:
 Example :
 Returns :
 Args    :

=cut

sub get_Ensembl_SeqFeatures_clone {
    my ($self, $contig) = @_;
    $self->throw("get_Ensembl_SeqFeatures_clone is unimplemented!");
    my @features = ();
    return (@features);
}

=head2 fetch_SeqFeature_by_contig_id

 Title   : fetch_SeqFeature_by_contig_id
 Usage   : $obj->fetch_SeqFeature_by_contig_id("Contig_X")
 Function: return DAS features for a contig
 Returns : 
 Args    : none


=cut

sub fetch_SeqFeature_by_contig_id {
    my ($self, $contig) = @_;
    $self->throw("fetch_SeqFeature_by_contig_id is unimplemented!");
    my @features = ();
    return (@features);
}

=head2 forwarded_for

 Title   : forwarded_for
 Usage   : $obj->forwarded_for($ENV{'HTTP_X_FORWARDED_FOR'})
 Function: store a DAS data source URL
 Returns : 
 Args    : none


=cut

sub forwarded_for {
    my ($self, $value) = @_;
    if (defined $value) {
        $self->{'_forwarded_for'} = $value;
    }
    return $self->{'_forwarded_for'};
}

=head2 _db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function:
 Example :
 Returns : value of _db_handle
 Args    : newvalue (optional)

=cut

sub _db_handle {
    my $caller = join(", ", (caller(0))[ 1 .. 2 ]);
    warn
        "\033[31m DEPRECATED use adaptor->_db_handle instead: \033[0m $caller";
    my $self = shift;
    return $self->adaptor->_db_handle(@_);
}

=head2 _types

 Title   :
 Usage   : DEPRECATED
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _types {
    my $caller = join(", ", (caller(0))[ 1 .. 2 ]);
    warn "\033[31m DEPRECATED use adaptor->types instead: \033[0m $caller";
    my $self = shift;
    return $self->adaptor->types(@_);
}

=head2 _dsn

 Title   :
 Usage   : DEPRECATED
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _dsn {
    my $caller = join(", ", (caller(0))[ 1 .. 2 ]);
    warn "\033[31m DEPRECATED use adaptor->dsn instead: \033[0m $caller";
    my $self = shift;
    return $self->adaptor->dsn(@_);
}

=head2 _url

 Title   : 
 Usage   : DEPRECATED
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _url {
    my $caller = join(", ", (caller(0))[ 1 .. 2 ]);
    warn "\033[31m DEPRECATED use adaptor->url instead: \033[0m $caller";
    my $self = shift;
    return $self->adaptor->url(@_);
}

=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub DESTROY {
    my ($obj) = @_;
    $obj->adaptor(undef());
    if ($obj->{'_db_handle'}) {
        $obj->{'_db_handle'} = undef;
    }
}

sub fetch_all_by_ID {
    my $self       = shift;
    my $parent_obj = shift;
    my $data_obj   = shift;

    my $id_type_base = $self->adaptor->type || 'swissprot';
    my $url          = $self->adaptor->url;
    my $dsn          = $self->adaptor->dsn;

    $parent_obj->can('get_all_DBLinks')
        || $self->throw(
        "Need a Bio::EnsEMBL obj (eg Translation) that can get_all_DBLinks");

    my $ensembl_id
        = $parent_obj->stable_id() ? $parent_obj->can('stable_id') : '';

    my %ids = ();

    my @id_types
        = $id_type_base eq 'mixed'
        ? @{ $self->adaptor->mapping }
        : ($id_type_base);
    foreach my $id_type (@id_types) {

        # If $id_type is prefixed with 'ensembl_', then ensembl id type
        if ($id_type =~ m/ensembl_(.+)/o) {
            my $type = $1;
            my @gene_ids;
            my @tscr_ids;
            my @tran_ids;
            if ($parent_obj->isa("Bio::EnsEMBL::Gene")) {
                push(@gene_ids, $parent_obj->stable_id);
                foreach my $tscr (@{ $parent_obj->get_all_Transcripts }) {
                    push(@tscr_ids, $tscr->stable_id);
                    my $tran = $tscr->translation || next;
                    push(@tran_ids, $tran->stable_id);
                }
            } elsif ($parent_obj->isa("Bio::EnsEMBL::Transcript")) {
                push(@tscr_ids, $parent_obj->stable_id);
                my $tran = $parent_obj->translation || next;
                push(@tran_ids, $tran->stable_id);
                push(@gene_ids, $data_obj->gene->stable_id);
            } elsif ($parent_obj->isa("Bio::EnsEMBL::Translation")) {
                push(@tran_ids, $parent_obj->stable_id);

                # if the source is ensembl_gene - get gene stable id
                if (defined(
                        my $gene
                            = $parent_obj->adaptor->db->get_GeneAdaptor
                            ->fetch_by_translation_stable_id(
                            $parent_obj->stable_id
                            )
                    )
                    )
                {
                    push(@gene_ids, $gene->stable_id);
                }
            } else {    # Assume protein
                warn("??? - ",
                    $parent_obj->transcript->translation->stable_id);
                push(@tran_ids,
                    $parent_obj->transcript->translation->stable_id);
            }
            if ($type eq 'gene') {
                map { $ids{$_} = 'gene' } @gene_ids;
            } elsif ($type eq 'transcript') {
                map { $ids{$_} = 'transcript' } @tscr_ids;
            } elsif ($type eq 'peptide') {
                map { $ids{$_} = 'peptide' } @tran_ids;
            }
        } elsif ($id_type eq 'mgi') {

            # MGI Accession IDs come from MarkerSymbol DB
            my $id_method = 'primary_id';
            foreach my $xref (grep { lc($_->dbname) eq 'markersymbol' }
                @{ $parent_obj->get_all_DBLinks })
            {
                my $id = $xref->primary_id || next;
                $id =~ s/\://g;
                $ids{$id} = $xref;
            }
        } else {

            # If no 'ensembl_' prefix, then DBLink ID
            # If $id_type is suffixed with '_acc', use primary_id call
            # rather than display_id
            my $id_method
                = $id_type =~ s/_acc$// ? 'primary_id' : 'display_id';
            foreach my $xref (@{ $parent_obj->get_all_DBLinks }) {
                lc($xref->dbname) ne lc($id_type) and next;
                my $id = $xref->$id_method || next;
                $ids{$id} = $xref;
            }
        }
    }

    # Return empty if no ids found
    return () if (!scalar keys(%ids));
    my $response;

    my @das_features = ();

    # Get features
    my @req;
    my $types = $self->adaptor->types || [];

    if (@$types) {
        foreach my $s (keys %ids) {
            my $rhash = { 'segment' => join(',', keys %ids) };

            if (my $maxbins = $self->adaptor->maxbins()) {
                $rhash->{'maxbins'} = $maxbins;
            }

            if (@$types) {
                $rhash->{'types'} = join ',', @$types;
            }
            push @req, $rhash;
        }

        $response = $self->adaptor->_db_handle->features(\@req);
    } else {
        $response = $self->adaptor->_db_handle->features([ keys %ids ]);
    }

    #   warn Data::Dumper::Dumper($response);
    foreach my $url (keys %$response) {
        foreach my $f (
            ref($response->{$url}) eq "ARRAY" ? @{ $response->{$url} } : ())
        {
            $self->_add_feature($f, $dsn, \@das_features);
        }
    }

    my @result_list = grep {
        $self->_map_DASSeqFeature_to_pep($ids{ $_->das_segment->ref }, $_)
            == 1
    } @das_features;

    my $STYLES   = $self->_get_stylesheet();
    my @segments = keys %ids;

    return (\@result_list, $STYLES, \@segments);
}

=head2 get_Ensembl_SeqFeatures_DAS

 Title   : get_Ensembl_SeqFeatures_DAS ()
 Usage   : get_Ensembl_SeqFeatures_DAS(['AL12345','13']);
 Function:
 Example :
 Returns :
 Args    :
 Notes   : This function sets the primary tag and source tag fields in the
           features so that higher level code can filter them by their type
           (das) and their data source name (dsn)

=cut

sub get_Ensembl_SeqFeatures_DAS {
    my $self         = shift;
    my $segments     = shift || [];
    my $dbh          = $self->adaptor->_db_handle();
    my $dsn          = $self->adaptor->dsn();
    my $types        = $self->adaptor->types() || [];
    my @das_features = ();
    @$segments || $self->throw("Need some segment IDs to query against");

    if (defined(my $error = $self->adaptor->verify)) {
        my $f = {
            'type'       => '__ERROR__',
            'type_id'    => '__ERROR__',
            'feature_id' => $error,
            'segment_id' => 1,
            'start'      => 1,
            'end'        => 1,
        };
        $self->_add_feature($f, $dsn, \@das_features);
        return (\@das_features, [], []);
    }

    # Get features
    my $response;
    my @req;
    foreach my $s (@$segments) {
        my $rhash = { 'segment' => $s };

        if (my $maxbins = $self->adaptor->maxbins()) {
            $rhash->{'maxbins'} = $maxbins;
        }

        if (@$types) {
            $rhash->{'types'} = join ',', @$types;
        }
        push @req, $rhash;
    }
    $response = $dbh->features(\@req);

  #$Data::Dumper::Indent = 3;
  #warn Data::Dumper::Dumper($response);
  #warn Data::Dumper::Dumper(\@req);
  #  if(@$types) {
  #    $response = $dbh->features({'segment' => $segments, 'type' => $types});
  #  } else {
  #    $response = $dbh->features($segments);
  #  }

# Parse the response. There is a problem using callbacks hence the explicit response handling
    foreach my $url (keys %$response) {
        foreach my $f (
            ref($response->{$url}) eq "ARRAY" ? @{ $response->{$url} } : ())
        {

            #warn Data::Dumper::Dumper($f);

            $self->_add_feature($f, $dsn, \@das_features);
        }
    }

    my $STYLES = $self->_get_stylesheet();
    return (\@das_features, $STYLES, $segments);
}

sub _get_stylesheet {
    my ($self) = @_;

    my $dbh = $self->adaptor->_db_handle();

    my $STYLES   = [];
    my $response = $dbh->stylesheet();

    foreach my $url (keys %$response) {
        foreach my $css (@{ $response->{$url} }) {
            my @categories = @{ $css->{category} || [] };
            foreach my $c (@categories) {
                my $c_id = $c->{category_id};
                my @types = @{ $c->{type} || [] };
                foreach my $t (@types) {
                    my $t_id   = $t->{type_id};
                    my @glyphs = $t->{glyph};
                    foreach my $g (@glyphs) {
                        my %ghash = %{ $g->[0] };
                        foreach my $gtype (keys %ghash) {
                            my $attr = $ghash{$gtype}->[0];
                            push @$STYLES,
                                {
                                'category' => $c_id,
                                'type'     => $t_id,
                                'glyph'    => $gtype,
                                'attrs'    => $attr,
                                };
                        }
                    }
                }
            }
        }
    }

    return $STYLES;
}

sub _add_feature {
    my ($self, $f, $dsn, $fa) = @_;

    # Filter out calls for non-features
    defined($f->{feature_id}) or next;

    my ($fstart, $fend) = ($f->{start}, $f->{end});
    if ($f->{type} =~ /^(INIT_MET|INIT_MET:)$/) {
        $fstart = $fend = 1;
    }

    my $das_sf = new Bio::EnsEMBL::ExternalData::DAS::DASSeqFeature;
    $das_sf->das_feature_id($f->{feature_id});
    $das_sf->das_feature_label($f->{feature_label});
    $das_sf->das_segment(
        Bio::Das::Segment->new(
            $f->{segment_id}, $f->{start}, $f->{end}, '1', $self, $dsn
        )
    );

    $das_sf->das_segment_label($f->{segment_id});
    $das_sf->das_id($f->{feature_id});
    $das_sf->das_dsn($dsn);
    $das_sf->source_tag($dsn);
    $das_sf->primary_tag('das');
    $das_sf->das_type($f->{type});
    $das_sf->das_type_id($f->{type_id});
    $das_sf->das_type_category($f->{type_category});
    $das_sf->das_type_reference($f->{type_reference});
    $das_sf->das_name($f->{feature_id});
    $das_sf->das_method($f->{method});
    $das_sf->das_method_id($f->{method_id});
    $das_sf->das_start($f->{start});
    $das_sf->das_end($f->{end});
    $das_sf->das_start($f->{start});
    $das_sf->das_end($f->{end});
    $das_sf->das_score($f->{score});
    $das_sf->das_orientation($f->{orientation} || 0);
    $das_sf->das_phase($f->{phase});
    $das_sf->das_target($f->{target});
    $das_sf->das_target_id($f->{target_id});
    $das_sf->das_target_label($f->{target_label});
    $das_sf->das_target_start($f->{target_start});
    $das_sf->das_target_stop($f->{target_stop});
    $das_sf->das_links(@{ $f->{link} }) if ($f->{link});

    if ($f->{group}) {
        $das_sf->das_groups(@{ $f->{group} });

        # For backward compatability
        $das_sf->das_group_id($f->{group}->[0]->{group_id});
        $das_sf->das_group_label($f->{group}->[0]->{group_label});
        $das_sf->das_group_type($f->{group}->[0]->{group_type});
    }

    my $note
        = ref($f->{note}) eq 'ARRAY'
        ? join('<br/>', @{ $f->{note} })
        : $f->{note};
    $note =~ s![\r\n]!!g;

    #    $das_sf->das_note($note);
    $das_sf->das_note($f->{note});
    push(@$fa, $das_sf);
    return;
}

1;
