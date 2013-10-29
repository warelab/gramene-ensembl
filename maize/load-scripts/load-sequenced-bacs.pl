#!/usr/local/bin/perl

=head1 NAME

load-sequenced-bacs.pl - Loads sequenced BAC data into a core Ensembl DB

=head1 SYNOPSIS

perl load-sequenced-bacs.pl [options]

Options:
 -h --help
 -m --man
 -r --registry_file
 -f --fpc_species
 -b --bac_species
    --dryrun
 -c --conf

=head1 OPTIONS

Migrates information about sequenced clones into the Ensembl DB.

B<--dryrun>
  Do not make any alterations to the database.

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-b --bac_species>
  Use this species entry from the registry file for the BAC database [REQUIRED].
  
B<-f --fpc_species>
  Use this species entry from the registry file for the FPC database [REQUIRED].
  
B<-c --conf>
  Reads all configuration information from file. Command-line arguments override.
  
B<--accession_number>
  Specify an accession number to update.

=head1 DESCRIPTION

B<This program> 

  Loads overgo-clone associations into a core Ensembl DB
  
  Maintained by Shiran Pasternak <shiran@cshl.edu>

=cut

use strict;
use warnings;

use Maize::Clonepath;

use Pod::Usage;

use English qw(-no_match_vars);
use Carp;
use DBI;
use DBI qw(:sql_types);
use FindBin qw($Bin);
use File::Basename qw( dirname );

use Log::Log4perl;

use vars qw($BASEDIR);

BEGIN {

    # Set the perl libraries
    $BASEDIR = dirname($Bin);
    unshift @INC, $BASEDIR . '/ensembl-live/ensembl/modules';
    Log::Log4perl::init($BASEDIR . '/conf/log4perl.conf');
}

use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;

use Bio::Seq;
use Bio::Seq::PrimaryQual;

use DBI;

use Readonly;
use English qw( -no-match-vars );
use List::Util;

Readonly my $DEFAULT_IMPROVED_REGION_LENGTH => 5000;

Readonly my $SEQUENCE_DNA_TYPE     => 'sequence';
Readonly my $SEQUENCE_QUALITY_TYPE => 'qualities';

Readonly my @QUALITY_WINDOW_SIZES => (1, 10, 100, 500);

my (%conf);
my %clone_feature_cache = ();
my $logger              = Log::Log4perl->get_logger('clone');

my $clonepath = undef;

MAIN: {
    $OUTPUT_AUTOFLUSH = 1;

    $clonepath = Maize::Clonepath->new();

    eval { $clonepath->load_configuration(\%conf); };
    if ($EVAL_ERROR) {
        warn "$EVAL_ERROR\n";
        pod2usage(2);
    }
    pod2usage(-verbose => 2) if $conf{'man'};
    pod2usage(1) if $conf{'help'};

    map {
        -e $_ || (warn("File $_ does not exist\n")    && pod2usage(1));
        -r $_ || (warn("Cannot read $_\n")            && pod2usage(1));
        -f $_ || (warn("File $_ is not plain-text\n") && pod2usage(1));
        -s $_ || (warn("File $_ is empty\n")          && pod2usage(1));
    } $conf{'registry_file'};

    my @synch_functions = ();
    if ($conf{'bac_species'}) {
        push(@synch_functions, \&synchronize_bac_clone);
        prepare_bac_database();
    }
    if ($conf{'fpc_species'}) {
        push(@synch_functions, \&synchronize_fpc_clone);
        prepare_fpc_database();
    }

    while (my $clone = $clonepath->next_sequenced_clone()) {
        $logger->info("Synchronizing clone $clone->{versioned_accession}");
        for my $func (@synch_functions) {
            $func->($clone);
        }
    }

    if ($conf{'bac_species'}) {
        postprocess_bac_database();
    }
    if ($conf{'fpc_species'}) {
        postprocess_fpc_database();
    }

    exit 0;
}

#======================================================================

=pod

=head2 prepare_bac_database
    Add meta information in BAC database before clones are added

=cut

sub prepare_bac_database {
    $logger->info("Preparing BAC database");
    ensure_clone_coord_system();
    ensure_contig_coord_system();
}

=pod

=head2 postprocess_fpc_database
    Post-process FPC database

=cut

sub postprocess_fpc_database {
    $logger->info("Post-processing FPC database");
}

=pod

=head2 postprocess_bac_database
    Post-process BAC database

=cut

sub postprocess_bac_database {
    $logger->info("Post-processing BAC database");
    $logger->info("Ensuring 'misc_feature' appears in 'meta_coord'");
    ensure_meta_coord_for_improved_regions();
}

=pod

=head2 synchronize_bac_clone
    Synchronizes sequenced clone information in BAC database
    
=cut

sub synchronize_bac_clone {
    my ($clone) = @_;

    our $database      ||= $clonepath->get_ensembl_db($BAC_DB_CODE);
    our $cs_adaptor    ||= $database->get_adaptor('CoordSystem');
    our $slice_adaptor ||= $database->get_SliceAdaptor;

    our $clone_cs  ||= $cs_adaptor->fetch_by_name('clone');
    our $contig_cs ||= $cs_adaptor->fetch_by_name('contig');
    our $counts    ||= {
        'total_clones'      => 0,
        'clone_updates'     => 0,
        'new_clone_slices'  => 0,
        'new_contig_slices' => 0,
    };

    $counts->{'total_clones'}++;
    my $clone_slice = $slice_adaptor->fetch_by_region('clone',
        $clone->{versioned_accession});

    if (defined $clone_slice) {
        if ($conf{'force'}) {
            $logger->debug(
                "Clone $clone->{versioned_accession} exists 'force' is on, so continuing."
            );
        } else {
            $logger->info(
                "Skipping existing clone $clone->{versioned_accession}.");

            return;
        }
    } else {
        $logger->info(
            "New versioned accession $clone->{versioned_accession}");
        $logger->debug("New seq_region for clone $clone->{accession_number}");
        $clone_slice = Bio::EnsEMBL::Slice->new(
            -coord_system    => $clone_cs,
            -start           => 1,
            -end             => $clone->{sequence_length},
            -strand          => 1,
            -seq_region_name => $clone->{versioned_accession},
            -adaptor         => $slice_adaptor,
        );
        if (!$conf{'dryrun'}) {
            $slice_adaptor->store($clone_slice);
            set_status_for_clone(
                {   'clone'   => $clone,
                    'slice'   => $clone_slice,
                    'db_code' => $BAC_DB_CODE
                }
            );
            set_clone_to_current(
                {   'clone'   => $clone,
                    'slice'   => $clone_slice,
                    'db_code' => $BAC_DB_CODE
                }
            );
        }
        $counts->{'new_clone_slices'}++;
    }

    remove_clone_contig_map($clone_slice, $BAC_DB_CODE);
    my $store_qualities = validate_qualities($clone);
    for my $contig (@{ $clone->{contigs} }) {
        my $num_improved_regions
            = defined($contig->{improved_regions})
            ? @{ $contig->{improved_regions} }
            : 0;
        $logger->debug(
            "$contig->{name} [$contig->{start_coord} $contig->{end_coord}] ",
            "for $clone->{versioned_accession} with $num_improved_regions ",
            "improved regions"
        );
        my $seq_region_name = seq_region_name_for($clone, $contig);
        if (   $contig->{end_coord} !~ m/^\d+$/
            || $contig->{start_coord} !~ m/^\d+$/)
        {
            qc_error(
                "Contig $seq_region_name with BAD coordinates ",
                "[contig->{start_coord} $contig->{end_coord}]! Skipping."
            );
            next;
        }
        my $contig_slice
            = $slice_adaptor->fetch_by_region('contig', $seq_region_name);
        my $create_new_contig_slice = 1;
        my $sequence = contig_sequence($clone, $contig);
        if (defined $contig_slice) {
            if (!$conf{'force'}) {
                qc_error(
                    "Contig '$seq_region_name' not expected in DB: @{[ $contig_slice->get_seq_region_id ]}"
                );
            }
            my $errors
                = assert_contigs_equal($contig_slice, $contig, \$sequence);
            if (scalar @$errors) {
                qc_error(
                    "[$seq_region_name] ",
                    "Contig conflicts between Clonepath and Ensembl: [",
                    join(', ', @$errors),
                    "]. REMOVING CONTIG."
                );
                remove_contig_slice($contig_slice);
                $create_new_contig_slice = 1;
            } else {
                $create_new_contig_slice = 0;
            }
        }

        if ($create_new_contig_slice) {
            $logger->debug("[$seq_region_name] Creating new contig slice");
            my $seq_region_start = 1;
            my $seq_region_end
                = $contig->{end_coord} - $contig->{start_coord} + 1;

            $contig_slice = create_new_contig_slice(
                {   'contig'       => $contig,
                    'name'         => $seq_region_name,
                    'adaptor'      => $slice_adaptor,
                    'sequence_ref' => \$sequence,
                    'start'        => $seq_region_start,
                    'end'          => $seq_region_end,
                    'contig_cs'    => $contig_cs,
                }
            );
            $counts->{'new_contig_slices'}++;
        } else {
            $logger->debug("[$seq_region_name] Skipping slice creation");
        }
        if ($create_new_contig_slice || $conf{force}) {
            if ($store_qualities) {
                my $qualities = contig_qualities($clone, $contig);
                set_contig_qualities($contig_slice, \$qualities,
                    \@QUALITY_WINDOW_SIZES);
            }
        }

        set_scaffold_for_contig(
            {   'region' => $contig,
                'slice'  => $contig_slice
            }
        );
        set_clone_end_for_contig(
            {   'region' => $contig,
                'slice'  => $contig_slice
            }
        );
        add_improved_regions_for_contig(
            {   'region' => $contig,
                'slice'  => $contig_slice
            }
        );
        $logger->debug(
            "Mapping contig $contig->{name} to ",
            "clone $clone->{accession_number}"
        );
        map_contig_to_clone(
            {   'clone_slice'  => $clone_slice,
                'contig_slice' => $contig_slice,
                'contig'       => $contig,
            }
        );
    }
}

=pod

=head2 synchronize_fpc_clone
    Synchronize FPC clone

=cut

sub synchronize_fpc_clone {
    my ($clone) = @_;
    update_accessioned_bac_set($clone);
    update_clone_feature($clone);
}

=pod

=head2 prepare_fpc_database
    Insert meta information before clones are added

=cut

sub prepare_fpc_database {
    $logger->info("Preparing FPC database");
}

sub ensure_meta_coord_for_improved_regions {
    my $max_length = $clonepath->improved_sequence_max_length()
        || $DEFAULT_IMPROVED_REGION_LENGTH;
    $logger->debug(
        "Using maximum length for improved sequences: $max_length");
    my $statement
        = $clonepath->prepared_statement($REPLACE_META_COORD_SQL_KEY,
        $BAC_DB_CODE);
    $statement->bind_param(1, 'misc_feature');
    $statement->bind_param(2, get_coord_system('contig', $BAC_DB_CODE)->dbID);
    $statement->bind_param(3, $max_length);
    $statement->execute();
}

=pod

=head2 get_coord_system
    Retrieves a coordinate system object

=cut

sub get_coord_system {
    my ($name, $db_code) = @_;
    our $cs_hash ||= {};
    $cs_hash->{$db_code}
        ||= $clonepath->get_ensembl_db($db_code)->get_adaptor('CoordSystem');

    return $cs_hash->{$db_code}->fetch_by_name($name);
}

sub ensure_clone_coord_system {
    ensure_coord_system(
        {   'name'           => 'clone',
            'rank'           => 1,
            'default'        => 1,
            'sequence_level' => 0,
        }
    );
}

sub ensure_contig_coord_system {
    ensure_coord_system(
        {   'name'           => 'contig',
            'rank'           => 2,
            'default'        => 1,
            'sequence_level' => 1
        }
    );
}

=pod

=head2 ensure_coord_system
    Ensures that a coordinate system is created

=cut

sub ensure_coord_system {
    my ($params)       = @_;
    my $name           = $params->{'name'};
    my $rank           = $params->{'rank'};
    my $default        = $params->{'default'};
    my $sequence_level = $params->{'sequence_level'};

    my $cs = get_coord_system($name, $BAC_DB_CODE);
    if (!defined($cs)) {
        $logger->info("Creating coord system '$name'");
        deactivate_clone_sequence_level();
        map_contig_coord_system_to_clone();
        $cs = Bio::EnsEMBL::CoordSystem->new(
            -NAME           => $name,
            -RANK           => $rank,
            -DEFAULT        => $default,
            -SEQUENCE_LEVEL => $sequence_level,
        );
        my $coord_system_adaptor = $clonepath->get_ensembl_db($BAC_DB_CODE)
            ->get_adaptor('CoordSystem');
        $coord_system_adaptor->store($cs) unless ($conf{'dryrun'});
    } else {
        $logger->debug("'$name' coord_system already in database");
    }
    return $cs;
}

sub map_contig_coord_system_to_clone {
    return if ($conf{'dryrun'});
    my $statement
        = $clonepath->prepared_statement($COORD_SYSTEM_CONTIG_SQL_KEY,
        $BAC_DB_CODE);
    $statement->execute();
}

sub create_new_contig_slice {
    my ($params)     = @_;
    my $contig       = $params->{'contig'};
    my $name         = $params->{'name'};
    my $adaptor      = $params->{'adaptor'};
    my $sequence_ref = $params->{'sequence_ref'};
    my $start        = $params->{'start'};
    my $end          = $params->{'end'};
    my $contig_cs    = $params->{'contig_cs'};
    $logger->debug("Adding contig $name");
    my $slice = Bio::EnsEMBL::Slice->new(
        -coord_system    => $contig_cs,
        -start           => $start,
        -end             => $end,
        -strand          => 1,
        -seq_region_name => $name,
        -adaptor         => $adaptor,
    );

    if (!$conf{'dryrun'}) {
        if (defined($$sequence_ref)) {
            $adaptor->store($slice, $sequence_ref);
        } else {
            return undef;
        }
    }
    return $slice;
}

sub remove_clone_contig_map {
    my ($clone_slice, $db_code) = @_;

    $logger->debug("Removing contig assembly for clone ",
        $clone_slice->display_id, " [id=", $clone_slice->get_seq_region_id,
        "]");
    my $statement
        = $clonepath->prepared_statement($REMOVE_CLONE_CONTIG_MAP_SQL_KEY,
        $db_code);
    $statement->bind_param(1, $clone_slice->get_seq_region_id);
    my $rows_modified = $statement->execute();

    $logger->debug("Removed contig assembly ($rows_modified rows)");
}

sub update_accession_version {
    my ($params) = @_;
    my $clone    = $params->{'clone'};
    my $slice    = $params->{'slice'};
    my $db_code  = $params->{'db_code'};
    set_region_attribute(
        {   'slice'       => $slice,
            'region'      => $clone,
            'db_code'     => $db_code,
            'db_field'    => 'version',
            'attrib_code' => 'acc-version',
            'attrib_name' => 'Accession version',
            'attrib_desc' => 'Accession version',
        }
    );
}

sub set_status_for_clone {
    my ($params) = @_;
    my $clone    = $params->{'clone'};
    my $slice    = $params->{'slice'};
    my $db_code  = $params->{'db_code'};
    set_region_attribute(
        {   'slice'       => $slice,
            'region'      => $clone,
            'db_code'     => $db_code,
            'db_field'    => 'status',
            'attrib_code' => 'seq-status',
            'attrib_name' => 'Sequence Status',
            'attrib_desc' => 'Sequence status',
        }
    );
}

=pod

=head2 set_clone_to_current
    Sets 'current-version' attribute and resets for previous clones

=cut

sub set_clone_to_current {
    my ($params) = @_;
    my $clone    = $params->{'clone'};
    my $slice    = $params->{'slice'};
    my $db_code  = $params->{'db_code'};
    remove_clone_current_attribute($clone->{accession_number}, $db_code);
    set_region_attribute(
        {   'slice'       => $slice,
            'region'      => $clone,
            'db_code'     => $db_code,
            'attrib_code' => 'current-version',
            'attrib_name' => 'Current Accession Version',
            'attrib_desc' =>
                'Identifies the most recent version of an accession',
            'value' => 'current',
        }
    );
}

=head2 contig_qualities

    Provides quality scores for given contig

=cut

sub contig_qualities {
    my ($clone, $contig) = @_;
    $logger->debug("Fetching contig qualities");
    return contig_subseq($clone, $contig, $SEQUENCE_QUALITY_TYPE);
}

sub contig_sequence {
    my ($clone, $contig) = @_;
    $logger->debug("Fetching contig sequence");
    return contig_subseq($clone, $contig, $SEQUENCE_DNA_TYPE);
}

=head2 contig_subseq

    Projects a clone sequence corresponding to the contig

=cut

sub contig_subseq {
    my ($clone, $contig, $type) = @_;
    my $sequence        = undef;
    my $subseq_function = undef;
    my $sequence_key    = undef;
    $logger->debug("Sequence type: $type");
    if ($type eq $SEQUENCE_DNA_TYPE) {
        $sequence_key = 'seq';
        $sequence     = Bio::Seq->new(
            -display_id => $clone->{accession_number},
            -seq        => $clone->{$sequence_key},
        );
        $subseq_function = sub {
            my ($sequence, $start, $end) = @_;
            return $sequence->subseq($start, $end);
        };
    } elsif ($type eq $SEQUENCE_QUALITY_TYPE) {
        $sequence_key = 'qualities';
        if (!defined($clone->{$sequence_key})) {
            qc_error(
                "[$clone->{accession_number}] ",
                "No qualities available for clone"
            );
            return undef;
        }
        $sequence = Bio::Seq::PrimaryQual->new(
            -display_id => $clone->{accession_number},
            -qual       => $clone->{$sequence_key},
        );
        $subseq_function = sub {
            my ($sequence, $start, $end) = @_;
            return join(" ", @{ $sequence->subqual($start, $end) });
        };
    } else {
        die "Unrecognized sequence type <$type>";
    }
    $logger->debug(
        "Sequence: $sequence with length @{[ $sequence->length ]}");
    my $subseq = undef;
    eval {
        $subseq = $subseq_function->(
            $sequence, $contig->{start_coord},
            $contig->{end_coord}
        );
    };
    if ($EVAL_ERROR) {
        qc_error(
            "[$clone->{accession_number}] $EVAL_ERROR",
            "   [Type=$type ",
            "Contig=($contig->{start_coord} $contig->{end_coord}) ",
            "length(seq)=@{[ length $clone->{$sequence_key} ]} ",
            "column(sequence_length)=$clone->{sequence_length}]"
        );
        return undef;
    }
    $logger->debug("Got @{[ length $subseq ]} bp for contig");
    return $subseq;

}

sub map_contig_to_clone {
    my ($params)     = @_;
    my $clone_slice  = $params->{clone_slice};
    my $contig_slice = $params->{contig_slice};
    my $contig       = $params->{contig};

    my $cp_start  = $contig->{start_coord};
    my $cp_end    = $contig->{end_coord};
    my $cp_length = $contig->{end_coord} - $contig->{start_coord} + 1;

    my $ens_start  = $contig_slice->start;
    my $ens_end    = $contig_slice->end;
    my $ens_length = $contig_slice->length;
    if ($logger->isDebugEnabled()) {
        $logger->debug(
            "ASSEMBLY: @{[$contig_slice->seq_region_name()]} ",
            "clonepath [$cp_start-$cp_end:$cp_length] ",
            "ensembl [$ens_start-$ens_end:$ens_length]"
        );
    }
    if ($ens_length != $cp_length) {
        qc_error(
            "[@{[$contig_slice->seq_region_name]}] ",
            "The contig length [$cp_length] does not agree ",
            "with the mapped length [$ens_length]"
        );
    }
    if ($conf{'dryrun'}) {
        return;
    }

    my $statement = $clonepath->prepared_statement($CLONE_CONTIG_MAP_SQL_KEY,
        $BAC_DB_CODE);
    $statement->bind_param(1, $clone_slice->get_seq_region_id);
    $statement->bind_param(2, $contig_slice->get_seq_region_id);
    $statement->bind_param(3, $contig->{start_coord});
    $statement->bind_param(4, $contig->{end_coord});
    $statement->bind_param(5, $contig_slice->start);
    $statement->bind_param(6, $contig_slice->end);
    $statement->bind_param(7, 1);
    my $rows_modified = $statement->execute();

    $logger->debug(
        "Mapped @{[$contig_slice->seq_region_name]} ",
        "to @{[$clone_slice->seq_region_name]} ",
        "at [$cp_start $cp_end] ($rows_modified rows)"
    );
}

sub set_scaffold_for_contig {
    my ($params) = @_;
    $params->{'db_field'}    = 'scaffold';
    $params->{'attrib_code'} = 'contig-scaffold';
    $params->{'attrib_name'} = 'Contig Scaffold';
    $params->{'attrib_desc'}
        = 'Scaffold that contains mutually ordered contigs';
    set_region_attribute($params);
}

sub set_clone_end_for_contig {
    my ($params) = @_;
    set_region_attribute(
        {   %$params,
            'db_field'    => 'clone_end',
            'attrib_code' => 'clone-end',
            'attrib_name' => 'Clone-End',
            'attrib_desc' =>
                'Side of the contig on which a vecgtor lies (enum:RIGHT, LEFT)'
        }
    );
    set_region_attribute(
        {   %$params,
            'db_field'    => 'vector',
            'attrib_code' => 'clone-vector',
            'attrib_name' => 'Vector sequence',
            'attrib_desc' =>
                'A clone-end vector associated with a contig (enum:SP6, T7)'
        }
    );
}

sub set_region_attribute {
    my ($params)    = @_;
    my $slice       = $params->{'slice'};
    my $region      = $params->{'region'};
    my $db_field    = $params->{'db_field'};
    my $attrib_code = $params->{'attrib_code'};
    my $attrib_name = $params->{'attrib_name'};
    my $attrib_desc = $params->{'attrib_desc'};
    my $value       = $params->{'value'};
    my $db_code = $params->{'db_code'} || $BAC_DB_CODE;

    $logger->debug("Setting '$attrib_code' attribute for region");

    my $attrib_value = $value || $region->{$db_field};
    if (!defined $attrib_value) {
        $logger->debug("No value defined. Attribute not stored");
        return;
    }

    my $attribute_adaptor
        = $clonepath->get_ensembl_db($db_code)->get_AttributeAdaptor();
    my $attrib = Bio::EnsEMBL::Attribute->new(
        -CODE        => $attrib_code,
        -NAME        => $attrib_name,
        -DESCRIPTION => $attrib_desc,
        -VALUE       => $attrib_value,
    );

    if (my $preexisting = $slice->get_all_Attributes($attrib_code)) {
        $logger->debug(
            "Removing @{[ scalar @$preexisting ]} preexisting '$attrib_code' attributes"
        );
        remove_seq_region_attrib($slice, $attrib_code, $db_code);
    }
    $attribute_adaptor->store_on_Slice($slice, [$attrib]);
}

sub add_improved_regions_for_contig {
    my ($params) = @_;
    my $slice    = $params->{'slice'};
    my $contig   = $params->{'region'};

    if (!defined $slice) {
        log_and_die(
            "Attempting to map improved sequences on an undefined slice! (Contig: $contig)"
        );
    }
    if (!defined $contig) {
        log_and_die(
            "Attempting to get improved sequences for an undefined contig! (Slice: $slice)"
        );
    }
    $logger->debug("Adding improved regions as misc_features");

    ensure_misc_set(
        {   'db_code'       => $BAC_DB_CODE,
            'misc_set_key'  => $IMPROVED_REGION_SET_KEY,
            'misc_set_name' => 'Improved Sequence',
            'misc_set_desc' => 'Sequence improved in finishing'
        }
    );

    for my $region (@{ $contig->{improved_regions} }) {
        $logger->debug(
            "Adding improved region ",
            "($region->{start_coord} $region->{end_coord}) for ",
            "$contig->{name} ($contig->{start_coord} $contig->{end_coord})"
        );
        my $feature = create_misc_feature(
            {   'slice' => $slice,
                'start' => $region->{start_coord} 
                    - $contig->{start_coord} + 1,
                'end' => $region->{end_coord} - $contig->{start_coord} + 1,
                'db_code' => $BAC_DB_CODE,
                'set_key' => $IMPROVED_REGION_SET_KEY,
            }
        );
    }
}

sub ensure_misc_set {
    my ($params)      = @_;
    my $db_code       = $params->{'db_code'};
    my $misc_set_key  = $params->{'misc_set_key'};
    my $misc_set_name = $params->{'misc_set_name'};
    my $misc_set_desc = $params->{'misc_set_desc'};
    my $length        = 1e6;

    my $misc_set_adaptor
        = $clonepath->get_ensembl_db($db_code)->get_adaptor('MiscSet');
    my $misc_set = Bio::EnsEMBL::MiscSet->new(
        -CODE            => $misc_set_key,
        -NAME            => $misc_set_name,
        -DESCRIPTION     => $misc_set_desc,
        -LONGEST_FEATURE => $length
    );
    $misc_set_adaptor->store($misc_set);
}

sub remove_seq_region_attrib {
    my ($slice, $attrib_code, $db_code) = @_;

    $logger->debug(
        "Removing attrib $attrib_code for @{[ $slice->get_seq_region_id ]}");
    my $statement
        = $clonepath->prepared_statement($REMOVE_CONTIG_ATTRIB_SQL_KEY,
        $db_code);
    $statement->bind_param(1, $slice->get_seq_region_id);
    $statement->bind_param(2, $attrib_code);

    $statement->execute() unless ($conf{'dryrun'});

}

=pod

=head2 remove_clone_current_attribute
    Remove current attribute for all clones with a given accession_number

=cut

sub remove_clone_current_attribute {
    my ($accession_number, $db_code) = @_;

    my $statement = $clonepath->prepared_statement(
        $REMOVE_CURRENT_ACCESSION_ATTRIB_SQL_KEY, $db_code);
    $statement->bind_param(1, $accession_number);
    $statement->execute() unless ($conf{'dryrun'});
}

sub deactivate_clone_sequence_level {

    unless ($conf{'dryrun'}) {
        my $statement
            = $clonepath->prepared_statement($COORD_SYSTEM_VERSION_SQL_KEY,
            $BAC_DB_CODE);
        $statement->execute();
    }

    # This needs to be done to sync adaptor with the update.
    $clonepath->get_ensembl_db($BAC_DB_CODE)->get_adaptor('CoordSystem')
        ->{'_is_sequence_level'} = {};
}

=pod

=head2 validate_clone
    Run a quality-control check to ensure the GenBank and Ensembl datasets jive

=cut

sub validate_clone {
    my ($feature, $genbank_value, $attribute_name) = @_;
    my $value = $feature->get_scalar_attribute($attribute_name);
    my $name  = $feature->get_scalar_attribute('name');
    if (defined $value && $value ne q{} && $value ne $genbank_value) {
        qc_error("[$name] Clone mismatch on field '$attribute_name' ",
            "(GenBank=$genbank_value Ensembl=$value)");
    }
}

=pod

=head2 update_clone_attribute
    Adds an attribute for the sequenced clone

=cut

sub update_clone_attribute {
    my ($feature, $attribute_name, $new_value) = @_;
    my $clone_name = $feature->get_scalar_attribute('name');
    my $old_value  = $feature->get_scalar_attribute($attribute_name);
    if (!defined $new_value) {
        qc_error("[$clone_name] Undefined new value for '$attribute_name'");
        return;
    }
    if (defined $old_value && $old_value ne q{}) {
        if ($old_value ne $new_value) {
            clone_debug_print($clone_name,
                "Overriding attribute '$attribute_name' ($old_value => $new_value)"
            );
            remove_attribute_from_feature($feature, $attribute_name,
                $FPC_DB_CODE);
        } else {
            clone_debug_print($clone_name,
                "Identical attribute '$attribute_name' => $old_value");
            return;
        }
    }
    my $attribute = Bio::EnsEMBL::Attribute->new(
        -CODE  => $attribute_name,
        -VALUE => $new_value,
    );
    my $log_message
        = "[$clone_name] Added '$attribute_name' attribute => $new_value\n";
    if ($conf{'dryrun'}) {
        dryrun_print($log_message);
        return;
    }
    $feature->add_Attribute($attribute);
    my $attribute_adaptor
        = $clonepath->get_ensembl_db($FPC_DB_CODE)->get_AttributeAdaptor();
    $attribute_adaptor->store_on_MiscFeature($feature, [$attribute]);

    $logger->debug($log_message);
}

=pod

=head2 remove_attribute_from_feature
    Since the API cannot do it for a single attribute

=cut

sub remove_attribute_from_feature {
    my ($feature, $attribute_name, $db_code) = @_;
    if ($conf{'dryrun'}) {
        dryrun_print("Removing attribute ",
            "[attribute=$attribute_name feature=@{[$feature->dbID]}]");
        return;
    }
    my $statement
        = $clonepath->prepared_statement($REMOVE_MISC_FEATURE_ATTRIB_SQL_KEY,
        $db_code);

    $statement->bind_param(1, $feature->dbID, SQL_INTEGER);
    $statement->bind_param(2, $attribute_name, SQL_VARCHAR);
    $statement->execute();
}

sub clone_debug_print {
    my ($clone, @messages) = @_;
    $logger->debug("[$clone] ", join('', @messages));
}

=pod

=head2 qc_error
    Print a QC error

=cut

sub qc_error {
    $logger->warn("[QC] ", join('', @_));
}

=pod

=head2 reset_previous_clones
    Sets initial attributes for historical clones

=cut

sub reset_previous_clones {
    my $adaptor = $clonepath->get_ensembl_db($FPC_DB_CODE)
        ->get_adaptor('MiscFeature');
    my @previous_clones
        = @{ $adaptor->fetch_all_by_attribute_type_value('embl_acc') };
    for my $clone (@previous_clones) {
        update_clone_attribute($clone, 'external',   'true');
        update_clone_attribute($clone, 'annotation', 'II');
    }
}

=pod

=head2 update_clone_feature
    Adds attributes for incoming sequenced clones

=cut

sub update_clone_feature {
    my ($seq_clone) = @_;
    $logger->debug(
        "Clone [name=$seq_clone->{clone_name} ",
        "accession=$seq_clone->{accession_number}]"
    );

    my ($clone_feature);
    eval { $clone_feature = fetch_feature_for($seq_clone); };
    if ($EVAL_ERROR) {
        qc_error($EVAL_ERROR);
        return;
    }

    validate_clone($clone_feature, $seq_clone->{clone_name},      'name');
    validate_clone($clone_feature, $seq_clone->{sequence_length}, 'seq_len');

    update_clone_attribute($clone_feature, 'embl_acc',
        $seq_clone->{versioned_accession});
    update_clone_attribute($clone_feature, 'seq_len',
        $seq_clone->{sequence_length});
    update_clone_attribute($clone_feature, 'seqstatus', $seq_clone->{status});
    update_clone_attribute($clone_feature, 'external',  'false');
    update_clone_attribute($clone_feature, 'annotation', 'II');
    update_clone_attribute($clone_feature, 'state',      '12:Accessioned');
}

=pod

=head2 update_accessioned_bac_set
    Updates the acc_bac_map misc_set to contain the sequenced BACs

=cut

sub update_accessioned_bac_set {
    my ($clone) = @_;
    my ($feature);
    eval { $feature = fetch_feature_for($clone); };
    if ($@) {
        qc_error($@);
        return;
    }
    add_misc_feature_to_misc_set($feature, $ACC_BAC_MAP_SET_KEY,
        $FPC_DB_CODE);
}

=pod

=head2 add_misc_feature_to_misc_set
    Store the correspondence between a feature and its misc_set. Required
    because no such method is provided in the API (must be a new misc_feature)

=cut

sub add_misc_feature_to_misc_set {
    my ($feature, $misc_set_key, $db_code) = @_;
    my $name = $feature->get_scalar_attribute('name');
    if (scalar @{ $feature->get_all_MiscSets($misc_set_key) } > 0) {
        clone_debug_print($name, "Already on $misc_set_key");
        return;
    }
    my $misc_set_adaptor
        = $clonepath->get_ensembl_db($db_code)->get_adaptor('MiscSet');
    my $misc_set = $misc_set_adaptor->fetch_by_code($misc_set_key);

    my $statement
        = $clonepath->prepared_statement($MISC_FEATURE_SET_SQL_KEY, $db_code);

    if (not defined $feature->dbID) {
        $logger->warn(
            "Cannot create correspondence for feature ",
            $feature->get_scalar_attribute('name'),
            ": No database ID"
        );
        return;
    }
    if (not defined $misc_set->dbID) {
        $logger->warn("Cannot create correspondence for set ",
            $misc_set->display_name, ": No database ID");
        return;
    }
    my @log_message = (
        "Adding feature correspondence: ",
        $feature->dbID, " -> ", $misc_set->dbID, "\n"
    );
    if ($conf{'dryrun'}) {
        dryrun_print(@log_message);
        return;
    }
    $logger->info(@log_message);
    $statement->bind_param(1, $feature->dbID,  SQL_INTEGER);
    $statement->bind_param(2, $misc_set->dbID, SQL_INTEGER);

    $statement->execute();
}

=pod

=head2 fetch_feature_for
    Fetches a feature for a given sequenced clone

=cut

sub fetch_feature_for {
    my ($sequenced_clone) = @_;
    my $clone_name = $sequenced_clone->{clone_name};
    if (exists $clone_feature_cache{$clone_name}) {
        return $clone_feature_cache{$clone_name};
    }
    my $misc_feature_adaptor = $clonepath->get_ensembl_db($FPC_DB_CODE)
        ->get_adaptor('MiscFeature');
    my $clone_features
        = $misc_feature_adaptor->fetch_all_by_attribute_type_value('name',
        $clone_name);
    if (scalar @$clone_features != 1) {
        if (scalar @$clone_features == 0) {
            die "[$clone_name] Clone not found in Ensembl";
        } else {
            die "[$clone_name] Unexpected error: ",
                "Multiple (@{[scalar @$clone_features]}) features";
        }
        $clone_feature_cache{$clone_name} = undef;
        return;
    }
    my $feature = $clone_features->[0];
    $clone_feature_cache{$clone_name} = $feature;
    return $feature;
}

sub dryrun_print {
    $logger->info("[DRYRUN] ", @_);
}

sub create_misc_feature {
    my ($params)     = @_;
    my $start        = $params->{'start'};
    my $end          = $params->{'end'};
    my $slice        = $params->{'slice'};
    my $db_code      = $params->{'db_code'};
    my $misc_set_key = $params->{'set_key'};

    my $strand = 1;

    my $misc_set = undef;

    $logger->debug("New MiscFeature ($start $end)");

    my $misc_feature_adaptor
        = $clonepath->get_ensembl_db($db_code)->get_adaptor('MiscFeature');
    my $misc_feature = Bio::EnsEMBL::MiscFeature->new(
        -START  => $start,
        -END    => $end,
        -STRAND => $strand,
        -SLICE  => $slice
    );

    if ($misc_set_key) {
        my $misc_set_adaptor
            = $clonepath->get_ensembl_db($db_code)->get_adaptor('MiscSet');
        $misc_set = $misc_set_adaptor->fetch_by_code($misc_set_key);
        $misc_feature->add_MiscSet($misc_set);
    }

    unless ($conf{'dryrun'}
        || misc_feature_exists($misc_feature, $BAC_DB_CODE))
    {
        $misc_feature_adaptor->store($misc_feature);
    }
    return $misc_feature;
}

sub misc_feature_exists {
    my ($misc_feature, $db_code) = @_;

    my $statement
        = $clonepath->prepared_statement($MISC_FEATURE_EXISTS_SQL_KEY,
        $db_code);
    $statement->bind_param(1, $misc_feature->slice->get_seq_region_id);
    $statement->bind_param(2, $misc_feature->start);
    $statement->bind_param(3, $misc_feature->end);
    $statement->bind_param(4, $misc_feature->strand);
    $statement->execute();
    my $count = $statement->fetchrow_arrayref()->[0];
    $logger->debug("MiscFeature exists? $count");
    return $count != 0;
}

sub clone_version {
    my ($slice) = @_;
    my @attributes = @{ $slice->get_all_Attributes('acc-version') };
    return defined($attributes[0]) ? $attributes[0]->value() : undef;
}

sub log_and_die {
    $logger->fatal(@_);
    die "FATAL! @_\n";
}

=pod

=head2 assert_contigs_equal
    Compares underlying information in contig slice and clonepath object

=cut

sub assert_contigs_equal {
    my ($slice, $contig, $sequence) = @_;

    my @errors = ();
    push @errors, "Sequences differ"
        if ($slice->seq() ne uc $$sequence);
    return \@errors;
}

=pod

=head2 remove_contig_slice
    Removes the contig from the Ensembl DB

=cut

sub remove_contig_slice {
    my ($slice) = @_;

    $logger->debug("[@{[$slice->seq_region_name]}] Removing contig slice");
    return if ($conf{'dryrun'});
    my $statement
        = $clonepath->prepared_statement($REMOVE_CONTIG_SLICE_SQL_KEY,
        $BAC_DB_CODE);
    $statement->bind_param(1, $slice->get_seq_region_id);
    $statement->bind_param(2, $slice->seq_region_name);
    my $rows_modified = $statement->execute();
    $logger->debug("Result: $rows_modified rows affected");

    # HACK: Remove the old slice ID from the adaptor's cache. This ensures
    # that when the contig gets recreated, the seq_region_id refreshed.
    # Before the cached value prevented the new contig from being mapped in
    # the 'assembly' table.
    my $slice_adaptor
        = $clonepath->get_ensembl_db($BAC_DB_CODE)->get_SliceAdaptor();
    my $key
        = join(":", $slice->seq_region_name(), $slice->coord_system->dbID());
    delete($slice_adaptor->{'sr_name_cache'}->{"$key"});

}

=head2 seq_region_name_for

    Constructs a name based on clone and contig

=cut

sub seq_region_name_for {
    my ($clone, $contig) = @_;
    return "$clone->{versioned_accession}-$contig->{name}";
}

=head2 set_contig_qualities

    Stores qualities at different window sizes for given contig

=cut

sub set_contig_qualities {
    my ($slice, $qualities, $sizes) = @_;

    our $database              ||= $clonepath->get_ensembl_db($BAC_DB_CODE);
    our $quality_score_adaptor ||= $database->get_adaptor('QualityScore');

    my $bioseq = Bio::Seq::PrimaryQual->new(-qual => $$qualities);
    my @qual_values = @{ $bioseq->qual() };
    for my $window_size (@$sizes) {
        my $extrapolated_values
            = extrapolate_qualities(\@qual_values, $window_size);
        my $extrapolated_bioseq
            = Bio::Seq::PrimaryQual->new(-qual => $extrapolated_values);
        my $quality_score = Bio::EnsEMBL::QualityScore->new(
            -POSITION    => 1,
            -SLICE       => $slice,
            -WINDOW_SIZE => $window_size,
            -SCORES      => join(' ', @{ $extrapolated_bioseq->qual }),
        );
        $quality_score_adaptor->store($quality_score);
    }
}

=head2 extrapolate_qualities

    Returns a list of qualities extrapolated over a given window size

=cut

sub extrapolate_qualities {
    my ($qualities, $window_size) = @_;
    if ($window_size == 1) {
        return $qualities;
    }
    my @extrapolated = ();
    for (my $i = 0; $i < scalar @$qualities; $i += $window_size) {
        my $last_index = $i + $window_size - 1;
        $last_index = $#$qualities if $last_index > $#$qualities;
        my @sub_array = @$qualities[ $i .. $last_index ];
        push @extrapolated,
            int(List::Util::sum(@sub_array) / scalar(@sub_array));
    }
    return \@extrapolated;
}

=head2 validate_qualities

    Validates that there are as many quality scores as there are nucleotides.

=cut

sub validate_qualities {
    my ($clone) = @_;
    my @errors = ();
    my $bioseq = Bio::Seq::PrimaryQual->new(-qual => $clone->{qualities});
    my $quality_length  = scalar @{ $bioseq->qual };
    my $sequence_length = length($clone->{seq});
    if ($quality_length == 0) {
        push @errors, "No qualities for clone";
    } elsif ($quality_length != $sequence_length) {
        push @errors,
            "Base qualities [length=$quality_length] do no match up with sequence [length=$sequence_length]";
    }
    if (scalar @errors) {
        qc_error("[$clone->{accession_number}] ", @errors);
        return 0;
    } else {
        return 1;
    }
}

1;
