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

=head1 DESCRIPTION

B<This program> 

  Loads overgo-clone associations into a core Ensembl DB
  
  Maintained by Shiran Pasternak <shiran@cshl.edu>

=cut

use strict;
use warnings;
use Getopt::Long;
use Config::General;
use Pod::Usage;
use Data::Dumper qw(Dumper);    # For debug

use English;
use Carp;
use DBI;
use DBI qw(:sql_types);
use FindBin qw($Bin);
use File::Basename qw( dirname );

use Log::Log4perl;

use vars qw($BASEDIR);

use Readonly;
use English qw( -no-match-vars );

use HTML::Template;
use POSIX;

Readonly my %DEFAULT_OPTIONS => (
    'clonepath_db'   => 'clonepath',
    'clonepath_host' => 'ascutney.cshl.edu',
    'clonepath_port' => 3306,
    'clonepath_user' => 'maize_rw',
    'clonepath_pass' => 'z3@m@ys',
    'registry_file'  => "${BASEDIR}/conf/SiteDefs.pm",
    'template_file'  => 'htdocs/sequence_status.html.tmpl',
    'output_file'    => 'htdocs/sequence_status.html',
);

Readonly my @REQUIRED_OPTIONS =>
    ('registry_file');

Readonly my @CLONEPATH_COLUMNS => qw(id accession_number version
    gi_number sequence_length
    chromosome status clone_name);

Readonly my $ORDERED_ORIENTED_ATTRIBUTE => 'order_orient';

my (%conf);
my ($misc_feature_adaptor, $attribute_adaptor, $misc_set_adaptor);
my ($clone_contig_map_statement);
my $registry_loaded     = 0;
my %clone_feature_cache = ();
my $logger              = Log::Log4perl->get_logger('bacstats');

BEGIN {

    # Set the perl libraries
    $BASEDIR = dirname($Bin);
    unshift @INC, $BASEDIR . '/ensembl-live/ensembl/modules';
    Log::Log4perl::init($BASEDIR . '/conf/log4perl.conf');
}

MAIN: {
    $OUTPUT_AUTOFLUSH = 1;
    load_configuration(\%conf);
    
    my $stats = get_bac_stats();
    update_bac_stats($stats);
    exit;
}

#======================================================================

sub get_bac_stats {
    my $dbh = get_clonepath_handle();
    
    my $results =
        $dbh->selectall_hashref(q(
select status, count(*) count from clone group by status with rollup
        ), 'status');
    
    my $stats = +{};
    for my $stat (keys %$results) {
        my $key = $stat || 'TOTAL';
        $stats->{$key} = $results->{$stat}->{'count'};
    }
    return $stats;    
}

sub update_bac_stats {
    my ($stats) = @_;
    
    $logger->info("Loading template file '$conf{'template_file'}'");
    my $template = HTML::Template->new(filename => $conf{'template_file'});
    
    for my $stat (keys %$stats) {
        $logger->debug("Templatizing [$stat -> $stats->{$stat}]");
        $template->param($stat => $stats->{$stat});
    }
    my $last_updated = formatted_date();
    $logger->debug("Templatizing [LAST_UPDATED -> $last_updated]");
    $template->param('LAST_UPDATED', $last_updated);
    
    open my $output, '>', $conf{'output_file'};
    print $output "Content-Type: text/html\n\n", $template->output;
    close $output;
}

sub formatted_date {
    return POSIX::strftime("%d %B %Y", localtime());
}


=pod

=head2 synchronize_fpc_database
    Synchronizes sequenced clone information in FPC database
    
=cut

sub synchronize_fpc_database {
    my ($clones) = @_;

    my $fpc_adaptor = load_db_adaptor($conf{'fpc_species'});

    $misc_feature_adaptor = $fpc_adaptor->get_adaptor('MiscFeature');
    $attribute_adaptor    = $fpc_adaptor->get_adaptor('Attribute');
    $misc_set_adaptor     = $fpc_adaptor->get_adaptor('MiscSet');

    reset_previous_clones();
    update_accessioned_bac_set($clones);
    update_clone_features($clones);
}

=pod

=head2 synchronize_bac_database
    Synchronizes sequenced clone information in BAC database
    
=cut

sub synchronize_bac_database {
    my ($clones) = @_;

    my $bac_adaptor = load_db_adaptor($conf{'bac_species'});

    my $coord_system_adaptor = $bac_adaptor->get_CoordSystemAdaptor;

    $logger->info("Ensuring 'contig' coord_system");
    ensure_contig_coord_system($coord_system_adaptor);
    add_seq_regions_for_contigs($bac_adaptor, $coord_system_adaptor, $clones);
}

sub ensure_contig_coord_system {
    my ($adaptor) = @_;

    my $contig_cs = $adaptor->fetch_by_name('contig');
    if (!defined($contig_cs)) {
        $logger->info('Creating contig coord system');
        deactivate_clone_sequence_level($adaptor);
        map_contig_coord_system_to_clone($adaptor);
        $contig_cs = Bio::EnsEMBL::CoordSystem->new(
            -NAME           => 'contig',
            -RANK           => 2,
            -DEFAULT        => 1,
            -SEQUENCE_LEVEL => 1
        );
        $adaptor->store($contig_cs);
    } else {
        $logger->debug("'contig' coord_system already in database");
    }
}

sub map_contig_coord_system_to_clone {
    my ($adaptor) = @_;
    return if ($conf{'dryrun'});
    my $statement = $adaptor->prepare(
        qq(
REPLACE INTO meta (meta_id, meta_key, meta_value)
      VALUES      (1, 'assembly.mapping', 'clone|contig')
    )
    );
    $statement->execute();
}

sub add_seq_regions_for_contigs {
    my ($adaptor, $cs_adaptor, $clones) = @_;

    my $clone_cs      = $cs_adaptor->fetch_by_name('clone');
    my $contig_cs     = $cs_adaptor->fetch_by_name('contig');
    my $slice_adaptor = $adaptor->get_SliceAdaptor();
    my %counts        = (
        'clone_updates'     => 0,
        'new_clone_slices'  => 0,
        'new_contig_slices' => 0,
    );
    for my $clone (grep { $_->{status} eq 'IMPROVED' } @$clones) {
        $counts{'clone_updates'}++;
        my $clone_slice = $slice_adaptor->fetch_by_region('clone',
            $clone->{accession_number});
        if (!defined($clone_slice)) {
            $logger->debug(
                "New seq_region for clone $clone->{accession_number}");
            $clone_slice = Bio::EnsEMBL::Slice->new(
                -coord_system    => $clone_cs,
                -start           => 1,
                -end             => $clone->{sequence_length},
                -strand          => 1,
                -seq_region_name => $clone->{accession_number},
                -adaptor         => $slice_adaptor,
            );
            if (!$conf{'dryrun'}) {
                $slice_adaptor->store($clone_slice);
            }
            $counts{'new_clone_slices'}++;
        } else {
            $logger->debug(
                "Found slice '@{[ $clone_slice->display_id ]}' in database");
        }
    CONTIG: for my $contig (@{ $clone->{contigs} }) {
            $logger->debug(
                "$contig->{name} [$contig->{start_coord} $contig->{end_coord}] ",
                "for $clone->{accession_number}"
            );
            my $seq_region_name
                = "$clone->{accession_number}-$contig->{name}";
            my $contig_slice
                = $slice_adaptor->fetch_by_region('contig', $seq_region_name);

            if (!defined($contig_slice)) {
                $logger->debug("Adding contig $seq_region_name");
                if (   $contig->{end_coord} !~ m/^\d+$/
                    || $contig->{start_coord} !~ m/^\d+$/)
                {
                    qc_error(
                        "Contig $seq_region_name with NULL coordinates ",
                        "[contig->{start_coord} $contig->{end_coord}]! Skipping."
                    );
                    next CONTIG;
                }
                my $seq_region_start = 1;
                my $seq_region_end
                    = $contig->{end_coord} - $contig->{start_coord} + 1;
                $contig_slice = Bio::EnsEMBL::Slice->new(
                    -coord_system    => $contig_cs,
                    -start           => $seq_region_start,
                    -end             => $seq_region_end,
                    -strand          => 1,
                    -seq_region_name => $seq_region_name,
                    -adaptor         => $slice_adaptor,
                );
                if (!$conf{'dryrun'}) {
                    my $sequence = contig_sequence($clone, $contig);
                    if (defined($sequence)) {
                        $slice_adaptor->store($contig_slice, \$sequence);
                    } else {
                        next CONTIG;
                    }
                }
                $counts{'new_contig_slices'}++;
            } else {
                $logger->debug("Found in DB: $contig_slice");
            }

            # TODO: Parse o-and-o from contig. For now just set to true
            $logger->debug(
                "Setting $ORDERED_ORIENTED_ATTRIBUTE attribute for contig");
            set_order_oriented_attribute_for_contig(
                {   'adaptor' => $adaptor,
                    'contig'  => $contig,
                    'slice'   => $contig_slice
                }
            );
            $logger->info(
                "Mapping contig $contig->{name} to ",
                "clone $clone->{accession_number}"
            );
            map_contig_to_clone(
                {   'adaptor'      => $slice_adaptor,
                    'clone_slice'  => $clone_slice,
                    'contig_slice' => $contig_slice,
                    'contig'       => $contig,
                }
            );
        }
    }
    $logger->info(<<__END_SUMMARY__);
SUMMARY:
    Clones with updates: $counts{'clone_updates'}
    New Clone Regions:   $counts{'new_clone_slices'}
    New Contig Regions:  $counts{'new_contig_slices'}
__END_SUMMARY__
}

sub contig_sequence {
    my ($clone, $contig) = @_;

    my $sequence = Bio::Seq->new(
        -display_id => $clone->{accession_number},
        -seq        => $clone->{seq},
    );
    my $subseq = undef;
    eval {
        $subseq
            = $sequence->subseq($contig->{start_coord}, $contig->{end_coord});
    };
    if ($EVAL_ERROR) {
        qc_error(
            "[$clone->{accession_number}] ",
            "Contig coordinates surpass sequence length ",
            "[Contig=($contig->{start_coord} $contig->{end_coord}) ",
            "length(seq)=@{[ length $clone->{seq} ]} ",
            "column(sequence_length)=$clone->{sequence_length}]"
        );
        return undef;
    }
    $logger->debug("Got @{[ length $subseq ]} bp for contig");
    return $subseq;
}

sub map_contig_to_clone {
    my ($params)     = @_;
    my $adaptor      = $params->{adaptor};
    my $clone_slice  = $params->{clone_slice};
    my $contig_slice = $params->{contig_slice};
    my $contig       = $params->{contig};
    if ($conf{'dryrun'}) {
        return;
    }
    if (!defined($clone_contig_map_statement)) {
        $clone_contig_map_statement = $adaptor->prepare(
            qq(
REPLACE INTO assembly
   (asm_seq_region_id, cmp_seq_region_id,
    asm_start, asm_end, cmp_start, cmp_end, ori) 
VALUES (?,?,?,?,?,?,?)
        )
        );
    }
    $clone_contig_map_statement->bind_param(1,
        $clone_slice->get_seq_region_id);
    $clone_contig_map_statement->bind_param(2,
        $contig_slice->get_seq_region_id);
    $clone_contig_map_statement->bind_param(3, $contig->{start_coord});
    $clone_contig_map_statement->bind_param(4, $contig->{end_coord});
    $clone_contig_map_statement->bind_param(5, $contig_slice->start);
    $clone_contig_map_statement->bind_param(6, $contig_slice->end);
    $clone_contig_map_statement->bind_param(7, 1);
    $clone_contig_map_statement->execute();
}

sub set_order_oriented_attribute_for_contig {
    my ($params) = @_;
    my $adaptor  = $params->{'adaptor'};
    my $slice    = $params->{'slice'};
    my $contig   = $params->{'contig'};

    my $attrib_value = $contig->{ordered} eq 'Y' ? 'true' : 'false';

    my $attribute_adaptor = $adaptor->get_AttributeAdaptor();
    my $o_and_o_attrib = Bio::EnsEMBL::Attribute->new(
        -CODE => $ORDERED_ORIENTED_ATTRIBUTE,
        -NAME => "Contig Order-Orientation",
        -DESCRIPTION =>
            "Whether or not a contig has order and orientation (boolean)",
        -VALUE => $attrib_value,
    );
    if (my $preexisting
        = $slice->get_all_Attributes($ORDERED_ORIENTED_ATTRIBUTE))
    {
        $logger->debug("Removing @{[ scalar @$preexisting ]} preexisting attributes");
        $attribute_adaptor->remove_from_Slice($slice);
    }
    $attribute_adaptor->store_on_Slice($slice, [$o_and_o_attrib]);
}

sub deactivate_clone_sequence_level {
    my ($adaptor) = @_;
    my $statement = $adaptor->prepare(
        qq(
UPDATE coord_system SET attrib = 'default_version' WHERE name = 'clone'
    )
    );
    $statement->execute();

    # This needs to be done to sync adaptor with the update.
    $adaptor->{'_is_sequence_level'} = {};
}

=pod

=head2 load_db_adaptor
    Loads a DB adaptor from the registry

=cut

sub load_db_adaptor {
    my ($species) = @_;

    if (!$registry_loaded) {

        # Load the ensembl file
        Bio::EnsEMBL::Registry->load_all($conf{'registry_file'});
        $registry_loaded = 1;
    }
    my $adaptor
        = ${ Bio::EnsEMBL::Registry->get_all_DBAdaptors(-species => $species)
        }[0];
    if (!defined $adaptor) {
        warn("No core DB for $species set in $conf{'registry_file'}: $@\n");
        pod2usage(1);
    }

    return $adaptor;
}

=pod

=head2 get_clonepath_handle
    Creates a connection to the clonepath database

=cut

sub get_clonepath_handle {
    $logger->debug("Connecting to $conf{'clonepath_host'}:$conf{'clonepath_user'}\@$conf{'clonepath_db'}");
    my $data_source = sprintf('DBI:mysql:database=%s;host=%s;port=%s',
        $conf{'clonepath_db'}, $conf{'clonepath_host'},
        $conf{'clonepath_port'});
    my $db_handle = undef;
    eval {
        $db_handle = DBI->connect($data_source, $conf{'clonepath_user'},
            $conf{'clonepath_pass'});
    };
    if (!defined($db_handle)) {
        die "Can't connect to database: @{[ $DBI::errstr ]}\n";
    }
    return $db_handle;
}

=pod

=head2 fetch_sequenced_clones
    Returns clones from the clonepath database

=cut

sub fetch_sequenced_clones {
    my $clonepath_handle = get_clonepath_handle();
    my $fields = join(', ', @CLONEPATH_COLUMNS);
    my $clones
        = $clonepath_handle->selectall_arrayref("SELECT $fields FROM clone",
        { Slice => {} });

    $logger->debug("Fetching sequence and contigs for clones");
    fetch_contigs_for_clones($clonepath_handle, $clones);
    return $clones;
}

sub fetch_contigs_for_clones {
    my ($handle, $clones) = @_;
    my $contig_statement = $handle->prepare(
        qq(
SELECT *
  FROM contig WHERE clone_id = ?
    )
    );
    my $sequence_statement = $handle->prepare(
        qq(
SELECT seq FROM clone WHERE id = ?
    )
    );
    for my $clone (@$clones) {
        my $clone_desc
            = "clone $clone->{accession_number} ($clone->{status})";
        if ($clone->{status} ne 'IMPROVED') {
            $logger->debug("Skipping contigs for clone $clone_desc");
            next;
        }
        $logger->debug("Fetching contigs for $clone_desc");
        $sequence_statement->execute($clone->{id});
        ($clone->{seq}) = $sequence_statement->fetchrow_array;
        $contig_statement->execute($clone->{id});
        $clone->{contigs} = $contig_statement->fetchall_arrayref({});

=pod        
        print <<__END_CLONE__;
CLONE [
    id      : $clone->{id}
    status  : $clone->{status}
    contigs : [
@{[ join("\n", map { my $ctg=$_; ' 'x8 . join(' ', map { "$_=$ctg->{$_}" } keys %$ctg) } @{ $clone->{contigs} }) ]}
    ]
]
__END_CLONE__
=cut

        $logger->debug(
            "Got @{[ scalar @{$clone->{contigs}} ]} contigs, ",
            "sequence with @{[ length $clone->{seq} ]} bp"
        );

    }
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
            qc_error("[$clone_name] Overriding attribute '$attribute_name' ",
                "($old_value => $new_value)");
            remove_attribute_from_feature($feature, $attribute_name);
        } else {
            clone_print($clone_name,
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
    $attribute_adaptor->store_on_MiscFeature($feature, [$attribute]);

    $logger->info($log_message);
}

=pod

=head2 remove_attribute_from_feature
    Since the API cannot do it for a single attribute

=cut

sub remove_attribute_from_feature {
    my ($feature, $attribute_name) = @_;
    if ($conf{'dryrun'}) {
        dryrun_print("Removing attribute ",
            "[attribute=$attribute_name feature=@{[$feature->dbID]}]");
        return;
    }
    my $statement = $misc_feature_adaptor->prepare(
        q(
DELETE
  FROM misc_attrib
 WHERE misc_feature_id = ?
   AND attrib_type_id = (SELECT attrib_type_id 
                           FROM attrib_type WHERE code = ?)
    )
    );

    $statement->bind_param(1, $feature->dbID, SQL_INTEGER);
    $statement->bind_param(2, $attribute_name, SQL_VARCHAR);
    $statement->execute();
}

sub clone_print {
    my ($clone, @messages) = @_;
    $logger->info("[$clone] ", join('', @messages));
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
    my @previous_clones
        = @{ $misc_feature_adaptor->fetch_all_by_attribute_type_value(
            'embl_acc') };
    for my $clone (@previous_clones) {
        update_clone_attribute($clone, 'external',   'true');
        update_clone_attribute($clone, 'annotation', 'II');
    }
}

=pod

=head2 update_clone_features
    Adds attributes for incoming sequenced clones

=cut

sub update_clone_features {
    my ($clones) = @_;
    for my $seq_clone (@$clones) {
        $logger->info(
            "Clone [name=$seq_clone->{clone_name} ",
            "accession=$seq_clone->{accession_number}]"
        );

        my ($clone_feature);
        eval { $clone_feature = fetch_feature_for($seq_clone); };
        if ($EVAL_ERROR) {
            qc_error($EVAL_ERROR);
            next;
        }

        validate_clone($clone_feature, $seq_clone->{clone_name}, 'name');
        validate_clone($clone_feature, $seq_clone->{sequence_length},
            'seq_len');

        update_clone_attribute($clone_feature, 'embl_acc',
            $seq_clone->{accession_number});
        update_clone_attribute($clone_feature, 'seq_len',
            $seq_clone->{sequence_length});
        update_clone_attribute($clone_feature, 'seqstatus',
            $seq_clone->{status});
        update_clone_attribute($clone_feature, 'external',   'false');
        update_clone_attribute($clone_feature, 'annotation', 'I');
        update_clone_attribute($clone_feature, 'state', '12:Accessioned');
    }
}

=pod

=head2 update_accessioned_bac_set
    Updates the acc_bac_map misc_set to contain the sequenced BACs

=cut

sub update_accessioned_bac_set {
    my ($clones) = @_;
    my $acc_bac_map = $misc_set_adaptor->fetch_by_code('acc_bac_map');
    for my $clone (@$clones) {
        my ($feature);
        eval { $feature = fetch_feature_for($clone); };
        if ($@) {
            qc_error($@);
            next;
        }
        add_feature_set($feature, $acc_bac_map);
    }
}

=pod

=head2 add_feature_set
    Store the correspondence between a feature and its misc_set. Required
    because no such method is provided in the API (must be a new misc_feature)

=cut

sub add_feature_set {
    my ($feature, $misc_set) = @_;
    my $name = $feature->get_scalar_attribute('name');
    if (scalar @{ $feature->get_all_MiscSets('acc_bac_map') } > 0) {
        clone_print($name, 'Already on acc_bac_map');
        return;
    }
    $feature->add_MiscSet($misc_set);    # Does nothing in the DB

    my $feature_set_sth = $misc_feature_adaptor->prepare(
        q(
INSERT IGNORE misc_feature_misc_set SET misc_feature_id = ?, misc_set_id = ?
    )
    );

    my @log_message = (
        "Adding feature correspondence: ",
        $feature->dbID, " -> ", $misc_set->dbID, "\n"
    );
    if ($conf{'dryrun'}) {
        dryrun_print(@log_message);
        return;
    }
    $logger->info(@log_message);
    $feature_set_sth->bind_param(1, $feature->dbID,  SQL_INTEGER);
    $feature_set_sth->bind_param(2, $misc_set->dbID, SQL_INTEGER);

    $feature_set_sth->execute();
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

sub load_configuration {
    my ($options) = @_;
    GetOptions(
        $options,          "help|h",
        "man|m",           "fpc_species|f=s",
        "bac_species|b=s", "registry_file|r=s",
        "dryrun",          "conf|c=s",
	"template_file=s", "output_file=s",
        "printconf"
    ) or pod2usage(2);
    pod2usage(-verbose => 2) if $options->{'man'};
    pod2usage(1) if $options->{'help'};

    if ($options->{'conf'}) {

        # Unable to pass Readonly %DEFAULT_OPTIONS, 'new' tries to modify it
        my %defaults      = %DEFAULT_OPTIONS;
        my $configuration = new Config::General(
            -ConfigFile            => $options->{'conf'},
            -DefaultConfig         => \%defaults,
            -MergeDuplicateOptions => 'true',
            -InterPolateVars       => 'true',
        );
        my %config_object = $configuration->getall;
        while (my ($key, $value) = each %config_object) {
            if (!defined($options->{$key})) {
                $options->{$key} = $value;
            }
        }
    }

    for my $required (@REQUIRED_OPTIONS) {
        if (!defined($options->{$required})) {
            warn("Option --$required is required\n");
            pod2usage(1);
        }
    }

    if ($options->{'printconf'}) {
        print "CONFIGURATION\n";
        for my $key (grep { $_ ne 'printconf' } sort keys %conf) {
            printf "   \%-20s : \%s\n", $key, $conf{$key};
        }
        exit 0;
    }
}

sub dryrun_print {
    $logger->info("[DRYRUN] ", @_);
}

1;
