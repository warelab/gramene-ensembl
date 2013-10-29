
package Maize::Clonepath;

=head1 NAME

Maize::Clonepath

=head1 SYNOPSIS

  use Maize::Clonepath

  
=head1 DESCRIPTION

A utility module for interfacing between the maize clonepath database and
the Ensembl Database

=cut

use strict;
use warnings;

use Getopt::Long;
use Config::General;
use Readonly;
use Log::Log4perl;
use FindBin qw($Bin);
use File::Basename qw( dirname );

use Bio::EnsEMBL::Registry;

BEGIN {
    use Exporter qw(import);

    our (@EXPORT);

    @EXPORT = map {"*${_}_SQL_KEY"}
        qw(IMPROVED_REGION CLONE_CONTIG_MAP REMOVE_MISC_FEATURE_ATTRIB
        REMOVE_CONTIG_ATTRIB COORD_SYSTEM_CONTIG COORD_SYSTEM_VERSION
        SEQUENCE CLONE_CONTIGS MISC_FEATURE_SET MISC_FEATURE_EXISTS
        REMOVE_CLONE_CONTIG_MAP REMOVE_SEQ_REGION_SEQ REMOVE_SEQ_REGION_ATTRIB
        REMOVE_SEQ_REGION UPDATE_SEQ_REGION_NAME REPLACE_META_COORD
        REMOVE_CURRENT_ACCESSION_ATTRIB CLONE MAX_IMPROVED_SEQUENCE_LENGTH
        REMOVE_CONTIG_SLICE);
    push @EXPORT, qw(*BAC_DB_CODE *FPC_DB_CODE
        *ACC_BAC_MAP_SET_KEY *IMPROVED_REGION_SET_KEY);

    Log::Log4perl::init("$Bin/../conf/log4perl.conf");
}

Readonly my @CLONEPATH_COLUMNS => qw(id accession_number version
    gi_number sequence_length chromosome status clone_name);

Readonly my %DEFAULT_OPTIONS => (
    'clonepath_db'   => 'clonepath',
    'clonepath_host' => 'ascutney.cshl.edu',
    'clonepath_port' => 3306,
    'clonepath_user' => 'maize_rw',
    'clonepath_pass' => 'z3@m@ys',
    'registry_file'  => '${basedir}/conf/SiteDefs.pm',
);

Readonly my @REQUIRED_OPTIONS => ('registry_file');

Readonly our $BAC_DB_CODE => 'bac';
Readonly our $FPC_DB_CODE => 'fpc';

Readonly our $ACC_BAC_MAP_SET_KEY     => 'acc_bac_map';
Readonly our $IMPROVED_REGION_SET_KEY => 'improved_sequence';

Readonly our $CLONE_SQL_KEY                      => 'clone';
Readonly our $SPECIFIC_CLONE_SQL_KEY             => 'specific_clone';
Readonly our $IMPROVED_REGION_SQL_KEY            => 'contig_improved_region';
Readonly our $CLONE_CONTIG_MAP_SQL_KEY           => 'clone_contig_map';
Readonly our $REMOVE_MISC_FEATURE_ATTRIB_SQL_KEY => 'rem_misc_feature_attrib';
Readonly our $REMOVE_CONTIG_ATTRIB_SQL_KEY       => 'rem_seq_region_attrib';
Readonly our $COORD_SYSTEM_CONTIG_SQL_KEY        => 'coord_system_contig';
Readonly our $COORD_SYSTEM_VERSION_SQL_KEY       => 'cs_def_ver';
Readonly our $SEQUENCE_SQL_KEY                   => 'clone_sequence';
Readonly our $QUALITIES_SQL_KEY                  => 'clone_qualities';
Readonly our $CLONE_CONTIGS_SQL_KEY              => 'clone_contigs';
Readonly our $MISC_FEATURE_SET_SQL_KEY           => 'insert_misc_feature_set';
Readonly our $MISC_FEATURE_EXISTS_SQL_KEY        => 'misc_feature_exists';
Readonly our $REMOVE_CLONE_CONTIG_MAP_SQL_KEY    => 'rem_clone_contig_map';
Readonly our $REMOVE_SEQ_REGION_SQL_KEY          => 'rem_seq_region';
Readonly our $REMOVE_SEQ_REGION_SEQ_SQL_KEY      => 'rem_seq_region_seq';
Readonly our $REMOVE_SEQ_REGION_ATTRIB_SQL_KEY =>
    'rem_all_seq_region_attribs';
Readonly our $UPDATE_SEQ_REGION_NAME_SQL_KEY => 'upd_seq_region_name';
Readonly our $REPLACE_META_COORD_SQL_KEY     => 'replace_meta_coord';
Readonly our $REMOVE_CURRENT_ACCESSION_ATTRIB_SQL_KEY =>
    'rem_current_acc_attribs';
Readonly our $MAX_IMPROVED_SEQUENCE_LENGTH_SQL_KEY =>
    'max_improved_seq_length';
Readonly our $REMOVE_CONTIG_SLICE_SQL_KEY => 'remove_contig_slice';

Readonly my %SQL_QUERY_FOR => (
    $CLONE_SQL_KEY => "select @{[join(',', @CLONEPATH_COLUMNS)]} from clone",
    $SPECIFIC_CLONE_SQL_KEY => <<SQL,
select @{[join(',', @CLONEPATH_COLUMNS)]}
  from clone where accession_number = ?
SQL
    $IMPROVED_REGION_SQL_KEY => <<SQL,
select * from improved_region where contig_id = ?
SQL
    $CLONE_CONTIG_MAP_SQL_KEY => <<SQL,
replace into assembly
   (asm_seq_region_id, cmp_seq_region_id,
    asm_start, asm_end, cmp_start, cmp_end, ori) 
values (?,?,?,?,?,?,?)
SQL
    $REMOVE_CLONE_CONTIG_MAP_SQL_KEY => <<SQL,
delete from assembly
 where asm_seq_region_id = ?
SQL
    $REMOVE_MISC_FEATURE_ATTRIB_SQL_KEY => <<SQL,
delete
  from misc_attrib
 where misc_feature_id = ?
   and attrib_type_id = (select attrib_type_id 
                           from attrib_type where code = ?)
SQL
    $COORD_SYSTEM_CONTIG_SQL_KEY => <<SQL,
replace into meta (meta_id, meta_key, meta_value)
      values      (1, 'assembly.mapping', 'clone|contig')
SQL
    $REMOVE_CONTIG_ATTRIB_SQL_KEY => <<SQL,
delete from seq_region_attrib
 where seq_region_id = ?
   and attrib_type_id = (select attrib_type_id
                           from attrib_type
                          where code = ?)
SQL
    $COORD_SYSTEM_VERSION_SQL_KEY => <<SQL,
update coord_system set attrib = 'default_version' where name = 'clone'
SQL
    $CLONE_CONTIGS_SQL_KEY => <<SQL,
select c.id, c.name, c.start_coord, c.end_coord,
       s.name as scaffold, cec.vector, cec.clone_end
  from contig c left join scaffold s on (c.scaffold_id = s.id)
                left join clone_end_contig cec on (c.id = cec.id)
 where c.clone_id = ?
SQL
    $SEQUENCE_SQL_KEY => <<SQL,
select seq from clone where id = ?
SQL
    $QUALITIES_SQL_KEY => <<SQL,
select raw_scores from clone_quality_scores where clone_id = ?
SQL
    $MISC_FEATURE_SET_SQL_KEY => <<SQL,
insert ignore misc_feature_misc_set set misc_feature_id = ?, misc_set_id = ?
SQL
    $MISC_FEATURE_EXISTS_SQL_KEY => <<SQL,
select count(*)
  from misc_feature
 where seq_region_id = ?
   and seq_region_start = ?
   and seq_region_end = ?
   and seq_region_strand = ?
SQL
    $REMOVE_SEQ_REGION_SQL_KEY => <<SQL,
delete from seq_region where seq_region_id = ?        
SQL
    $REMOVE_SEQ_REGION_SEQ_SQL_KEY => <<SQL,
delete from dna where seq_region_id = ?        
SQL
    $REMOVE_SEQ_REGION_ATTRIB_SQL_KEY => <<SQL,
delete from seq_region_attrib where seq_region_id = ?        
SQL
    $UPDATE_SEQ_REGION_NAME_SQL_KEY => <<SQL,
update seq_region set name = ? where seq_region_id = ?
SQL
    $REPLACE_META_COORD_SQL_KEY => <<SQL,
replace into meta_coord (table_name, coord_system_id, max_length) values (?, ?, ?)
SQL
    $REMOVE_CURRENT_ACCESSION_ATTRIB_SQL_KEY => <<SQL,
delete
  from seq_region_attrib
 using seq_region_attrib
  left join seq_region
         on seq_region_attrib.seq_region_id = seq_region.seq_region_id
 where seq_region.name like concat(?, '%')
SQL
    $MAX_IMPROVED_SEQUENCE_LENGTH_SQL_KEY => <<SQL,
select max(end_coord - start_coord + 1) from improved_region
SQL
    $REMOVE_CONTIG_SLICE_SQL_KEY => <<SQL,
delete seq_region, dna, assembly, seq_region_attrib
  from seq_region
       left join dna
              on dna.seq_region_id = seq_region.seq_region_id
       left join assembly
              on assembly.cmp_seq_region_id = seq_region.seq_region_id
       left join seq_region_attrib
              on seq_region_attrib.seq_region_id = seq_region.seq_region_id
 where seq_region.seq_region_id = ?
   and seq_region.name = ?
SQL
);

our ($logger);

=pod

=head2 new
    Constructor

=cut

sub new {
    my $class = shift;
    my ($params) = @_;
    $params ||= +{};
    my $self = bless {%$params}, $class;

    $self->{'basedir'} ||= dirname($Bin);

    $self->{'databases'} = +{
        'fpc' => +{
            'handle'     => undef,
            'statements' => +{},
        },
        'bac' => +{
            'handle'     => undef,
            'statements' => +{},
        },
        'clonepath' => +{
            'handle'     => undef,
            'statements' => +{},
        }
    };

    $logger = Log::Log4perl->get_logger('clonepath');

    return $self;
}

# Initializes an Ensembl database
sub _init_ensembl_database {
    my $self     = shift;
    my ($code)   = @_;
    my $database = $self->{'databases'}->{$code};
    $database->{'adaptor'} = $self->load_db_adaptor($code);
    $database->{'handle'}  = $database->{'adaptor'}->dbc;
    if ($self->{dryrun}) {
        $logger->debug("Setting Ensembl DB connection to read-only");
        $database->{'handle'}->{ReadOnly} = 1;
    }
}

=pod

=head2 load_db_adaptor
    Loads a DB adaptor from the registry

=cut

sub load_db_adaptor {
    my $self = shift;
    my ($code) = @_;

    if (!$self->{'_registry_loaded'}) {

        # Load the ensembl file
        Bio::EnsEMBL::Registry->load_all($self->{'registry_file'});
        $self->{'_registry_loaded'} = 1;
    }
    my $species = $self->{"${code}_species"};
    my $adaptor
        = ${ Bio::EnsEMBL::Registry->get_all_DBAdaptors(-species => $species)
        }[0];
    if (!defined $adaptor) {
        die("No core DB for $species set in $self->{'registry_file'}: $@\n");
    }

    return $adaptor;
}

=pod

=head2 get_clonepath_handle
    Creates a connection to the clonepath database

=cut

sub get_clonepath_handle {
    my $self     = shift;
    my $database = $self->{'databases'}->{'clonepath'};
    if (!defined $database->{'handle'}) {
        $database->{'handle'} = $self->_init_clonepath_handle;
    }
    return $database->{'handle'};
}

sub _init_clonepath_handle {
    my $self        = shift;
    my $data_source = sprintf(
        'DBI:mysql:database=%s;host=%s;port=%s',
        $self->{'clonepath_db'},
        $self->{'clonepath_host'},
        $self->{'clonepath_port'}
    );
    my $db_handle = undef;
    eval {
        $db_handle = DBI->connect(
            $data_source,
            $self->{'clonepath_user'},
            $self->{'clonepath_pass'}
        );
    };
    if (!defined($db_handle)) {
        die "Can't connect to database: @{[ $DBI::errstr ]}\n";
    }
    if ($self->{dryrun}) {
        $logger->debug("Setting database handle to read-only");
        $db_handle->{ReadOnly} = 1;
    }
    return $db_handle;
}

=pod

=head2 load_configuration
    Loads configuration from configuration file and command-line arguments.

=cut

sub load_configuration {
    my $self = shift;
    my ($options) = @_;
    GetOptions(
        $options,          "help|h",
        "man|m",           "fpc_species|f=s",
        "bac_species|b=s", "registry_file|r=s",
        "dryrun",          "conf|c=s",
        "printconf",       "accession_number=s",
        "force!",
    ) or die 'Bad options';

    my %defaults = %DEFAULT_OPTIONS;

    my $configuration = new Config::General(
        -DefaultConfig         => \%defaults,
        -ConfigFile            => $options->{'conf'},
        -MergeDuplicateOptions => 1,
        -InterPolateVars       => 0,
    );
    my %config_object = $configuration->getall;
    while (my ($key, $value) = each %config_object) {
        if (!defined($options->{$key})) {
            $options->{$key} = $value;
        }
    }

    for my $default_option (keys %defaults) {
        $options->{$default_option} ||= $defaults{$default_option};
    }

    for my $required (@REQUIRED_OPTIONS) {
        if (!defined($options->{$required})) {
            die "Option --$required is required";
        }
    }

    # Interpolate nested variables;
    for my $option (keys %$options) {
        my $value = $options->{$option};
        $value =~ s/\$\{(\w+)\}/$self->{$1}/g;
        $options->{$option} = $value;
    }

    if ($options->{'printconf'}) {
        print "CONFIGURATION\n";
        for my $key (grep { $_ ne 'printconf' } sort keys %$options) {
            printf "   \%-20s : \%s\n", $key, $options->{$key};
        }
        exit 0;
    }

    for my $option (keys %$options) {
        $self->{$option} = $options->{$option};
    }
}

=pod

=head2 prepared_statement
    Returns a prepared statement. If the statement is not already, it is
    prepared.

=cut

sub prepared_statement {
    my $self      = shift;
    my $key       = shift;
    my $db_code   = shift || 'clonepath';
    my $database  = $self->{'databases'}->{$db_code};
    my $statement = $database->{'statements'}->{$key};
    if (not defined $statement) {
        $statement = $database->{'statements'}->{$key}
            = $database->{'handle'}->prepare($SQL_QUERY_FOR{$key});
    }
    return $statement;
}

=pod

=head2 get_ensembl_db
    Retrieves the Ensembl root DB adaptor

=cut

sub get_ensembl_db {
    my $self     = shift;
    my ($code)   = @_;
    my $database = $self->{'databases'}->{$code};
    if (!defined $database->{'adaptor'}) {
        $self->_init_ensembl_database($code);
    }
    return $database->{'adaptor'};
}

=pod

=head2 next_sequenced_clone
    Returns the next clone from the clonepath database

=cut

sub next_sequenced_clone {
    my $self = shift;
    our $db_handle ||= $self->get_clonepath_handle();
    our $clone_statement;
    if (!defined $clone_statement) {
        if (defined $self->{'accession_number'}) {
            $clone_statement
                = $self->prepared_statement($SPECIFIC_CLONE_SQL_KEY);
            $clone_statement->bind_param(1, $self->{'accession_number'});
        } else {
            $clone_statement = $self->prepared_statement($CLONE_SQL_KEY);
        }
        $clone_statement->execute();
    }
    my $clone = $clone_statement->fetchrow_hashref();
    if (defined $clone) {
        $self->fetch_contigs_for_clone($clone);
        $clone->{versioned_accession}
            = join('.', map { $clone->{$_} } qw/accession_number version/);
        $logger->debug("Next clone: $clone->{'versioned_accession'}");
    } else {
        $logger->debug("No more clones!");
    }
    return $clone;
}

sub fetch_contigs_for_clone {
    my $self = shift;
    my ($clone) = @_;

    our $contig_statement
        ||= $self->prepared_statement($CLONE_CONTIGS_SQL_KEY);
    our $sequence_statement ||= $self->prepared_statement($SEQUENCE_SQL_KEY);
    our $qualities_statement
        ||= $self->prepared_statement($QUALITIES_SQL_KEY);

    my $clone_desc = "clone $clone->{accession_number} ($clone->{status})";
    $logger->debug("Fetching contigs for $clone_desc");
    $sequence_statement->execute($clone->{id});
    ($clone->{seq}) = $sequence_statement->fetchrow_array;
    $qualities_statement->execute($clone->{id});
    ($clone->{qualities}) = $qualities_statement->fetchrow_array;
    if ($clone->{qualities}) {
        $clone->{qualities} =~ s/^>.*?\n+//s; # Remove comment row
        $clone->{qualities} =~ s/^\s+//s; # Remove any leading spaces        
    }
    
    $contig_statement->execute($clone->{id});
    $clone->{contigs} = $contig_statement->fetchall_arrayref({});

    $logger->debug(
        "Got @{[ scalar @{$clone->{contigs}} ]} contigs, ",
        "sequence with @{[ length $clone->{seq} ]} bp"
    );

    $self->fetch_improved_sequences_for_contigs($clone->{contigs});
}

sub fetch_improved_sequences_for_contigs {
    my $self = shift;
    my ($contigs) = @_;

    $logger->debug(
        "Fetching improved regions for @{[ scalar @$contigs ]} contigs");
    our $statement ||= $self->prepared_statement($IMPROVED_REGION_SQL_KEY);
    for my $contig (@$contigs) {
        $statement->execute($contig->{id});
        $contig->{improved_regions} = $statement->fetchall_arrayref({});

        $logger->debug(
            "   Got @{[ scalar @{ $contig->{improved_regions} } ]} ",
            "improved regions for contig $contig->{name}"
        );
    }
}

sub improved_sequence_max_length {
    my $self = shift;
    my $statement
        = $self->prepared_statement($MAX_IMPROVED_SEQUENCE_LENGTH_SQL_KEY);
    $statement->execute();
    my $max_length = $statement->fetchrow_arrayref()->[0];
    $logger->debug("Improved sequence maximum length: $max_length");
    return $max_length;
}

1;

=pod

=head1 AUTHOR

Shiran Pasternak E<lt>shiran@cshl.eduE<gt>

=cut

=head1 COPYRIGHT

Copyright (c) 2007 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

