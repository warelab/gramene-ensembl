#!/usr/local/bin/perl

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

my $logger = Log::Log4perl->get_logger('bacstats');

my %QC = ();

BEGIN {

    # Set the perl libraries
    $BASEDIR = dirname($Bin);
    unshift @INC, $BASEDIR . '/ensembl-live/ensembl/modules';
    Log::Log4perl::init($BASEDIR . '/conf/log4perl.conf');
}

MAIN: {
    $OUTPUT_AUTOFLUSH = 1;

    my $dbh = get_db_handle();

    my %sql = ();
    map { $sql{$_} = sql_for_feature_with_seq_region($_) }
        qw(gene transcript exon);
    $sql{'translation'} = <<SQL;
select si.translation_id, si.stable_id
  from translation_stable_id si
       left join translation tx on si.translation_id = tx.translation_id
       left join transcript ts on ts.transcript_id = tx.transcript_id
       left join transcript_stable_id tsi on tsi.transcript_id = ts.transcript_id
 order by tsi.stable_id
SQL

    for my $feature (qw(gene transcript exon translation)) {
        update_stable_id_for_features(
            {   'feature'    => $feature,
                'dbh'        => $dbh,
                'select_sql' => $sql{$feature},
            }
        );
    }

    # print "QC\n";
    # for my $i (map { $_->[1] }
    #     sort { $a->[0] <=> $b->[0] } map { [ (split ':')[1], $_ ] } keys %QC)
    # {
    #     my $qc = $QC{$i};
    #     print "$i\t", join("\t", map {"$_=$qc->{$_}"} sort keys %$qc), "\n";
    # }

    $dbh->disconnect();
    exit;
}

#======================================================================

=pod

=head2 get_clonepath_handle
    Creates a connection to the clonepath database

=cut

sub get_db_handle {
    my $username    = 'maize_rw';
    my $password    = 'z3@m@ys';
    my $data_source = sprintf(
        'DBI:mysql:database=%s;host=%s;port=%s',
        'zea_mays_core_43_bac',
        'ascutney.cshl.edu', 3307
    );
    my $db_handle = undef;
    eval { $db_handle = DBI->connect($data_source, $username, $password); };
    if (!defined($db_handle)) {
        die "Can't connect to database: @{[ $DBI::errstr ]}\n";
    }
    return $db_handle;
}

sub update_stable_id_for_features {
    my ($params)   = @_;
    my $dbh        = $params->{'dbh'};
    my $feature    = $params->{'feature'};
    my $select_sql = $params->{'select_sql'};

    my $pk_field = "${feature}_id";

    # print '=' x 60, "\n";
    # print uc($feature), "\n";
    # print "SQL:\n$select_sql\n";
    my $rows = $dbh->selectall_arrayref($select_sql, { Slice => {} });

    my @index_updates = ();

    my $new_index   = 0;
    my $last_prefix = '';
    my $last_index  = 0;
    for my $row (@$rows) {
        my $pk        = $row->{$pk_field};
        my $stable_id = $row->{'stable_id'};
        my ($prefix, $index, $suffix)
            = ($stable_id =~ m/^(\w+\.\d+_[A-Z]+)(\d+)(.*)$/);

        if ($prefix ne $last_prefix) {
            $new_index = 0;
        }
        $new_index++;

        my $new_stable_id
            = sprintf('%s%0.3d%s', $prefix, $new_index, $suffix);

        $last_prefix = $prefix;
        push @index_updates,
            +{
            'stable_id'     => $new_stable_id,
            'pk'            => $pk,
            'old_stable_id' => $stable_id
            };

        my ($qc_key) = ($prefix =~ m/^(.*)_/);
        $qc_key = sprintf('%s:%0.3d', $qc_key, $new_index);

        $QC{$qc_key} ||= {};
        $QC{$qc_key}->{$feature} = $new_stable_id;
    }

    # In order to prevent duplicate stable_id errors
    # $dbh->do("UPDATE ${feature}_stable_id SET stable_id = NULL")
    # || die
    # "Cannot clean stable_ids in ${feature}_stable_id: $DBI::errstr\n";

    my $update_statement = $dbh->prepare(
        "UPDATE ${feature}_stable_id SET stable_id = ? WHERE $pk_field = ?");

    for my $update (@index_updates) {
        print join("\t",
            $feature, $update->{'pk'},
            $update->{'old_stable_id'},
            $update->{'stable_id'}),
            "\n";
        eval {
            $update_statement->bind_param(1, $update->{'stable_id'});
            $update_statement->bind_param(2, $update->{'pk'});
            $update_statement->execute();
        };
        if ($EVAL_ERROR) {
            die "Cannot perform update: @{[$dbh->errstr]}\n";
        }
    }
}

sub sql_for_feature_with_seq_region {
    my ($feature) = @_;
    my $sql = <<SQL;
select f.${feature}_id, si.stable_id, contig.name contig_name, clone.name clone_name,
       a.asm_start + f.seq_region_start as clone_start
  from ${feature} f
       left join ${feature}_stable_id si on f.${feature}_id = si.${feature}_id
       left join seq_region contig on f.seq_region_id = contig.seq_region_id
       left join assembly a on a.cmp_seq_region_id = contig.seq_region_id
       left join seq_region clone on clone.seq_region_id = a.asm_seq_region_id
 order by clone.name, clone_start
SQL
    return $sql;
}

1;
