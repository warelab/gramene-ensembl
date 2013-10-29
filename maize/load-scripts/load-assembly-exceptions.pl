#!/usr/local/bin/perl

=head1 NAME

load-assembly-exceptions.pl - Loads the assembly_exception table for analyzed BACs

=head1 SYNOPSIS

perl load-assembly-exceptions.pl [options]

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -v --verbose

=head1 OPTIONS


B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s --species>
  Use this species entry from the registry file [REQUIRED].
  
B<-v --verbose>
  Verbose output

=head1 DESCRIPTION

B<This program> 

    Creates rows in the assembly_exception table for each analyzed BACs in order to view BAC annotations inline.
  
  Maintained by Shiran Pasternak <shiran@cshl.edu>

=cut

use strict;
use warnings;

use FindBin qw($Bin);
use File::Basename qw( dirname );
use Pod::Usage;
use Getopt::Long;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::AssemblyExceptionFeature;
# use Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;

use vars qw($BASEDIR $db_adaptor $verbose);

BEGIN {
    $BASEDIR = dirname($Bin);
}

MAIN: {
    my $help    = 0;
    my $man     = 0;
    my ($species_name, $registry_file);

    GetOptions(
        "help|?"          => \$help,
        "man"             => \$man,
        "species=s"       => \$species_name,
        "registry_file=s" => \$registry_file,
        "verbose"         => \$verbose,
        )
        or pod2usage(2);

    pod2usage(-verbose => 2) if $man;
    pod2usage(1) if $help;

    # Validate file paths
    $registry_file ||= $BASEDIR . '/conf/SiteDefs.pm';

    map {
        -e $_ || (warn("File $_ does not exist\n")    && pod2usage(1));
        -r $_ || (warn("Cannot read $_\n")            && pod2usage(1));
        -f $_ || (warn("File $_ is not plain-text\n") && pod2usage(1));
        -s $_ || (warn("File $_ is empty\n")          && pod2usage(1));
    } $registry_file;

    $species_name || (warn("Need a --species\n") && pod2usage(1));
    Bio::EnsEMBL::Registry->load_all($registry_file);
    $db_adaptor = load_db_adaptor($species_name)
        || (warn("No core DB for $species_name set in $registry_file: $@\n")
        && pod2usage(1));

    my $bac_features = fetch_analyzed_bac_features();
    my $seq_regions  = create_simulated_seq_regions($bac_features);
    map_simulated_seq_regions_to_chromosomes($seq_regions);
    create_assembly_exceptions($seq_regions);
}

=pod

=head2 load_db_adaptor
    Loads a DB adaptor from the registry

=cut

sub load_db_adaptor {
    my ($species) = @_;

    ${ Bio::EnsEMBL::Registry->get_all_DBAdaptors(-species => $species) }[0];
}

=pod

=head2 fetch_analyzed_bac_features
    Returns misc_features for analyzed BACs.

=cut

sub fetch_analyzed_bac_features {
    my $misc_feature_adaptor = $db_adaptor->get_adaptor('MiscFeature');
    return $misc_feature_adaptor->fetch_all_by_attribute_type_value(
        'annotation', 'II');

}

=pod

=head2 create_simulated_seq_regions
    Create a new seq_region for each analyzed bac so that we can properly display assembly_exceptions

=cut

sub create_simulated_seq_regions {
    my ($bac_features) = @_;
    my $slice_adaptor  = $db_adaptor->get_adaptor('Slice');
    my $seq_regions    = {};
    for my $bac (@$bac_features) {
        my $accession      = $bac->get_scalar_attribute('embl_acc');
        my $chromosome     = $bac->slice;
        my $simulated_name = "chr@{[$chromosome->seq_region_name]}_$accession";
        my $slice          = Bio::EnsEMBL::Slice->new(
            -coord_system      => $chromosome->coord_system,
            -seq_region_name   => $simulated_name,
            -start             => $chromosome->start,
            -end               => $chromosome->end,
            -strand            => $chromosome->strand,
            -seq_region_length => $chromosome->seq_region_length,
            -adaptor           => $slice_adaptor,
        );

        # $Data::Dumper::Maxdepth = 2;
        # print Data::Dumper::Dumper($slice);
        my $seq_region_id = undef;
        eval {
            $seq_region_id = $slice_adaptor->get_seq_region_id($slice);
        };
        if (!$@) {
            if ($verbose) {
                print "Simulated slice $seq_region_id: @{[$slice->name]} already stored. Skipping.\n";
            }
        } else {
            $seq_region_id = $slice_adaptor->store($slice);
            if ($verbose) {
                print "Added simulated slice $seq_region_id : @{[$slice->name]}\n";
            }
        }
        $seq_regions->{$seq_region_id} = +{
            'slice' => $slice,
            'bac'   => $bac
        };
    }
    return $seq_regions;
}

=pod

=head2 map_simulated_seq_regions_to_chromosomes
    Map the simulated seq_regions to their respective chromosomes in the assembly table

=cut

sub map_simulated_seq_regions_to_chromosomes {
    my ($seq_regions) = @_;
    my $sql = <<__END_SQL__;
INSERT INTO assembly
   (asm_seq_region_id, cmp_seq_region_id,
    asm_start, asm_end, cmp_start, cmp_end, ori) 
VALUES (?,?,?,?,?,?,?)
__END_SQL__
    my $check_sql = <<__END_SQL__;
SELECT COUNT(*) FROM assembly
 WHERE asm_seq_region_id = ? AND cmp_seq_region_id = ?
__END_SQL__
    my $statement = $db_adaptor->dbc->prepare($sql);
    my $check_statement = $db_adaptor->dbc->prepare($check_sql);
    for my $id (keys %$seq_regions) {
        my $slice       = $seq_regions->{$id}->{'slice'};
        my $bac         = $seq_regions->{$id}->{'bac'};
        my $chrom       = $bac->slice;
        my $orientation = ($slice->strand == $chrom->strand) ? 1 : -1;
        my $index       = 0;
        $check_statement->execute($chrom->get_seq_region_id, $id);
        my $exists = $check_statement->fetchrow_arrayref->[0];
        print "EXISTS? $exists\n";
        next;
        if ($exists == 0) {
            $statement->bind_param(++$index, $chrom->get_seq_region_id);
            $statement->bind_param(++$index, $id);
            $statement->bind_param(++$index, $bac->start);
            $statement->bind_param(++$index, $bac->end);
            $statement->bind_param(++$index, $bac->start);
            $statement->bind_param(++$index, $bac->end);
            $statement->bind_param(++$index, $orientation);
            my $result = $statement->execute ||
                die "Cannot insert into DB: @{[$statement->errstr]}\n";

            if ($verbose) {
                print "Mapped slice ",
                    join(' to ',
                    map { $_->name, " [", $_->start, " ", $_->end, "]" }
                        ($slice, $chrom)),
                    "\n";
            }
        } else {
            if ($verbose) {
                print "Mapping already exists for @{[$chrom->name]}=>@{[$slice->name]}\n";
            }
        }
    }
}

=pod

=head2 create_assembly_exceptions
    Create assembly exceptions for each of the simulated seq_regions

=cut

sub create_assembly_exceptions {
    my ($seq_regions) = @_;
    my $adaptor = $db_adaptor->get_adaptor('AssemblyExceptionFeature');
    for my $id (keys %$seq_regions) {
        my $slice     = $seq_regions->{$id}->{'slice'};
        my $bac       = $seq_regions->{$id}->{'bac'};
        my $chrom     = $bac->slice;
        my $exception = Bio::EnsEMBL::AssemblyExceptionFeature->new(
            -start           => $bac->start,
            -end             => $bac->end,
            -type            => 'HAP',
            -slice           => $chrom,
            -alternate_slice => $slice,
        );
        $adaptor->store($exception);
        if ($verbose) {
            print
                "Added exception: Slice=@{[$exception->slice->name]} Alternate=@{[$exception->alternate_slice->name]} at [@{[$exception->start]} @{[$exception->end]}]\n";
        }
    }
}
