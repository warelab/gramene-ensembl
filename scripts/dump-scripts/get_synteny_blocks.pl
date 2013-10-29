#!/usr/local/bin/perl

# vim: tw=78: sw=4: ts=4: et: 

use strict;
use warnings;

BEGIN {
    $ENV{'EnsemblDir'} ||= '/usr/local/ensembl-live/';
}

use lib map { $ENV{'EnsemblDir'} . "/$_" }
    qw ( bioperl-live modules ensembl/modules ensembl-external/modules
    ensembl-variation/modules ensembl-draw/modules ensembl-compara/modules );

use Bio::EnsEMBL::Registry;
use English qw( -no_match_vars );
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Readonly;

my $registry_file = '/usr/local/gramene/conf/ensembl.registry';
my ( $help, $man_page );
GetOptions(
    'r|registry:s' => \$registry_file,
    'help'         => \$help,
    'man'          => \$man_page,
) or pod2usage(2);

if ( $help || $man_page ) {
    pod2usage({
        -exitval => 0,
        -verbose => $man_page ? 2 : 1
    });
}; 

if ( scalar @ARGV != 2 ) {
    die "Need two species\n";
}

if ( !-e $registry_file ) {
    pod2usage("Registry file '$registry_file' does not exist!");
}

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $compara_db = Bio::EnsEMBL::Registry->get_DBAdaptor('compara', 'compara');
my $gdba       = $compara_db->get_adaptor('GenomeDB');
my $ synrega   = $compara_db->get_adaptor('SyntenyRegion');

my %genome_db_by_name = map { $_->{'name'}, $_ } @{ $gdba->fetch_all }
    or die "No genome dbs found!\n";

my @genome_dbs;
for my $species ( @ARGV ) {
    $species =~ s/_/ /g;

    my $db = $genome_db_by_name{ $species } or die "No db named '$species'\n";

    push @genome_dbs, $db;
}

my $mlssa = $compara_db->get_adaptor( 'MethodLinkSpeciesSet' );
my $mlss  = $mlssa->fetch_by_method_link_type_GenomeDBs( 
    'SYNTENY', [@genome_dbs] 
);

print join("\t", qw[ 
    species1 chr1 start1 stop1 strand1 
    species2 chr2 start2 stop2 strand2 
]), "\n";

my $num_records = 0;
for my $synreg ( @{ $synrega->fetch_all_by_MethodLinkSpeciesSet( $mlss ) } ) {
    $num_records++;
    my @result;
    for my $dnafrag_region ( @{ $synreg->get_all_DnaFragRegions() } ) {
        push @result, 
            $dnafrag_region->genome_db->name,
            $dnafrag_region->dnafrag->name,
            $dnafrag_region->dnafrag_start,
            $dnafrag_region->dnafrag_end,
            $dnafrag_region->dnafrag_strand;
    }

    print join( "\t", @result ), "\n";
}

print STDERR "Done, exported $num_records records.\n";

__END__

# ----------------------------------------------------

=pod

=head1 NAME

get_synteny_blocks.pl - a script

=head1 SYNOPSIS

  get_synteny_blocks.pl [options] Oryza_sativa Sorghum_bicolor > os_sb.tab

Options:

  -r|--registry  Ensembl registry file, default is
                 /usr/local/gramene/conf/ensembl.registry
  --help         Show brief help and exit
  --man          Show full documentation

=head1 DESCRIPTION

Exports the syntenic blocks from Ensembl's Compara database for two 
given species in a tab-delimited format like so:

  ************ Record 1 ************
  species1: Oryza sativa
      chr1: 10
    start1: 10765474
     stop1: 10824327
   strand1: 0
  species2: Brachypodium distachyon
      chr2: 1
    start2: 3882483
     stop2: 3921995
   strand2: 1

=head1 SEE ALSO

Ensembl Compara.

=head1 AUTHORS

Jerry Lu E<lt>luj@cshl.eduE<gt>,
Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2010 Cold Spring Harbor Laboratory

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
