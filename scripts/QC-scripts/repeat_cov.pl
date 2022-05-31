#! /usr/bin/perl -w

=head1 NAME

repeats_cov.pl

=head1 SYNOPSIS

 repeats_cov.pl 

Options:

 --help      Show brief help and exi
 --registry  the ensembl registry file
 --species   the species in the registry file
 --logic_name the logic_name of the repeatMask analysis, could be more than one or none
--detail      print details of the computation
 
=head1 DESCRIPTION

=head1 SEE ALSO

perl.

=head1 AUTHOR

Bonnie Hurwitz E<lt>hurwitz@cshl.eduE<gt>,

=head1 COPYRIGHT

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::RepeatConsensus;
use Getopt::Long;
use Pod::Usage;

use Readonly;

Readonly my $CHROM_STAT_TMPL         => "%-15s %9d %9d %9d %6.2f%%\n";
Readonly my $TOTAL_DNA               => 'Total DNA';
Readonly my $MASKED_DNA              => 'Masked DNA';
Readonly my $NATIVE_NS              => 'Native Ns';
Readonly my $CHROMOSOME_COORD_SYSTEM => 'chromosome';

our $ENS_DBA;


sub count_Ns {
    my ($seqref) = shift;

    my $cnt = $seqref =~ s/N/N/g;
    
    return $cnt;
}


MAIN: {


    my ($help, $registry, $species, %sums, @logic_names, $detail);
   

    GetOptions (
		"help"       => \$help,
		"registry=s" => \$registry,
		"species=s"  => \$species,
		"logic_name=s" => \@logic_names,
		"detail"     => \$detail,
) or pod2usage;
    
    if( $help || !$registry || !$species){
	
	pod2usage(-verbose => 2);
    }
    
    Bio::EnsEMBL::Registry->load_all( $registry );
    $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
    $ENS_DBA || ( print( "\n[*DIE] No $species set in $registry\n\n" ) &&
		  pod2usage() );


    my $chromosomes = fetch_chromosomes();
    warn("Fetched ", scalar @$chromosomes, " toplevel seq regions");

    my %chrom_sums   = ();
    $sums{$TOTAL_DNA} = 0;
    for my $chrom (@$chromosomes) {
        my $seq     = $chrom->seq();
        my $n_count = count_Ns($seq);

	warn (
            "seq region @{[$chrom->seq_region_name]} has $n_count native Ns");
        $sums{$TOTAL_DNA} += $chrom->length;

	warn ("Getting all repeats");
	my $masked;
	if( scalar @logic_names > 0){
		$masked = $chrom->get_repeatmasked_seq(\@logic_names)->seq();
	}else{	
        	$masked = $chrom->get_repeatmasked_seq()->seq();
	}

	my $total_Ns = count_Ns($masked);
	print "total_Ns=$total_Ns, nativeNs=$n_count\n" if $detail;
	my $repeat_ns = $total_Ns - $n_count;
	$chrom_sums{ $chrom->seq_region_name }
	= { length => $chrom->length, masked => $repeat_ns, native_ns => $n_count };
	$sums{$MASKED_DNA} += $repeat_ns;
	$sum{$NATIVE_NS} += $n_count;
    
    }

    for my $chrom_name (sort keys %chrom_sums) {
	my $length = $chrom_sums{$chrom_name}->{length};
	my $masked = $chrom_sums{$chrom_name}->{masked};
	my $native_ns = $chrom_sums{$chrom_name}->{native_ns};
	print (
	       sprintf($CHROM_STAT_TMPL,
		       $chrom_name, $length, $masked, $native_ns, 100.0 * ($masked / ($length-$native_ns)))
	       ) if $detail;
    }
    print (
            sprintf($CHROM_STAT_TMPL,
                "$species-Total", $sums{$TOTAL_DNA}, $sums{$MASKED_DNA}, $sums{$NATIVE_NS},
                100.0 * ($sums{$MASKED_DNA} / ($sums{$TOTAL_DNA})-$sums{$NATIVE_NS}))
        );

}

sub fetch_chromosomes {

    my $slice_adaptor = $ENS_DBA->get_adaptor('Slice');

    warn("Getting all toplevel regions");
    my $slices = $slice_adaptor->fetch_all('toplevel');
    return $slices;
    
}
