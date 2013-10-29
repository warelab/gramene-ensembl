package EnsEMBL::Maize::Object::Rss;

use strict;
use warnings;

use EnsEMBL::Web::Object;
use XML::RSS;
use Data::Dumper qw(Dumper); # For debug

our @ISA = qw(EnsEMBL::Web::Object);

use constant FPC_SPECIES          => 'Zea_mays';
use constant BAC_SPECIES          => 'Zea_mays2';
use constant EXTERNAL_SPECIES     => 'Zea_mays_external';
use constant GRAMENE_MARKERS_LINK => 'http://www.gramene.org/db/markers';

sub rss_xml {
    my $self = shift;

    my $server = "http://".$SiteDefs::ENSEMBL_SERVERNAME;
    
    # get user data
    my $chrom = $self->param('chromosome');
    my $start = $self->param('start');
    my $end = $self->param('end');

    # slice adaptors
    my $fpc_slice_adaptor      = Bio::EnsEMBL::Registry->get_adaptor(FPC_SPECIES, 'core', 'Slice');
    my $bac_slice_adaptor      = Bio::EnsEMBL::Registry->get_adaptor(BAC_SPECIES, 'core', 'Slice');
    my $external_slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor(EXTERNAL_SPECIES, 'core', 'Slice');

    # get the fpc slice
    my $fpc_slice = $fpc_slice_adaptor->fetch_by_region('chromosome', $chrom, $start, $end);
    
    # get the clones in the slice and their db locations
    my ($clones, $clone_locations) = _clone ($fpc_slice);
    
    # get all fpc markers in the slice
    my @marker_features = @{ $fpc_slice->get_all_MarkerFeatures() };
    
    # create the RSS XML
    (my $date, my $time) = _time_stamp();
    my $rss = new XML::RSS (version => '2.0');
    
    # headers
    $rss->channel(title          => "Maize notification for Chromosome $chrom: $start-$end",
		  link           => $server,
		  description    => "Requested data found on Chromosome $chrom: $start-$end"
		  );
    
    # clones
    foreach my $clone (keys %$clones){
	my $clone_link;
	my $bac_slice;
	
	# see where the clone lives and link/slice accordingly

	# for the external case
	if ($clone_locations->{$clone} eq "external"){
	    $clone_link = $server."/".EXTERNAL_SPECIES;
	    $bac_slice = $external_slice_adaptor->fetch_by_region( 'clone', $clone );

	    $rss->add_item(title          => "Clone: $clone",
			   link           => "$clone_link"."/"."contigview?contig=$clone",
			   pubDate        => "$date $time",
			   description    => "$clone maps to chromsome $chrom from $clones->{$clone}[0] to $clones->{$clone}[1]"
			   );
	    my @genes = @{ $bac_slice->get_all_Genes };
	    
	    foreach my $gene (@genes){
		my $stable_id  = $gene->stable_id();
		my $start      = $gene->start();
		my $end        = $gene->end();
		$rss->add_item(title          => "Gene: $stable_id",
			       link           => "$clone_link"."/"."geneview?gene=$stable_id",
			       pubDate        => "$date $time",
			       description    => "Gene $stable_id reported in $clone"
			       );
	    }
	}

	# for the internal case
	else {
	    $clone_link = $server."/".BAC_SPECIES;
	    $bac_slice = $bac_slice_adaptor->fetch_by_region( 'clone', $clone );
	    
	    $rss->add_item(title          => "Clone: $clone",
			   link           => "$clone_link"."/"."contigview?contig=$clone",
			   pubDate        => "$date $time",
			   description    => "$clone maps to chromsome $chrom from $clones->{$clone}[0] to $clones->{$clone}[1]"
			   );
	    
	    # project to the contig coord system
	    my $contig_projection = $bac_slice->project('contig');
	    foreach my $segment ( @{$contig_projection} ){
		
		my $contig_slice = $segment->to_Slice();
		my @genes = @{ $contig_slice->get_all_Genes };
		
		foreach my $gene (@genes){
		    my $stable_id  = $gene->stable_id();
		    my $start      = $gene->start();
		    my $end        = $gene->end();
		    $rss->add_item(title          => "Gene: $stable_id",
				   link           => "$clone_link"."/"."geneview?gene=$stable_id",
				   pubDate        => "$date $time",
				   description    => "Gene $stable_id reported in $clone"
				   );
		}
	    }
	}
    }
    
    # markers
    foreach my $marker (@marker_features){
	my $marker_id    = $marker->display_id();
	my $marker_start = $marker->start();
	my $marker_end   = $marker->end();
	$rss->add_item(title          => "Marker: $marker_id",
                       link           => GRAMENE_MARKERS_LINK."/"."marker_view?marker_name=$marker_id",
                       pubDate        => "$date $time",
                       description    => "$marker_id maps to chromsome $chrom from $marker_start to $marker_end"
                       );
    }
    
    print $rss->as_string;
}

sub _time_stamp {
    my ($d,$t);
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

    $year += 1900;
    $mon++;
    $d = sprintf("%4d-%2.2d-%2.2d",$year,$mon,$mday);
    $t = sprintf("%2.2d:%2.2d:%2.2d",$hour,$min,$sec);
    return($d,$t);
}

sub _clone {
    my $slice = shift;
    my @clones = @{ $slice->get_all_MiscFeatures('bac_map') };

    my %clones;
    my %clone_location;
    foreach my $clone (@clones){
        my $accession = $clone->get_scalar_attribute('embl_acc');

        # get accessioned clones in contigview only                                 
        next unless ($accession ne '');

        my $name     = $clone->get_scalar_attribute('name');
	my $start    = $clone->seq_region_start;
	my $end      = $clone->seq_region_end;
	my $location = $clone->get_scalar_attribute('external');
	
        my @range = ($start, $end);

        $clones{$accession} = [@range];
	if ($location eq 'false'){
	    $clone_location{$accession} = "internal";
	}
	else {
	    $clone_location{$accession} = "external";
	}
    }

    return (\%clones, \%clone_location);
}
