package EnsEMBL::Maize::Component::Chromosome;

use EnsEMBL::Web::Component::Chromosome;
use EnsEMBL::Maize::Util::FPC;

our @ISA = qw(EnsEMBL::Web::Component::Chromosome);

# Returns a display of the feature counts for each track
sub stats {
    my $panel = shift;
    my $chr   = shift;

    my $chr_name = $chr->chr_name;
    my $label    = "Chromosome $chr_name";

    my $table_t = qq(
<table cellspacing="0" cellpadding="0" border="0" class="stats">
%s
</table>
    );
    my $row_t = qq(
<tr>
    <th>%s</th>
    <td align="right">%s</td>
    <td>%s</td>
</tr>
);

    my @stat_data = ([ 'Chromosome length', $chr->length, 'bps' ]);
    _add_fpc_stats($chr, \@stat_data);
    for my $attrib (@{ $chr->Obj->get_all_Attributes() }) {
        push @stat_data, [ $attrib->name, $attrib->value ];
    }

    my @sorted_stat_data =    #sort { $a->[0] cmp $b->[0] }
        @stat_data;

    my @rows;
    for my $stat (@sorted_stat_data) {
        my $label  = $stat->[0];
        my $number = $stat->[1];
        my $units  = $stat->[2];
        $number = $chr->thousandify($number);
        my $row_html = sprintf($row_t, $label, $number, ($units || '&nbsp;'));
        push @rows, $row_html;
    }
    my $html = sprintf($table_t, join("", @rows));

    $panel->add_row($label, $html);
    return 1;
}

sub _add_fpc_stats {
    my ($chromosome, $stats) = @_;
    my $slice   = $chromosome->Obj;
    my $fpcutil = EnsEMBL::Maize::Util::FPC->new;
    push @$stats,
        [ 'FPC Contigs', scalar @{ $fpcutil->fetch_contigs($chromosome) } ];
    push @$stats,
        [ 'FPC Clones', scalar @{ $fpcutil->fetch_all_fpc_clones($chromosome) } ];
    push @$stats,
        [
        'Accessioned BACs',
        scalar @{ $fpcutil->fetch_accessioned_bacs($chromosome) }
        ];
}

sub disclaimer {
    my $panel = shift;
    my $chr   = shift;
    if (my $disclaimer = $chr->species_defs->ASSEMBLY_DISCLAIMER) {
        $panel->add_row($disclaimer);
    }
    return 1;
}

sub chr_map {
    my ($panel, $object) = @_;
    my $config_name = 'Vmapview';
    my $species     = $object->species;
    my $chr_name    = $object->chr_name;

    # Maize specific change - Ideogram links to only cytoview at this point
    my $script = 'cytoview';

    #    my $script           = 'contigview';
    my $seq_region_width = '2500000';

    #    if (defined $object->species_defs->CONTIGVIEW_ENABLED
    #        and $object->species_defs->CONTIGVIEW_ENABLED == 0)
    #    {
    #        $script = 'cytoview';
    #        $seq_region_width *= 25;
    #    }

    my $config = $object->get_userconfig($config_name);

    #    warn Data::Dumper::Dumper($config);

    my $ideo_height = $config->{'_image_height'};
    my $top_margin  = $config->{'_top_margin'};
    my $hidden      = {
        'seq_region_name'  => $chr_name,
        'seq_region_width' => $seq_region_width,
        'seq_region_left'  => '1',
        'seq_region_right' => $object->length,
        'click_right'      => $ideo_height + $top_margin,
        'click_left'       => $top_margin,
    };

    # make a drawable container
    my $image = $object->new_karyotype_image();

    #    warn Data::Dumper::Dumper($config->{'general'}->{'Vmapview'});

    $image->imagemap   = 'no';
    $image->cacheable  = 'yes';
    $image->image_name = 'mapview-' . $species . '-' . $chr_name;
    $image->set_button(
        'form',
        'id'     => 'vclick',
        'URL'    => "/$species/$script",
        'hidden' => $hidden
    );
    $image->add_tracks($object, $config_name);
    $image->karyotype($object, '', $config_name);
    $image->caption
        = 'Click on the image above to zoom into that point on CytoView';

    #    warn Data::Dumper::Dumper($image->{'general'});

    $panel->add_image($image->render, $image->{'width'});
    return 1;

}

sub jump_to_contig {
    my $panel = shift;
    my ($chromosome) = @_;

    EnsEMBL::Maize::Component::GenomeEntryPoints::fpc_form($panel, $chromosome);
}

1;
