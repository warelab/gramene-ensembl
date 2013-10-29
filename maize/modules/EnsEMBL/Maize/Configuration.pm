package EnsEMBL::Maize::Configuration;

use strict;
use EnsEMBL::Web::Configuration;

our @ISA = qw( EnsEMBL::Web::Configuration );

sub context_location {
    my $self = shift;
    warn("New module?");
    my $obj = $self->{object};
    return unless $obj->can('location_string');
    my $species = $obj->species;
    my ($q_string, $header) = $obj->location_string;
    $header
        = "@{[$obj->seq_region_type_and_name]}<br />@{[$obj->thousandify(floor($obj->seq_region_start))]}";
    if (floor($obj->seq_region_start) != ceil($obj->seq_region_end)) {
        $header .= " - @{[$obj->thousandify(ceil($obj->seq_region_end))]}";
    }
    my $flag = "location$self->{flag}";
    return if $self->{page}->menu->block($flag);
warn("Overriding dynamic menus");
    if ($q_string) {
        my $flag = "location$self->{flag}";
        $self->add_block($flag, 'bulletted', $header, 'raw' => 1);   ##RAW HTML!
        $header =~ s/<br \/>/ /;
        if ($self->mapview_possible($obj->seq_region_name)) {
            $self->add_entry(
                $flag,
                'text'  => "View of @{[$obj->seq_region_type_and_name]}",
                'href'  => "/$species/mapview?chr=" . $obj->seq_region_name,
                'title' => 'MapView - show chromosome summary'
            );
        }

        # $self->add_entry( $flag, 'text' => 'Graphical view',
        #   'href'=> "/$species/contigview?l=$q_string",
        #   'title'=> "ContigView - detailed sequence display of $header" );
        $self->add_entry(
            $flag,
            'text'  => 'Graphical overview',
            'href'  => "/$species/cytoview?l=$q_string",
            'title' => "CytoView - sequence overview of $header"
        );
        $self->add_entry(
            $flag,
            'text'  => 'Export information about region',
            'title' => "ExportView - export information about $header",
            'href'  => "/$species/exportview?l=$q_string"
        );
        $self->add_entry(
            $flag,
            'text'  => 'Export sequence as FASTA',
            'title' => "ExportView - export sequence of $header as FASTA",
            'href'  =>
                "/$species/exportview?l=$q_string;format=fasta;action=format"
        );
        $self->add_entry(
            $flag,
            'text'  => 'Export EMBL file',
            'title' => "ExportView - export sequence of $header as EMBL",
            'href'  =>
                "/$species/exportview?l=$q_string;format=embl;action=format"
        );
        if ($obj->species_defs->multidb('ENSEMBL_MART_ENSEMBL')
            && !$obj->species_defs->ENSEMBL_NOMART)
        {
            $self->add_entry(
                $flag,
                'icon'  => '/img/biomarticon.gif',
                'text'  => 'Export Gene info in region',
                'title' => "BioMart - export Gene information in $header",
                'href'  => "/$species/martlink?l=$q_string;type=gene_region"
            );
        }
        if ($obj->species_defs->multidb('ENSEMBL_MART_SNP')) {
            $self->add_entry(
                $flag,
                'icon'  => '/img/biomarticon.gif',
                'text'  => 'Export SNP info in region',
                'title' => "BioMart - export SNP information in $header",
                'href'  => "/$species/martlink?l=$q_string;type=snp_region"
                )
                if $obj->species_defs->databases->{'ENSEMBL_VARIATION'};
        }
        if ($obj->species_defs->multidb('ENSEMBL_MART_VEGA')) {
            $self->add_entry(
                $flag,
                'icon'  => '/img/biomarticon.gif',
                'text'  => 'Export Vega info in region',
                'title' => "BioMart - export Vega gene features in $header",
                'href'  => "/$species/martlink?l=$q_string;type=vega_region"
                )
                if $obj->species_defs->databases->{'ENSEMBL_VEGA'};
        }
    }
}

sub really_delete_menu_block {
    my ($menu, $code) = @_;
    my $priority = $menu->{'blocks'}{$code}{'options'}{'priority'};
    $menu->delete_block($code);
    $menu->{'block_order'}[$priority]
        = [ grep { $_ ne $code } @{ $menu->{'block_order'}[$priority] } ];
}

