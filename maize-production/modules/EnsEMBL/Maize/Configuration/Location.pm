package EnsEMBL::Maize::Configuration::Location;

use strict;
use warnings;

use EnsEMBL::Web::Configuration;
use EnsEMBL::Web::Tools::Ajax;

use CGI;
use EnsEMBL::Maize::Configuration;
use EnsEMBL::Maize::Util::FPC;
use POSIX qw(floor ceil);

our @ISA = qw( EnsEMBL::Web::Configuration );

#----------------------------------------------------------------------
# Maize-specific ContigView
sub contigview {
    my $self     = shift;
    my $document = $self->{page};
    my $content  = $document->content;

    my $panel;
    foreach my $code ($content->panels) {
        if ($code =~ m/bottom/) {
            $panel = $content->panel($code);
            last;
        }
    }

    if ($panel) {

        # Replace menu with Maize's (EST and GSS options)
        #    $panel->remove_component('menu');
        $panel->replace_component('menu' =>
                'EnsEMBL::Maize::Component::Location::contigviewbottom_menu');
    }

    return 1;

}

sub is_fpc_database {
    return lc $ENV{'ENSEMBL_SPECIES'} eq 'zea_mays';
}

sub context_menu {
    my $self    = shift;
    my $obj     = $self->{object};
    my $species = $obj->real_species;

    return unless $self->{page}->can('menu');
    my $menu = $self->{page}->menu;
    return unless $menu;
    my $script
        = $ENV{'ENSEMBL_SCRIPT'} eq 'cytoview' ? 'CytoView' : 'ContigView';
    my $q_string = sprintf('%s:%s-%s',
        $obj->seq_region_name, $obj->seq_region_start, $obj->seq_region_end);
    my $flag = join(q{}, 'contig', ($self->{'flag'} || q{}));

    my $header = join('<br/>',
        $obj->seq_region_type_and_name,
        $obj->thousandify(floor($obj->seq_region_start)));

    if (floor($obj->seq_region_start) != ceil($obj->seq_region_end)) {
        $header .= " - @{[$obj->thousandify(ceil($obj->seq_region_end))]}";
    }

    EnsEMBL::Maize::Configuration::really_delete_menu_block($menu, $flag);
    $menu->add_block($flag, 'bulleted', $header, 'raw' => 1);
    if ($self->mapview_possible($obj->seq_region_name)) {
        $menu->add_entry(
            $flag,
            'code' => 'mv_link',
            'text' => "@{[$obj->seq_region_type_and_name]} in MapView",
            'title' =>
                "MapView - Overview of @{[$obj->seq_region_type_and_name]} including feature sumarries",
            'href' => "/$species/mapview?chr=" . $obj->seq_region_name
        );
    }
    $header =~ s/<br \/>/ /;

    if ($script eq 'ContigView') {
        my $cytoview_species = 'Zea_mays';
        my $accession        = $obj->seq_region_name;
        my $ctg_name = EnsEMBL::Maize::Util::FPC->fetch_contig_by_accession(
            $accession);
        $self->add_other_versions_menu_item($menu, $accession, $species, $flag);
        if ($ctg_name) {
            $menu->add_entry(
                $flag,
                'text'  => "View FPContig $ctg_name in CytoView",
                'title' => "CytoView - genome browser overview of $ctg_name",
                'href'  => "/$cytoview_species/cytoview?mapfrag=$ctg_name"
            );
        }
        my $clone_name
            = EnsEMBL::Maize::Util::FPC->fetch_clone_by_accession($accession);
        if ($clone_name) {
            $menu->add_entry(
                $flag,
                'text' => "View Clone $clone_name in CytoView",
                'title' =>
                    "CytoView - genome browser overview of $clone_name",
                'href' => "/$cytoview_species/cytoview?mapfrag=$clone_name"
            );
        }
        $menu->add_entry(
            $flag,
            'text'  => "View GenBank record for $accession",
            'title' => "Jump to GenBank record of $accession",
            'href' =>
                "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&val=${accession}&dopt=brief",
        );

    }
    $menu->add_entry(
        $flag,
        'text'  => 'Export information about region',
        'title' => "ExportView - export information about $header",
        'href'  => "/$species/exportview?l=$q_string"
    );
    unless ($obj->species_defs->ENSEMBL_NOMART) {
        $menu->add_entry(
            $flag,
            'icon'  => '/img/biomarticon.gif',
            'text'  => 'Export Gene info in region',
            'title' => "BioMart - export Gene information in $header",
            'href'  => "/$species/martlink?l=$q_string;type=gene_region"
        );
        $menu->add_entry(
            $flag,
            'icon'  => '/img/biomarticon.gif',
            'text'  => 'Export SNP info in region',
            'title' => "BioMart - export SNP information in $header",
            'href'  => "/$species/martlink?l=$q_string;type=snp_region"
        ) if $obj->species_defs->databases->{'ENSEMBL_VARIATION'};
        $menu->add_entry(
            $flag,
            'icon'  => '/img/biomarticon.gif',
            'text'  => 'Export Vega info in region',
            'title' => "BioMart - export Vega gene features in $header",
            'href'  => "/$species/martlink?l=$q_string;type=vega_region"
        ) if $obj->species_defs->databases->{'ENSEMBL_VEGA'};
    }
    my $URL = qq(/$species/urlsource?l=$q_string;script=) . lc($script);
    $menu->add_entry(
        $flag,
        'text' => "View URL based data on " . $script,
        'text' => 'View data stored on another webserver in ' . $script,
        'href' =>
            qq(javascript:X=window.open('$URL','urlsources','left=10,top=10,height=400,width=750,scrollbars=yes');X.focus()),
        ''
    );

    my @options_as = ();

    my %alignments = $obj->species_defs->multiX('ALIGNMENTS');

    for my $id (
        sort {
            10 * ($alignments{$b}->{'type'} cmp $alignments{$a}->{'type'})
                + ($a <=> $b)
        }
        grep {
            $alignments{$_}->{'species'}->{$species}
        }
        keys(%alignments)
        )
    {

        my $label = $alignments{$id}->{'name'};
        my $KEY   = "opt_align_${id}";

        my @species = grep { $_ ne $species }
            sort keys %{ $alignments{$id}->{'species'} };
        if (scalar(@species) == 1) {
            ($label = $species[0]) =~ s/_/ /g;
        }

        push @options_as,
            {
            'text' => "... <em>$label</em>",
            'raw'  => 1,
            'href' => sprintf(
                "/%s/alignsliceview?c=%s:%s;w=%s;align=%s",
                $species,          $obj->seq_region_name,
                $obj->centrepoint, $obj->length,
                $KEY
            )
            };

    }

    if (@options_as) {
        $menu->add_entry(
            $flag,
            'text'    => "View alignment with ...",
            'href'    => $options_as[0]{'href'},
            'options' => \@options_as,
            'title'   => "AlignSliceView - graphical view of alignment"
        );
    }

    my %species
        = (map { $obj->species_defs->multi($_, $species) }
            qw(BLASTZ_RAW BLASTZ_NET BLASTZ_RECIP_NET PHUSION_BLASTN TRANSLATED_BLAT BLASTZ_GROUP)
        );
    my @options = ();
    for (sort keys %species) {
        (my $HR = $_) =~ s/_/ /;
        push @options,
            {
            'text' => "... <em>$HR</em>",
            'raw'  => 1,
            'href' => sprintf(
                "/%s/multicontigview?s1=%s;c=%s:%s;w=%s",
                $species,              $_,
                $obj->seq_region_name, $obj->centrepoint,
                $obj->length
            )
            };
    }
    if (@options) {
        $menu->add_entry(
            $flag,
            'code'    => "mcv_link",
            'text'    => "View alongside ...",
            'href'    => $options[0]{'href'},
            'options' => \@options,
            'title' =>
                "MultiContigView - side by side view of genomic sequence"
        );
    }

    if (@{  $obj->species_defs->other_species(
                $species, 'ENSEMBL_CHROMOSOMES'
                )
                || []
        }
        )
    {
        my %species = ($obj->species_defs->multi('SYNTENY', $species));
        my @options = ();
        foreach (sort keys %species) {
            (my $HR = $_) =~ s/_/ /;
            push @options,
                {
                'text' => "... with <em>$HR</em>",
                'raw'  => 1,
                'href' => sprintf(
                    "/%s/syntenyview?otherspecies=%s;chr=%s;loc=%s",
                    $species,              $_,
                    $obj->seq_region_name, $obj->centrepoint
                )
                }
                if @{ $obj->species_defs->other_species($_,
                    'ENSEMBL_CHROMOSOMES')
                    || [] };
        }
        if (@options) {
            $menu->add_entry(
                $flag,
                'text'    => 'View Syntenic regions ...',
                'href'    => $options[0]{'href'},
                'options' => \@options
            );
        }
    }

    my %browsers = %{
        $obj->species_defs->other_species($species,
            'EXTERNAL_GENOME_BROWSERS')
            || {}
        };
    foreach (sort keys %browsers) {
        $menu->add_entry(
            $flag,
            'text' => "View region in $browsers{$_}",
            'href' => $obj->get_ExtURL(
                $_,
                {   'CHR'   => $obj->seq_region_name,
                    'START' => $obj->seq_region_start,
                    'END'   => $obj->seq_region_end
                }
            )
        );
    }
}

=pod

=head2 add_other_versions_menu_item
    Adds links to other versions of this BAC

=cut

sub add_other_versions_menu_item {
    my $self = shift;
    my ($menu, $accession, $species, $flag) = @_;
    my @archive
        = @{ EnsEMBL::Maize::Util::FPC->archive_for_accession($accession) };
    my @other_versions
        = map {
        +{  'text' => $_,
            'href' => "/$species/$ENV{'ENSEMBL_SCRIPT'}?region=$_"
            }
        } @archive;

    if (!scalar @other_versions) {
        push @other_versions, +{ 'text' => 'No other versions available' };
    }

    $menu->add_entry(
        $flag,
        'href'    => '',
        'title'   => 'View other versions of this clone',
        'options' => \@other_versions,
        'raw'     => 1,
        'text'    => <<HTML,
Other Versions&nbsp;<img src="/images/new.png" alt="[NEW]" style="vertical-align: bottom"/>
HTML
    );
}

sub export_step1 {
    ### Alternative context menu for step 1 of exportview
    my $self    = shift;
    my $obj     = $self->{object};
    my $species = $obj->real_species;
    return unless $self->{page}->can('menu');
    my $menu = $self->{page}->menu;
    return unless $menu;

    my $flag = 'species';
    $menu->add_block($flag, 'bulleted', 'Export a different species',
        'raw' => 1);

    my @species_inconf = @{ $obj->species_defs->ENSEMBL_SPECIES };
    my @group_order    = qw( Mammals Chordates Eukaryotes );
    my %spp_tree       = (
        'Mammals'    => { 'label' => 'Mammals',          'species' => [] },
        'Chordates'  => { 'label' => 'Other chordates',  'species' => [] },
        'Eukaryotes' => { 'label' => 'Other eukaryotes', 'species' => [] },
    );

    foreach my $sp (@species_inconf) {
        my $bio_name
            = $obj->species_defs->other_species($sp, "SPECIES_BIO_NAME");
        my $group = $obj->species_defs->other_species($sp, "SPECIES_GROUP")
            || 'default_group';
        unless ($spp_tree{$group}) {
            push @group_order, $group;
            $spp_tree{$group} = { 'label' => $group, 'species' => [] };
        }
        my $hash_ref = {
            'href' => "/$sp/exportview",
            'text' => "<i>$bio_name</i>",
            'raw'  => 1
        };
        push @{ $spp_tree{$group}{'species'} }, $hash_ref;
    }
    foreach my $group (@group_order) {
        next unless @{ $spp_tree{$group}{'species'} };
        my $text = $spp_tree{$group}{'label'};
        $menu->add_entry(
            'species',
            'href'    => '/',
            'text'    => $text,
            'options' => $spp_tree{$group}{'species'},
            'code'    => 'export_' . $group,
        );
    }
}

sub exportview {
    my $self = shift;
    my $obj  = $self->{object};
    if (is_fpc_database()) {
        $self->add_format(
            'fpc_features', 'FPC Features',
            'EnsEMBL::Maize::Component::Export::fpc_features_form',
            'EnsEMBL::Maize::Component::Export::fpc_features',
            'gff' => 'GFF format',
            'tab' => 'Tab separated values',
            'csv' => 'CSV (Comma Separated values)'
        );
    } else {
        $self->add_format(
            'bac_features', 'BAC Features',
            'EnsEMBL::Maize::Component::Export::bac_features_form',
            'EnsEMBL::Maize::Component::Export::bac_features',
            'gff' => 'GFF format',
            'tab' => 'Tab separated values',
            'csv' => 'CSV (Comma Separated values)'
        );
        $self->add_format(
            'flat', 'Flat File',
            'EnsEMBL::Maize::Component::Export::flat_form',
            'EnsEMBL::Maize::Component::Export::flat',
            'embl'    => 'EMBL',
            'genbank' => 'GenBank'
        );
        $self->add_format(
            'fasta',
            'FASTA File',
            'EnsEMBL::Maize::Component::Export::fasta_form',
            'EnsEMBL::Maize::Component::Export::fasta',
            'fasta' => 'FASTA format text file'
        );
        $self->add_format(
            'pipmaker', 'PIP (%age identity plot)',
            'EnsEMBL::Maize::Component::Export::pip_form', undef,
            'pipmaker' => 'Pipmaker / zPicture format',
            'vista'    => 'Vista Format'
        );
    }

    if ($obj->seq_region_name) {
        map { $obj->param($_) || $obj->param($_, '') }
            qw/type1 type2 anchor2/;
        if ($obj->param('type2') eq 'none' || !$obj->param('anchor2')) {
            if (   $obj->param('type1') eq 'transcript'
                || $obj->param('type1') eq 'peptide')
            {
                $self->{object}
                    ->alternative_object_from_factory('Transcript');
                if ((@{ $self->{object}->__data->{'objects'} || [] })
                    && !@{ $self->{object}->__data->{'transcript'} || [] })
                {
                    $self->{object}->param('db',
                        $self->{object}->__data->{'objects'}->[0]{'db'});
                    $self->{object}->param('transcript',
                        $self->{object}->__data->{'objects'}
                            ->[0]{'transcript'});
                    $self->{object}
                        ->alternative_object_from_factory('Transcript');
                }
            } elsif ($obj->param('type1') eq 'gene') {
                $self->{object}->alternative_object_from_factory('Gene');
                if ((@{ $self->{object}->__data->{'objects'} || [] })
                    && !@{ $self->{object}->__data->{'gene'} || [] })
                {
                    $self->{object}->param('db',
                        $self->{object}->__data->{'objects'}->[0]{'db'});
                    $self->{object}->param('gene',
                        $self->{object}->__data->{'objects'}->[0]{'gene'});
                    $self->{object}->alternative_object_from_factory('Gene');
                }
            }
        }
        $self->{object}->clear_problems();
        if ($obj->param('action')) {
            my $format = $self->get_format($obj->param('format'));
            if ($format) {
                if ($obj->param('action') eq 'export') {
                    my $panel3 = $self->new_panel(
                        '',
                        'code'    => 'stage3',
                        'caption' => 'Results',
                    );
                    $panel3->add_components(
                        'results' => $format->{'superdisplay'});
                    $self->add_panel($panel3);
                    return;
                } else {
                    my $panel2 = $self->new_panel(
                        '',
                        'code' => 'stage2_form',
                        'caption' =>
                            qq(Configuring $format->{'supername'} output for $format->{'name'})
                    );
                    $self->add_form($panel2,
                        'stage2_form' => $format->{'superform'});
                    $panel2->add_components(
                        qw(select EnsEMBL::Web::Component::Export::stage2));
                    $self->add_panel($panel2);
                    return;
                }
            }
        }
    } else {
        if ($obj->param('format')) {
            ## We have an error here... so we will need to pass it through to the webform...
        }
    }
    ## Display the form...
    my $panel1 = $self->new_panel(
        '',
        'code'    => 'stage1_form',
        'caption' => qq(Select region/feature to Export)
    );
    $self->add_form($panel1,
        qw(stage1_form EnsEMBL::Maize::Component::Export::stage1_form));
    $panel1->add_components(
        qw(stage1 EnsEMBL::Maize::Component::Export::stage1));
    $self->add_panel($panel1);
}

sub cytoview {
    my $self = shift;
}

##############################################################################
## Helper functions....
##############################################################################
## add_format, get_format are helper functions for configuring ExportView #####
##############################################################################

sub add_format {
    my ($self, $code, $name, $form, $display, %options) = @_;
    unless ($self->{object}->__data->{'formats'}{$code}) {
        $self->{object}->__data->{'formats'}{$code} = {
            'name'    => $name,
            'form'    => $form,
            'display' => $display,
            'sub'     => {}
        };
        foreach (keys %options) {
            $self->{object}->__data->{'formats'}{$code}{'sub'}{$_}
                = $options{$_};
        }
    }
}

sub get_format {
    my ($self, $code) = @_;
    my $formats = $self->{object}->__data->{'formats'};
    foreach my $super (keys %$formats) {
        foreach (keys %{ $formats->{$super}{'sub'} }) {
            return {
                'super'        => $super,
                'supername'    => $formats->{$super}{'name'},
                'superform'    => $formats->{$super}{'form'},
                'superdisplay' => $formats->{$super}{'display'},
                'code'         => $_,
                'name'         => $formats->{$super}{'sub'}{$_}
                }
                if $code eq $_;
        }
    }
}
