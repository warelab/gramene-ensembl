package EnsEMBL::Web::Renderer::Chromosome::HTML;

=head1 NAME

EnsEMBL::Web::Renderer::Chromosome::HTML.pm 

=head1 SYNOPSIS

This object creates HTML to be output to the HTML output object

=head1 DESCRIPTION

    $Chromosome_renderer = $Chromosome_data->renderer;
    $Chromosome_renderer->outputGenericChromosomeTable();         
        
 This object contains wrappers for common display 'groups' and also more granular calls that
 can be reused to create different page layouts/structure

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 CONTACT

Brian Gibbins - bg2@sanger.ac.uk

=cut

use strict;
use warnings;
no warnings "uninitialized";
use CGI;
use vars qw( @ISA );
use EnsEMBL::Web::Renderer::HTML;
use EnsEMBL::Web::Renderer::Chromosome;

use Data::Dumper;

@ISA = qw( EnsEMBL::Web::Renderer::HTML EnsEMBL::Web::Renderer::Chromosome );

=head2 ExtURL

 Arg[1]      : none
 Example     : $self->ExtURL
 Description : Lazy loader for ExtURL... used to create external URLs
 Return type : A ExtURL object

=cut

sub ExtURL {
    my $self = shift;
    $self->{'_ext_url_'} ||= ExtURL->new();
}

=head2 EC_URL

 Arg[1]      : string
                 url string
 Example     : $self->EC_URL
 Description : Creates EC URLs
 Return type : string - HTML

=cut

sub EC_URL {
    my $self       = shift;
    my $string     = shift;
    my $urls       = $self->ExtURL();
    my $URL_string = $string;
    $URL_string =~ s/-/\?/g;
    return qq(<a href=")
        . ($urls->get_url('EC_PATHWAY', $URL_string))
        . qq(">EC $string</a>);
}

#-------

=head2 ensembl_help_link

 Arg[1]      : None
 Example     : $self->printHTML($self->ensembl_help_link)
 Description : Creates ensembl eh? help link
 Return type : string - HTML

=cut

sub ensembl_help_link {
    my $self = shift;
    my ($kw, $se) = @_;
    my $LINK = (defined $se ? "se=$se&" : "") . "kw=$kw";
    return
        qq(javascript:X=window.open('/$ENV{'ENSEMBL_SPECIES'}/helpview?$LINK','helpview','height=400,width=500,left=100,screenX=100,top=100,screenY=100,resizable,scrollbars=yes');X.focus();void(0));
}

#--------

=head2 print

 Arg[1]         : String
 Example     : $self->print($string)
 Description : wrapper for the printHTML call on the output object, This allows output straight to the browser
 Return type : none

=cut

sub print {
    my $self = shift;
    $self->Output->printHTML(@_);
}

=head2 two_col_table

 Arg[1]         : List
 Example     : $self->two_col_table(list .....)
 Description : wrapper for the Two_col_table call on the output object,
 Return type : none

=cut

sub two_col_table {
    my $self = shift;
    $self->Output->print_two_col_table(@_);
}

#-----------------

sub do_chromoview {
    my $self = shift;
    my ($imageHTML, $feature_count) = $self->outputChromoviewMapImage;
    $self->print(
        qq(<table border="0"><tr>
            <td valign="top">$imageHTML</td>
	        <td>&nbsp;</td>	
            <td valign="top">)
            . $self->output_chromo_forms($feature_count)
            . qq(</td><td>&nbsp;</td>
            </tr></table>)
    );
}

sub doMapview {
    my $self      = shift;
    my $imageHTML = $self->outputMapviewMapImage;

    if (my $disclaimer = $self->DataObj->species_defs->ASSEMBLY_DISCLAIMER) {
        $self->print(
            qq(
<p align="left" class="datawarning">$disclaimer</p> )
        );
    }

    $self->print(
        qq(
        <table border="0" width="100%%"><tr>         
            <td valign="top">$imageHTML</td>
            <td valign="top" width="100%%"><br />) . $self->mapviewForm . qq(
            </td>
            </tr></table>)
    );
}

=head2 outputChromoviewMapImage

 Arg[1]      : none
 Example     : $transdata->renderer->outputChromoviewMapImage
 Description : Wrapper for the table for chromoview
 Return type : none

=cut

sub outputChromoviewMapImage {
    my $self = shift;
    my $data = $self->DataObj;
    my $dc   = $data->chromoviewMapImage;
    return $self->mapImage($dc);
}

sub outputMapviewMapImage {
    my $self = shift;
    my $data = $self->DataObj;
    my $dc   = $data->mapviewMapImage;

    # check if image is cached; in this case, $dc does not contain a
    # Bio::EnsEMBL::VDrawableContainer object, but rather a string containing
    # the name of the cached file
    my $cached_file;
    unless (ref $dc eq "Bio::EnsEMBL::VDrawableContainer") {
        $cached_file = $dc;
    }

    my ($html) = $self->mapImage($dc, $cached_file, "mapview");
    return $html;
}

sub mapImage {
    my ($self, $dc, $ideo_filename, $prefix) = @_;
    $prefix ||= "chromoview";
    my $data    = $self->DataObj;
    my $species = $ENV{'ENSEMBL_SPECIES'};
    my $chr     = $data->chr_name;
    my $chr_length
        = ($data->chr_name eq 'ALL') ? $data->max_chr_length : $data->length;

    if ($data->has_a_problem) {
        $self->Output->error_page($data->problem->[0]);
        $self->Output->exit_script;
    }

    my ($feature_count, $ideo_height, $top_margin);
    if ($ideo_filename) {

        # parse ideo_height and top_margin from cached file name
        ($ideo_height, $top_margin) = (split(/\./, $ideo_filename))[ -3, -2 ];
    } else {

        # if no cached image is supplied, render the image
        $ideo_height   = $dc->{'config'}->{'_image_height'};
        $top_margin    = $dc->{'config'}->{'_top_margin'};
        $feature_count = $dc->{'_feature_count'};
        ($ideo_filename)
            = $self->render_dc($prefix, $dc,
            "$species.$chr.$ideo_height.$top_margin", 1);
    }

    # if we still have no image, something went wrong ...
    return unless $ideo_filename;

    my $HTML;
    if ($data->chr_name eq 'ALL') {
        $HTML = qq(<br/><img border="0" src="$ideo_filename"><br><br>);
    } else {
        my $cv_width = 100000;
        my $scaling = $self->DataObj->species_defs->ENSEMBL_GENOME_SIZE || 1;
        $cv_width = sprintf("%d", sprintf("%.1g", $cv_width * $scaling));
        my $script = 'contigview';

        # Maize hack for Zea_mays2 AGI FPC database
        if ($species eq 'Zea_mays') {
            $script   = 'cytoview';
            $cv_width = 1000000;
        }
        my $form = $self->Output->startForm("/$species/$script", 'GET');
        $HTML = qq(<br/>$form
            <input type="hidden" name="seq_region_name" value="$chr">
            <input type="hidden" name="seq_region_width" value="$cv_width">
            <input type="hidden" name="seq_region_left" value="1">
            <input type="hidden" name="seq_region_right" value="$chr_length">
            <input type="hidden" name="click_right" value="@{[$ideo_height+$top_margin]}">
            <input type="hidden" name="click_left" value="$top_margin">
            <input type="IMAGE" name="vclick" src="$ideo_filename" BORDER="0" alt="Click to jump to contigview at this point"><br><br>);
        $HTML .= $self->Output->endForm;
    }
    return ($HTML, $feature_count);

}

sub mapviewForm {
    my $self     = shift;
    my $chr_name = $self->DataObj->chr_name;
    my $data     = $self->DataObj;

    my $html = $self->Output->stacked_table(
        $self->change_chromosome_form,
        $self->mapview_chr_stats($data->getChrStats),
        $self->jump2anchorviewform,
        $self->syntenyForm,
        $self->diseaseview_link,
        $self->chromoview_link,
    );
    return $html;
}

=head2 output_chromo_forms

 Arg[1]		: (optional) String $option
=cut

sub output_chromo_forms {

    my $self          = shift;
    my $feature_count = shift;
    my $chr_name      = $self->DataObj->chr_name;

    my $species = $self->Input->species;
    my $script  = $self->Input->script;

    my $html = '';
    $html .= $self->chromoview_chr_stats($feature_count);
    $html .= $self->data_input_form('', 1);

    return $html;

}

sub chromoview_chr_stats {
    my $self          = shift;
    my $feature_count = shift;
    my $chr_name      = $self->DataObj->chr_name;
    my $title         = "Chromosome $chr_name";
    my $html = qq(<br /><table cellspacing="0" cellpadding="0" border="0">);
    my $stats;
    for (keys %$feature_count) {
        next unless ($feature_count->{$_});
        $html .= qq(<tr><td align="right"><b>) . ucfirst($_) . qq(</b>: </td>
                        <td>&nbsp;</td>
                        <td align="right"> $feature_count->{$_} </td>
                    </tr>);
        $stats = 1;
    }
    unless ($stats) {
        $html
            .= qq(<tr><td><p class="normal"><b>No data has been uploaded for this chromosome</b>
        <br />If you wish to display your own data, please select a chromosome before uploading.</p></td></tr>);
    }
    my $table = $self->Output->stacked_table($title => $html);
    return $table;
}

sub mapview_chr_stats {
    my $self          = shift;
    my $feature_count = shift;
    my $chr_name      = $self->DataObj->chr_name;
    my $title         = qq(Chromosome $chr_name);
    my @orderlist     = (
        'Length',
        'Gene Count',
        'Known Gene Count',
        'PseudoGene Count',
        'SNP Count'
    );

    my $html = qq(<table cellspacing="0" cellpadding="0" border="0">);

#GRAMENE - Add link to SeqTable
#  $html .= qq(
# <tr>
#  <td><a href="/Oryza_sativa/SeqTable?chr=$chr_name"><b>View Chromome $chr_name Clones</b></a></td>
#  <td colspan="2"><a href="/Oryza_sativa/SeqTable?chr=$chr_name"><img src="/gfx/buttons/go.gif" border="0" alt="lookup" title="lookup"></td></td>
# </tr> );

    #GRAMENE - All links to CMAP
    # TODO: Move hyperlink into DEFAULTS.ini
    my $cmap_link = ("/db/cmap/viewer"
            . "?changeMenu=1&compMenu=1"
            . "&aggregate=1&min_correspondences=10"
            . "&data_source=Build17&ref_map_aids=gt0205-$chr_name"
            . "&comparative_map_right=-1"
            . "&comp_map_set_right=%s");

    my %cmap_maps = (
        'Maize'   => 'cmf1104',
        'Sorghum' => 'patt2003',
        'Wheat'   => 'west'
    );
    foreach my $sp (sort keys %cmap_maps) {
        my $link = sprintf($cmap_link, $cmap_maps{$sp});

#    $html .= qq(
# <tr>
#  <td><a href="$link"><b>Jump to CMap cf. $sp</b></a></td>
#  <td colspan="2"><a href="$link"><img src="/gfx/buttons/go.gif" border="0" alt="lookup" title="lookup"></td></td>
#</tr> );
    }

    my $stats;

    my %captions = (

        # Feature Tracks
        'SWISSPROT_TREMBL_PROTEINS' => 'Rice_Protein_SpTrEMBL',
        'TIGR_GENE'                 => 'Rice_GeneModel_TIGR',
        'GENEMODEL_TIGR'            => 'Rice_GeneModel_TIGR',
        'SUBMITTERGENEANNOTATION'   => 'Rice_GeneModel_Submitted',
        'FGENESH'                   => 'Rice_GeneModel_FGENESH',

        # EST tracks
        'RICE_EST'                    => 'Rice_EST',
        'RICE_GI'                     => 'Rice_ESTCluster_TGI',
        'RICE_TUG'                    => 'Rice_ESTCluster_TUG',
        'RICE_IND_CLUSTER'            => 'RiceIndica_ESTCluster_BGI',
        'RICE_IND_EST'                => 'RiceIndica_EST_BGI',
        'RICE_JAP_CDNA_KOME'          => 'RiceJaponica_cDNA_KOME',
        'RICE_MARKER'                 => 'Rice_Marker_RFLP',
        'RICE_RFLP_MARKER'            => 'Rice_Marker_RFLP',
        'TOS17'                       => 'Rice_FST_Tos17',
        'RICE_TOS17_INSERT'           => 'Rice_FST_Tos17',
        'RICE_DS_INSERT'              => 'Rice_FST_Ds',
        'RICE_T_DNA_INSERT'           => 'Rice_FST_T-DNA',
        'RICE_TRANSPOSON_INSERT_SITE' => 'Rice_FST_IS',
        'RICE_SSR'                    => 'Rice_Marker_SSR',
        'BARLEY1_GENECHIP_EXEMPLARS'  => 'Barley_Exemplar_GeneChip',
        'BARLEY_EST'                  => 'Barley_EST',
        'BARLEY_GI'                   => 'Barley_ESTCluster_TGI',
        'BARLEY_TUG'                  => 'Barley_ESTCluster_TUG',
        'MAIZE_CDS'                   => 'Maize_CDS',
        'MAIZE_EST'                   => 'Maize_EST',
        'MAIZE_GI'                    => 'Maize_ESTCluster_TGI',
        'MAIZE_TUG'                   => 'Maize_ESTCluster_TUG',
        'MAIZE_CORNSENSUS'            => 'Maize_ESTCluster_MMPcornsensus',
        'MAIZE_MU_INSERT'             => 'Maize_FST_Mu',
        'MILLET_EST'                  => 'Millet_EST',
        'RYEGRASS_EST'                => 'Ryegrass_EST_Vialactia',
        'RYEGRASS_CLUSTER'            => 'Ryegrass_ESTCluster_Vialactia',
        'SORGHUM_EST'                 => 'Sorghum_EST',
        'SORGHUM_GI'                  => 'Sorghum_ESTCluster_TGI',
        'SORGHUM_TUG'                 => 'Sorghum_ESTCluster_TUG',
        'SORGHUM_CLUSTER_PRATT'       => 'Sorghum_ESTCluster_Pratt',
        'SUGARCANE_EST'               => 'Sugarcane_EST',
        'WHEAT_EST'                   => 'Wheat_EST',
        'WHEAT_GI'                    => 'Wheat_ESTCluster_TGI',
        'WHEAT_TUG'                   => 'Wheat_ESTCluster_TUG',

        #GSS Tracks
        'RICE_JAPONICA_BACEND'           => 'RiceJaponica_BACend_IRGSP',
        'OR_BBA'                         => 'RiceNivara_BACend_OMAP',
        'OR_CBA'                         => 'RiceRufipogon_BACend_OMAP',
        'RICE_BRACHYANTHA_BACEND'        => 'RiceBrachyantha_BACend_OMAP',
        'RICE_ALTA_BACEND'               => 'RiceAlta_BACend_OMAP',
        'RICE_AUSTRALIENSIS_BACEND'      => 'RiceAustraliensis_BACend_OMAP',
        'RICE_GLABERRIMA_BACEND'         => 'RiceGlaberrima_BACend_OMAP',
        'RICE_NIVARA_BACEND'             => 'RiceNivara_BACend_OMAP',
        'RICE_PUNCTATA_BACEND'           => 'RicePunctata_BACend_OMAP',
        'RICE_RUFIPOGON_BACEND'          => 'RiceRufipogon_BACend_OMAP',
        'MAIZE_HI_COT_BENNETZEN'         => 'Maize_HiCot_Bennetzen',
        'MAIZE_METH_FILT_TIGR'           => 'Maize_MethylFilter_Orion',
        'MAIZE_METH_FILT_CSHL/MCCOMBIE'  => 'Maize_MethylFilter_CSHL',
        'MAIZE_METH_FILT_CSHL_MCCOMBIE'  => 'Maize_MethylFilter_CSHL',
        'MAIZE_HI_COT_TIGR'              => 'Maize_HiCotCluster_TIGR',
        'MAIZE_METH_FILT_HI_COT_CLUSTER' =>
            'Maize_HiCotMethylFilterCluster_TIGR',
        'RYEGRASS_SEQUENCE'         => 'Ryegrass_MethylFilter_Orion',
        'RYEGRASS_ASSEMBLY'         => 'Ryegrass_MethylFilterCluster_Orion',
        'SORGHUM_RYEGRASS_ASSEMBLY' => 'Ryegrass_MethylFilterCluster_Orion',
        'SORGHUM_RYEGRASS_SEQ'      => 'Ryegrass_MethylFilter_Orion',
        'SORGHUM_GSS-READ_KLEIN'    => 'Sorghum_GSS_Klein',
        'SORGHUM_ORION'             => 'Sorghum_MethylFilter_Orion',

        # BAC tracks
        'BAC_MAP'             => 'BACs',
        'ACCESSIONED_BAC_MAP' => 'Acc_BACs',
        'FPC_CONTIG'          => 'FPContigs',
    );

    my %stats;
    delete($feature_count->{'Top Level'});
    foreach my $name (keys %{ $feature_count || {} }) {
        my $value = $feature_count->{$name};
        $name =~ s/\s+/_/g;
        my $key = $name;
        $key =~ s/_Count//i;
        $name = $captions{ uc($key) } || $key;
        $stats{$name} = $value;
    }
    my @sorted = (
        map {
            $_->[1] =~ s/ /_/g;
            $_->[1]
            }
            sort {
            ($a->[0] <=> $b->[0])
                || ($a->[1] cmp $b->[1])
            }
            map {
            my $i = 4;
            $_ =~ s/_/ /g;
            if    (/sptrembl/i)  { $i = 2 }
            elsif (/genemodel/i) { $i = 2 }
            elsif (/^length/i)   { $i = 1 }
            elsif (/^rice/i)     { $i = 3 }
            [ $i, $_ ]
            } keys %stats
    );

    foreach my $name (@sorted) {

        #for (@orderlist){
        next if ($name eq 'TE');
        my $bps_label;
        if ($name eq "Length") {
            $bps_label = 'bps';
            $feature_count->{$name}
                = $self->DataObj->commify($feature_count->{$name});
        }
        $html .= qq(<tr><td align="left"><b><small>)
            . ucfirst($name)
            . qq(</small></b>: </td>
                        <td>&nbsp;</td>
                        <td align="right"><small> $stats{$name} <small></td>
                        <td> &nbsp;$bps_label</td>
                    </tr>);
        $stats = 1;
    }
    unless ($stats) {
        $html .= qq(<tr><td><b>Could not load chromosome stats</b><td></tr>);
    }
    $html .= qq(  </table>  );
    return ($title, $html);
}

sub chromoview_chr_stats_annot {
    my $self          = shift;
    my $feature_count = $self->DataObj->getChrStats;
    my $chr           = $self->DataObj->chr_name;
    my $species       = $ENV{'ENSEMBL_SPECIES'};
    my $vega_helplink = $self->DataObj->species_defs->SPECIES_SHORT_NAME
        . "_gene_classification";
    my $title = $feature_count->{'Name'} || qq(Chromosome $chr);
    my $ack = {
        "Homo_sapiens" => {
            6       => "Havana",
            "6-COX" => "Havana",
            7       => "Wash. U.",
            9       => "Havana",
            10      => "Havana",
            13      => "Havana",
            14      => "Genoscope",
            20      => "Havana",
            22      => "Collins et al.",
        },
        "Mus_musculus" => { 13 => "Havana", },
    };

    if ($ack->{$species}->{$chr}) {
        $title
            .= qq(&nbsp;&nbsp;<a href="acknowledgements.html?#chr$chr"><span class = "small">$ack->{$species}->{$chr} Group</span></a>);
    }
    my $html = qq(<br /><table cellspacing="0" cellpadding="0" border="0">);

    $feature_count->{'01:Length'} = $feature_count->{'Length'};
    delete $feature_count->{'Length'};
    delete $feature_count->{'Top Level'};
    delete $feature_count->{'Name'};
    my $stats;
    foreach (sort keys %$feature_count) {
        next unless $feature_count->{$_};
        (my $label = $_) =~ s/\d+://;
        $label = ucfirst($label);

        my $bps_label
            = ($label eq "Length" or $label eq "Annotated sequence length")
            ? 'bps'
            : '';
        my $value = $self->DataObj->commify($feature_count->{$_});

        if (/^snp/) {
            $label =~ s/^snp/SNP/ig;
        }
        $html .= qq(<tr><td align="right"><b>) . $label . qq(</b>:</td>
							 <td>&nbsp;</td>
							 <td align="right">$value</td>
							 <td>&nbsp;$bps_label</td>
							 </tr>);
        $stats = 1;
    }
    $html
        .= qq(<tr><td><b>No data has been uploaded for this chromosome</b><td></tr>)
        unless $stats;
    $html .= qq(
		<tr>
		<td align="center" colspan="4">
		<br /><a href="javascript:X=hw\('$species', '$vega_helplink', ''\)">Definitions of indices</a><br />
		</td>
		</tr>
		<tr>);
    if (@{ $self->DataObj->chromosome->get_all_MiscFeatures('NoAnnotation') }) {
        $html .= qq(<td align="center" colspan="4">
		    <br />Shaded regions have not been manually annotated.<br />
		    </td>);
    }
    $html .= qq (</tr></table>);

    return ($title, $html);
}

sub haplotypes_info {
    my $self     = shift;
    my $chr_name = $self->DataObj->chr_name;
    my $species  = $ENV{'ENSEMBL_SPECIES'};
    my $script   = $ENV{'ENSEMBL_SCRIPT'};

    my $has_haplotypes;
    my $haplo_chr;
    my $haplotypes = $self->DataObj->species_defs
        ->ENSEMBL_HAPLOTYPES;    #returns an array ref of haplotypes?
    unless ($haplotypes)
    { # checks that the species has haplotypes, if not returns an anon array containing just the chromosome ID
        $has_haplotypes = 0;
        $haplo_chr      = [$chr_name];
    }

    foreach my $haploset (@{$haplotypes}) {
        my @haplo_count = split(",", $haploset);
        foreach (@haplo_count) {
            if ($_ eq $chr_name) {
                $has_haplotypes = 1;
                $haplo_chr      = \@haplo_count;
            }
        }
    }

    unless ($has_haplotypes) {return}

    my $num   = scalar(@{$haplo_chr});
    my $first = shift @{$haplo_chr};
    my $title = qq(Haplotypes);
    my $html  = qq(<br /><table cellspacing="0" cellpadding="0" border="0">

          <td>
            <br />);
    if ($chr_name eq $first) {
        $html
            .= qq(There are $num annotated haplotypes for chromosome $chr_name:);
    }
    $html .= qq(<table cellpadding="0" cellspacing="0" border="0" width="100%">
             <tr>
              <td width=1%><img src="/gfx/blank.gif" height="10" width="15" alt=""></td>
              <td><img src="/gfx/blank.gif" height="10" width="200" alt=""></td>
             </tr>

             <tr class="background1">
              <td valign="top"><img width="8" height="8" src="/gfx/bullet.blue.gif" alt="o"></td>
              <td valign="top">
                <a href=/$species/karyotypes/chr$first.html>Overview</a> - more information on the haplotypes
              </td> 
             </tr>
                
             <!--
             <tr class="background1">
              <td valign="top"><img width="8" height="8" src="/gfx/bullet.blue.gif" alt="o"></td>
              <td valign="top">
                <a href=/$species/haplomapview?chr=$first>Compare haplotypes</a> - density plots for all haplotypes
              </td> 
             </tr>
             -->
                
             <tr class="background1">
              <td valign="top"><img width="8" height="8" src="/gfx/bullet.blue.gif" alt="o"></td>
              <td valign="top">
                Individual haplotypes:
                <ul>
    );
    foreach (@{$haplo_chr}) {
        $html .= "<li><a href='/$species/$script?chr=$_'>$_</a>\n";
    }
    $html .= qq(
                  <li><a href='/$species/$script?chr=$first'>$first (reference sequence)</a>
                </ul>
              </td> 
             </tr>
            </table>
            <br />
	  </td>
          <td><img src="/gfx/blank.gif" width="15" alt=""></td>
        </tr> );

    $html .= qq(</table>);
    return ($title, $html);
}

sub acknowledgements_annot {
    my $self = shift;
    return if ($ENV{'ENSEMBL_SPECIES'} eq "Danio_rerio");
    my $chr_name = $self->DataObj->chr_name;
    my $title    = qq(Acknowledgements);
    my $html     = qq(<br /><table cellspacing="0" cellpadding="0" border="0">
        <tr class="background1">
          <td><img src="/gfx/blank.gif" width="15" alt=""></td>
          <td align="center" colspan="4">
          <br /><a href="acknowledgements.html?#chr$chr_name">Chromosome $chr_name Project Acknowledgments</a><br /><br />
      </td>
          <td><img src="/gfx/blank.gif" width="15" alt=""></td>
          </tr>
      </table>
    );
    return ($title, $html);
}

=head2 change_chromosome_form

  Arg[1]      : (optional) String $extra
                additional chromosomes to add the the pulldown list
  Arg[2]      : (optional) Boolean $http_get
                flag to force HTTP GET (default: POST, multipart/form-data)
  Example     : $renderer->change_chromosome_form(undef, 1)
  Description : prints a form to change chromosomes
  Return type : title/html pair for use in stacked table output
  Exceptions  : none
  Caller      : general

=cut

sub change_chromosome_form {
    my ($self, $extra, $http_get) = @_;
    my @form_args = ("/$ENV{'ENSEMBL_SPECIES'}/$ENV{'ENSEMBL_SCRIPT'}");

    if ($http_get) {
        push @form_args, 'get';
    } else {
        push @form_args, ('post', 'multipart/form-data');
    }
    my $html = $self->Output->startForm(@form_args) . "\n";
    my ($title, $widgets) = $self->change_chromosome_widgets($extra);
    $html .= $widgets;
    $html .= $self->Output->endForm();
    return ($title, $html);
}

sub change_chromosome_widgets {

    my ($self, $extra) = @_;

    my $chr_name = $self->DataObj->chr_name || 'ALL';
    my $title = 'Change Chromosome';
    my $html_options;
    my @chromosomes = @{ $self->DataObj->all_chromosomes };
    push @chromosomes, $extra if ($extra);
    foreach (@chromosomes) {
        $html_options .= qq(<option value="$_")
            . ($_ eq $chr_name ? ' selected' : '')
            . qq(>$_</option>\n);
    }
    my $html = qq(<br />            
  <table cellspacing="0" cellpadding="0" border="0">
    <tr>
      <td>Chromosome: </td>
      <td> <select name="chr">$html_options</select></td>
      <td> &nbsp;<INPUT TYPE="image" VALUE="go" src="/gfx/buttons/go.gif" border="0" valign="middle"></td>
    </tr>
  </table>);

    return ($title, $html);
}

sub jumptoChromosome {
    my $self     = shift;
    my $extra    = shift;
    my $chr_name = $self->DataObj->chr_name || 'ALL';
    my $title    = qq(Change Chromosome);
    my $html_options;
    my @chromosomes = @{ $self->DataObj->all_chromosomes };
    push @chromosomes, $extra if ($extra);
    my $OTHER = $self->param('otherspecies');

    foreach (@chromosomes) {
        $html_options .= qq(<option value="$_")
            . ($_ eq $chr_name ? ' selected' : '')
            . qq(>$_</option>\n);
    }

    my ($html, $form);
    unless ($self->Output->in_form) {
        $html = $self->Output->startForm(
            "/$ENV{'ENSEMBL_SPECIES'}/$ENV{'ENSEMBL_SCRIPT'}",
            "post", "multipart/form-data");
        $form = 1;
    }
    $html .= qq(<br />            
  <table cellspacing="0" cellpadding="0" border="0">
    <tr>
      <td><b>Jump to Chromosome: </b></td>
      <td><select name="chr">$html_options</select></td>
      <td>&nbsp;<INPUT TYPE="image" VALUE="go" src="/gfx/buttons/go.gif" border="0" valign="middle"></td>
    </tr>
    <tr>
      <td colspan="3"><a href="/$ENV{'ENSEMBL_SPECIES'}/mapview?chr=$chr_name">Jump to MapView</a> for chromosome statistics.</td>
    </tr>
  </table>
  <input type="hidden" name="otherspecies" value="$OTHER">
  );
    $html .= $self->Output->endForm if ($form);
    return ($title, $html);
}

sub upload_data_form {
    my $self  = shift;
    my $title = qq(Upload Data Set:);

    my $html
        = qq(Accepted <a href="/$ENV{'ENSEMBL_SPECIES'}/helpview?se=1&kw=chromoview">file formats</a><br /><br />
    <table cellpadding="0" cellspacing="0" border="0">
            <tr>
                <td>Paste file: </td><td><textarea name="paste_file" rows="10" cols="40">)
        . $self->param('paste_file')
        . qq(</textarea><br /><br /></td>
            </tr><tr valign="top">
                <td>Upload file: </td><td><input type="file" name="upload_file" size="40" value=")
        . $self->param('upload_file')
        . qq("><br /><br /></td>
            </tr><tr valign="top">
                <td>File URL: </td><td><input type="text" name="url_file" size="40" value=")
        . $self->param('url_file')
        . qq("><br /></td>
               </tr><tr valign="top">
                <td>&nbsp;</td><td>
                <table><tr>
                    <td><input type="checkbox" name="maxmin"  value="checked")
        . $self->param('maxmin')
        . qq(>&nbsp;Show Max/Min lines<br /><br /></td>
                    <td><img src="/gfx/blank.gif" width="50" height="1" alt=""></td>
                    <td><input type="checkbox" name="snpfrq"  value="checked" )
        . $self->param('snpfrq')
        . qq(>&nbsp;Show SNP frequency<br /><br /></td>
                </tr>
                <tr>    
                    <td><input type="checkbox" name="genefrq"  value="checked" )
        . $self->param('genefrq')
        . qq(>&nbsp;Show gene frequency<br /><br /></td>
                    <td><img src="/gfx/blank.gif" width="50" height="1" alt=""></td>
                    <td><input type="checkbox" name="pctGC"  value="checked" )
        . $self->param('pctGC')
        . qq(>&nbsp;Show GC content frequency<br /><br /></td>
                </tr></table>
                
                
                </td>
            </tr>            
            </table>
            <INPUT TYPE="image" VALUE="go" src="/gfx/buttons/go.gif" border="0" valign="middle">);
    return ($title, $html);
}

sub jump2anchorviewform {
    my $self        = shift;
    my $data        = $self->DataObj;
    my $chr         = $data->chr_name;
    my $ensembl_map = $data->species_defs->DB_FEATURES->{'MARKERS'};
    my $marker_from = $self->param('marker_start');
    my $marker_to   = $self->param('marker_end');
    my $title       = "Jump to Contigview";
    my $html
        = qq(<br>Click anywhere on the chromosome ideogram or one of the feature distribution plots to jump to a contig-level view of features at that point. 
);

    if (0) {    #if($ensembl_map){

        $html
            .= "Alternatively, you can jump to contigview between any two markers on this chromosome:";
        $html .= $self->Output->startForm("/$ENV{'ENSEMBL_SPECIES'}/anchorview",
            "GET", "multipart/form-data");
        $html .= qq( 
            <table align="center" cellpadding="0" cellspacing="5" border="0">
            <tr>
                <td><input type="hidden" name="chr" value="$chr">
                    <input type="hidden" name="type1" value="marker">
                    <input type="hidden" name="type2" value="marker">Between: </td>
                <td><input type="text" name="marker_start" value="$marker_from"></td>
                <td rowspan="2">&nbsp;</td>
                <td rowspan="2" align="center" valign="middle">
                    <INPUT TYPE="image" VALUE="go" src="/gfx/buttons/go.gif" border="0"></td>
            </tr><tr>
                <td>and: </td>
                <td><input type="text" name="marker_end" value="$marker_to"></td>
            </tr> 
            </table>
         );
        $html .= $self->Output->endForm;
    }

    # Maize hack
    if (my $slice = $self->DataObj->EnsemblObj) {
        my $table_t = qq(
<table align="center" cellpadding="0" cellspacing="5" border="0">
 <tr>
  <td>Jump to %s:</td>
  <td>
   <select name="contig">%s
   </select>
  </td>
  <td>&nbsp;</td>
  <td align="center" valign="middle">
   <INPUT TYPE="image" VALUE="go" src="/gfx/buttons/go.gif" border="0">
  </td>
 </tr>
<tr>
</table>);
        my $option_t = qq(
<option value="%s">%s</option>);

        my @seq_levels;
        my $seq_level_name;
        foreach (@{ $slice->project('seqlevel') }) {
            push @seq_levels, $_->to_Slice->seq_region_name;
            $seq_level_name ||= $_->to_Slice->coord_system_name;
        }
        my $script = $seq_level_name eq 'fpc' ? 'cytoview' : 'contigview';
        $html .= $self->Output->startForm("/$ENV{'ENSEMBL_SPECIES'}/$script",
            "GET", "multipart/form-data");
        $html .= sprintf($table_t,
            $seq_level_name,
            join('', map { sprintf($option_t, $_, $_) } @seq_levels));
        $html .= $self->Output->endForm;
    }

    $html .= qq( 
        <p align="center"><a href="/$ENV{'ENSEMBL_SPECIES'}/anchorview?chr=$chr">Display contig-level view between any two features.</a>
        <br />&nbsp;
    );

    return ($title, $html);
}

sub syntenyForm {
    my $self    = shift;
    my $data    = $self->DataObj;
    my $chr     = $data->chr_name;
    my $label   = 'Synteny';
    my %hash    = $data->species_defs->multi('SYNTENY');
    my @SPECIES = grep {
        $data->species_defs->other_species($_, 'ASSEMBLY_STATUS') eq 'FULL'
    } keys(%hash);
    my ($species_options, $common_name);

    if (@SPECIES) {
        my $html
            = $self->Output->startForm("/$ENV{'ENSEMBL_SPECIES'}/syntenyview",
            "GET");
        $common_name = $data->species_defs->SPECIES_COMMON_NAME;
        $species_options = join '', map {
            qq(<option value="$_">)
                . (
                $data->species_defs->other_species($_, "SPECIES_COMMON_NAME"))
                . qq(</option>)
        } @SPECIES;

        $html .= qq(<br />
        <table cellspacing="0" cellpadding="0" border="0">
          <tr>
            <td>View $common_name Chr $chr vs &nbsp;</td>
            <td><select name="otherspecies">$species_options</select></td>
            <input type="hidden" name="chr" value="$chr">
            <td>&nbsp; <INPUT TYPE="image" VALUE="go" src="/gfx/buttons/go.gif" border="0" valign="middle"></td>
          </tr>
            </table>
        );
        $html .= $self->Output->endForm;
        return ($label, $html);
    } else {
        return;
    }
}

sub diseaseview_link {
    my $self = shift;
    my $data = $self->DataObj;
    return unless ($data->species_defs->databases->{'ENSEMBL_DISEASE'});

    my $chr   = $data->chr_name;
    my $title = 'OMIM Diseases';
    my $html
        = qq(<br><a href="/$ENV{'ENSEMBL_SPECIES'}/diseaseview?chr=$chr">Browse OMIM Diseases</a> on this chromosome.<br>);
    return ($title, $html);
}

sub chromoview_link {
    my $self = shift;
    my $data = $self->DataObj;

    my $chr   = $data->chr_name;
    my $title = 'Map your data';
    my $html
        .= qq(<br><a href="/$ENV{'ENSEMBL_SPECIES'}/karyoview">Map your own data on the karyotype</a> using KaryoView.<br>);
    return ($title, $html);
}

sub doSyntenyview {
    my $self      = shift;
    my $imageHTML = $self->renderSyntenyviewImage;

    $self->print(
        qq(
        <table border="0" width="80%%">
           <tr>
            <td valign="top"><br />$imageHTML</td>
            <td valign="top" width="100%%"><br />)
            . $self->syntenyviewTable
            . qq(</td>
            </tr></table>)
    );
}

sub renderSyntenyviewImage {
    my $self = shift;
    my $data = $self->DataObj;
    my $dc   = $data->syntenyviewImage;
    if ($data->has_a_problem) {
        $self->Output->error_page($data->problem->[0]);
        exit;
    }

    #my $TMP_URL = $SiteDefs::ENSEMBL_TMP_URL_IMG;
    my ($filename, $map) = $self->render_image_imagemap("syntenyview", $dc);
    return $filename
        ? (
        qq(<IMG SRC="$filename" BORDER="0" usemap="#karyo"><map name="karyo">\n$map</map>)
        )
        : '';
}

sub syntenyviewTable {
    my $self = shift;
    my $data = $self->DataObj;
    my $html = $self->Output->stacked_table(

        #$self->syntenyMatches,
        #$self->navigateSynteny,
        $self->jumptoChromosome,
        $self->syntenyForm,
    );

    return $html;
}

sub syntenyMatches {
    my $self  = shift;
    my $data  = $self->DataObj;
    my $title = ' Homology Matches';
    my $OTHER = $self->param('otherspecies')
        || ($ENV{'ENSEMBL_SPECIES'} eq 'Homo_sapiens'
        ? 'Mus_musculus'
        : 'Homo_sapiens');
    my $urls         = $self->ExtURL;
    my $table_header = [
        {   'key'   => '_species_gene',
            'title' => "<i>" . $ENV{'ENSEMBL_SPECIES'} . "</i> Genes",
            'width' => '40%',
            'align' => 'left'
        },
        {   'key'    => '_join',
            'title'  => ' ',
            'width'  => '10%',
            'align'  => 'center',
            'valign' => 'middle'
        },
        {   'key'   => '_synteny_gene',
            'title' => "<i>$OTHER</i> Homologues",
            'width' => '40%',
            'align' => 'left'
        },
    ];
    my $tabledata = $data->getSyntenyMatches;
    my $html      = $self->Output->spreadsheet_table($table_header, $tabledata,
        { 'rows' => [qw(background1 background3)] });
    return ($title, '<br />' . $html);
}

sub navigateSynteny {
    my $self  = shift;
    my $data  = $self->DataObj;
    my $chr   = $data->chr_name;
    my $OTHER = $self->param('otherspecies')
        || ($ENV{'ENSEMBL_SPECIES'} eq 'Homo_sapiens'
        ? 'Mus_musculus'
        : 'Homo_sapiens');
    my $title       = ' Navigate Homology';
    my @localgenes  = @{ $data->get_synteny_local_genes };
    my $first_start = @localgenes ? $localgenes[0]->start : 0;
    my $last_end    = @localgenes ? $localgenes[-1]->end : 0;

    my $html = qq(<table><tr>
                    <td>
<a href= "/$ENV{'ENSEMBL_SPECIES'}/syntenyview?otherspecies=$OTHER&chr=$chr&loc=$first_start&pre=1">Upstream</a> \(&lt;)
        . $data->bp_to_nearest_unit($first_start)
        . qq(\)&nbsp;&nbsp;&nbsp;</td>
                    <td>
<a href = "/$ENV{'ENSEMBL_SPECIES'}/syntenyview?otherspecies=$OTHER&chr=$chr&loc=$last_end" >Downstream</a> \(&gt;)
        . $data->bp_to_nearest_unit($last_end)
        . qq(\)&nbsp;&nbsp;&nbsp;</td>
                  </tr></table>);
    return ($title, $html);
}

#------------------

sub do_karyoview {
    my $self = shift;
    $self->print("<br />");
    (          $self->param('paste_file')
            || $self->param('upload_file')
            || $self->param('url_file'))
        ? $self->output_karyo_image()
        : $self->output_karyo_intro();
    $self->output_karyo_form();
}

sub output_karyo_intro {
    my $self  = shift;
    my $title = "KaryoView";
    my $html
        = qq(<p>This page enables you to display your own data on a customisable karyotype image.</p>
    <p>Use the form below to enter or upload your data and to configure the display.  Hit the "Display Karyotype" button to view the results.</p>
    <p>You can continue to modify the configuration or data after displaying the karyotype.</p>
    <p>If you have any <a href="/helpdesk/index.html">suggestions or feedback</a> about KaryoView, we'd like to hear from you.</p>)
        ;    #'

    $self->Output->print_stacked_table_formatted($title, $html,
        { 'width' => '55%' });
}

sub output_karyo_image {
    my $self     = shift;
    my $karyo_dc = $self->DataObj->render_karyotype_image;
    my ($image, $image_map)
        = $self->render_image_imagemap('karyoview', $karyo_dc);
    $self->print(
        qq(
        <table border="0" align="center"><tr>
            <td valign="top"><img src="$image" border="0" usemap="ideo"><br><br>\n<map name="ideo">$image_map</map></td>
            </tr></table>)
    );
}

sub output_karyo_form {

    my $self = shift;

    my $species = $self->Input->species;
    my $script  = $self->Input->script;

    my $dataHTML = $self->data_input_form('85%', 1);

    $self->print($dataHTML);

#$self->Output->print_stacked_table_formatted('Input Data' => $dataHTML, {'width'=>'85%'});

}

=head2 data_input_form

 Arg[1]		: String $width
 		 Width of table
 Arg[2]		: (optional) Boolean
 		  Show appropriate display options for this view

=cut

sub data_input_form {

    my $self         = shift;
    my $width        = shift;
    my $show_options = shift;
    my $script       = $ENV{'ENSEMBL_SCRIPT'};

    my $data   = CGI::escapeHTML($self->param('paste_file'));
    my $upload = $self->param('upload_file');
    my $formats
        = qq(Accepted <a href="javascript:window.open('/$ENV{'ENSEMBL_SPECIES'}/helpview?se=1&kw=chromoview#FileFormats', 'helpview', 'width=400,height=500,resizable,scrollbars'); void(0);">file formats</a>);
    my $blurb;

    if ($script eq 'karyoview') {
        $blurb = qq(
    <p>Type or paste into the box to the right a list of features you want to
    display on a karyotype, with each feature on a separate line.<br>
    Feature data must be separated by spaces or tabs, and must contain, in
    order:
        <ul>
        <li><b>chromosome</b>
        <li><b>start base</b> (in chromosome coordinates)
        <li><b>end base</b> (in chromosome coordinates)
        <li><b>feature name or ID</b>
        </ul>    
    e.g. <b>3 12000000 15000000 my_feature</b></p>
    );
    }

    my ($html, @sections);

    $html .= $self->Output->startForm("$script", 'post', 'multipart/form-data');

    if ($script eq 'chromoview') {
        my ($chr_title, $chr_table) = $self->change_chromosome_widgets('ALL');
        push @sections, ($chr_title => $chr_table);
    }

    my $upload_title = 'Upload Data Set';
    my $chr_name     = $self->DataObj->chr_name;

    my $upload_table = qq(
      $formats<br /><br />
      <table class="hidden"><tr>
      <td>Paste file:
      $blurb
       <td align="left"><br>
    <textarea name="paste_file" rows="12" cols="40" >$data</textarea>
       </td>
       </tr>
       <tr>
       <td><img src="/gfx/blank.gif" width="10" height="10"></td>
     </tr>
       <tr>
       <td align="left">Upload file:</td>
       <td><input type="file" name="upload_file"></td>
       </td>
       </tr>
       <tr><td>File URL:</td>
       <td><input type="text" name="url_file" size="40" value=")
        . $self->param('url_file')
        . qq("><br /></td></tr>
       <tr>
       <td colspan="2" align="center">
       <input type="hidden" name="chr" value="$chr_name">
       <input type="image" src="/gfx/buttons/go.gif" name="submit" value="Display Karyotype"></td>
       </tr>
       </table>);

    push @sections, ($upload_title => $upload_table);

    if ($show_options) {
        my $options_title = 'Configure Display Options';
        my $method_name   = "_image_options_$script";
        my $options_table = $self->$method_name();
        push @sections, ($options_title => $options_table);
    }

    $html .= $self->Output->stacked_table_formatted(@sections,
        { 'width' => $width });

    $html .= $self->Output->endForm();

    return $html;
}

#-----------------

=head2 image_chromosome_number

 Example     : $self->image_chromosome_number()
 Description : Generates form widgets for user selection of 
                which chromosome(s) to display
 Return type : String (HTML) 

=cut

sub image_chromosome_number {

    my $self  = shift;
    my $extra = shift;

    my $display = $self->param('display');
    my $default = 'ALL';

    my $title = 'Chromosome(s) to display';
    my $html_options;
    my @chromosomes = @{ $self->DataObj->all_chromosomes };
    push @chromosomes, $extra if ($extra);
    foreach (@chromosomes) {
        $html_options .= qq(<option value="$_");
        if ($_ eq $self->param('chr')
            || ($_ eq $default && !$self->param('chr')))
        {
            $html_options .= ' selected="selected"';
        }
        $html_options .= qq(>$_</option>\n);
    }
    my $html = qq(<br />            
  <table cellspacing="5" cellpadding="0" border="0">
    <tr>
      <td>Chromosome:</td>
      <td> <select name="chr">$html_options</select></td>
    </tr>
  </table>);

    return ($title, $html);
}

#-----------------

=head2 image_options_location

 Example     : $self->image_options_location()
 Description : Generates form widgets for user configuration of 
                feature location maps
 Return type : String (HTML) 

=cut

sub image_options_location {

    my $self    = shift;
    my $extra   = shift;
    my @colours = $self->get_fg_cols();

    my $zmenu = lc($self->param('zmenu')) || 'on';
    my $col   = lc($self->param('col'))   || 'red';
    my $style = lc($self->param('style')) || 'rharrow';

    # Available drawing styles
    my %style = (
        'box'           => "Filled box",
        'filledwidebox' => "Filled wide box",
        'widebox'       => "Outline wide box",
        'outbox'        => "Oversize outline box",
        'wideline'      => "Line",
        'lharrow'       => "Arrow left side",
        'rharrow'       => "Arrow right side",
        'bowtie'        => "Arrows both sides",
        'text'          => "Text label (+ wide box)",
    );

    my @stylelist
        = qw(box filledwidebox widebox outbox wideline lharrow rharrow bowtie);

    push @stylelist, 'text' if $extra eq 'text';

    my $stylehtml;
    foreach (@stylelist) {
        $stylehtml .= qq(<option value="$_")
            . ($_ eq $style ? ' selected' : '') . qq(>)
            . $style{$_}
            . qq(</option>\n);
    }

    my $colourhtml;
    foreach (@colours) {
        $colourhtml .= qq(<option value="$_")
            . ($_ eq $col ? ' selected' : '')
            . qq(>$_</option>\n);
    }

    my $zmenuhtml;
    foreach ("on", "off") {
        $zmenuhtml .= qq(<option value="$_")
            . ($_ eq $zmenu ? ' selected' : '')
            . qq(>$_</option>\n);
    }

    my $html = qq(
      <table>
    <tr class="background1">
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="10"></td>
     </tr>
    <tr>
      <td align="left">Display data points in the following style:</td>
      <td align="left">&nbsp;&nbsp;</td>
      <td><select name="style">$stylehtml</select></td>
    </tr>
    <tr class="background1">
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
    <tr>
      <td align="left">Display data points in the following colour:</td>
      <td align="left">&nbsp;&nbsp;</td>
      <td><select name="col">$colourhtml</select></td>
    </tr>
    <tr class="background1">
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
    <tr>
      <td align="left">Show menus on mouseover:</td>
      <td align="left">&nbsp;&nbsp;</td>
      <td><select name="zmenu">$zmenuhtml</select></td>
    </tr>
    <tr>
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
    </table>);

    # Add chromosome layout widgets
    $html .= _image_options_karyotype($self);

    return $html;

}

#-----------------

=head2 image_options_density

 Example     : $self->image_options_density()
 Description : Generates form widgets for user configuration of 
                feature density maps
 Return type : String (HTML) 

=cut

sub image_options_density {

    my $self    = shift;
    my @colours = $self->get_fg_cols();

    my $html = '<table>';

    my %checkboxes = (
        'maxmin'    => 'Show Max/Min lines',
        'Vgenes'    => 'Show gene frequency',
        'Vpercents' => 'Show GC content frequency'
    );

    # check if this species has SNPs and add option if it has
    my $species_defs = $self->DataObj->species_defs;
    my $has_snps     = $species_defs->get_table_size(
        { -db => 'ENSEMBL_VARIATION', -table => 'source' });
    if ($has_snps) {
        $checkboxes{'Vsnps'} = 'Show SNP frequency';
    }

    my $count;
    foreach my $box (keys %checkboxes) {
        my $newrow = $count % 2;

        # if even number, start a new row
        if (!$newrow) {
            $html .= '<tr>';
        }
        my $tick = $self->param($box) eq 'on' ? ' checked="checked"' : '';
        $html
            .= qq(<td><input type="checkbox" name="$box" value="on"$tick>&nbsp;$checkboxes{$box}<br /><br /></td>\n);

        # if even number, insert a spacer cell
        if (!$newrow) {
            $html
                .= qq(<td><img src="/gfx/blank.gif" width="50" height="1" alt=""></td>\n);
        }

        # otherwise end the row
        else {
            $html .= "</tr>\n";
        }
        $count++;
    }
    $html
        .= qq(<tr><td colspan="4" align="center">Track colour:&nbsp;<select name="col">);
    my $col = lc($self->param('col')) || 'purple';
    foreach (@colours) {
        $html .= qq(<option value="$_")
            . ($_ eq $col ? ' selected' : '')
            . qq(>$_</option>\n);
    }
    $html .= qq(</select></td></tr></table>);

    # Add chromosome layout widgets
    $html .= $self->_image_options_karyotype();

    return $html;
}

#-----------------

=head2 _image_options_karyotype

 Example     : $self->image_options_karyotype()
 Description : Private function to add form widgets for chromosome layout
 Return type : String (HTML) 

=cut

sub _image_options_karyotype {

    my $self = shift;

    my $chr_length = CGI::escapeHTML($self->param('chr_length')) || 200;
    my $v_padding  = CGI::escapeHTML($self->param('v_padding'))  || 50;
    $v_padding = 20 unless $v_padding >= 20;

    my $h_padding = CGI::escapeHTML($self->param('h_padding')) || 4;
    $h_padding = 1 unless $h_padding >= 1;

    my $h_spacing = CGI::escapeHTML($self->param('h_spacing')) || 6;
    $h_spacing = 1 unless $h_spacing >= 1;

    my $rows = $self->param('rows') > 0 ? $self->param('rows') : 1;
    my $rowhtml;
    foreach (1 .. 4) {
        $rowhtml .= qq(<option value="$_")
            . ($_ eq $rows ? ' selected' : '')
            . qq(>$_</option>\n);
    }

    my $html .= qq(
      <table class="hidden">
      <tr>
       <td colspan="7"><img src="/gfx/blank.gif" width="10" height="10"></td>
     </tr>
    <tr >
       <td align="left">Number of rows of chromosomes:</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><select name="rows">$rowhtml</select></td>
       <td><img src="/gfx/blank.gif" width="20" alt=""></td>
       <td align="left">Height of the longest chromosome (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size="4" name="chr_length" value="$chr_length"></td>
     </tr>
     <tr>
       <td colspan="7"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
     <tr>
       <td align="left">Padding around chromosome (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size="3" name="h_padding" value="$h_padding"></td>
       <td><img src="/gfx/blank.gif" width="10" alt=""></td>
       <td align="left">Spacing between chromosomes (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size ="3" name="h_spacing" value="$h_spacing"></td>
     </tr>
     <tr>
       <td colspan="7"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
     <tr>
       <td align="left">Spacing between rows (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size="3" name="v_padding" value="$v_padding"></td>
       <td><img src="/gfx/blank.gif" width="10" alt=""></td>
     </tr>
     </table>);

    return $html;

}

#-----------------

sub _image_options_karyoview {

    my $self = shift;

    my $chr_length = CGI::escapeHTML($self->param('chr_length')) || 300;
    my $v_padding  = CGI::escapeHTML($self->param('v_padding'))  || 50;
    $v_padding = 20 unless $v_padding >= 20;

    my $h_padding = CGI::escapeHTML($self->param('h_padding')) || 4;
    $h_padding = 1 unless $h_padding >= 1;

    my $h_spacing = CGI::escapeHTML($self->param('h_spacing')) || 6;
    $h_spacing = 1 unless $h_spacing >= 1;

    my $rows = $self->param('rows') > 0 ? $self->param('rows') : 1;
    my $zmenu = lc($self->param('zmenu')) || 'on';
    my $col   = lc($self->param('col'))   || 'red';
    my $style = lc($self->param('style')) || 'widebox';

    # Available drawing styles
    my %style = (
        'box'           => "Filled box",
        'filledwidebox' => "Filled wide box",
        'widebox'       => "Outline wide box",
        'outbox'        => "Oversize outline box",
        'wideline'      => "Line",
        'lharrow'       => "Arrow left side",
        'rharrow'       => "Arrow right side",
        'bowtie'        => "Arrows both sides",
    );

    my $stylehtml;
    foreach (
        qw(box filledwidebox widebox outbox wideline lharrow rharrow bowtie ))
    {
        $stylehtml .= qq(<option value="$_")
            . ($_ eq $style ? ' selected' : '') . qq(>)
            . $style{$_}
            . qq(</option>\n);
    }

    my $colourhtml;
    foreach (
        qw( black purple magenta red orange brown green darkgreen blue darkblue violet grey darkgrey)
        )
    {
        $colourhtml .= qq(<option value="$_")
            . ($_ eq $col ? ' selected' : '')
            . qq(>$_</option>\n);
    }

    my $zmenuhtml;
    foreach ("on", "off") {
        $zmenuhtml .= qq(<option value="$_")
            . ($_ eq $zmenu ? ' selected' : '')
            . qq(>$_</option>\n);
    }

    my $rowhtml;
    foreach (1 .. 4) {
        $rowhtml .= qq(<option value="$_")
            . ($_ eq $rows ? ' selected' : '')
            . qq(>$_</option>\n);
    }

    # data display section
    my $html = qq(
      <table>
    <tr class="background1">
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="10"></td>
     </tr>
    <tr>
      <td align="left">Display data points in the following style:</td>
      <td align="left">&nbsp;&nbsp;</td>
      <td><select name="style"/>$stylehtml</select></td>
    </tr>
    <tr class="background1">
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
    <tr>
      <td align="left">Display data points in the following colour:</td>
      <td align="left">&nbsp;&nbsp;</td>
      <td><select name="col"/>$colourhtml</select></td>
    </tr>
    <tr class="background1">
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
    <tr>
      <td align="left">Show menus on mouseover:</td>
      <td align="left">&nbsp;&nbsp;</td>
      <td><select name="zmenu"/>$zmenuhtml</select></td>
    </tr>
    <tr>
       <td colspan="3"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
    </table>);

    # image options section
    $html .= qq(
      <table class="hidden">
      <tr>
       <td colspan="7"><img src="/gfx/blank.gif" width="10" height="10"></td>
     </tr>
    <tr >
       <td align="left">Number of rows of chromosomes:</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><select name="rows"/>$rowhtml</select></td>
       <td><img src="/gfx/blank.gif" width="20" alt=""></td>
       <td align="left">Height of the longest chromosome (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size="4" name="chr_length" value="$chr_length"></td>
     </tr>
     <tr>
       <td colspan="7"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
     <tr>
       <td align="left">Padding around chromosome (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size="3" name="h_padding" value="$h_padding"></td>
       <td><img src="/gfx/blank.gif" width="10" alt=""></td>
       <td align="left">Spacing between chromosomes (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size ="3" name="h_spacing" value="$h_spacing"></td>
     </tr>
     <tr>
       <td colspan="7"><img src="/gfx/blank.gif" width="10" height="5"></td>
     </tr>
     <tr>
       <td align="left">Spacing between rows (pixels):</td>
       <td align="left">&nbsp;&nbsp;</td>
       <td><input type="text" size="3" name="v_padding" value="$v_padding"></td>
       <td><img src="/gfx/blank.gif" width="10" alt=""></td>
     </tr>
    </table>);

    return $html;

}

sub _image_options_chromoview {

    my $self = shift;
    my $html = '<table>';

    my %checkboxes = (
        'maxmin'    => 'Show Max/Min lines',
        'Vsnps'     => 'Show SNP frequency',
        'Vgenes'    => 'Show gene frequency',
        'Vpercents' => 'Show GC content frequency'
    );
    my $count;
    foreach my $box (keys %checkboxes) {
        my $newrow = $count % 2;

        # if even number, start a new row
        if (!$newrow) {
            $html .= '<tr>';
        }
        my $tick = $self->param($box) eq 'on' ? ' checked="checked"' : '';
        $html
            .= qq(<td><input type="checkbox" name="$box" value="on"$tick>&nbsp;$checkboxes{$box}<br /><br /></td>\n);

        # if even number, insert a spacer cell
        if (!$newrow) {
            $html
                .= qq(<td><img src="/gfx/blank.gif" width="50" height="1" alt=""></td>\n);
        }

        # otherwise end the row
        else {
            $html .= "</tr>\n";
        }
        $count++;
    }
    $html .= '</table>';

    return $html;
}

1;
