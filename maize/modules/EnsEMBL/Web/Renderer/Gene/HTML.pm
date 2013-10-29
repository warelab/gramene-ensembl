package EnsEMBL::Web::Renderer::Gene::HTML;

=head1 NAME

EnsEMBL::Web::Renderer::Gene::HTML.pm 

=head1 SYNOPSIS

This object creates HTML to be output to the HTML output object

=head1 DESCRIPTION

    $gene_renderer = $gene_data->renderer;
    $gene_renderer->outputGenericGeneTable();        
    $gene_renderer->gene_transcript_table_title;
    
 This object contains wrappers for common display 'groups' and also more
 granular calls that can be reused to create different page layouts/structure

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 CONTACT

Jim Stalker- jws@sanger.ac.uk

=cut

use strict;
use warnings;
no warnings "uninitialized";

use vars qw( @ISA );
use EnsEMBL::Web::Renderer::HTML;
use EnsEMBL::Web::Renderer::Gene;
use EnsEMBL::Web::Output::HTML::DropDown::MenuContainer::GeneSNPView;
use Data::Dumper;

@ISA = qw( EnsEMBL::Web::Renderer::HTML EnsEMBL::Web::Renderer::Gene );

=head2 outputGenericGeneTable

 Arg[1]      : none
 Example     : $genedata->renderer->outputGenericGeneTable
 Description : Wrapper for the table at top of geneview
 Return type : none

=cut

sub outputGenericGeneTable() {
  my $self = shift;
  $self->gene_table_title;
  $self->Output->print_two_col_table
    (
     $self->gene_display_ID,
     $self->stable_id,
     $self->genomic_location,
     $self->alternative_alleles,
     $self->gene_description,
     $self->prediction_method,
     $self->genesequence_link,
     $self->exportview_link,
     $self->genesnpview_link,
     $self->transcript_neighbourhood,
     $self->orthologue_matches,
     $self->paralogue_matches,
     $self->disease_matches,
     #$self->das_configurator,
     $self->das_annotation,
    );

}


=head2 outputGeneSequenceTable

  Arg [1]   : 
  Function  : Produces table containing gene sequence
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub outputGeneSequenceTable {
  my $self = shift;
  $self->gene_table_title;
  $self->Output->print_two_col_table
    (
     $self->gene_display_ID,
     $self->stable_id,
     $self->genomic_location,
     #$self->transcript_neighbourhood,
     $self->markup_options,
     $self->markedup_geneseq,
    );

}


=head2 gene_table_title

 Arg[1]      : none
 Example     : $genedata->renderer->gene_table_title
 Description : prints the table title for the gene table
 Return type : none

=cut

# may move this onto the table call somehow
sub gene_table_title{
    my $self = shift;
    my $type = $self->DataObj->gene_type;
    $self->Output->generic_table_title("<h3>$type Gene Report</h3>" );
}

=head2 gene_transcript_table_title

 Arg[1]      : none
 Example     : $genedata->renderer->gene_transcript_table_title
 Description : prints the table title for the gene=transcript tables
 Return type : string - html

=cut

sub gene_transcript_table_title{
    my $self = shift;
    $self->Output->generic_table_title( '<h3>Transcript/Translation Summary</h3>' );
}

=head2 gene_display_ID

 Arg[1]      : (optional) String
               table header text
 Example     : $genedata->renderer->gene_display_ID($string)
 Description : renders the xref_display_id and type in two_col_table format
 Return type : key-value pair, label and html

=cut

sub gene_display_ID {
    my $self = shift;
    my $gene = $self->DataObj();
    my $label = shift || 'Gene';
    return unless $gene->display_xref();    

    my ($display_name, $dbname) = $gene->display_xref();

    # link to external database
    my $ext_id = $gene->get_external_id($dbname);
    $ext_id =~ s/\.\d+//; # GRAMENE Strip accession. TODO: Fix in DB
    my $linked_display_name = $display_name;
    if ($ext_id) {
        $linked_display_name = $self->url(
            $display_name, $self->ExtURL->get_url($dbname, $ext_id)
        );
    }

    my $site_type = ucfirst(lc($SiteDefs::ENSEMBL_SITETYPE));
    my $html = qq {<b>$linked_display_name </b> <small>($dbname ID)</small>&nbsp;<small>};
    # GRAMENE - the following is not needed
    # (to view all $site_type genes linked to the name <a href="/$ENV{ENSEMBL_SPECIES}/featureview?type=Gene&id=$display_name">click here</a>)<small>};
    return ($label, $html);
}

=head2 stable_id

 Arg[1]      : (optional) String
               table header text
 Example     : $genedata->renderer->stable_id($string)
 Description : Renders the ensembl_stable_id in two_col_table format
 Return type : key-value pair, label and html

=cut

sub stable_id {
  my $self = shift;
  my $gene =  $self->DataObj();
  my $db_type = $gene->db_type ;
  my $db = $self->param('db');
  my $label = shift || "$db_type Gene ID";
  my $exturl = $self->ExtURL;
  my $geneid = $gene->stable_id;

  my $type = '';
  if( $db_type eq 'Vega' ){
    my $tmpl = qq(<small>(%s)</small> %s);
    my $link =  $self->url("[View in Vega]",
                           $exturl->get_url('VEGA_GENE', $gene->stable_id),
                           '_blank');
    $type = sprintf( $tmpl, $gene->feature_type, $link );
  }
  if( $self->param('_gene_sequence') ){
    $geneid = $self->url($geneid, "geneview?gene=$geneid&db=$db");
  }

  $geneid = sprintf( "<b>%s</b>",$geneid ). $type;

  return ($label, $geneid);
}

=head2 vega_stable_id

 Arg[1]      : (optional) String
               table header text
 Example     : $genedata->renderer->vega_stable_id($string)
 Description : Renders the stable_id in two_col_table format
               Vega version of stable_id, with link to ensembl
 Return type : key-value pair, label and html

=cut

sub vega_stable_id {
    my $self = shift;
    my $label = shift || "Locus ID";
    my $gene =  $self->DataObj();
    my $db = $self->param('db');
    my $geneid = $gene->stable_id;
    my $ensembl_url = $self->get_ExtURL('ENS_GENEVIEW', $geneid);
    if ($self->param('_gene_sequence')) {
        $geneid = $self->url($geneid, "geneview?gene=$geneid&db=$db");
    }
    my $html = qq(<b>$geneid</b>);

    ## add crosslink to ensembl only if enabled in species.ini and not on
    ## a haplotype
    if ($self->DataObj->species_defs->ENSEMBL_CROSSLINKS && !($gene->gene->slice->seq_region_name =~ /_/)) {
        $html .= qq( [<a href="$ensembl_url">View in Ensembl</a>]);
    }
    
    return ($label, $html);
}

=head2 version

 Arg[1]      : (optional) String
               table header text
 Example     : $genedata->renderer->version($string)
 Description : Renders the gene version in two_col_table format
 Return type : key-value pair, label and html

=cut

sub version {
    my $self = shift;
    my $gene =  $self->DataObj();
    my $label = shift || "Version";
    my $html = $gene->version;
    return ($label, $html);
}

=head2 mod_date

 Arg[1]      : none
 Example     : $genedata->renderer->mod_date
 Description : Renders the dates the gene was modified and created in two_col_table format
 Return type : key-value pair, label and html

=cut

sub mod_date {
    my $self = shift;
    my $gene =  $self->DataObj();  

    my $label = "Date";
    my ($html1,$mod_date) = $gene->mod_date;
    my ($html2,$created_date) = $gene->created_date;

    my $html = qq[Gene last modified on $html1 (<small>Created on $html2</small>)];

#compare with SPECIES.INI entry
#    my $Species_defs = SpeciesDefs->new();
#    my $old_db_date = $Species_defs->LAST_UPDATE;
    
#    if ($mod_date gt $old_db_date ) {
#	$html .= qq[ (<b>modified since last release</b>)];
#    }

#    if ($created_date gt $old_db_date ) {
#	$html .= qq( <b>New Gene<b>);
#    }

    return ($label, $html);
}

=head2 author

 Arg[1]      : (optional) String
               Label
 Example     : $genedata->renderer->author
 Description : returns the annotation author for use in two_col_table
 Return type : key-value pair of label/html

=cut

sub author {
    my $self = shift;
    my $label = shift || "Author";
    my $gene =  $self->DataObj();
    my $author = $gene->get_author_name;
    my $html;
    if ($author) {
        $html .= "This locus was annotated by " . $author . " ";
        $html .= $self->email_URL($gene->get_author_email);
    } else {
        $html = "unknown";
    }
    return ($label, $html);
}

=head2 remarks

 Arg[1]      : (optional) String
               Label
 Example     : $genedata->renderer->remarks
 Description : returns annotation remarks for use in two_col_table
 Return type : key-value pair of label/html

=cut

sub remarks {
    my $self = shift;
    my $label = shift || "Remarks";
    my $gene =  $self->DataObj();
    my $remarks = $gene->get_remarks;
    my $html;
    if (@{$remarks} && $remarks->[0]) {
        $html = join ('<br />', @{$remarks});
    } else {
        $html = "No remarks";
    }
    return ($label, $html);
}

=head2 synonyms

 Arg[1]      : (optional) String
               Label
 Example     : $genedata->renderer->synonyms
 Description : returns locus synonyms for use in two_col_table
 Return type : key-value pair of label/html

=cut

sub synonyms {
    my $self = shift;
    my $label = shift || "Alternative Symbols";
    my $gene =  $self->DataObj();
    my $synonyms = $gene->get_synonyms;
    my $db = $self->param('db') || 'core';
    my $html;
    foreach (keys %{$synonyms}) {
        if ($synonyms->{$_}) {
            $html .= qq(<a href='geneview?gene=$synonyms->{$_}&db=$db'>$_</a><br />);
        } else {
            $html .= "$_<br />";
        }
    }
    return unless $html;
    return ($label, $html);
}

=head2 

  Arg[1]      : (optional) String
                Label
  Example     : $genedata->renderer->gene_type
  Description : returns gene gene_type for use in two_col_table
  Return type : key-value pair of label/html

=cut

sub gene_type {
    my $self = shift;
    my $label = shift || "Type";
    my $gene =  $self->DataObj();
    my $gene_type = $gene->get_gene_type;
    my $species = $gene->Input->species;
    my $helplink = $gene->species_defs->SPECIES_SHORT_NAME . "_gene_classification";
    my $html = qq($gene_type [<a href="javascript:X=hw\('$species', '$helplink', ''\)">Definition</a>]);
    return ($label, $html);
}

=head2 genomic_location

 Arg[1]      : none
 Example     : $genedata->renderer->genomic_location
 Description : Renders the genomic location of the gene in two_col_table format
 Return type : key-value pair, label and html

=cut

sub genomic_location{
    my $self = shift;
    my $gene =  $self->DataObj;
    my $geneid = $gene->stable_id;
    my ($chr, $start, $end, $contig, $contig_start) = $gene->get_genomic_location();
    my @alt_locs = @{$gene->get_alternative_locations};
    my $label  = 'Genomic Location';        
    my $html = '';
    
    if (!$start && !$end && !$chr){
        $html .=  qq(<i>This gene cannot be located on the current assembly</i>\n);
    }
    else {
        my $mbase = $gene->bp_to_nearest_unit($start,1);
	my $seqregion = $gene->coord_system_name;
        my $contigview_link = "contigview?region=$contig&fpos_start=$contig_start&fpos_end=".($contig_start + 1)."&fpos_context=200000&highlight=$geneid";
	
	# Genomic Location
        $html .= qq(
            <b>View gene in genomic location: </b>
            <a href="/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$start&vc_end=$end&highlight=$geneid">$start - $end bp ($mbase)</a> 
                on $seqregion $chr<br />
            <b>This gene is located in sequence: </b>
            <a href="$contigview_link">$contig</a>
	);

	# Haplotype/PAR locations
	foreach my $loc (@alt_locs){
	    my ($altchr, $altstart, $altend, $altseqregion) = @$loc;
	    my $altmbase = $gene->bp_to_nearest_unit($altstart,1);
	    $html .= qq(
		<br><b>Haplotype/PAR location: </b>
            <a href="/$ENV{'ENSEMBL_SPECIES'}/contigview?l=$altchr:$altstart-$altend&highlight=$geneid">$altstart - $altend bp ($altmbase)</a> 
                on $altseqregion $altchr
	    );
	}
    }
    $html .=  qq(<br />);
    return ($label, $html)
}


=head2 alternative_alleles

 Arg[1]      : none
 Example     : $genedata->renderer->alternative_alleles
 Description : Renders the haplotype alleles of the gene in two_col_table format
 Return type : key-value pair, label and html

=cut

sub alternative_alleles{
    my $self = shift;
    my $gene =  $self->DataObj;
    my $geneid = $gene->stable_id;
    my @alleles = @{$gene->get_alternative_alleles};
    return unless @alleles;
    my $label  = 'Haplotype Alleles';        
    my $html = '';
    
    # Genomic Location
    foreach my $allele(@alleles){
	my ($allele_name, $haplotype) = @$allele;
	$html .= qq(
	    <b>Allele: </b>
	    <a href="/$ENV{'ENSEMBL_SPECIES'}/geneview?gene=$allele_name">$allele_name</a> 
		on haplotype $haplotype<br />
	);
    }
    return ($label, $html)
}



=head2 gene_description

 Arg[1]      : none
 Example     : $genedata->renderer->gene_description
 Description : Renders the gene description of the gene in two_col_table format
 Return type : key-value pair, label and html

=cut

sub gene_description {
    my $self = shift ;
    my $urls = $self->ExtURL;
    my $description = $self->DataObj->gene_description();
    $description =~ s/EC\s+([-*\d]+\.[-*\d]+\.[-*\d]+\.[-*\d]+)/$self->EC_URL($1)/e;
    $description =~ s/\[\w+:(\w+)\;\w+:(\w+)\]//g;
    my ($edb, $acc) = ($1, $2);
    my $url         = "[<a href='".$urls->get_url($edb, $acc)."'>Source: $1 ($2)</a>]" unless !$acc ;
    my $label = 'Description';
    my $html = qq($description <small>$url</small>&nbsp;<br />);
    return ($label, $html);
}

=head2 prediction_method

 Arg[1]      : none
 Example     : $genedata->renderer->prediction_method
 Description : Renders the prediction_method of the gene in two_col_table format
 Return type : key-value pair, label and html

=cut

sub prediction_method {
    my $self = shift ;
    my $genedata = $self->DataObj;
    my $db = $genedata->feature_type ;  
    my $label = ( $db eq 'vega' ? 'Curation' : 'Prediction' ).' Method';
    my $gene_prediction_text = $genedata->get_prediction_method ."&nbsp;";      
    return ($label, $gene_prediction_text );
}

=head2 transcript_neighbourhood

 Arg[1]      : (optional) String
               table header text
 Example     : $genedata->renderer->transcript_neighbourhood($string)
 Description : Renders the transcript_neighbourhood in two_col_table format
 Return type : key-value pair, label and html

=cut

sub transcript_neighbourhood {
    my $self = shift ;
    my $genedata = $self->DataObj;
    my $db = $genedata->get_db() ;
    #my $label = shift || ( $db eq 'vega' ? 'Curation' : 'Prediction' ).' Transcript';
    my $label = shift || " Transcript Structure";
    my $transcripts = $genedata->get_all_transcripts;
    my $message = qq(<p class="gv_warning">A large number of transcripts have been returned for this gene. To reduce render time for this page the protein and transcript  information
                     will not be displayed. To view this information please follow the transview and protview links
                     below. <br /></p> ) if(scalar(@{$transcripts}) > 17);              

# add table call here
    my $html = qq($message <table class="hidden">);
    my $i = 1 ;         
    foreach my $trans (sort{$a->stable_id cmp $b->stable_id} @{$transcripts}){
        $html .= qq(<tr><td valign="top" align="right"><b>$i:</b>&nbsp;</td>);
        my $trans_stable_id = $trans->stable_id();      

        if ($trans->display_xref()){
            my ($trans_display_id, $db_name) = $trans->display_xref();
            my $display_id = $message ? $trans_display_id : qq(<a href="#$trans_stable_id">$trans_display_id</a>) ;
            $html .=  qq(<td nowrap> $display_id <small> ($trans_stable_id) </small> </td>);
        }else{
            my $display_id = $message ? $trans_stable_id : qq(<a href="#$trans_stable_id">$trans_stable_id </a>) ;
            $html .= qq(<td nowrap> $trans_stable_id &nbsp; </td>);
        }   
        $html .= qq(<td valign="top" nowrap> [<A href="/$ENV{'ENSEMBL_SPECIES'}/transview?transcript=$trans_stable_id&db=$db">Transcript information</a>] </td>);

    if ($self->param('show_vega_evidence_link')) {
       $html .= qq(<td valign="top" nowrap> [<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$trans_stable_id&db=$db">Exon information & supporting evidence</a>] </td>);
   } else {
        $html .= qq(<td valign="top" nowrap> [<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$trans_stable_id&db=$db">Exon information</a>] </td>);
   }
        if (my $pep_id = $trans->translation->stable_id){
            $html .= qq(<td valign="top" nowrap> [<a href="/$ENV{'ENSEMBL_SPECIES'}/protview?peptide=$pep_id&db=$db">Protein information</a>] </td>);
        }else{
            $html .= qq( <td valign="top" nowrap>[No translation]</td> );
        }
        $html .=  qq(</tr>);    
        $i++;
    }
    $html .= qq( <tr><td colspan = "5">&nbsp;</td></tr>);   
    my $image = $self->transcript_neighbourhood_image;       
    $html .= qq(</table>$image);     
    return ($label, $html);
}

=head2 transcript_neighbourhood_image

 Arg[1]      : none
 Example     : $genedata->renderer->transcript_neighbourhood_image
               $image = $self->transcript_neighbourhood_image
 Description : Renders the transcript_neighbourhood_image and clickable 
               imagemap
 Return type : AHREF image url

=cut

sub transcript_neighbourhood_image {
  my $self = shift;
  my $wu_conf = $self->DataObj->configure_neighbourhood_image();
  my $mc = $self->DataObj->get_neighbourhood_menu(  $wu_conf );
  my $dc = $self->DataObj->get_neighbourhood_image( $wu_conf );
  my ($image, $imagemap) = $self->render_image_imagemap( "geneneigh", $dc );
  my $htmpl = qq(
<table align="center" class="hidden">%s
<tr><td><img border="0" src="%s" usemap="#geneview_transnhood"></td></tr>
</table>
<map name="geneview_transnhood">%s</map>);

  my $src = "$image";

  # for now, disable menu in Vega
  my $menu_html;
  unless ($SiteDefs::ENSEMBL_SITETYPE eq 'Vega') {
    $menu_html = $mc->render_html;
  }
  
  my $html = sprintf( $htmpl, $menu_html, $src, $imagemap );
  $self->Output->appendHTML( $mc->render_js ); # Save JS to end
  return $html;
}

sub genesnpview_image {
  my $self = shift;
  my( $dc,$conf) = $self->DataObj->render_genesnpview_image( );
  my $mc = EnsEMBL::Web::Output::HTML::DropDown::MenuContainer::GeneSNPView->new(
    'species' => $self->Input->species,
    'script'  => 'genesnpview',
    'config'  => $conf,
    'panel'   => 'genesnpview',
    'db'      => $self->param('db'),
    'gene'    => $self->DataObj->stable_id
  );
  $mc->add_left_menu(  'SNPClasses' );
  $mc->add_left_menu(  'SNPValid'   );
  $mc->add_left_menu(  'SNPTypes'   );
  $mc->add_left_menu(  'SNPContext' );
  $mc->add_left_menu(  'ImageSize'  );
  $mc->add_left_menu(  'THExport'  );
  $mc->add_right_menu( 'SNPHelp'    );

  $self->Output->printHTML( qq(<table border="0" cellspacing="0" cellpadding="0" align="center">) );
  $self->Output->printHTML( $mc->render_html );
  $self->{'additional_image_types'} = [];
  foreach my $format (qw(pdf postscript svg)) {
    push @{ $self->{'additional_image_types'} }, $format if $conf->get('_settings',"opt_$format");
  }

  my( $filename, $imagemap, $image_links_hash ) = $self->render_image_imagemap( "gsv", $dc );
  my $extra_html;
  if( keys %$image_links_hash ) {
    $extra_html = '<div align="right">'. join( '; ',
      map { qq(<a href="$image_links_hash->{$_}" target="_blank">Render as ).uc($_).'</a>' } keys %$image_links_hash
    ).'.</div>';
  }
  $self->Output->printHTML( qq(<tr><td><img border="0" src="$filename" usemap="#genesnpview">$extra_html</td></tr></table>) );
  $self->Output->printHTML( qq(<map name="genesnpview">$imagemap</map>) );
  return $mc->render_js;
}

#----------------------------------------------------------------------


my $martview_href = qq( /Multi/martview?species=%s&focus=snp&stage=output&stage_initialised=start&stage_initialised=filter&named_snp=1&named_snp_filter=FS_RefSNP_ID&named_snp_list=%s );

my $martview_link = "<A href=%s>[Export SNPs]</A>";

sub genesnpview_table {
  my $self = shift;
  my $genesnp_data = $self->DataObj->genesnp_table_data();
  my @snp_ids = map{ $_->{ID} =~ />\s*(\d+)\s*</; $1 } @{$genesnp_data}; #Geesh!
 # if( scalar @snp_ids ){
 #   my $link = sprintf( $martview_link, 
 #                       (sprintf( $martview_href, 
 #                                 $ENV{'ENSEMBL_SPECIES'},
 #                                 join( ",", @snp_ids ) ) ) );
 #   push @$genesnp_data, {ID=>"$link"  };
 # }
  return $self->Output->print_spreadsheet_table(
    [ { 'key' => 'ID', 'align' => 'center' },
      { 'key' => 'class', 'align' => 'center' },
      { 'key' => 'alleles', 'align' => 'center' },
      { 'key' => 'ambiguity', 'align' => 'center' },
      { 'key' => 'status', 'align' => 'center' },
      { 'key' => 'chr' , 'align' => 'center' },
      { 'key' => 'pos' , 'align' => 'center' },
    ],
    $genesnp_data,
    { 'align' => 'center', 'width'=> 900 }
  );
}


sub genesnpview_table_transcript {
  my $self = shift;
  my $genesnp_data = $self->DataObj->genesnp_table_data();
  my @snp_ids = map{ $_->{ID} =~ />\s*(\w+)\s*</; $1 } @{$genesnp_data};#Geesh!
#  if( scalar @snp_ids ){
#    my $link = sprintf( $martview_link,
#                        (sprintf( $martview_href,
#                                  $ENV{'ENSEMBL_SPECIES'},
#                                  join( ",", @snp_ids ) ) ) );
#    push @$genesnp_data, {ID=>"$link" , 'skip' => 1 };
#  }
  for my $transcript_data ( @{$self->DataObj->get_all_transcripts} ) {
    my $peptide = $transcript_data->translation;
    $self->print( "<h3>SNPs for Transcript @{[$transcript_data->stable_id]}",
      ($peptide ? " (Peptide @{[$peptide->stable_id]})" : '' ), "</h3>" );
    my $transcript_snps = $self->DataObj->{'genesnp_transcripts'}->{$transcript_data->stable_id};
    my @transsnp_data = ();
    foreach my $row ( @$genesnp_data ) { 
      # warn Data::Dumper::Dumper( $row );
      if( $row->{'skip'} == 1 ) {
        push @transsnp_data, $row;
      } else {
        (my $T = $row->{'pos'} ) =~ s/-\d+//;
        my $transsnp = $transcript_snps->{'snps'}->{ $row->{'raw_id'} };
        if( $transsnp && $row->{'end'} >= $transcript_snps->{'start'}-$transcript_snps->{'extent'} &&
                         $row->{'start'} <= $transcript_snps->{'end'}  +$transcript_snps->{'extent'} ) {
          $row->{'snptype'}  = $transsnp->consequence_type;
          if( $transsnp->translation_start ) { 
            $row->{'aachange'} = $transsnp->pep_allele_string;
            $row->{'aacoord'}  = $transsnp->translation_start.' ('.(($transsnp->cdna_start-$transcript_data->transcript->cdna_coding_start)%3+1).')';
          } else {
            $row->{'aachange'} = '-'; $row->{'aacoord'} = '-';
          }
          push @transsnp_data, $row;
        }
      } 
    }
    $self->Output->print_spreadsheet_table(
      [ { 'key' => 'ID', 'align' => 'center' },
        { 'key' => 'class', 'align' => 'center' },
        { 'key' => 'alleles', 'align' => 'center' },
        { 'key' => 'ambiguity', 'align' => 'center' },
        { 'key' => 'status', 'align' => 'center' },
        { 'key' => 'chr' , 'align' => 'center' },
        { 'key' => 'pos' , 'align' => 'center' },
        { 'key' => 'snptype', 'title' => 'SNP type', 'align' => 'center' },
        { 'key' => 'aachange', 'title' => 'AA change', 'align' => 'center' },
        { 'key' => 'aacoord',  'title' => 'AA co-ordinate', 'align' => 'center' } 
      ],
      \@transsnp_data,
      { 'align' => 'center', 'width'=> 900 }
    );
  }
}

=head2 orthologue_matches

 Arg[1]      : none
 Example     : $genedata->renderer->orthologue_matches
 Description : Renders the orthologue_matches in two_col_table format
 Return type : key-value pair, label and html

=cut

sub orthologue_matches {
  my $self =shift;
  my $dataobj = $self->DataObj;
  my $orthologue = $dataobj->get_homology_matches('ENSEMBL_ORTHOLOGUES');
  return unless $orthologue;

## call table method here
  my $urls = $self->ExtURL;   
  my $db = $dataobj->get_db() ; 
  my %orthologue_list = %{$orthologue};
  my $label = 'Orthologue Prediction';
  my $html = qq(The following gene(s) have been identified as putative orthologues by reciprocal BLAST analysis:<br />
             <table class="hidden">);
  $html .= qq(<tr valign='top'><td><b>Species</b></td><td><b>Type</b></td><td nowarp><b>dN/dS</b></td><td nowrap><b>Gene identifier</b></td></tr>);
  my %orthologue_map = qw(SEED BRH PIP RHS);

  my %SPECIES;
  my $STABLE_ID = $dataobj->stable_id; my $C = 1;
  my $FULL_URL = qq(/$ENV{'ENSEMBL_SPECIES'}/multicontigview?gene=$STABLE_ID);
  my $ALIGNVIEW = 0;
  foreach my $stable_id (sort keys %orthologue_list){
    $html .= qq(<tr><td nowrap valign='top'><b><i>$orthologue_list{$stable_id}{'spp'}</i></b></td>);
    my $description = $orthologue_list{$stable_id}{'description'};
       $description = "No description" if $description eq "NULL";
    my $orthologue_desc = $orthologue_map{ $orthologue_list{$stable_id}{'homology_desc'} } || $orthologue_list{$stable_id}{'homology_desc'};
    my $orthologue_dnds_ratio = $orthologue_list{$stable_id}{'homology_dnds_ratio'};
       $orthologue_dnds_ratio = "--" unless (defined $orthologue_dnds_ratio);
       $html .= qq(<td valign='top'>$orthologue_desc</td><td valign='top'><small class='normal'>$orthologue_dnds_ratio</small></td>);
    if($orthologue_list{$stable_id}{'display_id'}) {
      (my $spp = $orthologue_list{$stable_id}{'spp'}) =~ tr/ /_/ ;
      $SPECIES{ $spp } = 1;
      my $EXTRA = qq(<small class="normal">[<a href="/$ENV{'ENSEMBL_SPECIES'}/multicontigview?gene=$STABLE_ID&s1=$spp&g1=$stable_id">MultiContigView</a>]</small>);
      if(  $orthologue_desc ne 'DWGA' ) {
        $EXTRA .= qq(&nbsp;<small class="normal">[<a href="/$ENV{'ENSEMBL_SPECIES'}/alignview?class=Homology&gene=$STABLE_ID&g1=$stable_id">Align</a>]</small>);
        $ALIGNVIEW = 1;
      }
      $FULL_URL .= "&s$C=$spp&g$C=$stable_id";$C++;
      my $link = qq(/$spp/geneview?gene=$stable_id&db=$db);
      if( $description =~ s/\[\w+:(\w+)\;\w+:(\w+)\]//g ) {
        my ($edb, $acc) = ($1, $2);
        $description .= "[<a href='".$urls->get_url($edb, $acc)."'>Source: $edb ($acc)</a>]" unless !$acc;
      }        
      $html .= qq(<td><a href="$link">$stable_id</a> )."( ".$orthologue_list{$stable_id}{'display_id'}." )".qq( $EXTRA<br />
                <small class='normal'>$description</small></td></tr>);
    } else {
      $html .= qq(<td>$stable_id<br /><small class='normal'>$description</small></td></tr>);
    }
  }
  $html .= qq(</table>);
  if( keys %SPECIES ) {
    # $html .= qq(<p><a href="$FULL_URL">View all genes in MultiContigView</a>;);
    $html .= qq( <a href="/$ENV{'ENSEMBL_SPECIES'}/alignview?class=Homology&gene=$STABLE_ID">View alignments of homologies</a>.</p>) if $ALIGNVIEW;
  }
  $html .= qq(<small>UBRH - (U)nique (B)est (R)eciprocal (H)it<br />
                MBRH - one of (M)any (B)est (R)eciprocal (H)its<br />
              RHS   = Reciprocal Hit based on Synteny around BRH<br />
              DWGA  = Derived from Whole Genome Alignment</small>);
  return($label, $html);
}


=head2 paralogue_matches

 Arg[1]      : none
 Example     : $genedata->renderer->paralogue_matches
 Description : Renders the paralogue_matches in two_col_table format
 Return type : key-value pair, label and html

=cut

sub paralogue_matches {
# Yes, this is very similar to orthologue_matches, but with a number of minor diffs
  my $self =shift;
  my $dataobj = $self->DataObj;
  my $paralogue = $dataobj->get_homology_matches('ENSEMBL_PARALOGUES');
  return unless $paralogue;

## call table method here
  my $urls = $self->ExtURL;   
  my $db = $dataobj->get_db() ; 
  my %paralogue_list = %{$paralogue};
  my $label = 'Paralogue Prediction';
  my $html = qq(The following gene(s) have been identified as putative paralogues:<br />
             <table class="hidden">);
  $html .= qq(<tr valign='top'><td nowarp><b>dN/dS</b></td><td nowrap><b>Gene identifier</b></td></tr>);

  my $STABLE_ID = $dataobj->stable_id; my $C = 1;
  foreach my $stable_id (sort keys %paralogue_list){
    my $description = $paralogue_list{$stable_id}{'description'};
       $description = "No description" if $description eq "NULL";
    my $paralogue_dnds_ratio = $paralogue_list{$stable_id}{'homology_dnds_ratio'};
       $paralogue_dnds_ratio = "--" unless (defined $paralogue_dnds_ratio);
       $html .= qq(<td valign='top'><small class='normal'>$paralogue_dnds_ratio</small></td>);
    if($paralogue_list{$stable_id}{'display_id'}) {
      (my $spp = $paralogue_list{$stable_id}{'spp'}) =~ tr/ /_/ ;
       my $link = qq(/$spp/geneview?gene=$stable_id&db=$db);
       if( $description =~ s/\[\w+:(\w+)\;\w+:(\w+)\]//g ) {
         my ($edb, $acc) = ($1, $2);
         $description .= "[<a href='".$urls->get_url($edb, $acc)."'>Source: $edb ($acc)</a>]" unless !$acc;
       }        
       $html .= qq(<td><a href="$link">$stable_id</a> )."( ".$paralogue_list{$stable_id}{'display_id'}." )".
       qq(<br /><small class='normal'>$description</small></td></tr>);
    } 
    else {
      $html .= qq(<td>$stable_id<br /><small class='normal'>$description</small></td></tr>);
    }
  }
  $html .= qq(</table>);
  return($label, $html);
}


=head2 database_matches

  Arg[1]      : (optional) String
                Label
  Example     : $genedata->renderer->database_matches
  Description : Renders database matches for gene in two_col_table format
                This is a tempory hack for Vega to get xrefs on genes which
                should be mapped to transcripts/translation (but aren't)
  Return type : Key / value pair - label and HTML

=cut

sub database_matches {
    my $self = shift;
    my $label = shift || 'Database Matches';
    my $data = $self->DataObj;
    my $gene = $data->gene;
    my $database_matches = $data->get_database_matches($gene);
    my $db = $data->get_db;
    my %links = %{$self->_sort_database_matches($database_matches)};
    return unless (%links);

    my $html = qq(<table class="hidden">);
    foreach my $key (sort keys %links){
        if (scalar (@{$links{$key}}) > 0) {
            my @sorted_links = sort @{$links{$key}};
            $html .= qq(<tr valign="top"><td nowrap="nowrap"><b>$key:</b></td><td>);
          
            if ($sorted_links[0] =~ /<br/i) {
                $html .= join(' ', @sorted_links );
            } else {
                # Need a BR each 5 entries
                $html .= qq(<table class="hidden"><tr>);
                my @sorted_lines;
                for (my $i=0; $i<@sorted_links; $i++) {
                    my $line_num = int($i/4);
                    if( ref(  $sorted_lines[$line_num] ) ne 'ARRAY' ){$sorted_lines[$line_num] = [];}
                    push( @{$sorted_lines[$line_num]}, "<td>".$sorted_links[$i]."</td>" );
                }
                $html .= join( '</tr><tr>', map{ join( ' ', @$_ ) } @sorted_lines );
                $html .= qq(</tr></table>);
            }
            $html .= qq(</td></tr>);
        }
    }
    $html .= qq(</table>); 
    return ($label, $html);
}

=head2 _sort_database_matches

 Arg[1]      : Arrayref
               A listref of database matches
 Example     : $genedata->renderer->_sort_datbase_matches
 Description : sorts the database matches
 Return type : hashref of database matches

=cut

sub _sort_database_matches {
    my $self = shift;
    my $database_matches = shift;
    my $data = $self->DataObj;
    my $database = $data->database;
    my $db = $data->get_db() ;
    my $urls = $self->ExtURL;
    my %links ;
    # Nice names    
    my %nice_names = (  
            'protein_id'            => 'Protein ID', 
            'drosophila_gene_id'    => 'Drosophila Gene',
            'flybase_gene'          => 'Flybase Gene',
            'flybase_symbol'        => 'Flybase Symbol',
            'affy_hg_u133'          => 'Affymx Microarray U133',
            'affy_hg_u95'           => 'Affymx Microarray U95',
            'anopheles_symbol'      => 'Anopheles symbol',
            'sanger_probe'          => 'Sanger Probe',
            'wormbase_gene'         => 'Wormbase Gene',
            'wormbase_transcript'   => 'Wormbase Transcript',
            'wormpep_id'            => 'Wormpep ID',
            'briggsae_hybrid'       => 'Briggsae Hybrid',
            'sptrembl'              => 'SpTrEMBL',
        );
                       
    foreach my $type (sort @{$database_matches}) { 
        my $link = "";
        my $join_links = 0;
        my $externalDB = $type->database();
        my $display_id = $type->display_id();
        my $primary_id = $type->primary_id();
    
        # remove all orthologs  
        next if ($type->status() eq 'ORTH');
        
        # Ditch celera genes from FlyBase
        next if ($externalDB =~ /^flybase/i && $display_id =~ /^CG/ );

        # Genes shouldn't really have protein associations of their own
        next if uc($externalDB) eq "REFSEQ";
        next if uc($externalDB) eq "SWISSPROT";

        # remove internal links to self and transcripts
        next if $externalDB eq "Vega_gene";
        next if $externalDB eq "Vega_transcript";
        next if $externalDB eq "Vega_translation";
        
        if ($externalDB eq "GO") {
            if ($data->database('go')) {
                $link = '<a href="'.
                $urls->get_url($externalDB, $display_id).
                '">'. $display_id. '</a>';
            }
            else {
                $link = '<a href="'.
                $urls->get_url('AMIGO', $display_id).
                '">'. $display_id. '</a>';
            }
        } elsif ($externalDB eq "protein_id") {
            # Can't link to srs if there is an Version - so strip it off
            $primary_id =~ s/(.*)\.\d+$/$1/o;
        }
        if ($urls->is_linked($externalDB)) {
            $link = '<a href="'.$urls->get_url($externalDB, $primary_id).'">'. $display_id. '</a>';
            if ($externalDB =~ /^AFFY_HG_U\d*/i) {
                $link = '<a href="' .$urls->get_url($externalDB, $display_id) .'">'. $display_id. '</a>';
            }
            if( $type->isa('Bio::EnsEMBL::IdentityXref') ) {
                $link .=' <small> [Target %id: '.$type->target_identity().'; Query %id: '.$type->query_identity().']</small>';            
                $join_links = 1;    
            }
            if ($externalDB =~/^(SWISS|SPTREMBL)/i) { # add Search GO link            
                $link .= ' [<a href="'.$urls->get_url('GOSEARCH',$primary_id).'" class="small">Search GO</A>]';
            }
            if( $join_links  ) {
                $link .= '<br>';
            }
        } else {
            $link = " $display_id ";
        }
        my $display_name = $nice_names{lc($externalDB)} || ($externalDB =~ s/_/ /g, $externalDB)  ;
        push (@{$links{$display_name}}, $link);         
    }
    return \%links;
}

=head2 disease_matches

 Arg[1]      : none
 Example     : $genedata->renderer->disease_matches
 Description : Renders the disease_matches in two_col_table format
 Return type : key-value pair, label and html

=cut

sub disease_matches{
    my $self = shift;
    my $genedata = $self->DataObj;
    my $omim_list = $genedata->get_disease_matches;
    if( ! ref($omim_list) or ! scalar(%$omim_list) ){ return }
    
    my $urls = $self->ExtURL;
    my $label = 'Disease Matches';
    my $html  = qq(
<b>This Ensembl entry corresponds to the following OMIM disease identifiers:</b><br />);

    my $table_tmpl = qq(
<TABLE class="hidden">%s
</TABLE>);  
    my $row_tmpl = qq(
<TR><TD>%s</TD><TD><SMALL>[%s] - [%s]</SMALL></TD></TR>);
    my $link_tmpl = qq(<A href="%s" target="%s">%s</A>);

#call table method here
   for my $description (sort keys %{$omim_list}){ 
       for my $omim (sort @{$omim_list->{$description}}){

       my $omim_href = sprintf($link_tmpl, 
                   $urls->get_url('OMIM',$omim), 
                   'new', 
                   "Omim ID: $omim" );
       
       my $diseaseview_href = sprintf($link_tmpl,
                      "diseaseview?omimid=$omim",
                      '',
                      "View disease info" );

       $html .= sprintf( $row_tmpl, $description, 
                 $omim_href, $diseaseview_href );
       $description="&nbsp;";
       }
   }  
   $html = sprintf( $table_tmpl, $html );
   return ($label, $html);
}

=head2 genesequence_link

 Arg[1]      : none
 Example     : 
 Description : Renders the genesequence link in two_col_table format
 Return type : key-value pair, label and html

=cut

sub genesequence_link {
    my $self = shift;
    my $geneid = $self->DataObj->stable_id();
    my $db = $self->param('db');
    my $label = 'Sequence Markup';
    my $html = qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/geneview?gene=$geneid&db=$db&_gene_sequence=1">View genomic sequence for this gene with exons highlighted</a>);
    return ($label, $html );
    #return ($label, $html);
}


=head2 exportview_link

 Arg[1]      : none
 Example     : $genedata->renderer->exportview_link
 Description : Renders the export_view link in two_col_table format
 Return type : key-value pair, label and html

=cut

sub exportview_link {
    my $self = shift;
    my $geneid = $self->DataObj->stable_id();
    my $label = 'Export Data';
    my $html = qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/exportview?type=feature&ftype=gene&id=$geneid">Export gene data in EMBL, GenBank or FASTA</a>);
    return ($label, $html);
}

=head2 genesnpview_link

 Arg[1]      : none
 Example     : $genedata->renderer->genesnpview_link
 Description : Renders the genesnpview_view link in two_col_table format
 Return type : key-value pair, label and html

=cut

sub genesnpview_link {
  my $self = shift;
  return  if $self->param('db') and $self->param('db') ne 'core';
  my $T = $self->DataObj->database('VARIATION');
  if( $self->DataObj->species_defs->databases->{ENSEMBL_VARIATION} ) {
    return (
      'SNP information' ,
      qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/genesnpview?db=@{[$self->param('db')]}&gene=@{[$self->DataObj->stable_id]}">View information about variations on this gene.</a>)
    );
  } else {
    return ;
  }
}

#----------------------------------------------------------------------

=head2 markup_options

  Arg [1]   : none
  Function  : Creates HTML form used to specify markup options
  Returntype: array
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub markup_options {
  my $self = shift;
  my $TABLE_TMPL = qq(
<form>
<input type="hidden" name="_gene_sequence" value="1">
<input type="hidden" name="gene" value="%s">
<input type="hidden" name="db" value="%s">
<table cellspacing=2 cellpadding=2 border=0>%s
</table>
<input type="submit" value="update">
</form>);

  my $ROW_TMPL = qq(
 <tr>
  <td>%s</td>
  <td>%s</td>
 </tr>);

  my @rows;

  my @f5form = $self->DataObj->get_slice->renderer->flank5_display_form;
  @f5form && push( @rows, sprintf( $ROW_TMPL, @f5form ) );
 
  my @f3form = $self->DataObj->get_slice->renderer->flank3_display_form;
  @f3form && push( @rows, sprintf( $ROW_TMPL, @f3form ) );

  my @eform  = $self->DataObj->get_slice->renderer->exon_display_form;
  my @eoform = $self->DataObj->get_slice->renderer->exon_ori_form;
  if( @eform ){
    my $form = $eform[1].$eoform[1];
    @eform && push( @rows, sprintf( $ROW_TMPL, $eform[0], $form) );
  }

  my @sform = $self->DataObj->get_slice->renderer->snp_display_form;
  @sform && push( @rows, sprintf( $ROW_TMPL, @sform ) );

  my @lform = $self->DataObj->get_slice->renderer->line_numbering_form;
  @lform && push( @rows, sprintf( $ROW_TMPL, @lform ) );

  my $geneid = $self->DataObj->stable_id;
  my $db = $self->param('db') || 'core';

  my $label = "Markup Options";
  my $data  = sprintf( $TABLE_TMPL, $geneid, $db, join( '', @rows ) );
  return( $label, $data );
}

#----------------------------------------------------------------------

=head2 markedup_geneseq

  Arg [1]   : none
  Function  : Creates marked-up FASTA-like string coeeesponding to gene slice
  Returntype: key-value pair, label and html
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub markedup_geneseq {
  my $self = shift;
  return( "Gene Sequence", $self->DataObj->get_slice->renderer->markedup_seq );
}




=head2 geneview_link

 Arg[1]      : none
 Example     : $genedata->renderer->geneview_link
 Description : Renders the geneview link in two_col_table format
 Return type : key-value pair, label and html

=cut

sub geneview_link {
    my $self = shift;
    return () if $self->param('db') and $self->param('db') ne 'core';
    return unless $self->DataObj->database('SNP');
    return (
      'Gene information' ,
      qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/geneview?db=@{[$self->param('db')]}&gene=@{[$self->DataObj->stable_id]}">View information about this gene.</a>)
    );
}

sub genesnpview_legend{
    my $self = shift;
    $self->print(qq(
<div align="center">
  <img src="/gfx/genesnpview-key.gif" height="160" width="800" border="0" alt="" /><br />&nbsp;
</div> ));

}


#----------------------------------------------------------------------

=head2 das_configurator

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub das_configurator {
  my $self = shift;

  my $gene = $self->param( "gene" );
  my $db      = $self->param( "db" );
  my $param_str = "gene=$gene&db=$db";

  my $conf_submit = qq(
<FORM name="dasConfigForm" action="dasconfview" method="post">
  <INPUT type="hidden" name="conf_script" value="geneview">
  <INPUT type="hidden" name="conf_script_params" value="$param_str">
  <INPUT type="submit" value="Manage Sources">
</FORM>);

  my $html_tmpl = qq(
<FORM name="dasForm" id="dasForm" method="get"> %s 
</FORM>);
  my $check_tmpl = qq(
  <INPUT type="checkbox" name="%s" value="1" %s
         onClick="javascript:document.dasForm.submit()"> %s 
  <INPUT type="hidden" name=":%s" value="0"> );
  my $hidden_tmpl = qq(
  <INPUT type="hidden" name="%s" value="%s"> );

  my $a_tmpl = qq(<A href="%s" target="new"> %s </A>);

  my $label = sprintf( $a_tmpl, "/Docs/gene_das.html",,"GeneDAS" );
  $label .= " Sources";

  my $hidden = '';
  foreach my $param( "gene", "db" ){
    $hidden .= sprintf( $hidden_tmpl, $param, $self->param($param) );
  }

  my @checks = ();
  my $data = $self->DataObj;
  my $das_attribute_data = $data->get_das_attributes( 
					"name", "authority", "label", "active" );
  foreach my $source( @$das_attribute_data ){
    my $name = $source->{name} ||
      ( warn "DAS source found with no name attribute" ) && next;
    my $source_desc = $source->{label};

    my $param_name = join( '!!', 'ENSEMBL_GENE_DAS_SOURCES', $name, 'active' );

    my $selected = $self->param($param_name) ? 'CHECKED' : '';
    my $href = $source->{authority} || undef();
    my $label = $href ? sprintf( $a_tmpl, $href, $name ) : $name;
    $label .= " ($source_desc)" if $source_desc;

    push @checks, sprintf( $check_tmpl, $param_name,
               $selected, $label, $param_name);
  }
  my $html = sprintf( $html_tmpl, join( "<BR>", @checks ).$hidden );
  $html .= "$conf_submit";

  return( $label, $html );
}

#----------------------------------------------------------------------

=head2 das_annotation

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub das_annotation {
  my $self = shift;
  my $data = $self->DataObj;

  my $table_tmpl = qq(
<TABLE class="hidden">%s
</TABLE>);
  my $head_row = qq(
 <TR><TD><I> Source </I></TD><TD><I> ID    </I></TD>
     <TD><I> Type   </I></TD><TD><I> Notes </I></TD></TR> );
  my $space_row = qq(
 <TR><TD colspan=4 height=5><IMG src="/gfx/blank.gif" height=5></TD></TR>);
  my $row_tmpl = qq(
 <TR><TD valign="top"><B> %s </B></TD><TD valign="top"> %s </TD>
     <TD valign="top">    %s     </TD><TD><SMALL> %s </SMALL></TD></TR> );
  my $link_tmpl = qq(<A href="%s" target="%s">%s</A>);

  my @table_data = ();
  my $das_attribute_data = $data->get_das_attributes
    ("name", "authority","active", "type" );

  my %FLabels = (
		 0 => 'Feature contained by gene', 
		 1 => 'Feature overlaps gene',
		 2 => 'Feature overlaps gene',
		 3 => 'Feature contains gene',
		 4 => 'Feature colocated with gene'
		);

  foreach my $source( @$das_attribute_data ){
    $source->{active} || next;
    my $source_nm = $source->{name};
    my $source_lab = $source_nm;
    if( my $ln = $source->{authority} ){
      $source_lab = sprintf( $link_tmpl, $ln, $source_nm, $source_nm );
    }
    my $label = "GeneDAS";
    push( @table_data, $label.": ". $source_lab,'' );

    my @rhs_rows;

    if ($source->{type} eq 'ensembl_gene_location') {
      my $slice = $data->get_gene_slice();
      my @features = $data->get_das_annotation_by_slice($source_nm, $slice);
      my $slice_length = $slice->end - $slice->start;
      my (%uhash) = (); # to filter out duplicates

      my (@filtered_features) = ();

      foreach my $feature (@features) {
#	print Dumper($feature), "<hr>";	

	next if ($feature->das_method_id eq 'Component');

	my $id = $feature->das_feature_id;
	if (defined($uhash{$id})) {
	  $uhash{$id}->{start} = $feature->das_start if ($uhash{$id}->{start} > $feature->das_start);
	  $uhash{$id}->{end} = $feature->das_end if ($uhash{$id}->{end} < $feature->das_end);
	  $uhash{$id}->{merged} = 1;
	} else {
	  $uhash{$id}->{start} = $feature->das_start;
	  $uhash{$id}->{end} = $feature->das_end;
	  $uhash{$id}->{type} = $feature->das_type;

	  my $segment = $feature->das_segment->ref;
	  if( my $href = $feature->das_link ){
	    $uhash{$id}->{label} = sprintf( $link_tmpl, $href, $segment, $feature->das_feature_label )
	  } else {
	    $uhash{$id}->{label} =  $feature->das_feature_label;
	  }

	  if( my $note = $feature->das_note ){
	    $note=~s|((\S+?):(http://\S+))|
	      <A href="$3" target="$segment">[$2]</A>|ig;
	    $note=~s|([^"])(http://\S+)([^"])|
	      $1<A href="$2" target="$segment">$2</A>$3|ig;
    
	    my $script = $self->param('script');
	    $note=~s|((\S+?):navigation://(\S+))|
	      <A href="$script?gene=$3" >[$2]</A>|ig;
	    $uhash{$id}->{note} = $note;
	  }

	}

      }

      foreach my $id ( sort { $uhash{$a}->{type} cmp $uhash{$b}->{type} || $a cmp $b || $uhash{$a}->{note} cmp $uhash{$b}->{note} } keys(%uhash )) {
# Build up the type of feature location	: see FLabels hash few lines above for location types
	my $ftype = 0;

	if ($uhash{$id}->{start} == $slice->start) {
	  if ($uhash{$id}->{end} == $slice_length) {
	    # special case - feature fits the gene exactly
	    $ftype = 4;
	  }
	} else {
	  if ($uhash{$id}->{start} < 0) {
	    # feature starts before gene starts
	    $ftype |= 2;
	  }
	  if ($uhash{$id}->{end} > $slice_length) {
	    # feature ends after gene ends
	    $ftype |= 1;
	  }
	}
	
#	print "$id>", $uhash{$id}->{start}, ":", $uhash{$id}->{end}, ": ", (defined($uhash{$id}->{merged})) ? "Merged " : " ", "[", $FLabels{$ftype}, "]<br>";
#	my $fnote = sprintf("%s%s [%d: %ld-%ld : %ld-%ld(%ld)]", (defined($uhash{$id}->{merged})) ? "Merged " : "", $FLabels{$ftype}, $ftype, $uhash{$id}->{start}, $uhash{$id}->{end}, $slice->start, $slice->end, $slice_length);
	my $fnote = sprintf("%s%s", (defined($uhash{$id}->{merged})) ? "Merged " : "", $FLabels{$ftype});
	push( @rhs_rows, sprintf( $row_tmpl, 
				  $uhash{$id}->{type} || '&nbsp;',
				  $uhash{$id}->{label}  || "&nbsp",
				  $fnote || "&nbsp",
				  $uhash{$id}->{note} || '&nbsp;' ) );
      }

      if( scalar( @rhs_rows ) == 0 ){
	push( @rhs_rows, "No annotation" );
      }

    } else {
      my @features = $data->get_das_annotation_by_name($source_nm,'global');
      if( ! scalar( @features ) ){
	push( @rhs_rows, "No annotation" );
      }

      foreach my $feature( sort{ 
	$a->das_type_id    cmp $b->das_type_id ||
	  $a->das_feature_id cmp $b->das_feature_id ||
	    $a->das_note       cmp $b->das_note
	  } @features ){
	
	my $segment = $feature->das_segment->ref;
	my $id = $feature->das_feature_id;
	if( my $href = $feature->das_link ){
	  $id = sprintf( $link_tmpl, $href, $segment, $id )
	}
	my $note;
	if( $note = $feature->das_note ){
	  $note=~s|((\S+?):(http://\S+))|
	    <A href="$3" target="$segment">[$2]</A>|ig;
	  $note=~s|([^"])(http://\S+)([^"])|
	    $1<A href="$2" target="$segment">$2</A>$3|ig;
	  
	  my $script = $self->param('script');
	  $note=~s|((\S+?):navigation://(\S+))|
	    <A href="$script?gene=$3" >[$2]</A>|ig;
	}

	push( @rhs_rows, sprintf( $row_tmpl, 
                #$feature->{-source} || '&nbsp;',
                $feature->das_type || '&nbsp;',
                $id                || '&nbsp;',
                $note              || '&nbsp;' ) );
      }
    }
    $table_data[-1] = sprintf( $table_tmpl, join($space_row, @rhs_rows ) );
  }

  return (@table_data);   
}


sub mouse_note{
    my $self = shift;
    my $gene =  $self->DataObj();
    my $label = qq(<font color="red">Assembly m32</font>);
    my $params = join ('&', (map {$_."=".$self->param($_)} $self->Input->param));
    my $link = qq(<a href="http://mouse30.ensembl.org/Mus_musculus/$ENV{'ENSEMBL_SCRIPT'}/?$params">here</a>);

    my $M32_GENES = $gene->species_defs->M32_GENES || {};
    my $html;
    if (exists $M32_GENES->{$gene->stable_id}){
        $html = qq (
        There are some long-range mapping issues between mouse assemblies
        30 and 32, and <b>this gene has been identified as problematic</b>. For
        more information see <a href="/Mus_musculus/m32analysis.html">our
        analysis</a>.  We would encourage you to assess the evidence for
        positioning if you are designing experiments sensitive to long
        range mapping.
       <br><br> 
        To view the corresponding feature on mouse NCBIm30 please click $link .
        );
    }
    else {
        $html = qq(
        
        There are some long range mapping issues between mouse assemblies
        30 and 32, but this gene has not been identified as problematic.
        For more information see <a
        href="/Mus_musculus/m32analysis.html">our analysis</a>, though we
        would encourage you to always assess the evidence for positioning
        if you are designing experiments sensitive to long range mapping.
       <br><br> 
        To view the corresponding feature on mouse NCBIm30 please click $link .
    );
    
    }   
    return ($label, $html);
}


1;  
