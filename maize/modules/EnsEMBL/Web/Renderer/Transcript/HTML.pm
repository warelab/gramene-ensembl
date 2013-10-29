package EnsEMBL::Web::Renderer::Transcript::HTML;

=head1 NAME

EnsEMBL::Web::Renderer::Transcript::HTML.pm 

=head1 SYNOPSIS

This object creates HTML to be output to the HTML output object

=head1 DESCRIPTION

    $transcript_renderer = $trans_data->renderer;
    $transcript_renderer->outputGenericTranscriptTable();        
        
 This object contains wrappers for common display 'groups' and also more
 granular calls that can be reused to create different page layouts/structure

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 CONTACT

Brian Gibbins - bg2@sanger.ac.uk

=cut

use strict;
use warnings;
no warnings "uninitialized";

use vars qw( @ISA );
use EnsEMBL::Web::Renderer::HTML;
use EnsEMBL::Web::Renderer::Transcript;
use Data::Dumper;

@ISA = qw( EnsEMBL::Web::Renderer::HTML EnsEMBL::Web::Renderer::Transcript );

=head2 outputGenericTranscriptTable

 Arg[1]      : none
 Example     : $transdata->renderer->outputGenericTranscriptTable
 Description : Wrapper for the table at top of transview
 Return type : none

=cut

sub outputGenericTranscriptTable {
    my $self = shift;
    $self->transcript_table_title;
    $self->Output->print_two_col_table(
        $self->transcript_display_ID,
        $self->stable_id,
        $self->other_view_links,
        $self->genomic_location,
        $self->description,
        $self->prediction_method,
        $self->similarity_matches,
        $self->GKB_links,
        $self->go_links,
        $self->interpro_links,
        $self->exportview_link,
    );
}

=head2 outputTransviewPageBottom

 Arg[1]      : none
 Example     : $transdata->renderer->outputTransviewPageBottom
 Description : Wrapper for the bottom of transview, including transcript sequence markup, transcript
                 neigbourhood and structure images
 Return type : none

=cut

sub outputTransviewPageBottom{
    my $self = shift;
    $self->Output->print_hidden_table(
            $self->print_fasta,
           [$self->transview_transcript_image, $self->transcript_neighbourhood_image, $self->markup_key_image],
    );
}

=head2 outputGeneTranscriptTable

 Arg[1]      : none
 Example     : $transdata->renderer->outputGeneTranscriptTable
 Description : Wrapper for the standard geneview transcript tables. This method is in the transcript renderer
                object because it uses transcript data objects to build this part of the display
 Return type : none

=cut
    
sub outputGeneTranscriptTable {
    my $self = shift;
    $self->Output->print_two_col_table(
        $self->transcript_summary,
        $self->similarity_matches,
        $self->GKB_links,
        $self->go_links,
        $self->interpro_links,
        $self->family_links,
        $self->geneview_transcript_image,
        $self->geneview_peptide_image,
    );
}

=head2 outputGenericExonTable

 Arg[1]      : none
 Example     : $transdata->renderer->outputGenericExonTable
 Description : Wrapper for the standard exonview table at the top of the page. This method is in the transcript renderer
                object because it uses transcript data objects to build this part of the display
 Return type : none

=cut

sub outputGenericExonTable {
    my $self = shift;
    $self->exon_table_title;
    $self->Output->print_two_col_table(
        $self->transcript_display_ID,
        $self->stable_id,
        $self->other_view_links_exon,
        $self->genomic_location,
        $self->description,
    );
}

#----

=head2 transcript_table_title

 Arg[1]      : none
 Example     : $transdata->renderer->transcript_table_title
 Description : prints the tile for the transcript table at top of transview
 Return type : string - html

=cut

# may move this onto the table call somehow
sub transcript_table_title{
    my $self = shift;
    my $type = $self->DataObj->gene_type;
    $self->Output->generic_table_title("<h3>$type Transcript Report</h3>" );
}

=head2 exon_table_title

 Arg[1]      : none
 Example     : $transdata->renderer->exon_table_title
 Description : prints the tile for the exon table at top of exonview
 Return type : string - html

=cut

# may move this onto the table call somehow
sub exon_table_title{
    my $self = shift;
    my $type = $self->DataObj->gene_type;
    $self->Output->generic_table_title("<h3>$type Exon Report</h3>" );
}

=head2 transcript_display_ID

 Arg[1]      : none
 Example     : $transdata->renderer->transcript_display_ID
 Description : renders the xref_display_id and type in two_col_table format
 Return type : key-value pair, label and html

=cut

sub transcript_display_ID{
    my $self = shift;
    my $trans = $self->DataObj();
    my $label = 'Transcript';
    return unless $trans->display_xref();    
    my ($display_name, $dbname) = $trans->display_xref();
    my $html = "<b>$display_name </b> <small> ($dbname ID)</small>";
    return ($label, $html);     
}

=head2 stable_id

 Arg[1]      : (optional) String
                Label
 Example     : $transdata->renderer->stable_id
 Description : renders the ensembl_stable_id in two_col_table format
 Return type : key-value pair, label and html

=cut

sub stable_id {
    my $self = shift;
    my $trans =  $self->DataObj();
    my $type = $trans->gene_type ;
    my $label = shift || "$type Transcript ID";
       $type = " <small>(". $trans->feature_type .")</small>" if ($type eq 'Vega');
    my $transid = "<b>".$trans->stable_id."</b> $type";
    return ($label, $transid);
}

=head2 vega_stable_id

 Arg[1]      : (optional) String
                Label
 Example     : $transdata->renderer->vega_stable_id
 Description : renders the stable_id for Vega in two_col_table format
 Return type : key-value pair, label and html

=cut

sub vega_stable_id {
    my $self = shift;
    my $label = shift || "Vega Transcript ID";
    my $html = "<b>".$self->DataObj->stable_id."</b>";
    return ($label, $html);
}

=head2 version

 Arg[1]      : (optional) String
               Label
 Example     : $transdata->renderer->version($string)
 Description : Renders the transcript version in two_col_table format
 Return type : key-value pair, label and html

=cut

sub version {
    my $self = shift;
    my $trans =  $self->DataObj();
    my $label = shift || "Version";
    my $html = $trans->version;
    return ($label, $html);
}

=head2 transcript_class

 Arg[1]      : (optional) String
               Label
 Example     : $transdata->renderer->transcript_class("Class")
 Description : Renders the transcript class in two_col_table format
 Return type : key-value pair, label and html

=cut

sub transcript_class {
    my $self = shift;
    my $label = shift || "Class";
    my $trans = $self->DataObj;
    my $trans_class = $trans->get_transcript_class;
    my $species = $trans->Input->species;
    my $html = qq($trans_class [<a href="javascript:X=hw\('$species', 'Vega_transcript_classification', ''\)">Definition</a>]);
    return ($label, $html);
}

=head2 genomic_location

 Arg[1]      : none
 Example     : $transdata->renderer->genomic_location
 Description : renders the transcript location in two_col_table format
 Return type : key-value pair, label and html

=cut

sub genomic_location{
    my $self = shift;   
    my $trans =  $self->DataObj;
    my $transid = $trans->stable_id;
    my ($chr, $start, $end, $contig, $contig_start) = $trans->get_genomic_location();
    my $label = 'Genomic Location';     
    my $html = '';
    
    if (!$start && !$end && !$chr){
    $html .=  qq(<i>This transcript cannot be located on the current assembly</i>\n);
    }else{
    my $mbase = $trans->bp_to_nearest_unit($start,1);
    my $chr_display = $trans->is_valid_chr($chr) ? "chromosome $chr" : "chromosome fragment $chr" ;
    my $contigview_link = "contigview?contig=$contig&fpos_start=$contig_start&fpos_end=".($contig_start + 1)."&fpos_context=200000&highlight=$transid";
            
    $html .= qq(
            <b>View transcript in genomic location: </b>
            <a href="/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$start&vc_end=$end&highlight=$transid">$start - $end bp ($mbase)</a> 
                on $chr_display<br />
            <b>This transcript is located in sequence: </b>
            <a href="$contigview_link">$contig</a> );
    }   
    $html .=  qq(<br />);   
    return ($label, $html)
}

=head2 description

 Arg[1]      : none
 Example     : $transdata->renderer->description
 Description : renders the GENE description in two_col_table format (no description on transcript)
 Return type : key-value pair, label and html

=cut

sub description {
    my $self = shift ;
    my $urls = $self->ExtURL;
    my $description = $self->DataObj->trans_description();    
       $description =~ s/EC\s+([-*\d]+\.[-*\d]+\.[-*\d]+\.[-*\d]+)/$self->EC_URL($1)/e;
       $description =~ s/\[\w+:(\w+)\;\w+:(\w+)\]//g;
    my ($edb, $acc) = ($1, $2);
    my $url         = "[<a href='".$urls->get_url($edb, $acc)."'>Source: $1 ($2)</a>]" unless !$acc;
    
    my $label = 'Description';
    my $html = qq($description <small>$url</small>&nbsp;<br />);       
    return ($label, $html);
}

=head2 remarks

 Arg[1]      : (optional) String
               Label
 Example     : $transcriptdata->renderer->remarks
 Description : returns annotation remarks for use in two_col_table
 Return type : key-value pair of label/html

=cut

sub remarks {
    my $self = shift;
    my $label = shift || "Remarks";
    my $transcript =  $self->DataObj();
    my $remarks = $transcript->get_remarks;
    my $html;
    if (@{$remarks}) {
        $html = join ('<br />', @{$remarks});
    } else {
        $html = "No remarks";
    }
    return ($label, $html);
}

=head2 prediction_method

 Arg[1]      : none
 Example     : $transdata->renderer->prediction_method
 Description : renders the prediction method for a transcript in two_col_table format 
 Return type : key-value pair, label and html

=cut

sub prediction_method {
    my $self = shift ;
    my $transdata = $self->DataObj;
    my $db = $transdata->feature_type ; 
    my $label = ( $db eq 'vega' ? 'Curation' : 'Prediction' ).' Method';
    my $trans_prediction_text = $transdata->get_prediction_method ."&nbsp;";        
    return ($label, $trans_prediction_text );
}

=head2 other_view_links

 Arg[1]      : (optional) String
                Label
 Example     : $pepdata->renderer->other_view_links
 Description : renders the ensembl_links to other views in two_col_table format
 Return type : key-value pair, label and html

=cut

sub other_view_links{
    my $self  = shift;
    my $label = shift;

    my $web_tran = $self->DataObj; 
    $label ||= $web_tran->gene_type . " Transcript";

    my $ens_tran = $web_tran->transcript;
    my $ens_gene = $web_tran->gene;
    my $ens_pep  = $ens_tran->translation;

    my $stable_id  = $ens_tran->stable_id;
    my $gene_id    = $ens_gene ? $ens_gene->stable_id : undef;

    my $db = $web_tran->get_db;
    my $display_name =( $ens_tran->display_xref() ? 
                        $ens_tran->display_xref()->display_id : $stable_id );

    my $num_exons = scalar(@{$ens_tran->get_all_Exons});
    my $trans_length = $ens_tran->seq->length . "bp";

    my $html = qq(
<table class="hidden">
 <tr> );

    $html .= qq(
  <td><b> Exons: </b> $num_exons </td> );
     $html .= qq(
  <td><b> Transcript length: </b> $trans_length </td>);

    if( $ens_pep ){
      my $peptide_length = $ens_pep->length;
      $html .= qq( 
<td><b>Translation length:</b> $peptide_length residues </td>);
    } else {
      $html .= qq( 
<td><b>Transcript does not translate</b> $trans_length </td>);
    }

    $html .= qq(
</tr></table> );

    if( $ens_gene ){
      $html .= qq(
<table class="hidden">
 <tr>
  <td>This transcript is a product of gene: <a href="/$ENV{'ENSEMBL_SPECIES'}/geneview?gene=$gene_id&db=$db">$gene_id</a></td>
 </tr>
</table> );
    }

    if ($self->param('show_vega_evidence_link')) {
    $html .= qq(
<table class="hidden">
 <tr>
  <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$stable_id&db=$db">Exon information & supporting evidence</a>]</td>  );

} else {
    $html .= qq(
<table class="hidden">
 <tr>
  <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$stable_id&db=$db">Exon information</a>]</td>  );
}

    if( $ens_pep ){
      $html .= qq(
  <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/protview?transcript=$stable_id&db=$db">Protein information</a>]</td> );
    } else {
      $html .=  qq(
  <td>[No translation]</td> );
    }
    $html .=  qq(
 </tr>
</table> );
    return ($label, $html);
}

=head2 other_view_links_exon

 Arg[1]      : (optional) String
                Label
 Example     : $pepdata->renderer->other_view_links_exon
 Description : renders the ensembl_links to other views for exonview in two_col_table format
 Return type : key-value pair, label and html

=cut

sub other_view_links_exon{
    my $self = shift;
    my $label = shift;

    my $web_tran =  $self->DataObj();
    my $db = $web_tran->get_db;
    $label ||= $web_tran->gene_type . " Transcript";

    my $web_pep  = $web_tran->translation;

    my $ens_tran = $web_tran->transcript;
    my $ens_gene = $web_tran->gene;
    my $ens_pep  = $ens_tran->translation;

    my $tran_id = $ens_tran->stable_id;
    my $gene_id = $ens_gene ? $ens_gene->stable_id : '';
    my $pep_id  = $ens_pep  ? $ens_pep->stable_id  : '';

    my $html = '';

    if( $ens_gene ){
      $html .= qq(
This transcript is a product of Ensembl gene <a href="/$ENV{'ENSEMBL_SPECIES'}/geneview?gene=$gene_id&db=$db">$gene_id</a> <br /> );
    }

    $html .= qq(
[ <a href="/$ENV{'ENSEMBL_SPECIES'}/transview?transcript=$tran_id&db=$db">Transcript Information</a> ] );

    if( $web_tran->get_supporting_evidence ){
       # Have evidence, therefore need link
      $html .= qq(
[ <a href="#evidence">Supporting Evidence</a> ]  );
    }

    if( $ens_pep ){
      $html .= qq(
[ <a href="/$ENV{'ENSEMBL_SPECIES'}/protview?transcript=$tran_id&db=$db">Peptide Information</a> ] );
    } else{
      $html .= qq( [ No Translation ] );
    }
    return ($label, $html);
}

=head2 author

 Arg[1]      : (optional) String
               Label
 Example     : $transcriptdata->renderer->author
 Description : returns the annotation author for use in two_col_table
 Return type : key-value pair of label/html

=cut

sub author {
    my $self = shift;
    my $label = shift || "Author";
    my $transcript =  $self->DataObj();
    my $author = $transcript->get_author_name;
    my $html;
    if ($author) {
        $html .= "This locus was annotated by " . $author . " ";
        $html .= $self->email_URL($transcript->get_author_email);
    } else {
        $html = "unknown";
    }
    return ($label, $html);
}

=head2 transcript_summary

 Arg[1]      : none
 Example     : $transdata->renderer->transcript_summary
 Description : renders the transcript_summary for a transcript in
               two_col_table format for geneview transcript tables
 Return type : key-value pair, label and html

=cut

sub transcript_summary {
    my $self = shift;
    my $transdata = $self->DataObj;
    my $trans = $transdata->transcript;
    my $db = $transdata->get_db;
    my $stable_id = $trans->stable_id;
    my $display_name = $trans->display_xref() ? $trans->display_xref()->display_id : $stable_id ;
    my $label = qq(<a name="$stable_id">$display_name</a>);     
    my $html = qq(<table class="hidden"><tr>);
    
    $html .= qq(<td><b>Stable ID:</b> $stable_id </td>) if $trans->display_xref();  
    $html .= qq(<td><b>Exons:</b> ). scalar(@{$trans->get_all_Exons}) .qq(</td> 
                <td><b>Transcript length:</b> ). $trans->seq->length .qq( bp </td>);
    my $peptide    = $transdata->translation();
    my $protein    = $peptide->translation;
    my $peptide_id = $peptide->stable_id;
    if ($protein){
        my $peptide_length = $protein->length;
        $html .= qq( <td><b>Translation length:</b> $peptide_length residues </td>);
    }else {
        $html .= qq( <td><b>Transcript does not translate</b> </td>);
    }   
    
    $html .= qq(</tr></table>); 
    $html .= qq(<table class="hidden"><tr>
        <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/transview?transcript=$stable_id&db=$db">Transcript information</a>]</td>  
        <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$stable_id&db=$db">Exon information</a>]</td>  );
    if( $protein ){
      if( $peptide_id ){
        $html .= qq(
        <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/protview?peptide=$peptide_id&db=$db">Protein information</a>]</td></tr></table>);
      }
      else{
        $html .= qq(
        <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/protview?transcript=$stable_id&db=$db">Protein information</a>]</td></tr></table>);
        
      }
    }else{
        $html .=  qq(<td>[No translation]</td></tr></table> );
    }
    return ($label, $html);
}

=head2 vega_transcript_summary

 Arg[1]      : none
 Example     : $transdata->renderer->vega_transcript_summary
 Description : renders the Vega transcript summary for a transcript in
               two_col_table format for geneview transcript tables
 Return type : key-value pair, label and html

=cut

sub vega_transcript_summary {
    my $self = shift;
    my $transdata = $self->DataObj;
    my $trans = $transdata->transcript;
    my $version = $transdata->version;
    my $trans_class = $transdata->get_transcript_class;
    my $db = $transdata->get_db;
    my $stable_id = $trans->stable_id;
    my $display_name = $trans->display_xref ? $trans->display_xref->display_id : $stable_id ;
    my $label = qq(<a name="$stable_id">$display_name</a>);     
    my $html = qq(<table class="hidden"><tr>);
    $html .= qq(<td><b>Stable ID:</b> $stable_id </td>) if $trans->display_xref();  
    $html .= qq(<td><b>Version:</b> $version</td>);
    $html .= qq(<td><b>Class:</b> $trans_class</td>);
    $html .= qq(</tr></table><table class="hidden"><tr>);
    $html .= qq(<td><b>Exons:</b> ). scalar(@{$trans->get_all_Exons}) .qq(</td> 
                <td><b>Transcript length:</b> ). $trans->seq->length .qq( bp </td>);
    my $peptide    = $transdata->translation();
    my $protein    = $peptide->translation;
    my $peptide_id = $peptide->stable_id;
    if ($protein){
        my $peptide_length = $protein->length;
        $html .= qq( <td><b>Translation length:</b> $peptide_length residues </td>);
    }else {
        $html .= qq( <td><b>Transcript does not translate</b> </td>);
    }
    if ($self->param('show_vega_evidence_link')) {   
	$html .= qq(</tr></table><table class="hidden"><tr>
		    <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/transview?transcript=$stable_id&db=$db">Transcript information</a>]</td>  
		    <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$stable_id&db=$db">Exon information & supporting evidence</a>]</td>  );
    } else {
	$html .= qq(</tr></table><table class="hidden"><tr>
		    <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/transview?transcript=$stable_id&db=$db">Transcript information</a>]</td> );
    }
    
    
    if( $protein ){
      if( $peptide_id ){
        $html .= qq(
        <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/protview?peptide=$peptide_id&db=$db">Protein information</a>]</td></tr></table>);
      }
      else{
        $html .= qq(
        <td>[<a href="/$ENV{'ENSEMBL_SPECIES'}/protview?transcript=$stable_id&db=$db">Protein information</a>]</td></tr></table>);
        
      }
    }else{
        $html .=  qq(<td>[No translation]</td></tr></table> );
    }
    return ($label, $html);
}

=head2 transview_transcript_image

 Arg[1]      : none
 Example     : $transdata->renderer->transview_transcript_image
 Description : renders the transcript structure image for transview
 Return type : String - HTML 

=cut

sub transview_transcript_image{
    my $self = shift;
    my $transData = $self->DataObj;
    my $db = $transData->get_db ;
    my $transid = $transData->stable_id;    
    my $image = $self->transcript_image('transview');
    my $html = qq(
        <div align="center"><b>Transcript Structure</b><br /><br />
        <a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$transid&db=$db" title="View Exon Information">
        $image</div></p>
        </a>);  
    return ($html);
}

=head2 transview_transcript_structure

 Arg[1]      : (optional) String
               Label
 Example     : $transdata->renderer->transview_transcript_structure
 Description : renders the transcript structure image for transview. Similar
               to transview_transcript_image, but for use in two_col_table
 Return type : String - HTML 

=cut

sub transview_transcript_structure {
    my $self = shift;
    my $label = shift || 'Transcript Structure';
    my $transData = $self->DataObj;
    my $db = $transData->get_db ;
    my $transid = $transData->stable_id;    
    my $image = $self->transcript_image('transview');
    my $html = qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$transid&db=$db" title="View Exon Information">$image</a>); 
    return ($label, $html);
}

=head2 geneview_transcript_image

 Arg[1]      : none
 Example     : $transdata->renderer->geneview_transcript_image
 Description : renders the transcript structure image for geneview in two_col_table format
 Return type : Key /value pair - label and HTML 

=cut

sub geneview_transcript_image{
    my $self = shift;
    my $label = 'Transcript Structure' ;    
    my $transData = $self->DataObj;
    my $transid = $transData->stable_id;
    my $db = $transData->get_db ;
    my $image = $self->transcript_image('geneview');
    my $html = qq(<br><a href="/$ENV{'ENSEMBL_SPECIES'}/transview?transcript=$transid&db=$db" title="View Transcript Information">
                    <div align="center">$image</div></p></a>);
    return ($label, $html);
}

=head2 transcript_image

 Arg[1]      : String
                name of webuser-config to use (geneview, altsplice, etc)
 Example     : $transdata->renderer->transcript_image
 Description : renders the transcript structure image only
 Return type : String - Full file name for image

=cut

sub transcript_image{
    my $self = shift;
    my $script = shift;    
    my $transData = $self->DataObj;
    my $dc  = $transData->render_transcript_structure($script);
    return unless $dc;
    my $filename =    $self->render_dc( "transimg", $dc );  
    return qq(<img border="0" src="$filename">);   
}

=head2 geneview_peptide_image

 Arg[1]      : none
 Example     : $transdata->renderer->geneview_peptide_image
 Description : wrapper to print peptide image in two_col_table format
 Return type : Key / value pair - label and HTML

=cut

sub geneview_peptide_image{
    my $self = shift;   
    my $label = 'Protein Features ';
    my $html =  $self->peptide_image;
    return unless $html;
    return ($label, $html);
}

=head2 peptide_image

 Arg[1]      : none
 Example     : $transdata->renderer->peptide_image
 Description : renders peptide image for geneview (transcript tables in geneview)
 Return type : String - HTML for image and link

=cut

sub peptide_image{
    my $self = shift;   
    my $pepdata = $self->DataObj->translation;
    my $peptideid = $pepdata->stable_id;
    my $db = $pepdata->get_db ;
    my $dc = $pepdata->render_peptide_image( 'protview', 1 );
    return unless $dc;
    my ($filename, $map) = $self->render_image_imagemap('protimg', $dc);
    return qq(
        <br><div align="center">
        <a href="/$ENV{'ENSEMBL_SPECIES'}/protview?peptide=$peptideid&db=$db" title="View Peptide Information">
        <img border="0" src="$filename" usemap="$filename"></a><br/></div><br/>) if $filename;
}

=head2 similarity_matches

 Arg[1]      : (optional) String
               Label
 Example     : $transdata->renderer->similarity_matches
 Description : Renders similarity matches for transcript in two_col_table format
 Return type : Key / value pair - label and HTML

=cut

sub similarity_matches{
    my $self = shift;
    my $label = shift || 'Similarity Matches';
    my $trans = $self->DataObj->transcript;     
    my $data = $self->DataObj();
    return $trans->{'similarity_links'} if $trans->{'similarity_links'};
    my @similarity_links = @{$data->get_similarity_hash($trans)};    
    return unless (@similarity_links);

    my $db = $data->get_db();
    my $entry = $data->gene_type || 'Ensembl';
    #if ($db eq 'vega' or $data->species_defs->ENSEMBL_SITETYPE eq 'Vega') {
    #    $entry = 'Vega curated';
    #} else {
    #    $entry = 'Ensembl';
    #}

    # _sort_similarity_links formats the links
    my %links = %{$self->_sort_similarity_links(@similarity_links) ||{}};       
    return unless %links;

# add table call here
    my $html = qq(
        <b>This $entry entry corresponds to the following database identifiers:</b>
        <br /><table class="hidden">);
    foreach my $key (sort keys %links){
        if (scalar (@{$links{$key}}) > 0){
            my @sorted_links = sort @{$links{$key}};
            $html .= qq(<tr valign="top"><td nowrap="nowrap"><b>$key:</b></td><td>);
          
            if( $sorted_links[0] =~ /<br/i ){
                $html .= join(' ', @sorted_links );
            } else { # Need a BR each 5 entries
                $html .= qq(<table class="hidden"><tr>);
                my @sorted_lines;
                for( my $i=0; $i<@sorted_links; $i++ ){
                    my $line_num = int($i/4);
                    if( ref(  $sorted_lines[$line_num] ) ne 'ARRAY' ){$sorted_lines[$line_num] = [];}
                    push( @{$sorted_lines[$line_num]}, "<td>".$sorted_links[$i]."</td>" );
                }
                $html .= join( '</tr><tr>', map{ join( ' ', @$_ ) } @sorted_lines );
                $html .= qq(</tr></table>);
            }
            #if ($key eq "GO" && !$data->database('go')){
            #    $html .= qq(<small>GO mapping is inherited from swissprot/sptrembl</small>);
            #}
            $html .= qq(</td></tr>);
        }
    }   
    $html .= qq(</table>); 
    return ($label , $html);
}

=head2 _sort_similarity_links

 Arg[1]      : none
 Example     : $transdata->renderer->_sort_similarity_links
 Description : sorts the similarity matches
 Return type : hashref of similarity matches

=cut

sub _sort_similarity_links{
    my $self = shift;
    my $trans = $self->DataObj->transcript;
    my @similarity_links = @_;
    my $data = $self->DataObj();
    my $database = $data->database;
    my $db = $data->get_db() ;
    my $urls = $self->ExtURL;
    my %links ;
    my $ALIGN_LINK = qq( [<a href="/$ENV{'ENSEMBL_SPECIES'}/alignview?transcript=%s&sequence=%s&db=%s" class = "small" target="palignview">align</a>] );
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
            'ens_hs_transcript'	    => 'Ensembl Human Transcript',
            'ens_hs_translation'    => 'Ensembl Human Translation',
        );
                       
    foreach my $type (sort @similarity_links) { 
        my $link = "";
        my $join_links = 0;
        my $externalDB = $type->database();
        my $display_id = $type->display_id();
        my $primary_id = $type->primary_id();

        # remove all orthologs  
        next if ($type->status() eq 'ORTH');
    
        # Ditch celera genes from FlyBase
        next if ($externalDB =~ /^flybase/i && $display_id =~ /^CG/ );

        # remove internal links to self and transcripts
        next if $externalDB eq "Vega_gene";
        next if $externalDB eq "Vega_transcript";
        next if $externalDB eq "Vega_translation";

        if( $externalDB eq "GO" ){ 
            #&& $data->database('go')){
            push @{$trans->{'go_links'}}, $display_id;
            next;   
        } 
	elsif ($externalDB eq "GKB") {
            my ($key, $primary_id) = split ':', $display_id;
            push @{$trans->{'GKB_links'}{$key}} , $type ;
            next;
        } 
	elsif ($externalDB eq "REFSEQ") { 
	    # strip off version
	    $display_id =~ s/(.*)\.\d+$/$1/o;
	}  
        elsif ($externalDB eq "protein_id") { 
	    # Can't link to srs if there is an Version - so strip it off
	    $primary_id =~ s/(.*)\.\d+$/$1/o;
	}
         
        # Build external links
        if ($urls and $urls->is_linked($externalDB)) {
            $link = '<a href="'.$urls->get_url($externalDB, $primary_id).'">'. $display_id. '</a>';
            if ( uc( $externalDB ) eq "REFSEQ" and $display_id =~ /^NP/) {
                $link = '<a href="'.$urls->get_url('REFSEQPROTEIN',$primary_id).'">'. $display_id. '</a>';
            } elsif ($externalDB eq "HUGO") {
                $link = '<a href="' .$urls->get_url('GENECARD',$display_id) .'">Search GeneCards for '. $display_id. '</A>';
            } elsif ($externalDB eq "MarkerSymbol") { # hack for mouse MGI IDs
                $link = '<a href="' .$urls->get_url('MARKERSYMBOL',$primary_id) .'">'."$display_id ($primary_id)".'</A>';
            } 
            if( $type->isa('Bio::EnsEMBL::IdentityXref') ) {
                $link .=' <small> [Target %id: '.$type->target_identity().'; Query %id: '.$type->query_identity().']</small>';            
                $join_links = 1;    
            }
            if (( $data->species_defs->ENSEMBL_PFETCH_SERVER ) && 
             ( $externalDB =~/^(SWISS|SPTREMBL|LocusLink|protein_id|RefSeq|EMBL|Gene-name|Uniprot)/i ) ) {  
                my $seq_arg = $display_id;
                $seq_arg = "LL_$seq_arg" if $externalDB eq "LocusLink";
                $link .= sprintf( $ALIGN_LINK,
                $trans->stable_id,
                $seq_arg,
                $db );
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
        # override for Affys - we don't want to have to configure each type, and
        # this is an internal link anyway.
        if ($externalDB =~ /^AFFY_/i) {
                $link = '<a href="' .$urls->get_url('AFFY_FASTAVIEW', $display_id) .'">'. $display_id. '</A>';
        }
        my $display_name = $nice_names{lc($externalDB)} || ($externalDB =~ s/_/ /g, $externalDB)  ;
        push (@{$links{$display_name}}, $link);         
    }
    $trans->{'similarity_links'} = \%links ;
    return $trans->{'similarity_links'};
}

=head2 go_links

 Arg[1]      : none
 Example     : $transdata->renderer->go_links
 Description : renders the go links in two_col_table if the go
               database exists
               This is copy of the Translation::HTML::go_links
               method, but there is no way to get hold of a
               Translation::HTML object from here. Site-wide updates need
               making in both places.
 Return type : Key / value pair - label and HTML

=cut

sub go_links { 
    my $self = shift;
    my $trans = $self->DataObj->transcript;
    my $data = $self->DataObj();   

    my $databases = $data->DBConnection;    
    my $goview = $data->database('go') ? 1 : 0;

    unless ($trans->{'go_links'}){
      $self->similarity_matches();
      return unless $trans->{'go_links'};
    }
    my $go_hash = $data->get_go_list();
    my $GOIDURL  = "/$ENV{'ENSEMBL_SPECIES'}/goview?acc=";
    my $QUERYURL = "/$ENV{'ENSEMBL_SPECIES'}/goview?depth=2&query=";
    my $URLS     = $self->ExtURL;

    return unless ( $go_hash);
    my $label = 'GO';
    my $html =  qq( <b>The following GO terms have been mapped to this entry via <a href="http://www.ebi.ac.uk/interpro/README1.html">InterProScan</a>:</b><br />);

    foreach my $go (sort keys %{$go_hash}){
      my @go_data = @{$go_hash->{$go}||[]};
      my( $evidence, $description ) = @go_data;
      my $link_name = $description;
      $link_name =~ s/ /\+/g;

      my $goidurl  = "$GOIDURL$go";
      my $queryurl = "$QUERYURL$link_name";
      unless( $goview ){
        $goidurl  = $URLS->get_url('GO',$go);
        $queryurl = $URLS->get_url('GOTERMNAME', $link_name);
      }
      $html .= qq(<A HREF="$goidurl">$go</A> [<A HREF="$queryurl">$description</A>] <tt>$evidence</tt> <br />\n);           
    } 
    return ($label, $html);     
}


=head2 GKB_links

 Arg[1]      : none
 Example     : $transdata->renderer->GKB_links
 Description : renders the GKB links in it's own two_col_table 
 Return type : Key / value pair - label and HTML

=cut

sub GKB_links { 
    my $self = shift;
    my $trans = $self->DataObj->transcript;
    my $data = $self->DataObj();   

    unless ($trans->{'GKB_links'}){
      $self->similarity_matches();
    }
    my $GKB_hash = $trans->{'GKB_links'};
    return unless ( $GKB_hash);

    my $label = 'Genome KnowledgeBase';
    my $html =  qq( <b>The following identifiers have been mapped to this entry via Genome KnowledgeBase:</b><br />);

    my $urls = $self->ExtURL;
    $html .= qq(<table class="hidden">);
    foreach my $db (sort keys %{$GKB_hash}){
        $html .= qq(<tr><th>$db</th><td><table class="hidden">);
        foreach my $GKB (@{$GKB_hash->{$db}}){
            my $primary_id = $GKB->primary_id;
            my ($t, $display_id) = split ':', $primary_id ;
            my $description = $GKB->description;
            $html .= '<tr><td><a href="' .$urls->get_url('GKB', $primary_id) .'">'. $display_id. '</A>&nbsp;</td>
                <td>'.$description.'</td>
            </tr>';
        }
        $html .= qq(</table></td></tr>)
    }
    $html .= qq(</table>);
    return ($label, $html);     
}

=head2 interpro_links

 Arg[1]      : none
 Example     : $transdata->renderer->interpro_links
 Description : renders the Interpro links two_col_table format
               This is copy of the Translation::HTML::interpro_links
               method, but there is no way to get hold of a
               Translation::HTML object from here. Site-wide updates need
               making in both places.
 Return type : Key / value pair - label and HTML

=cut

sub interpro_links {
    my $self = shift;
    my $trans = $self->DataObj->transcript;
    my $pepdata = $self->DataObj->translation;
    my $interpro_hash = $pepdata->get_interpro_links($trans);
    my $urls = $self->ExtURL;
    my $site_type = $self->DataObj->species_defs->ENSEMBL_SITETYPE;
    return unless (%$interpro_hash);

    my $label = 'InterPro';

# add table call here
    my $html = qq(<table class="hidden">);
    for my $accession (keys %$interpro_hash){
        my $interpro_url = $urls->get_url('INTERPRO',$accession);
        my $desc = $interpro_hash->{$accession};
        $html .= qq(<tr>
                <td><A HREF="$interpro_url">$accession</A> </td><td> $desc 
             - [<A HREF="domainview?domainentry=$accession">View other $site_type genes with this domain</A>]</td></tr>);
    }
    $html .= qq( </table> );
    return ($label, $html);
}

=head2 family_links

 Arg[1]      : none
 Example     : $transdata->renderer->family_links
 Description : renders the family links two_col_table format
 Return type : Key / value pair - label and HTML

=cut

sub family_links {
    my $self = shift;
    my $pepdata = $self->DataObj->translation;
    my $families = $pepdata->get_family_links($pepdata);    
    return unless %$families;
   
    my $label = 'Protein Family';
    my $html;
    foreach my $family_id (keys %$families) {
        my $family_url = "/$ENV{'ENSEMBL_SPECIES'}/familyview?family=$family_id";
        my $family_count = $families->{$family_id}{'count'};
        my $family_desc = $families->{$family_id}{'description'};
        $html .= qq(
           <a href="$family_url">$family_id</A> : $family_desc<br />  
            This cluster contains $family_count Ensembl gene member(s)<br />);
    }
    return ($label, $html);
}


=head2 exportview_link

 Arg[1]      : none
 Example     : $transdata->renderer->exportview_link
 Description : renders the exportview link two_col_table format
 Return type : Key / value pair - label and HTML

=cut

sub exportview_link {
    my $self = shift;
    my $transid = $self->DataObj->stable_id();
    my $label = 'Export Data';
    my $html = qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/exportview?type=feature&ftype=transcript&id=$transid">Export
                    transcript data in EMBL, GenBank or FASTA</a>);
    return ($label, $html);
}


=head2 transcript_neighbourhood_image

 Arg[1]      : none
 Example     : $transdata->renderer->transcript_neighbourhood_image
 Description : renders the transcript neigbourhood image for transview
 Return type : String - HTML

=cut

sub transcript_neighbourhood_image {
    my $self = shift;
    my $dc = $self->DataObj->configure_neighbourhood_image();
    my ($image, $imagemap) = $self->render_image_imagemap( "transneigh", $dc );
    my $html = qq(<br /><div align='center'><b>Transcript Neighbourhood</b><br /><br /><img border="0" src="$image" usemap="#geneview_transnhood"><br><br>
                <map name="geneview_transnhood">$imagemap</map></div>);
    return $html;
}

=head2 transcript_neighbourhood

 Arg[1]      : (optional) String
               Label
 Example     : $transdata->renderer->transcript_neighbourhood
 Description : renders the transcript neigbourhood image for transview. Similar
               to transview_neighbourhood_image, but for use in two_col_table
 Return type : String - HTML

=cut

sub transcript_neighbourhood {
    my $self = shift;
    my $label = shift || 'Transcript Neighbourhood';
    my $dc = $self->DataObj->configure_neighbourhood_image();
    my ($image, $imagemap) = $self->render_image_imagemap( "transneigh", $dc );
    my $html = qq(<img border="0" src="$image" usemap="#geneview_transnhood"><map name="geneview_transnhood">$imagemap</map>);
    return ($label, $html);
}

=head2 print_fasta

 Arg[1]      : none
 Example     : $transdata->renderer->print_fasta
 Description : renders the transcript cdna fasta markup for transview (includes <form>)
 Return type : String - HTML

=cut

sub print_fasta {
  my $self = shift;
  my ($sel_number, $sel_no) =  ('','');
  my ( $sel_snps, $sel_codons, $sel_peptide, $sel_plain) = ('','','','');
  my $Data = $self->DataObj;
  my $db = $Data->get_db() ;
  my $stable_id = $Data->stable_id;
  my $number = $Data->param('number'); 
  my $show = $Data->param('show') || 'plain';
  
  if($show eq'snps') { $sel_snps = ' selected'; } 
  elsif($show eq 'codons') {$sel_codons=' selected'; } 
  elsif($show eq 'peptide') { $sel_peptide =' selected'; } 
  else{ $sel_plain =' selected'; }  
  
  my $fasta ;
  if ($show eq 'plain'){
    $fasta = $Data->get_trans_seq ;
    $fasta =~ s/([acgtn\*]+)/'<font color="blue">'.uc($1).'<\/font>'/eg;
  } else {
    $fasta = $self->markup_transcript_fasta( $Data->get_markedup_trans_seq );
  } 
  
  if($number eq 'on') { $sel_number = ' selected'; } else {$sel_no=' selected'; }
    
  my $SNP_LINE = "";
  if (exists($Data->species_defs->databases->{'ENSEMBL_VARIATION'}) ||
      exists($Data->species_defs->databases->{'ENSEMBL_GLOVAR'})) {
    $SNP_LINE = "<option value='snps' $sel_snps>Codons/peptide/SNPs</option>";
  }
  
  my $html = qq(
    <b>Transcript cDNA Sequence</b>
        <form action="/$ENV{'ENSEMBL_SPECIES'}/$ENV{'ENSEMBL_SCRIPT'}" id="ddForm" name="ddForm" method="get">
        <input type="hidden" name="transcript" value="$stable_id" />
        <input type="hidden" name="db" value="$db" />
        <select name="show" onChange="document.ddForm.submit()">
          <option value="plain" $sel_plain>No markup</option>
          <option value="codons" $sel_codons>Codons highlighted</option>
          <option value="peptide" $sel_peptide>Codons/peptide sequence</option>
          $SNP_LINE
        </select>
        <select name="number" onChange="document.ddForm.submit()">
          <option value="on" $sel_number>Number bps.</option>
          <option value="off" $sel_no>No numbers</option>
        </select>
    </form>
        <pre class="snpmarkup">$fasta</pre>
        <blockquote><b>[<i>Exons are shown in alternating blue/black</i>] &nbsp;- [<a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$stable_id&db=$db">View exon information</a>] &nbsp;</blockquote>
  );
}

=head2 markup_transcript_fasta

 Arg[1]      : Int - coding_start
 Arg[2]      : Int - coding_end
 Arg[3]      : Int - trans_strand
 Arg[4]      : array_ref of markup position
 Example     : $trans_seq = $transdata->markup_transcript_fasta
 Description : returns the the transcript sequence along with positions for markup
 Return type : a string

=cut

sub markup_transcript_fasta {
    my $self = shift;
    my ($cd_start, $cd_end, $trans_strand, $bps) = @_;
    my $Data = $self->DataObj;
    my $trans = $Data->transcript;
    my $number = $self->param('number');
    my $show = $self->param('show');    
    my $wrap = 60;
    my $count = 0;
    my ($pep_previous, $ambiguities, $previous, $output, $fasta, $peptide)  = '';
    my $pos = 1;
    my $SPACER = $number eq 'on' ? '       ' : '';
    my %bg_color = (  # move to constant MARKUP_COLOUR
        'utr' => $Data->species_defs->ENSEMBL_COLOURS->{'background0'},
        'c0'  => $Data->species_defs->ENSEMBL_COLOURS->{'white'},
        'c1'  => $Data->species_defs->ENSEMBL_COLOURS->{'background3'},
        'c99'  => 'ffcc99',
        'synutr' => '00cc00',
        'sync0'  => '99ff99',
        'sync1'  => '99ff96',
        'indelutr' => '9999ff',
        'indelc0'  => '99ccff',
        'indelc1'  => '99ccff',
        'snputr' => '00cc00',
        'snpc0'  => 'ff9999',
        'snpc1'  => 'ff9999',
    );
    foreach(@$bps) {
      if($count == $wrap) {
        my( $NUMBER, $PEPNUM ) = ('','');
        if($number eq 'on') {
          $NUMBER = sprintf("%6d ",$pos);
          $PEPNUM = ( $pos>=$cd_start && $pos<=$cd_end ) ? sprintf("%6d ",int( ($pos-$cd_start+3)/3) ) : $SPACER ;
          $pos += $wrap;
        }
        $output .= ($show eq 'snps' ? "$SPACER$ambiguities\n" : '' ).
              $NUMBER.$fasta. ($previous eq '' ? '':'</span>')."\n".
              ( ( $show eq 'snps' || $show eq 'peptide' ) ?
          "$PEPNUM$peptide". ($pep_previous eq ''?'':'</span>')."\n\n" : '' );
          $previous=''; $pep_previous=''; $count=0; $peptide = ''; $ambiguities = ''; $fasta ='';
      }
      my $bg = $bg_color{"$_->{'snp'}$_->{'bg'}"};
      my $style = qq(style="color:$_->{'fg'};). ( $bg ? qq( background-color:$bg;) : '' ) .qq(");
      my $pep_style = '';
      if( $show eq 'snps') {
        if($_->{'snp'} ne '') {
           if( $trans_strand == -1 ) {
             $_->{'alleles'}=~tr/acgthvmrdbkynwsACGTDBKYHVMRNWS\//tgcadbkyhvmrnwsTGCAHVMRDBKYNWS\//;
             $_->{'ambigcode'} =~ tr/acgthvmrdbkynwsACGTDBKYHVMRNWS\//tgcadbkyhvmrnwsTGCAHVMRDBKYNWS\//;
           }
            $style .= qq( title="Alleles: $_->{'alleles'}");
        }
        if($_->{'aminoacids'} ne '') {
            $pep_style = qq(style="color: #ff0000" title="$_->{'aminoacids'}");
        }
        $ambiguities.=$_->{'ambigcode'};
      }
      if($style ne $previous) {
        $fasta.=qq(</span>) unless $previous eq '';
        $fasta.=qq(<span $style>) unless $style eq '';
        $previous = $style;
      }
      if($pep_style ne $pep_previous) {
        $peptide.=qq(</span>) unless $pep_previous eq '';
        $peptide.=qq(<span $pep_style>) unless $pep_style eq '';
        $pep_previous = $pep_style;
      }
      $count++;
      $fasta.=$_->{'letter'};
      $peptide.=$_->{'peptide'};
    }
    my( $NUMBER, $PEPNUM ) = ('','');
    if($number eq 'on') {
       $NUMBER = sprintf("%6d ",$pos);
       $PEPNUM = ( $pos>=$cd_start && $pos<=$cd_end ) ? sprintf("%6d ",int( ($pos-$cd_start-1)/3 +1) ) : $SPACER ;
#       warn( $NUMBER, " - ", $PEPNUM );
       $pos += $wrap;
    }
    $output .= ($show eq 'snps' ? "$SPACER$ambiguities\n" : '' ).
               $NUMBER.$fasta. ($previous eq '' ? '':'</span>')."\n".
               ( ( $show eq 'snps' || $show eq 'peptide' ) ?
             "$PEPNUM$peptide". ($pep_previous eq ''?'':'</span>')."\n\n" : '' );
        
    return $output;
}


=head2 markup_help

 Arg[1]      : none
 Example     : $transdata->renderer->markup_help
 Description : renders the transcript cdna fasta markup help image for transview - distinguishes between vega and ensembl on the basis of input parameter set in vega/transview
 Return type : String - HTML

=cut

sub markup_key_image {
    my $self = shift;
    my $show = $self->param('show') ;
    my $image_key ;
	if ($self->param('show_vega_markup')) {
		if     ($show eq 'codons'){$image_key = qq(<img src="/gfx/helpview/transview-key1.png" height="200" width="200" alt="[Key]" border="0">);}
		elsif ($show eq 'peptide'){$image_key = qq(<img src="/gfx/helpview/transview-key2.png" height="200" width="200" alt="[Key]" border="0">); }
		elsif ($show eq 'snps'){$image_key = qq(<img src="/gfx/helpview/transview-key3.png" height="350" width="300" alt="[Key]" border="0">); }
	}
	else {
		if     ($show eq 'codons'){$image_key = qq(<img src="/gfx/helpview/transview-key1.gif" height="200" width="200" alt="[Key]" border="0">);}
		elsif ($show eq 'peptide'){$image_key = qq(<img src="/gfx/helpview/transview-key2.gif" height="200" width="200" alt="[Key]" border="0">); }
		elsif ($show eq 'snps'){$image_key = qq(<img src="/gfx/helpview/transview-key3.gif" height="350" width="300" alt="[Key]" border="0">); }
	}
    return qq(<p>&nbsp;</p><div align='center'>$image_key <p>* - for HTML4 complient browsers \(NS6+, IE5+\)</p><br /></div><br /><br />) if $image_key; 
}

#--------------

=head2 supportingEvidenceImage

 Arg[1]      : none
 Example     : $transdata->supportingEvidenceImage
 Description : gets supporting evidence image
 Return type : string (html)

=cut

sub supportingEvidenceImage{
    my $self = shift;
    my $data = $self->DataObj;
    my $supporting_evidence = $data->get_supporting_evidence;
    return unless $supporting_evidence;

    my $transid = $data->stable_id;
    my $db = $data->get_db;
    my $show = $self->param('showall')   ;
    my $exon_count = $supporting_evidence->{ 'transcript' }{'exon_count'};

    my $hits = scalar(keys %{$supporting_evidence->{ 'hits' }});
    my $html = qq(
        <a name="evidence"><br/><h4>Supporting evidence</h4></a>
        <p>The supporting evidence below consists of the sequence matches
        on which the exon predictions were based and are sorted by alignment score.</p>
    );  
    if ($exon_count > 100 && !$show){
        $html.= qq(<p class="red1">
            NOTE: This transcript has a large number of exons.</p>
            <p>The supporting evidence image may take a while to load, please
            click <a
            href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$transid&db=$db&showall=1">here<a/>
            to view supporting evidence.<br /></p>);
        return $self->print($html);
    }   
    
    if ($hits > 10 && !$show){
        $html.= qq(<P class="red1">
            There are a large number of supporting evidence hits for this transcript. 
            </P><P>Only the top ten 10 hits have been shown. 
            <a href="/$ENV{'ENSEMBL_SPECIES'}/exonview?transcript=$transid&db=$db&showall=1">Click to view all $hits<a/> 
             supporting evidence hits.<br /></P>); 
    }       
    my $dc  = $data->renderSupportingEvidenceImage($supporting_evidence);
    return unless $dc ;
    
    my ( $filename, $imagemap ) =   $self->render_image_imagemap( "suppevid", $dc ) ;   
    $html.= qq(<img border="0" src="$filename" usemap="#$filename"> <map name="$filename">$imagemap </map> <br /><br />);
  
    $self->print($html);
}

sub exonListTable{
    my $self = shift;
    my $full_seq  = $self->param('fullseq');            
    my $only_exon = $self->param('oexon');
    my $sscon     = $self->param('sscon') ;             
    my $flanking  = $self->param('flanking') || 50;
    my $check_full = 'checked' unless !$full_seq ;
    my $check_exon = 'checked' unless !$only_exon;
    $sscon = 25 unless $sscon >= 1;   
    
    $self->print("<h3>Exon Information</h3>");
    $self->Output->print_spreadsheet_table($self->outputExonListTable);
     
    my $HTML = qq(
<table width="100%">
    <tr valign="top" align="center" CLASS="background1" >
        <td valign='top' align='right' width='5%'><font color="#5A85D6">
            <b>&nbsp;&nbsp;<br>KEY:</b></font></td>
        <td align="left" width='20%'>
            <span style="color:green">&nbsp;&nbsp;<br>-Up/downstream region </span>&nbsp;&nbsp;<br>
            <span style="color:blue">-Intron sequence </span>&nbsp;&nbsp;<br>
            <span style="color:9400d3">-UTR region</span></td>
        <td><form method="GET" action="/$ENV{'ENSEMBL_SPECIES'}/exonview">
            <table CELLPADDING="3" CELLSPACING="0" BORDER="0" CLASS="background1" align="center">
                <tr valign="top" >
                    <td><p>Display <input name="sscon" value="$sscon" size="4" maxlength="4">bases either side of intron.&nbsp;&nbsp;</td>     
                    <td><p>View full intron sequence </td>
                    <td><input type="checkbox" name="fullseq" $check_full></td>
                    </tr>
                <tr valign="center" >
                    <td><p>Display <input name="flanking" value="$flanking" size="4" maxlength="4">bases of 5' 3' flanking sequence.&nbsp;&nbsp;</td>
                    <td><p>View only exon sequence </td>
                    <td><input type="checkbox" name="oexon" $check_exon></td>
                    </tr>
                <tr><td colspan='3' align="center"><input type="image" src="/gfx/redraw2.gif" name="submit"
      border="0"></td>);
          
      foreach (qw(gene transcript peptide db)){
            $HTML.=qq(<INPUT TYPE="hidden" NAME="$_" VALUE=") . &CGI::param($_) . qq(">) if &CGI::param($_);
          }
          
      $HTML.= qq(
        </tr></table></form>
        </td>
        <td align="right" width='25%'>&nbsp;<br /></td></tr>
        </table> );
          
  $self->print($HTML);

}

sub outputExonListTable {
    my $self = shift;
    my $transdata = $self->DataObj;
    my @exon_list;
    my $exon_information = $self->get_exon_information;
    return unless @$exon_information;

    my $counter = 0;
    my $table_header = [{'key' => 'Number', 'title' => 'No.', 'width' => '5%', 'align' => 'center' },
                        {'key' => 'exint',  'title' => 'Exon / Intron', 'width' => '20%', 'align' => 'center' }, 
                        {'key' => 'Chr', 'title' => 'Chr', 'width' => '10%', 'align' => 'center' },
                        {'key' => 'Strand',     'title' => 'Strand', 'width' => '10%', 'align' => 'center' },
                        {'key' => 'Start', 'title' => 'Start', 'width' => '15%', 'align' => 'center' },
                        {'key' => 'End', 'title' => 'End', 'width' => '15%', 'align' => 'center' },
                        {'key' => 'Length', 'title' => 'Length', 'width' => '10%', 'align' => 'center' },
                        {'key' => 'Sequence', 'title' => 'Sequence', 'width' => '20%', 'align' => 'left' },
                        ];  
    
    return ( $table_header, $exon_information, {'rows' => [qw(background1 background3)]} );     
}

sub get_exon_information {
    my $self = shift ;
    my $trans = $self->DataObj->transcript;

    # no of bp to show either side of a splice site..
    my $sscon     = $self->param('sscon') ;             

    # no of bp up/down stream of transcript
    my $flanking  = $self->param('flanking') || 50;    

    # flag to display full sequence (introns and exons)
    my $full_seq  = $self->param('fullseq');            

    # display only exons flag
    my $only_exon = $self->param('oexon');              
    my $entry_exon = $self->param('exon');   
    my $coding_start = $trans->coding_region_start; 
    my $coding_end = $trans->coding_region_end;    
    my @el = @{$trans->get_all_Exons};    
    my $strand  = $el[0]->strand;
    my $chr_name = $el[0]->slice->seq_region_name;     
    my @exon_col = qw(blue black);
    my @back_col = qw(background1 background3);     
    my $background = 'background1'; 
    my @exon_info_list;
    my ($exonA, $exonB, $j, $upstream, $exon_info,$intron_info);

    $sscon = 25 unless $sscon >= 1;    
# works out length needed to join intron ends with dots
    my $sscon_dot_length = 60-2*($sscon %30);
    my $flanking_dot_length = 60-($flanking%60);   
# upstream flanking seq
    if ($flanking && !$only_exon){  
        my $exon = $el[0];
        if ($strand == 1){
            $upstream = $exon->slice()->subseq( ($exon->start)-($flanking),   ($exon->start)-1 , $strand);
        }else {
            $upstream = $exon->slice()->subseq( ($exon->end)+1,   ($exon->end)+($flanking),  $strand);
        }   
        $upstream =  lc(('.'x $flanking_dot_length).$upstream);
        $upstream =~ s/([\.\w]{60})/$1<BR>/g; 
        $exon_info = { 'exint'    => qq(5' upstream sequence&nbsp), #'
                       'Sequence' => qq(<font face="courier" color="green">$upstream</font>)};

        push @exon_info_list, $exon_info;
    }


       
    # Loop over each exon
    for ($j=1; $j<= scalar(@el); $j++) {
      my ($intron_start, $intron_end, $intron_len, $intron_seq) ; 
      my $col = $exon_col[$j%2];                    #choose exon text colour
      $exonA = $el[$j-1];
      $exonB = $el[$j];
      
      my $intron_id = "Intron $j-".($j+1)  ;
      my $dots = '.'x $sscon_dot_length;                 
      my $seq       = uc($exonA->seq()->seq());
      my $seqlen    = length($seq);                 
      my $exonA_ID  = $exonA->stable_id;
      my $exonA_start   = $exonA->start;
      my $exonA_end     = $exonA->end;
      my $exonB_start   = $exonB->start if $exonB ;
      my $exonB_end     = $exonB->end if $exonB ;
      my $utrspan_start = qq(<span style="color:9400d3">);  ##set colour of UTR
      my $count = 0;
      my $k = 0;


      # Is this exon entirely UTR?
      if ($coding_end < $exonA_start || $coding_start > $exonA_end){
            $seq   =~ s/([\.\w]{60})/$1<\/span><BR>$utrspan_start/g ;
            $seq   .= qq(</span>);
            $seq = "$utrspan_start"."$seq";
      }
      # Handle reverse strand transcripts.  Yes, this means we have a bunch of 
      # duplicated code to handle forward strand.
      elsif ($strand eq '-1'){
        my @exon_nt  = split '', $seq;
        my $coding_len =  ($exonA_end) - $coding_start + 1 ;
        my $utr_len =  $exonA_end - $coding_end   ;

	# CDS is within this exon, and we have UTR start and end
        if ($coding_start > $exonA_start &&  $coding_end < $exonA_end){
            $seq = qq($utrspan_start);
            for (@exon_nt){
                if ($count == 60 && ($k < $coding_len && $k > $utr_len)){
                    $seq .= "<br>";
                    $count =0;
                }elsif ($count == 60 && ($k > $coding_len || $k < $utr_len)){
                    $seq .= "</span><br>$utrspan_start";
                    $count =0;
                }elsif ($k == $utr_len || $k == $coding_len){
                    $seq .= "</span>";
		  if ($count == 60) {
		    $seq .= "<br>";
		    $count = 0;
		  }
                }   
                $seq .= $_ ;
                $count++;
                $k++;
            }       
        }
	# exon starts with UTR
	elsif ($coding_start > $exonA_start  ){
            $seq = "";
            for (@exon_nt){
                if ($count == 60 && ($k > $coding_len)){
                    $seq .= "</span><br>$utrspan_start";
                    $count =0;
                }elsif ($count == 60 && $k < $coding_len){
                    $seq .= "<br>";
                    $count =0;
                }elsif ($k == $coding_len){
		  if ($count == 60) {
		    $seq .= "<br>";
		    $count = 0;
		  }
                    $seq .= qq($utrspan_start);
                }   
                $seq .= $_ ;
                $count++;
                $k++;
            }       
        }
	# exon ends with UTR
	elsif($coding_end < $exonA_end ){
            $seq = $utrspan_start;
            for (@exon_nt){ 
                if ($count == 60 && $utr_len > $k ){
                    $seq .= "</span><br>$utrspan_start";
                    $count =0;
                }elsif ($count == 60 && $k > $utr_len){
                    $seq .= "<br>";
                    $count =0;  
                }elsif ($k == $utr_len){
                    $seq .= qq(</span>);
		  if ($count == 60) {
		    $seq .= "<br>";
		    $count = 0;
		  }
                }    
                $seq .= $_ ;
                $count++;
                $k++;
            }
        $seq .= "</span>";
        }
	else{ # entirely coding exon		
            $seq =~ s/([\.\w]{60})/$1<BR>/g ;
        }
      }
      # Now we handle the forward strand transcripts.  Same as before, except
      # the exons come in reverse order, and we have to flip utr & coding 
      # length
      else{
        my @exon_nt  = split '', $seq;
        my $utr_len =  $coding_start - $exonA_start ;
        my $coding_len =  $seqlen - ($exonA_end - $coding_end)  ;

	# CDS is within this exon, and we have UTR start and end
        if ($coding_start > $exonA_start &&  $coding_end < $exonA_end){
            $seq = qq($utrspan_start);
            for (@exon_nt){
                if ($count == 60 && ($k > $utr_len && $k < $coding_len)){
                    $seq .= "<br>";
                    $count =0;
                }elsif ($count == 60 && ($k < $utr_len || $k > $coding_len)){
                    $seq .= "</span><br>$utrspan_start";
                    $count =0;
                }elsif ($k == $utr_len || $k == $coding_len){
                    $seq .= "</span>";
		  if ($count == 60) {
		    $seq .= "<br>";
		    $count = 0;
		  }

                }   
                $seq .= $_ ;
                $count++;
                $k++;
            }       
        }

	# exon starts with UTR
	elsif ($coding_start > $exonA_start ){                         
            $seq = qq($utrspan_start);
            for (@exon_nt){
                if ($count == 60 && ($k > $utr_len)){
                    $seq .= "<br>";
                    $count =0;
                }elsif ($count == 60 && $k < $utr_len){
                    $seq .= "</span><br>$utrspan_start";
                    $count =0;
                }elsif ($k == $utr_len){
		  $seq .= "</span>";
		  if ($count == 60) {
		    $seq .= "<br>";
		    $count = 0;
		  }

                }   
                $seq .= $_ ;
                $count++;
                $k++;
            }       
        }
	# exon ends with UTR
	elsif($coding_end < $exonA_end ){
            $seq = '';      
            for (@exon_nt){ 

                if ($count == 60 && $coding_len > $k ){
                    $seq .= "<br>";
                    $count =0;
                }elsif ($count == 60 && $k > $coding_len){
                    $seq .= "</span><br>$utrspan_start";
                    $count =0;
                }elsif ($k == $coding_len){
		  if ($count == 60) {
		    $seq .= "<br>";
		    $count = 0;
		  }
		  $seq .= qq($utrspan_start);
                }    
                $seq .= $_ ;
                $count++;
                $k++;
            }
        $seq .= "</span>";
        }
	# Entirely coding exon.
	else{ 
            $seq =~ s/([\.\w]{60})/$1<BR>/g ;
        }             
      }
     
      if ($entry_exon && $entry_exon eq $exonA_ID){
        $exonA_ID = "<b>$exonA_ID</b>" ;
      }
     $exon_info = {     'Number'    => $j,
                        'exint'     => qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr_name&chr_start=$exonA_start&chr_end=$exonA_end">$exonA_ID</a>), 
                        'Chr'       => $chr_name,
                        'Strand'    => $strand,
                        'Start'     => $exonA_start,
                        'End'       => $exonA_end,
                        'Length'    => "$seqlen bp",
                        'Sequence'  => qq(<font face="courier" color="black">$seq</font>)};
        push @exon_info_list, $exon_info;                     
    if (!$only_exon && $exonB){      
        eval{
          if($strand == 1 ) { # ...on the forward strand
                $intron_start = $exonA_end+1; 
                $intron_end = $exonB_start-1;
                $intron_len = ($intron_end - $intron_start) +1;
                if (!$full_seq && $intron_len > ($sscon *2)){
                    my $seq_start_sscon = $exonA->slice()->subseq( ($intron_start),   ($intron_start)+($sscon-1),  $strand);
                    my $seq_end_sscon = $exonB->slice()->subseq( ($intron_end)-($sscon-1), ($intron_end), $strand);
                    $intron_seq = "$seq_start_sscon$dots$seq_end_sscon";
                }
                else {$intron_seq = $exonA->slice()->subseq( ($intron_start),   ($intron_end),   $strand);}
            } else { # ...on the reverse strand
                $intron_start = $exonB_end+1; 
                $intron_end = $exonA_start-1;
                $intron_len = ($intron_end - $intron_start) +1;
                if (!$full_seq && $intron_len > ($sscon *2)){
                    my $seq_end_sscon = $exonA->slice()->subseq( ($intron_start), ($intron_start)+($sscon-1), $strand);
                    my $seq_start_sscon = $exonB->slice()->subseq( ($intron_end)-($sscon-1), ($intron_end), $strand);
                    $intron_seq = "$seq_start_sscon$dots$seq_end_sscon";
                }
                else {$intron_seq = $exonA->slice()->subseq( ($intron_start),   ($intron_end),   $strand);}
            }       
        }; # end of eval     
        $intron_seq =  lc($intron_seq);
        $intron_seq =~ s/([\.\w]{60})/$1<BR>/g;
      
         $intron_info = {   'Number'    => "&nbsp",
                            'exint'     => $intron_id, 
                            'Chr'       => $chr_name,
                            'Strand'    => $strand,
                            'Start'     => $intron_start,
                            'End'       => $intron_end,
                            'Length'    => "$intron_len bp",
                            'Sequence'  => qq(<font face="courier" color="blue">$intron_seq</font>)};
        push @exon_info_list, $intron_info;
      }     
    }     #finished foreach loop

# Print last exon in transcript (change)
    if ($exonB){
    my $exonB_seq   = uc($exonB->seq()->seq());
    my $seqlen      = length($exonB_seq);
       $exonB_seq   =~ s/([\*\w]{60})/$1<BR>/g;
    my $exonB_ID    = $exonB->stable_id;
    my $exonB_start = $exonB->start;
    my $exonB_end   = $exonB->end;
    $exon_info = {      'Number'    => $j,
                        'exint'     => qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr_name&chr_start=$exonB_start&chr_end=$exonB_end">$exonB_ID</a>),  
                        'Chr'       => $chr_name,
                        'Strand'    => $strand,
                        'Start'     => $exonB_start,
                        'End'       => $exonB_end,
                        'Length'    => "$seqlen bp",
                        'Sequence'  => qq(<font face="courier" color="black">$exonB_seq</font>)};
    push @exon_info_list, $exon_info;
}   
# downstream HTML
    if ($flanking && !$only_exon){
        my $exon = $exonB ? $exonB : $exonA ;           # if single exon use exonA
        my $downstream;
        if ($strand == 1){
            $downstream = $exon->slice()->subseq( ($exon->end)+1,($exon->end)+($flanking),  $strand );
        }else {
            $downstream = $exon->slice()->subseq( ($exon->start)-($flanking), ($exon->start)-1 ,  $strand ) ;
        }
        $downstream = lc($downstream.('.'x $flanking_dot_length));
        $downstream =~ s/([\.\w]{60})/$1<BR>/g;  
        $exon_info = { 'exint'    => qq(3' downstream sequence&nbsp), #'
                       'Sequence' => qq(<font face="courier" color="green">$downstream</font>)};
        push @exon_info_list, $exon_info;
    }


return (\@exon_info_list);
} # end of sub


sub mouse_note{
	my $self = shift;
	my $transdata =  $self->DataObj();
	my $label = qq(<font color="red">Assembly m32</font>);
	my $params = join ('&', (map {$_."=".$self->param($_)} $self->Input->param));
	my $link = qq(<a href="http://mouse30.ensembl.org/Mus_musculus/$ENV{'ENSEMBL_SCRIPT'}/?$params">here</a>);

	my $M32_GENES = $transdata->species_defs->M32_GENES || {};
	my $html;
	my $gene = $transdata->gene;
	my $geneid= $gene->stable_id if $gene;
	if ( exists $M32_GENES->{$geneid}){
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
