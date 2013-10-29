package EnsEMBL::Web::Renderer::Translation::HTML;

=head1 NAME

EnsEMBL::Web::Renderer::Translation::HTML.pm 

=head1 SYNOPSIS

This object creates HTML to be output to the HTML output object

=head1 DESCRIPTION

    $translation_renderer = $pep_data->renderer;
	$translation_renderer->outputGenericPeptideTable();		 
   	
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
use Bio::Perl;
use EnsEMBL::Web::Renderer::HTML;
use EnsEMBL::Web::Renderer::Translation;

@ISA = qw( EnsEMBL::Web::Renderer::HTML EnsEMBL::Web::Renderer::Translation );

#----------------------------------------------------------------------

=head2 outputGenericPeptideTable

 Arg[1]      : none
 Example     : $pepdata->renderer->outputGenericPeptideTable
 Description : Wrapper for the table at top of protview
 Return type : none

=cut

sub outputGenericPeptideTable {
    my $self = shift;
    $self->prot_table_title;
    $self->Output->print_two_col_table
      ( 
       $self->xref_display_ID,
       $self->stable_id,
       $self->other_view_links,
       $self->description,
       $self->prediction_method,
       $self->go_links,
       $self->interpro_links,
       $self->family_links,
#       $self->das_configurator,
	$self->das_sources_selector,
	$self->protview_peptide_image,
#       $self->das_annotation,
	$self->das_sources_annotation,
       $self->exportview_link, 
      );
}

#----------------------------------------------------------------------

=head2 outputprotviewSeqInfo

 Arg[1]      : none
 Example     : $pepdata->renderer->outputprotviewSeqInfo
 Description : Wrapper for the markered up fasta seq and pepstats on protview
 Return type : none

=cut

sub outputprotviewSeqInfo {
    my $self = shift;	
    $self->Output->print_hidden_table
      ( $self->print_fasta_seq,
	[$self->pepstats , $self->markup_key_image], );
}

#----------------------------------------------------------------------

=head2 outputprotviewBottomTables

 Arg[1]      : none
 Example     : $pepdata->renderer->outputprotviewBottomTables
 Description : Wrapper for the list tables at the bottom of protview
 Return type : none

=cut

sub outputprotviewBottomTables {
  my $self = shift;
  $self->Output->print_spreadsheet_table($self->domain_list);
  $self->Output->print_spreadsheet_table($self->other_protein_feature_list);
  $self->Output->print_spreadsheet_table($self->snp_list);
}

#----------------------------------------------------------------------

=head2 prot_table_title

 Arg[1]      : none
 Example     : $pepdata->renderer->prot_table_title
 Description : prints the tile for the peptide table at top of protview
 Return type : string - html

=cut

sub prot_table_title{
    my $self = shift;
    my $type = $self->DataObj->gene_type;				
    $self->Output->generic_table_title("<h3>$type Protein Report</h3>");
}

#----------------------------------------------------------------------

=head2 xref_display_ID

 Arg[1]      : none
 Example     : $pepdata->renderer->xref_display_ID
 Description : renders the xref_display_id and type in two_col_table format
 Return type : key-value pair, label and html

=cut

sub xref_display_ID{
    my $self = shift;
    my $pep = $self->DataObj();
    my $label = 'Peptide';
    return unless $pep->display_xref();
    my ($display_name, $dbname) = $pep->display_xref();
    my $html = "<b>$display_name </b> <small> ($dbname ID)</small>";
    if (grep { $_->dbname eq 'CCDS' } @{$pep->translation->get_all_DBLinks}) {
        my $ccds_url = $self->ExtURL->get_url('CCDS') || '/Homo_sapiens/ccds.html';
        $html .= '<br />Member of Human '.$self->url('CCDS', $ccds_url).' set'
    }
    return ($label, $html);   	
}

#----------------------------------------------------------------------

=head2 stable_id

 Arg[1]      : (optional) String
                Label
 Example     : $pepdata->renderer->stable_id
 Description : renders the ensembl_stable_id in two_col_table format
 Return type : key-value pair, label and html

=cut

sub stable_id{
    my $self = shift;
    my $transl =  $self->DataObj();
    my $type = $transl->gene_type;

    my( $display_id, $id_type );
    if( $display_id = $transl->stable_id ){
      $id_type = 'Translation ID';
    }
    else{ 
      $display_id = $transl->transcript->stable_id;
      $id_type = 'ID';
    }

    my $label = shift || "$type $id_type";
    my $html  = "<b>$display_id</b>";

    return ($label, $html);
}

=head2 version

 Arg[1]      : (optional) String
               table header text
 Example     : $pepdata->renderer->version($string)
 Description : Renders the translation version and modification time in
               two_col_table format
 Return type : key-value pair, label and html

=cut

sub version {
    my $self = shift;
    my $trans =  $self->DataObj();
    my $label = shift || "Version";
    my $html = $trans->version;
    return ($label, $html);
}

#----------------------------------------------------------------------

=head2 other_view_links

 Arg[1]      : (optional) String
                Label
 Example     : $pepdata->renderer->other_view_links
 Description : renders the ensembl_links to other views in two_col_table format
 Return type : key-value pair, label and html

=cut

sub other_view_links{
    my $self = shift;	
    my $data =  $self->DataObj();
    my $label = shift || $data->gene_type ." Translation";
    my $db = $data->get_db() ;
    my $geneid  = $data->gene ? $data->gene->stable_id : '';
    my $transid = $data->transcript->stable_id;
    my $pepid   = $data->stable_id;

    my $sp = $ENV{'ENSEMBL_SPECIES'};

    my $gene_href = "/$sp/geneview?gene=$geneid&db=$db";
    my $tran_href = "/$sp/transview?transcript=$transid&db=$db";
    my $exon_href = "/$sp/exonview?transcript=$transid&db=$db";

    my $html = '';

    if( $geneid ){ 
	$html .= qq(
		    This peptide is a product of gene <a href="$gene_href">$geneid</a> <br /> );
    }
    if ($self->param('show_vega_evidence_link')) {
	$html .= qq(
		    [ <a href="$tran_href"> Transcript Information </a> ]
		    [ <a href="$exon_href"> Exon information and supporting evidence</a> ] );
    } else {
	$html .= qq(
		    [ <a href="$tran_href"> Transcript Information </a> ]
		    [ <a href="$exon_href"> Exon Information </a> ] );	
    }
    return ($label, $html);
}

#----------------------------------------------------------------------

=head2 description

 Arg[1]      : none
 Example     : $pepdata->renderer->description
 Description : renders the gene description in two_col_table format
 Return type : key-value pair, label and html

=cut

sub description {
    my $self = shift ;
    my $urls = $self->ExtURL;
    my $description = $self->DataObj->description();    
       $description =~ s/EC\s+([-*\d]+\.[-*\d]+\.[-*\d]+\.[-*\d]+)/$self->EC_URL($1)/e;
       $description =~ s/\[\w+:(\w+)\;\w+:(\w+)\]//g;
    my ($edb, $acc) = ($1, $2);
    my $url 	    = "[<a href='".$urls->get_url($edb, $acc)."'>Source: $1 ($2)</a>]" unless !$acc;    
    my $label = 'Description';
    my $html = qq($description <small>$url</small>&nbsp;<br />);	   
    return ($label, $html);
}

#----------------------------------------------------------------------

=head2 author

 Arg[1]	     : (optional) String
               Label
 Example     : $pepdata->renderer->author
 Description : returns the annotation author for use in two_col_table
 Return type : key-value pair of label/html

=cut

sub author {
    my $self = shift;
    my $label = shift || "Author";
    my $pep = $self->DataObj();
    my $author = $pep->get_author_name;
    my $html;
    if ($author) {
        $html .= "This locus was annotated by " . $author . " ";
        $html .= $self->email_URL($pep->get_author_email);
    } else {
        $html = "unknown";
    }
    return ($label, $html);
}

#----------------------------------------------------------------------
=head2 prediction_method

 Arg[1]      : none
 Example     : $pepdata->renderer->prediction_method
 Description : renders the prediction method text in two_col_table format
 Return type : key-value pair, label and html

=cut

sub prediction_method {
    my $self = shift ;
    my $data = $self->DataObj;
    my $db = $data->get_db() ;	
    my $label = ( $db eq 'vega' ? 'Curation' : 'Prediction' ).' Method';
    my $prediction_text = $data->get_prediction_method ."&nbsp;";		
    return ($label, $prediction_text );
}

#----------------------------------------------------------------------

=head2 exportview_link

 Arg[1]      : none
 Example     : $pepdata->renderer->exportview_link
 Description : renders the exportview link in two_col_table format
 Return type : key-value pair, label and html

=cut

sub exportview_link {
    my $self = shift;
    my $pepid = $self->DataObj->stable_id();
    my $label = 'Export Data';
    my $sp  = $ENV{'ENSEMBL_SPECIES'};
    my $href = "/$sp/exportview?type=feature&ftype=peptide&id=$pepid";
    my $txt  = "Export peptide data in EMBL, GenBank or FASTA";
    my $html = qq( <a href="$href"> $txt </a> );
    return ($label, $html);
}

#----------------------------------------------------------------------

=head2 interpro_links

 Arg[1]      : none
 Example     : $pepdata->renderer->interpro_links
 Description : renders the interpro links in two_col_table format
 Return type : key-value pair, label and html

=cut

sub interpro_links {
    my $self = shift;
    my $pep = $self->DataObj;
    my $interpro_hash = $pep->get_interpro_links();
    my $urls = $self->ExtURL;
    
    return unless (%$interpro_hash);
    my $label = 'InterPro';

# add table call here
    my $html = qq(<table class="hidden">);
    for my $accession (keys %$interpro_hash){
	my $interpro_url = $urls->get_url('INTERPRO',$accession);
	my $desc = $interpro_hash->{$accession};
	$html .= qq(<tr>
	    	    <td><A HREF="$interpro_url">$accession</A> </td><td> $desc 
		       - [<A HREF="domainview?domainentry=$accession">View other genes with this domain</A>]</td></tr>);
    }
    $html .= qq( </table> );
    return ($label, $html);
}


#----------------------------------------------------------------------
=head2 go_links

 Arg[1]      : none
 Example     : $transdata->renderer->go_links
 Description : renders the go links in two_col_table if the go
               database exists
               This is copy of the Transcript::HTML::go_links
               method, but there is no way to get hold of a
               Transcript::HTML object from here. Site-wide updates need
               making in both places.
 Return type : Key / value pair - label and HTML

=cut


use EnsEMBL::Web::Renderer::Transcript::HTML;
use EnsEMBL::Web::Data::Transcript;
sub get_transcript_renderer{
  # The grimest of hacks to allow Transcript data to de rendered directly!
  my $self = shift;

  # Return cached
  if( my $r = $self->{_transcript_renderer} ){ return $r }

  # Create fresh
  my $data = EnsEMBL::Web::Data::Transcript->new
      ( $self->DataObj->transcript,
        $self->DataObj->DBConnection,
        $self->Input, $self->Output );
  
  $self->{_transcript_renderer} = EnsEMBL::Web::Renderer::Transcript::HTML->new
      ( $data, $self->Input, $self->Output );

  return $self->{_transcript_renderer};
}


sub go_links { 
    my $self = shift;
    return $self->get_transcript_renderer->go_links;
}

#----------------------------------------------------------------------

=head2 family_links

 Arg[1]      : none
 Example     : $pepdata->renderer->family_links
 Description : renders the family links in two_col_table format
 Return type : key-value pair, label and html

=cut

sub family_links {
    my $self = shift;
    my $pep = $self->DataObj;
    my $families = $pep->get_family_links();
	
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

#----------------------------------------------------------------------

=head2 protview_peptide_image

 Arg[1]      : none
 Example     : $pepdata->renderer->protview_peptide_image
 Description : wrapper to print peptide image in two_col_table format for protview
 Return type : Key / value pair - label and HTML

=cut

sub protview_peptide_image{
    my $self = shift;   
    my $label = 'Protein Features ';
    my $html =  $self->peptide_image;
    return ($label, $html) if $html;
}

#----------------------------------------------------------------------

=head2 peptide_image

 Arg[1]      : none
 Example     : $pepdata->renderer->peptide_image
 Description : wrapper to print peptide image only for protview with snps
 Return type : string - HTML

=cut

sub peptide_image{
    my $self = shift;
    my $pepdata = $self->DataObj;
    my $peptideid = $pepdata->stable_id;
    my $db = $pepdata->get_db ;
    my $dc = $pepdata->render_peptide_image( 'protview' );
    my ($filename, $map) = $self->render_image_imagemap( "pepimg", $dc ) ;

    my $html_tmpl = qq(
<br />
<div align="center">
 <img border="0" src="%s" usemap="%s" >
 <map name="%s">%s</map>
</div>
<br/> );

    my $html = sprintf( $html_tmpl,
			"$filename", $filename, #src
			$filename, $map, );              #map

    return $html;
}

#----------------------------------------------------------------------

=head2 print_fasta_seq

 Arg[1]      : none
 Example     : $pepdata->renderer->print_fasta_seq
 Description : prints markedup peptide fasta sequence
 Return type : String - HTML

=cut

sub print_fasta_seq {
  my $self = shift;
  my ($sel_number, $sel_no) =  ('','');
  my( $sel_snps, $sel_exons, $sel_plain) = ('','','','');
  my $Data = $self->DataObj;
  my $db = $Data->get_db() ;
  my $stable_id = $Data->stable_id;
  my $trans_id = $Data->transcript->stable_id;
  my $number = $self->param('number');
  my $show = $self->param('show') || 'plain';

  if($show eq'snps') { $sel_snps = ' selected'; }
  elsif($show eq 'exons') {$sel_exons=' selected'; }
  else{ $sel_plain =' selected'; }

  my $fasta = ($show eq 'plain') ? $Data->get_pep_seq : $self->do_markedup_pep_seq;
  if($number eq 'on') { $sel_number = ' selected'; } else {$sel_no=' selected'; }
  my $SNP_LINE = "";
  if (exists($Data->species_defs->databases->{'ENSEMBL_VARIATION'}) ||
      exists($Data->species_defs->databases->{'ENSEMBL_GLOVAR'})) {
    $SNP_LINE = "<option value='snps' $sel_snps>Exons/SNPs</option>";
  }
  my $html = qq(
  	<h5>Peptide Sequence</h5>
    	<form action="/$ENV{'ENSEMBL_SPECIES'}/$ENV{'ENSEMBL_SCRIPT'}" id="ddForm" name="ddForm" method="get">
    	<input type="hidden" name="peptide" value="$stable_id" />
      <input type="hidden" name="transcript" value="$trans_id" />
    	<input type="hidden" name="db" value="$db" />
        <select name="show" onChange="document.ddForm.submit()">
          <option value="plain" $sel_plain>No markup</option>
          <option value="exons" $sel_exons>Exons highlighted</option>
          $SNP_LINE
        </select>
        <select name="number" onChange="document.ddForm.submit()">
          <option value="on" $sel_number>Number residues</option>
          <option value="off" $sel_no>No numbers</option>
        </select>
  	</form>
        <pre class="snpmarkup">$fasta</pre><br /><br />  );
}

#----------------------------------------------------------------------

=head2 do_markedup_pep_seq

 Arg[1]		   : none
 Example     : $pep_seq = $pepdata->do_markedup_pep_seq
 Description : returns the the peptide sequence  with markup
 Return type : a string

=cut


sub do_markedup_pep_seq {
  my $self       = shift;
  my $Data       = $self->DataObj;
  my $number     = $self->param('number');
  my $show       = $self->param('show');
  my $peptide    = $Data->translation;
  my $trans      = $Data->transcript;
  my $protein    = $Data->translation;
  my $pep_splice = $Data->pep_splice_site($protein);
  my $pep_snps   = $Data->pep_snps;
  my $wrap       = 60;
  my $db         = $Data->get_db;
  my $pep_id     = $Data->stable_id;
  my $pep_seq    = $protein->seq;
  my @exon_colours = qw(black blue red);
  my %bg_color = (
                  'c0'      => $Data->species_defs->ENSEMBL_COLOURS->{'white'},
                  'syn'     => '99ff99',
                  'insert'  => '99ccff',
                  'delete'  => '99ccff',
                  'snp'     => 'ff9999',
                 );
  my @aas = map {{'aa' => $_ }} split //, uc($pep_seq) ; # store peptide seq in hash
  my ($output, $fasta, $previous) = '';
  my ($count, $flip, $i) = 0;
  my $pos = 1;
  my $SPACER = $number eq 'on' ? '       ' : '';
  foreach (@aas) {						# build markup
    if($count == $wrap) {
      my $NUMBER = '';
      if($number eq 'on') {
        $NUMBER = sprintf("%6d ",$pos);
        $pos += $wrap;
      }
      $output .= ($show eq 'snps' ? "\n$SPACER" : '' ).
        $NUMBER.$fasta. ($previous eq '' ? '':'</span>')."\n" ;
      $previous=''; $count=0; $fasta ='';
    }
    if ( $pep_splice->{$i}{'exon'} ){ $flip = 1 - $flip }
   	my $fg = $pep_splice->{$i}{'overlap'} ? $exon_colours[2] : $exon_colours[$flip];
    my $bg = $bg_color{$pep_snps->[$i]{'type'}};
    my $style = qq(style="color:$fg;");
    my $type = $pep_snps->[$i]{'type'};
    if( $show eq 'snps') {
	    $style = qq(style="color:$fg;). ( $bg ? qq( background-color:$bg;) : '' ) .qq(");
      if ($type eq 'snp'){
        $style .= qq(title="Residues: $pep_snps->[$i]{'pep_snp'} ");
	    }
	    if ($type eq 'syn'){
        my $string = '';
        for my $letter ( 0..2 ){
		    	$string .= $pep_snps->[$i]{'ambigcode'}[$letter]  ? '('.$pep_snps->[$i]{'ambigcode'}[$letter].')' : $pep_snps->[$i]{'nt'}[$letter];
        }
        $style .= qq(title="Codon: $string ");
	    }
	    if($type eq 'insert') {
        $pep_snps->[$i]{'alleles'} = join '', @{$pep_snps->[$i]{'nt'}};
        $pep_snps->[$i]{'alleles'} = Bio::Perl::translate_as_string($pep_snps->[$i]{'alleles'});   # translate insertion.. bio::perl call
        $style .= qq(title="Insert: $pep_snps->[$i]{'allele'} ");
	    }
	    if($type eq 'delete') {
        $style .= qq(title="Deletion: $pep_snps->[$i]{'allele'} ");
      }
	    if($type eq 'frameshift') {
        $style .= qq(title="Frame-shift ");
      }
    }		# end if snp
    
    if($style ne $previous) {
      $fasta.=qq(</span>) unless $previous eq '';
      $fasta.=qq(<span $style>) unless $style eq '';
      $previous = $style;
    }
    $count++; $i++;
    $fasta .= $_->{'aa'};    
  }
  
  my $NUMBER = '';
  if($number eq 'on') {
    $NUMBER = sprintf("%6d ",$pos); $pos += $wrap;
  }
  $output .= ($show eq 'snps' ? "\n$SPACER" : '' ).$NUMBER.$fasta. ($previous eq '' ? '':'</span>')."\n";
  
  my( $sel_snps, $sel_exons,$sel_peptide)=('','','');
  if($show eq'snps') { $sel_snps = ' selected'; }
  elsif($show eq 'exons') {$sel_exons=' selected'; } 
  else { ($sel_snps, $sel_exons ) = ''; }
  
  my ( $sel_numbers, $sel_no)=('','');
  if($number eq'on') { $sel_numbers = ' selected'; }
  else {$sel_no=' selected'; }
  
  my $SNP_LINE = exists($Data->species_defs->databases->{'ENSEMBL_VARIATION'}) ? qq(<option value="snps" $sel_snps>Exons/SNPs</option>) : '' ;
  return ($output);
}

#----------------------------------------------------------------------

=head2 pepstats

 Arg[1]      : none
 Example     : $pepdata->renderer->pepstats
 Description : Prints pepstat table
 Return type : String - HTML

=cut

sub pepstats {
    my $self = shift;
    my $pepstats = $self->DataObj->get_pepstats() ;

    return ("<b>Pepstats could not be <br /> calculated for this peptide.</b>") unless (%{$pepstats||{}}) ;
    my $html .= qq(<h5>Peptide Stats</h5> <table class="hidden">);
# call table method
    foreach my $label (keys %$pepstats){
    	$html .= qq(<tr align="left">
	    	<td><b>$label</b> </td>
		    <td>= </td>
		    <td> $pepstats->{$label} </td>
		 	</tr>);
    }
    $html .= "</table><br />&nbsp;<br />" ;
    return $html;
}

#----------------------------------------------------------------------

=head2 markup_help

 Arg[1]      : none
 Example     : $pepdata->renderer->markup_key
 Description : renders the peptide fasta markup help image for protview, distinguishes between vega and ensembl images
 Return type : String - HTML

=cut

sub markup_key_image {
    my $self = shift;
    my $show = $self->param('show') ;
    my $image_key;
	if ($self->param('show_vega_markup')) {
		if ($show eq 'exons') {
			$image_key = qq(<img src="/gfx/helpview/protview_key1.png" alt="[Key]" border="0">);
		}
		if ($show eq 'snps') {
			$image_key = qq(<img src="/gfx/helpview/protview_key2.png" alt="[Key]" border="0">);
		}
	}
	else { 
		if ($show eq 'exons' || $show eq 'snps') {
			$image_key = qq(<img src="/gfx/helpview/protview_key1.gif" alt="[Key]" border="0">);
		}
    }
    return qq($image_key <br />* - for HTML4 complient browsers \(NS6+, IE5+\)<br /><br />) if $image_key;
}

#----------------------------------------------------------------------

=head2 domain_list

 Arg[1]      : none
 Example     : $pepdata->renderer->domain_list
 Description : Sorts domains into correct format for spreadsheet table, also set various
 				table configuration parameters
 Return type : list of array refs

=cut

sub domain_list{
    my $self = shift;
    my $pepdata = $self->DataObj;
    my @domains = @{$pepdata->get_protein_domains()};
    return unless @domains ;
	my @domain_list;
    my $urls = $self->ExtURL;
	my $table_header = [{'key' => 'Domain type', 'title' => 'Domain type', 'width' => '25%', 'align' => 'center' },
						{'key' => 'Accession number', 	'title' => 'Accession number', 'width' => '15%', 'align' => 'center' },
						{'key' => 'Description', 	'title' => 'Description', 'width' => '30%', 'align' => 'center' },
						{'key' => 'Start', 'title' => 'Start', 'width' => '15%', 'align' => 'center' },
						{'key' => 'End', 	'title' => 'End', 'width' => '15%', 'align' => 'center' },
						];

# may do a code reference to url call else clean up url creation on domain type
	foreach my $domain (@domains){
            my $db = uc($domain->analysis->db);
            my $id = $domain->hseqname;
            my $url = $urls->get_url($db,$id);
            if ($url){
                $url = qq(<A href="$url">$id</A>);
            } else {
                $url = $id;
	    }
    	    my $tmp_array = {
				'Domain type' 		=> $domain->analysis->db,
	    	                'Accession number' 	=> $url,
			  	'Start'   	     	=> $domain->start,
			  	'End'     	     	=> $domain->end,
			  	'Description'     	=> $domain->idesc}, ;

		push @domain_list, $tmp_array;
    }
    return ($table_header, \@domain_list);
}

#----------------------------------------------------------------------

=head2 other_protein_feature_list

 Arg[1]      : none
 Example     : $pepdata->renderer->other_protein_feature_list
 Description : Sorts domains into correct format for spreadsheet table, also set various
 				table configuration parameters
 Return type : list of array refs

=cut

sub other_protein_feature_list {
  my $self = shift ;
  my $pepdata = $self->DataObj;     
  my @other_list;
  
	my $other = [
               @{$pepdata->get_all_ProteinFeatures('tmhmm')},
               @{$pepdata->get_all_ProteinFeatures('Signalp')},
               @{$pepdata->get_all_ProteinFeatures('ncoils')},
               @{$pepdata->get_all_ProteinFeatures('Seg')},
              ];	
  
	return unless @$other;
	
	my $table_header = [{'key' => 'Domain type', 
                       'title' => 'Domain type', 
                       'width' => '40%', 
                       'align' => 'center' },
                      {'key' => 'Start', 	
                       'title' => 'Start', 
                       'width' => '25%', 
                       'align' => 'center' },
                      {'key' => 'End', 	
                       'title' => 'End', 
                       'width' => '25%', 
                       'align' => 'center' },];
	
  foreach my $oth (@$other) {
    my $analysis = $oth->analysis;
    my $domain_type = $analysis->db || $analysis->logic_name || 'unknown';
    $domain_type =~ s/_/ /g;
    my $tmp_array = { 'Domain type' => ucfirst($domain_type),
                      'Start'   	   => $oth->start,
                      'End'     	   => $oth->end},;	
    push @other_list, $tmp_array;
  }
  return ($table_header, \@other_list);		
}

#----------------------------------------------------------------------

=head2 snp_list

 Arg[1]      : none
 Example     : $pepdata->renderer->snp_list
 Description : Sorts snp list into correct format for spreadsheet table, also set various
 				table configuration parameters
 Return type : list of array refs

=cut

sub snp_list {	    
  my $self = shift ;
  my $pepdata = $self->DataObj;
  my @snp_list;
  my $snps = $pepdata->pep_snps();
  return unless @$snps;
	
  my $counter = 0;
  my $table_header = [{'key' => 'Residue', 'title' => 'Residue', 
		       'width' => '10%', 'align' => 'center' },
		      {'key' => 'SNP ID', 	'title' => 'SNP ID', 
		       'width' => '15%', 'align' => 'center' }, 
		      {'key' => 'SNP type', 'title' => 'SNP type', 
		       'width' => '20%', 'align' => 'center' },
		      {'key' => 'Alleles', 	'title' => 'Alleles', 
		       'width' => '20%', 'align' => 'center' },
		      {'key' => 'Ambiguity code', 'title' => 'Ambiguity code', 
		       'width' => '15%', 'align' => 'center' },
		      {'key' => 'Alternate residues', 
		       'title' => 'Alternate residues', 
		       'width' => '20%', 'align' => 'center' },
		     ];

  foreach my $residue (@$snps){	
    my $type = $residue->{'type'} eq 'snp' ? "Non-synonymous" : ($residue->{'type'} eq 'syn' ? 'Synonymous': ucfirst($residue->{'type'}));
    $counter++;
    next if !$residue->{'allele'};
    my $snp_id = $residue->{'snp_id'};
    my $source = $residue->{'snp_source'} ? "&source=".$residue->{'snp_source'} : "";
    my $tmp_array = {
		     'Residue' 	=> $counter,	    	    	  
		     'SNP type' => $type,
		     'SNP ID'   => qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/snpview?snp=$snp_id$source">$snp_id</a>),
		     'LDview'   => qq(<a href="/$ENV{'ENSEMBL_SPECIES'}/ldview?snp=$snp_id$source">$snp_id</a>),
		     'Alleles'     => $residue->{'allele'},
		     'Ambiguity code'     => join('', @{$residue->{'ambigcode'}||[]}),
		     'Alternate residues'     => $residue->{'pep_snp'} ? $residue->{'pep_snp'} : '-'};
    push @snp_list, $tmp_array;
  }
  return ( $table_header, \@snp_list);		
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

  my $transcript = $self->param( "transcript" );
  my $peptide = $self->param( "peptide" );
  my $db      = $self->param( "db" );
  my $param_str = "db=$db";
  $transcript and $param_str .= "&transcript=$transcript";
  $peptide    and $param_str .= "&peptide=$peptide";

  my $conf_submit = qq(
<FORM name="dasConfigForm" action="dasconfview" method="post">
  <INPUT type="hidden" name="conf_script" value="protview">
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

  my $label = sprintf( $a_tmpl, "/Docs/gene_das.html",,"ProteinDAS" );
  $label .= " Sources";

  my $hidden = '';
  foreach my $param( "peptide", "translation", "db" ){
    my $val = $self->param($param) || next;
    $hidden .= sprintf( $hidden_tmpl, $param, $val );
  }

  my @checks = ();
  my $data = $self->DataObj;
  my $das_attribute_data = $data->get_das_attributes
				( "name", "authority", "label", "active" ) || return ();
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
    ("name", "authority","active" );

  foreach my $source( @$das_attribute_data ){
    $source->{active} || next;
    my $source_nm = $source->{name};
    my $source_lab = $source_nm;
    if( my $ln = $source->{authority} ){
      $source_lab = sprintf( $link_tmpl, $ln, $source_nm, $source_nm );
    }
    my $label = "ProteinDAS";
    push( @table_data, $label.": ". $source_lab,'' );

    my @features = $data->get_das_annotation_by_name($source_nm,'global');

    my @rhs_rows;

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
	
	$note=~s|((\S+?):navigation://(\S+))|
	  <A href="protview?gene=$3" >[$2]</A>|ig;
	#$note=~s|([^"])(navigation://\S+)([^"])|
	#  $1<A href="$2">$2</A>$3|ig;
      }

      push( @rhs_rows, sprintf( $row_tmpl, 
				#$feature->{-source} || '&nbsp;',
				$feature->das_type || '&nbsp;',
				$id                || '&nbsp;',
				$note              || '&nbsp;' ) );
    }
    $table_data[-1] = sprintf( $table_tmpl, 
			       join($space_row, @rhs_rows ) );
  }

  return (@table_data);   
}

#======================================================================

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


=head2 similarity_matches

 Arg[1]      : (optional) String
               Label
 Example     : $pepdata->renderer->similarity_matches
 Description : Renders similarity matches for transcript in two_col_table format
 Return type : Key / value pair - label and HTML

=cut

sub similarity_matches {
    my $self = shift;
    my $label = shift || 'Similarity Matches';
    my $transl = $self->DataObj->translation;     
    my $data = $self->DataObj();
    # Check cache
    unless ($transl->{'similarity_links'}) {
        my @similarity_links = @{$data->get_similarity_hash($transl)};   

        return unless (@similarity_links);

        # sort links
        $self->_sort_similarity_links(@similarity_links);
    }

    my %links = %{$transl->{'similarity_links'}};
    return unless %links;

    my $db = $data->get_db();
    my $entry = $data->gene_type || 'Ensembl';
    # add table call here
    my $html;
    unless ($data->species_defs->ENSEMBL_SITETYPE eq 'Vega') {
        $html = qq(
                <b>This $entry entry corresponds to the following database identifiers:</b><br />);
    }
    $html .= qq(<table class="hidden">);
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
            $html .= qq(</td></tr>);
        }
    }   
    $html .= qq(</table>); 
    return ($label , $html);
}

=head2 _sort_similarity_links

 Arg[1]      : none
 Example     : $pepdata->renderer->_sort_similarity_links
 Description : sorts the similarity matches
 Return type : hashref of similarity matches

=cut

sub _sort_similarity_links{
    my $self = shift;
    my $transl = $self->DataObj->translation;
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
            'uniprot/sptrembl'      => 'UniProt/TrEMBL',
            'uniprot/swissprot'     => 'UniProt/Swiss-Prot',
            'pubmed'                => 'Sequence Publications',
        );
                       
    foreach my $type (sort @similarity_links) { 
        my $link = "";
        my $join_links = 0;
        my $externalDB = $type->database();
        my $display_id = $type->display_id();
        my $primary_id = $type->primary_id();

        # remove all orthologs  
        next if ($type->status() eq 'ORTH');
 
        # ditch medline entries - redundant as we also have pubmed
        next if lc($externalDB) eq "medline";
    
        # Ditch celera genes from FlyBase
        next if ($externalDB =~ /^flybase/i && $display_id =~ /^CG/ );

        # remove internal links to self and transcripts
        next if $externalDB eq "Vega_gene";
        next if $externalDB eq "Vega_transcript";
        next if $externalDB eq "Vega_translation";

        if( $externalDB eq "GO" ){ #&& $data->database('go')){
            push @{$transl->{'go_links'}} , $display_id;
            next;   
        } 
	elsif ($externalDB eq "GKB") {
            my ($key, $primary_id) = split ':', $display_id;
            push @{$transl->{'GKB_links'}{$key}} , $type ;
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
                $transl->stable_id,
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
        my $display_name = $nice_names{lc($externalDB)} || ($externalDB =~ s/_/ /g, $externalDB)  ;
        push (@{$links{$display_name}}, $link);         
    }
    $transl->{'similarity_links'} = \%links ;
    return $transl->{'similarity_links'};
}



1; 
