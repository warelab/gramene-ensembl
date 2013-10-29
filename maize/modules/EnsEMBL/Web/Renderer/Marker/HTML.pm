package EnsEMBL::Web::Renderer::Marker::HTML;

=head1 NAME

EnsEMBL::Web::Renderer::Marker::HTML.pm 

=head1 SYNOPSIS

This object creates HTML to be output to the HTML output object

=head1 DESCRIPTION

    $marker_renderer = $marker_data->renderer;
	$marker_renderer->outputGenericmarkerTable();		 
    $marker_renderer->marker_table_title;
	
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
use EnsEMBL::Web::Renderer::Marker;

@ISA = qw( EnsEMBL::Web::Renderer::HTML EnsEMBL::Web::Renderer::Marker );

=head2 outputGenericMarkerview

 Arg[1]      : none
 Example     : $markerdata->renderer->outputGenericMarkerview
 Description : Wrapper for the standard markerview page
 Return type : none

=cut

sub outputGenericMarkerview() {
  my $self = shift;
  my $marker = $self->param('marker');
  $self->Output->generic_table_title("Chromosome Map Marker $marker");
  $self->Output->print_two_col_table(
     $self->markers,
	 $self->marker_location,
         $self->marker_type, #Maize-specific
	 $self->marker_synonyms,
	 $self->marker_primers,
    );
	$self->Output->printHTML($self->Output->spreadsheet_table($self->map_locations));
}

sub markers{
	my $self = shift;
	my $Data =$self->DataObj;
	my $label = 'Marker Source';
	my $html ;
	my $important_synonyms = $Data->main_MarkerSynonyms;
	
	return unless $important_synonyms;
	
	foreach my $synonym (@$important_synonyms){
		my $db = $synonym->source; 
		my $id = $synonym->name;
		my $url = $self->get_ExtURL($db, $id) ;
		$id = $self->url($id, $url) if ($url)	;
		$html .= "$id &nbsp;&nbsp;( <b>Database:</b> $db ) <br />";
	}
	return unless ($html =~ /\w+/g) ;
	return ($label, $html);
}

sub marker_type{
  my $self = shift;
  my $type = $self->DataObj->markerType || return;
  return( "Marker Type","<b>$type</b>" );
}

sub marker_location{
	my $self = shift;
	my $Data =$self->DataObj;
	my $label = 'Marker Location';	
	my $marker = $self->param('marker');
	my $marker_feats = $Data->MarkerFeatures;
	my $count = scalar(@$marker_feats);
	my $MAX_MAP_WEIGHT = 100;
	my $html ;
	
	foreach my $feature (@$marker_feats){
          my $name  = $feature->contig->chr_name;
          my $start = $feature->start;
          my $end   = $feature->end;
          if ($count > $MAX_MAP_WEIGHT) {
            $html = "<p><b>$marker is currently mapped to $count different EnsEMBL locations</b><p>";
            last;
          } elsif (@$marker_feats){
            $html .= qq(<b>Location</b>  <a href="\/). $self->Input->species .qq(/contigview?&chr=$name&vc_start=). ($start-10000) .qq(&vc_end=). ($end+10000) .qq(">$start - $end bps</a> <b> on chromosome</b> <a href="\/). $self->Input->species .qq(/mapview?&chr=$name">$name</a> &nbsp;&nbsp;&nbsp;[<a href="\/). $self->Input->species .qq(/exportview?type=basepairs&chr=$name&bp_start=$start&bp_end=$end">Export Data</a>]<br />);
          } else {
            $html = "<p><B>Marker $marker is not mapped to the assembly in the current Ensembl database.</B> <p>";
          }		
	}
	$html = qq(<p><B>Marker $marker is not mapped to the assembly in the current Ensembl database.</B> <p>) unless $html;
	return ($label, $html);
}

sub marker_synonyms{
	my $self = shift;
	my $Data =$self->DataObj;
	my $label = 'Marker Synonyms';	
	my $other_synonyms = $Data->other_MarkerSynonyms;
	
	return unless $other_synonyms;
	
	my $html = qq(<table border="0" class="hidden" >);
	my %source;
	my %counter;
	my $max_count = 0 ;
	foreach my $synonym (@$other_synonyms){
		my $src_name = ucfirst($synonym->source) || 'Other';
		$counter{$src_name}++;
		push @{$source{$src_name}}, $synonym->name ;	
		$max_count = $counter{$src_name} if ($max_count < $counter{$src_name} );	
	    }

	my $cols = ($max_count / 5) ;
	$cols++ if $max_count % 5;

	foreach my $key (sort keys %source){
	    $html .= qq(<tr><td valign="top"><b>$key : </b></td>);
	    for (1..$cols) {
		my @list = splice (@{$source{$key}} ,0 ,5);
		$html .= qq(<td>);
		foreach my $id (@list) {
		    my $url = $self->get_ExtURL($key, $id);
		    if ($url) {
			$html .= qq(<a href=).$url.qq(>).$id.qq(</a>&nbsp;<br />);
		    }
		    else {
			$html .= $id.qq(&nbsp;<br />);
		    }
		}
		$html .= qq(</td>);
	    }
	    $html .= qq(</tr>);
	}	
	$html .= "</table>";
	return ($label, $html);
    }


sub marker_primers{
	my $self = shift;
	my $Data =$self->DataObj;
	my $label = 'Marker Primers';	
	my $marker = $self->param('marker');
	my $marker_obj = $Data->marker;
	my $l = $marker_obj->left_primer;
	my $r = $marker_obj->right_primer;
	my $min_psize = $marker_obj->min_primer_dist;
	my $max_psize = $marker_obj->max_primer_dist;
	my $product_size;
	my $html;

	if(!$min_psize) {
  		$product_size = "&nbsp";
	} elsif($min_psize == $max_psize) {
  		$product_size = "$min_psize";
	} else {
  		$product_size = "$min_psize - $max_psize";
	}

	if (! $l){
  		$html =
			 "<p><b>Marker $marker primers are not in the database</b></p>\n";
	} else {
		$l =~ s/([\.\w]{30})/$1<BR>/g ;
		$r =~ s/([\.\w]{30})/$1<BR>/g ;
  		$html = 
		 qq(<table class="hidden" width="100%"><tr align="left">
		     <td width="20%"><b>Expected Product Size</b></td>
			 <td><b>Left Primer</b></td>
			 <td><b>Right Primer</b></td>
			 </tr><tr>
			 <td align="center">$product_size</td>
			 <td align="left">$l</td>
			 <td align="left">$r</td>
			 </tr></table>);
  	}

	return ($label, $html);
}

sub map_locations{
	my $self = shift;
	my $Data =$self->DataObj;
	my $marker_map_locations = $Data->get_map_locations;
	return undef unless @$marker_map_locations;
	
	my $marker = $self->param('marker');
	$self->Output->generic_table_title("Marker $marker Map Locations");
	my $table_header = [{'key' => 'map', 'title' => 'Map Name', 'width' => '25%', 'align' => 'center' },
						{'key' => 'syn', 	'title' => 'Synonym', 'width' => '25%', 'align' => 'center' }, 
						{'key' => 'Chr', 'title' => 'Chromosome', 'width' => '15%', 'align' => 'center' },
						{'key' => 'pos', 	'title' => 'Position', 'width' => '25%', 'align' => 'center' },
						{'key' => 'lod', 'title' => 'LOD Score', 'width' => '35%', 'align' => 'center' },
						];	
						
	return ($table_header, $marker_map_locations, {'rows' => [qw(background1 background3)]});
}

1;	
