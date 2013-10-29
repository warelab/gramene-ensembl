package Bio::EnsEMBL::GlyphSet::Vtigr;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet;
@ISA = qw(Bio::EnsEMBL::GlyphSet);
use Sanger::Graphics::Glyph::Rect;
use Sanger::Graphics::Glyph::Poly;
use Sanger::Graphics::Glyph::Text;
use Sanger::Graphics::Glyph::Line;

sub init_label {
  my ($self) = @_;
  my $lab_method = undef;
  foreach my $lname( split( /\s+/, $self->my_config('logicname'),2 )){
    my $text = ${$self->my_config('track_labels')||{}}{$lname} || $lname;
    my $col  = ${$self->my_config('colours')||{}}{$lname} || 'black';
    $lab_method = $lab_method ? 'label2' : 'label';
    $self->$lab_method(new Sanger::Graphics::Glyph::Text({
      'text'      => $text,
      'font'      => 'Small',
      'colour'    => $col,
      'absolutey' => 1,
    }));
  }
  return 1;
}

sub _init {
  my ($self) = @_;
  my $chr    = $self->{'extras'}->{'chr'} || $self->{'container'}->{'chr'};
  
  my @logic_names = split( '\s', $self->my_config('logicname') );
  
  my $genes_col = $self->my_config( 'col_genes');
  my $pctgc_col = $self->my_config( 'col_gc' );
  
  my $slice_adapt   = $self->{'container'}->{'sa'};
  my $density_adapt = $self->{'container'}->{'da'};
  
  my $chr_slice = $slice_adapt->fetch_by_region('chromosome', $chr);
  
  my $num_bins = 150;
  my $pctgc = $density_adapt->fetch_Featureset_by_Slice
      ($chr_slice, 'PercentGC', $num_bins, 1);
  my $genes = $density_adapt->fetch_Featureset_by_Slice
      ($chr_slice, 'Density_tigr_gene', $num_bins, 1);
  
  #return unless $pctgc->size() && $genes->size();
  if( ! $pctgc->size and ! $genes->size ){ return }
  
  my $WIDTH =  $self->my_config( 'width' );

  $genes->scale_to_fit($WIDTH);
  $pctgc->scale_to_fit($WIDTH); 
  $genes->stretch(0);
  $pctgc->stretch(0);
  
  
  my @genes = @{ $genes->get_all_binvalues() || [] };
  my @pctgc = @{ $pctgc->get_all_binvalues() || [] };

  my $old_x = undef;
  my $old_y = undef;

  for( my $i=0; $i<$num_bins; $i++ ){
    
    if( my $gene = $genes[$i] ){
      my $g_x ;
      my( $gstart, $gend ) = ( $gene->start, $gene->end );
      $g_x = new Sanger::Graphics::Glyph::Rect({
        'x'      => $gstart,
        'y'      => 0,
        'width'  => $gend - $gstart,
        'height' => $gene->scaledvalue,
        'bordercolour' => $genes_col,
        'absolutey' => 1,
        'href'   => "/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$gstart&vc_end=$gend"
          });
      $self->push($g_x);
    }

    if( my $pctgc = $pctgc[$i] ){					
      my $new_x = ($pctgc->end + $pctgc->start)/2;
      my $new_y = $pctgc->scaledvalue;
      if(defined $old_x) {      
        my $g_x = new Sanger::Graphics::Glyph::Line({
          'x'      => $old_x,
          'y'      => $old_y,
          'width'  => $new_x-$old_x,
          'height' => $new_y-$old_y,
          'colour' => $pctgc_col,
          'absolutey' => 1,
        });			
        $self->push($g_x);
      }
      $old_x = $new_x;
      $old_y = $new_y;
    }
  }
}

1;
