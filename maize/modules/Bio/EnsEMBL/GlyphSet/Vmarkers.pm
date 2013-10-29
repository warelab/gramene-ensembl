package Bio::EnsEMBL::GlyphSet::Vmarkers;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet;
@ISA = qw(Bio::EnsEMBL::GlyphSet);
use Sanger::Graphics::Glyph::Rect;
use Sanger::Graphics::Glyph::Poly;
use Sanger::Graphics::Glyph::Text;
use Sanger::Graphics::Glyph::Line;

#warn "passat\n";

sub init_label {

  my ($self) = @_;
  my $lab_method;
  foreach my $lname( split( /\s+/, $self->my_config('logicname'), 2 ) ){

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

  my $self = shift;
  my $chr  = $self->{'extras'}->{'chr'} || $self->{'container'}->{'chr'};

  my $num_bins = 150;

  my $slice_adapt   = $self->{'container'}->{'sa'};
  my $density_adapt = $self->{'container'}->{'da'};
  my $slice = $slice_adapt->fetch_by_region('chromosome', $chr);
  my @lnames = split( /\s+/, $self->my_config('logicname') );
  my $width =  $self->my_config( 'width' );

  my %features;
  my $found = 0;
  foreach my $lname( @lnames ){
    my $fs = $density_adapt->fetch_Featureset_by_Slice
        ($slice,$lname,$num_bins,1);
    $fs->size || ( $features{$lname}=[] and next );
    $found++;
    $fs->scale_to_fit($width);
    $fs->stretch(1);
    $features{$lname} = $fs->get_all_binvalues() || [];
  }
  $found || return;

  my %colours = %{$self->my_config('colours')||{}};
  my $nudge  = -1;
  my $bin_size = int( $slice->length / $num_bins );
  for( my $i=0; $i<$num_bins; $i++ ){
    foreach my $lname( @lnames ){
      my $dfeat = $features{$lname}->[$i] || next;
      my $g_x = new Sanger::Graphics::Glyph::Line({
        'x'      =>  ( $i * $bin_size ) + ( $nudge>0 ? $bin_size/2 : 0 ),
        'Y'      => 0,
        'width'  => 0,
        'height' => $dfeat->scaledvalue, 
        'colour' => $colours{$lname} || 'black',
        'absolutey' => 1,
      });
      $self->push($g_x);
      $nudge *=-1;
    }
  }
  return 1;
}

1;
