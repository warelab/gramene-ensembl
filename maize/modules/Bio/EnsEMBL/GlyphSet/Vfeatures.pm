package Bio::EnsEMBL::GlyphSet::Vfeatures;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet;
@ISA = qw(Bio::EnsEMBL::GlyphSet);
use Sanger::Graphics::Glyph::Rect;
use Sanger::Graphics::Glyph::Poly;
use Sanger::Graphics::Glyph::Text;
use Sanger::Graphics::Glyph::Line;

use Data::Dumper;

sub init_label {
  my ($self) = @_;
  my @logic_names = split( '\s', $self->my_config('logicname') );
  my %labels  = %{$self->my_config('labels')  || {}};
  my %colours = %{$self->my_config('colours') || {}};
  my $col     = $self->my_config('col');
  foreach my $method( 'label', 'label2' ){
    my $lname = shift @logic_names || last;
    $self->$method( $self->label(new Sanger::Graphics::Glyph::Text({
      'text'      => $labels{$lname}  || $lname,
      'colour'    => $colours{$lname} || $col || 'black',
      'font'      => 'Small',
      'absolutey' => 1,
    })));
  }
}

sub _init {
  my ($self) = @_;
  my $Config = $self->{'config'};
  my $chr    = $self->{'extras'}->{'chr'} || $self->{'container'}->{'chr'};

  my @logic_names = split( '\s', $self->my_config('logicname') );
  my %colours = %{$self->my_config('colours') || {}};
  my $col     = $self->my_config('col');

#  my $genes_col = $Config->get( 'Vgenes','col_genes' );
#  my $known_col = $Config->get( 'Vgenes','col_known' );

  my $slice_adapt   = $self->{'container'}->{'sa'};
  my $density_adapt = $self->{'container'}->{'da'};

  my $chr_slice = $slice_adapt->fetch_by_region('chromosome', $chr);

  my %features;
  foreach my $lname( @logic_names ){
    $features{$lname} = $density_adapt->fetch_Featureset_by_Slice
        ($chr_slice, $lname, 150, 1);
  }

  my $v_offset = $Config->container_width() - ($chr_slice->length() || 1);

#  return unless scalar(grep{$_->size} values %features);
  return unless $features{$logic_names[0]}->size;

  my $width = $self->my_config( 'width' );
  map{ $_->scale_to_fit($width); $_->stretch(0) } values %features;

  my @feats1 = @{$features{$logic_names[0]}->get_all_binvalues()};
  #my @feats2 = @{$features{$logic_names[1]}->get_all_binvalues()};

  foreach (@feats1){
    #my $known_gene = shift @known_genes;  
    #my $g_x = new Sanger::Graphics::Glyph::Rect({
    #  'x'      => $v_offset + $known_gene->start,
    #  'y'      => 0,
    #  'width'  => $known_gene->end - $known_gene->start,
    #  'height' => $known_gene->scaledvalue,
    #  'colour' => $known_col,
    #  'absolutey' => 1,
    #});
    #$self->push($g_x);
    my $g_x = new Sanger::Graphics::Glyph::Rect({
      'x'      => $v_offset + $_->start,
      'y'      => 0,
      'width'  => $_->end - $_->start,
      'height' => $_->scaledvalue,
      'bordercolour' => $col,
      'absolutey' => 1,
      'href'   => "/@{[$self->{container}{_config_file_name_}]}/contigview?chr=$chr&vc_start=$_->start&vc_end=$_->end"
    });
    $self->push($g_x);
  }
  
}

1;
