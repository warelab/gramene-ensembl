package Bio::EnsEMBL::GlyphSet::Vsorghum1;
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
    my $Config = $self->{'config'};	
    $self->label(new Sanger::Graphics::Glyph::Text({
	        'text'      => 'Klein GSS',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vsorghum1','col_klein'),
		'absolutey' => 1,
    }));
    $self->label2(new Sanger::Graphics::Glyph::Text({
		'text'      => 'Sorghum ESTs',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vsorghum1','col_est'),		
		'absolutey' => 1,
    }));
}

sub _init {
    my ($self) = @_;
    my $Config = $self->{'config'};
    my $chr      = $self->{'container'}->{'chr'};
   	my $klein_col = $Config->get( 'Vsorghum1','col_klein' );
   	my $est_col = $Config->get( 'Vsorghum1','col_est' );
	
	
    my $ests = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'sorghum_est');
    my $kleins = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'klein_gss_reads');

    return unless $ests->size() && $kleins->size();

	my $biggest = $ests->{'_biggest_value'} > $kleins->{'_biggest_value'}
	            ? $ests->{'_biggest_value'} : $kleins->{'_biggest_value'};
   	$kleins->scale_to_fit( $Config->get( 'Vmaize1', 'width' )*$kleins->{'_biggest_value'}/$biggest);
	$kleins->stretch(0);
   	$ests->scale_to_fit( $Config->get( 'Vmaize1', 'width' )*$ests->{'_biggest_value'}/$biggest);
	$ests->stretch(0);

	my @kleins = $kleins->get_binvalues();
	my @ests = $ests->get_binvalues();	

    foreach (@kleins){
       my $est = shift @ests;	
	    my $g_x = new Sanger::Graphics::Glyph::Line({
			'x'      => ($_->{'chromosomeend'}+$_->{'chromosomestart'})/2,
			'Y'      => 0,
			'width'  => 0,
			'height' => $est->{'scaledvalue'},
			'colour' => $est_col,
			'absolutey' => 1,
		});
	    $self->push($g_x);
		$g_x = new Sanger::Graphics::Glyph::Rect({
			'x'      => $_->{'chromosomestart'},
			'y'      => 0,
			'width'  => $_->{'chromosomeend'}-$_->{'chromosomestart'},
			'height' => $_->{'scaledvalue'},
			'bordercolour' => $klein_col,
			'absolutey' => 1,
			'href'   => "/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$_->{'chromosomestart'}&vc_end=$_->{'chromosomeend'}"
		});
	    $self->push($g_x);
	}
}

1;
