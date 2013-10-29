package Bio::EnsEMBL::GlyphSet::Vmaize1;
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
	        'text'      => 'Maize BACends',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vmaize1','col_bacend'),
		'absolutey' => 1,
    }));
    $self->label2(new Sanger::Graphics::Glyph::Text({
		'text'      => 'Maize Hicot Methyl Clusters',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vmaize1','col_hicot'),		
		'absolutey' => 1,
    }));
}

sub _init {
    my ($self) = @_;
    my $Config = $self->{'config'};
    my $chr      = $self->{'container'}->{'chr'};
   	my $bacend_col = $Config->get( 'Vmaize1','col_bacend' );
   	my $hicot_col = $Config->get( 'Vmaize1','col_hicot' );
	
	
    my $hicots = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'tigr_methyl_hicot_clusters');
    my $bacends = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'maize_bacends');

    return unless $hicots->size() && $bacends->size();


	my $biggest = $hicots->{'_biggest_value'} > $bacends->{'_biggest_value'}
	            ? $hicots->{'_biggest_value'} : $bacends->{'_biggest_value'};
   	$bacends->scale_to_fit( $Config->get( 'Vmaize1', 'width' )*$bacends->{'_biggest_value'}/$biggest);
	$bacends->stretch(0);
   	$hicots->scale_to_fit( $Config->get( 'Vmaize1', 'width' )*$hicots->{'_biggest_value'}/$biggest);
	$hicots->stretch(0);
		

	my @bacends = $bacends->get_binvalues();
	my @hicots = $hicots->get_binvalues();	

    foreach (@bacends){
       my $hicot = shift @hicots;	
	    my $g_x = new Sanger::Graphics::Glyph::Line({
			'x'      => ($_->{'chromosomeend'}+$_->{'chromosomestart'})/2,
			'Y'      => 0,
			'width'  => 0,
			'height' => $hicot->{'scaledvalue'},
			'colour' => $hicot_col,
			'absolutey' => 1,
		});
	    $self->push($g_x);
		$g_x = new Sanger::Graphics::Glyph::Rect({
			'x'      => $_->{'chromosomestart'},
			'y'      => 0,
			'width'  => $_->{'chromosomeend'}-$_->{'chromosomestart'},
			'height' => $_->{'scaledvalue'},
			'bordercolour' => $bacend_col,
			'absolutey' => 1,
			'href'   => "/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$_->{'chromosomestart'}&vc_end=$_->{'chromosomeend'}"
		});
	    $self->push($g_x);
	}
}

1;
