package Bio::EnsEMBL::GlyphSet::Vmaize2;
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
	        'text'      => 'Maizedb Cornsensus',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vmaize2','col_cornsensus'),
		'absolutey' => 1,
    }));
    $self->label2(new Sanger::Graphics::Glyph::Text({
		'text'      => 'Maize ESTs',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vmaize2','col_est'),		
		'absolutey' => 1,
    }));
}

sub _init {
    my ($self) = @_;
    my $Config = $self->{'config'};
    my $chr      = $self->{'container'}->{'chr'};
   	my $cornsensus_col = $Config->get( 'Vmaize2','col_cornsensus' );
   	my $est_col = $Config->get( 'Vmaize2','col_est' );
	
	
    my $ests = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'maize_est');
    my $cornsensuss = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'maizedb_cornsensus');

    return unless $ests->size() && $cornsensuss->size();

	my $biggest = $cornsensuss->{'_biggest_value'} > $ests->{'_biggest_value'}
	            ? $cornsensuss->{'_biggest_value'} : $ests->{'_biggest_value'};
   	$ests->scale_to_fit( $Config->get( 'Vmaize1', 'width' )*$ests->{'_biggest_value'}/$biggest);
	$ests->stretch(0);
   	$cornsensuss->scale_to_fit( $Config->get( 'Vmaize1', 'width' )*$cornsensuss->{'_biggest_value'}/$biggest);
	$cornsensuss->stretch(0);


	my @cornsensuss = $cornsensuss->get_binvalues();
	my @ests = $ests->get_binvalues();	

    foreach (@cornsensuss){
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
			'bordercolour' => $cornsensus_col,
			'absolutey' => 1,
			'href'   => "/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$_->{'chromosomestart'}&vc_end=$_->{'chromosomeend'}"
		});
	    $self->push($g_x);
	}
}

1;
