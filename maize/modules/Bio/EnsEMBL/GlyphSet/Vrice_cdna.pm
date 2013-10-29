package Bio::EnsEMBL::GlyphSet::Vrice_cdna;
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
	        'text'      => 'Clustered ESTs',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vrice_cdna','col_cluster'),
		'absolutey' => 1,
    }));
    $self->label2(new Sanger::Graphics::Glyph::Text({
		'text'      => 'Rice cDNAs',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vrice_cdna','col_cdna'),		
		'absolutey' => 1,
    }));
}

sub _init {
    my ($self) = @_;
    my $Config = $self->{'config'};
    my $chr      = $self->{'container'}->{'chr'};
   	my $cluster_col = $Config->get( 'Vrice_cdna','col_cluster' );
   	my $cdna_col = $Config->get( 'Vrice_cdna','col_cdna' );
	
	
    my $cdnas = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'rice_cdna_clones');
    my $clusters = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'OStug.trim.crop');

    return unless $cdnas->size() && $clusters->size();


   	$clusters->scale_to_fit( $Config->get( 'Vrice_cdna', 'width' ) );
	$clusters->stretch(0);
	my $Hscale_factor = $cdnas->{'_biggest_value'} / $clusters->{'_biggest_value'};
   	$cdnas->scale_to_fit( $Config->get( 'Vrice_cdna', 'width' ) * $Hscale_factor );	
	$cdnas->stretch(0);
		

	my @clusters = $clusters->get_binvalues();
	my @cdnas = $cdnas->get_binvalues();	

    foreach (@clusters){
       my $cdna = shift @cdnas;	
	    my $g_x = new Sanger::Graphics::Glyph::Line({
			'x'      => ($_->{'chromosomeend'}+$_->{'chromosomestart'})/2,
			'Y'      => 0,
			'width'  => 0,
			'height' => $cdna->{'scaledvalue'},
			'colour' => $cdna_col,
			'absolutey' => 1,
		});
	    $self->push($g_x);
		$g_x = new Sanger::Graphics::Glyph::Rect({
			'x'      => $_->{'chromosomestart'},
			'y'      => 0,
			'width'  => $_->{'chromosomeend'}-$_->{'chromosomestart'},
			'height' => $_->{'scaledvalue'},
			'bordercolour' => $cluster_col,
			'absolutey' => 1,
			'href'   => "/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$_->{'chromosomestart'}&vc_end=$_->{'chromosomeend'}"
		});
	    $self->push($g_x);
	}
}

1;
