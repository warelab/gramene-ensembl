package Bio::EnsEMBL::GlyphSet::Vfgenesh;
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
#    $self->label(new Sanger::Graphics::Glyph::Text({
#	        'text'      => 'FgenesH(Known)',
#		'font'      => 'Small',
#		'colour'	=> $Config->get('Vfgenesh','col_known'),
#		'absolutey' => 1,
#    }));
    my $label = new Sanger::Graphics::Glyph::Text({
		'text'      => '% GC',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vfgenesh','col_gc'),
		'absolutey' => 1,
    });
    $self->label($label);
    $self->label2(new Sanger::Graphics::Glyph::Text({
		'text'      => 'FgenesH',
		'font'      => 'Small',
		'colour'	=> $Config->get('Vfgenesh','col_genes'),		
		'absolutey' => 1,
    }));
}

sub _init {
    my ($self) = @_;
warn("Fgenesh track");
    my $Config = $self->{'config'};
    my $chr      = $self->{'container'}->{'chr'};
   	my $genes_col = $Config->get( 'Vfgenesh','col_genes' );
#   	my $known_col = $Config->get( 'Vfgenesh','col_known' );
   	my $gc_col    = $Config->get( 'Vfgenesh','col_gc' );
	
	
#    my $known_genes = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'knfgenesh');
    my $gc 		    = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'gc');
    my $genes = $self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,'fgenesh');

#    return unless $known_genes->size() && $genes->size();
    return unless $gc->size() && $genes->size();


#   	$genes->scale_to_fit( $Config->get( 'Vfgenesh', 'width' ) );
#	$genes->stretch(0);
#	my $Hscale_factor = $known_genes->{'_biggest_value'} / $genes->{'_biggest_value'};
#   	$known_genes->scale_to_fit( $Config->get( 'Vfgenesh', 'width' ) * $Hscale_factor );	
#	$known_genes->stretch(0);
        
	my $max_genes=$genes->{'_biggest_value'};
	my $max_gc=$gc->{'_biggest_value'};
	my $MAX         = $max_genes > $max_gc ? $max_genes : $max_gc;
	$MAX ||= 1;
	my $WIDTH =  $Config->get( 'Vfgenesh', 'width' );
	$genes->scale_to_fit($max_genes/$MAX*$WIDTH);
	$genes->stretch(0);
	$gc->scale_to_fit($max_gc/$MAX*$WIDTH);
	$gc->stretch(0);

		

	my @genes = $genes->get_binvalues();
	my @gc = $gc->get_binvalues();
#	my @known_genes = $known_genes->get_binvalues();	

	my $old_x = undef;
	my $old_y = undef;
    foreach (@genes){
#       my $known_gene = shift @known_genes;	
       my $g_x ;
#	  $g_x = new Sanger::Graphics::Glyph::Rect({
#			'x'      => $known_gene->{'chromosomestart'},
#			'y'      => 0,
#			'width'  => $known_gene->{'chromosomeend'}-$_->{'chromosomestart'},
#			'height' => $known_gene->{'scaledvalue'},
#			'colour' => $known_col,
#			'absolutey' => 1,
#		});
#	    $self->push($g_x);
		$g_x = new Sanger::Graphics::Glyph::Rect({
			'x'      => $_->{'chromosomestart'},
			'y'      => 0,
			'width'  => $_->{'chromosomeend'}-$_->{'chromosomestart'},
			'height' => $_->{'scaledvalue'},
			'bordercolour' => $genes_col,
			'absolutey' => 1,
			'href'   => "/$ENV{'ENSEMBL_SPECIES'}/contigview?chr=$chr&vc_start=$_->{'chromosomestart'}&vc_end=$_->{'chromosomeend'}"
		});
	    $self->push($g_x);

		my $gcvalue = shift @gc;					
		my $new_x = ($gcvalue->{'chromosomeend'}+$gcvalue->{'chromosomestart'})/2;
		my $new_y = $gcvalue->{'scaledvalue'};
		if(defined $old_x) {

		    my $g_x = new Sanger::Graphics::Glyph::Line({
				'x'      => $old_x,
				'y'      => $old_y,
				'width'  => $new_x-$old_x,
				'height' => $new_y-$old_y,
				'colour' => $gc_col,
				'absolutey' => 1,
			});			
			$self->push($g_x);
		}
		$old_x = $new_x;
		$old_y = $new_y;
	}
}

1;
