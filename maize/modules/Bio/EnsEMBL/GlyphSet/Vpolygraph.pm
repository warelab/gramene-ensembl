package Bio::EnsEMBL::GlyphSet::Vpolygraph;
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
    #How to get >1 color in a label? can it be a Composite?
    my $label = new Sanger::Graphics::Glyph::Text({
		'text'      =>  $Config->get('Vpolygraph','what')->[0]{'label'},
		'font'      => 'Small',
		'colour'	=> $Config->get('Vpolygraph','what')->[0]{'col'},
		'absolutey' => 1,
    });
    my $label2 = new Sanger::Graphics::Glyph::Text({
		'text'      =>  $Config->get('Vpolygraph','what')->[1]{'label'},
		'font'      => 'Small',
		'colour'	=> $Config->get('Vpolygraph','what')->[1]{'col'},
		'absolutey' => 1,
    });
		
    $self->label(  $label  );
    $self->label2( $label2 );
}

sub _init {
    my ($self) 		= @_;
    my $Config 		= $self->{'config'};
    my $chr      	= $self->{'container'}->{'chr'};
    my $what		= $Config->get( 'Vpolygraph','what' );
    my @colors		= map  { $_->{'col'} } @$what;
	
    my @densities 	= map {
	$self->{'container'}->{'da'}->get_density_per_chromosome_type($chr,$_->{'type'}) } @$what; 

    print STDERR "missing a size\n" and
    return if grep { ! $_->size() } @densities;
    my @max= map { $_->{'_biggest_value'} } @densities;
    my $MAX=1;
       for(@max) { $MAX=$_ if $MAX<$_ }
    
    my $WIDTH=$Config->get( 'Vpolygraph', 'width' ) /$MAX ;
    my @binvals;
    for my $i (0..$#densities ) {
	$densities[$i]->scale_to_fit( $max[$i] * $WIDTH);
	$densities[$i]->stretch(0);
	$binvals[$i] = [ $densities[$i]->get_binvalues()];
    }

	my @old_x = (undef) x scalar(@densities);
	my @old_y = (undef) x scalar(@densities);
    foreach my $pt (0..$#{$binvals[0]} ) {
	for my $i (0..$#densities ) {
	    my $value = $binvals[$i][$pt];
		my $new_x = ($value->{'chromosomeend'}+$value->{'chromosomestart'})/2;
		my $new_y = $value->{'scaledvalue'};
		if(defined $old_x[$i]) {

		    my $g_x = new Sanger::Graphics::Glyph::Line({
				'x'      => $old_x[$i],
				'y'      => $old_y[$i],
				'width'  => $new_x-$old_x[$i],
				'height' => $new_y-$old_y[$i],
				'colour' => $colors[$i],
				'absolutey' => 1,
			});			
			$self->push($g_x);
		}
		$old_x[$i] = $new_x;
		$old_y[$i] = $new_y;
	}
   }
}

1;
