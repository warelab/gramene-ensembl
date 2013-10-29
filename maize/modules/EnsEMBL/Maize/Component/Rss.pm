package EnsEMBL::Maize::Component::Rss;

use strict;
use EnsEMBL::Web::Component;
use EnsEMBL::Web::Form;

use Readonly;

my @chromosomes = map { +{ 'name' => "Chromosome $_", 'value' => $_ } } (1 .. 10);
#my @data_types = { "value" =>"bacs", "name"=>"Sequenced BACs", 'checked' => 1
#	       };
Readonly my $FORM_NAME     => 'rss_form';
Readonly my %FORM_ELEMENTS => (
    'chromosome' => {
        'type'  => 'DropDown',
        'required' => 'yes',
        'label' => 'Chromosome:',
        'name'  => 'chromosome',
	'values' => \@chromosomes,
        'order' => 1,
	'select' => 'select',
#	'spanning' => 'inline',
    },
    'start' => {
        'type'     => 'String',
        'required' => 'yes',
        'label'    => "FPC Start Coordinate (base pairs):",
        'name'     => 'start',
        'order'    => 2,
    },
    'end' => {
        'type'     => 'String',
        'required' => 'yes',
        'label'    => "FPC End Coordinate (base pairs):",
        'name'     => 'end',
        'order'    => 3,
    },
#    'data' => {
#	'type'     => 'MultiSelect',
#	'name'     => 'data',
#	'label'    => "Data Types for RSS Notification",
#	'values'   => \@data_types,
#	'order'    => 4,
#    },
    'action' => {
        'type'  => 'Hidden',
        'name'  => 'action',
        'value' => 'submit',
        'order' => 5,
        'required' => 'no',
    },
);

sub show_form {
    my $panel    = shift;
    my ($object) = @_;
    my $script   = $object->script;
    my $form
        = EnsEMBL::Web::Form->new($FORM_NAME, "/@{[$object->species]}/$script",
        'get');

    $form->add_element(
		       'type'  => 'Information',
		       'value' => '<br></br><p>Note: BAC Notification is provided as RSS (Really Simple Syndication), an XML-based format for information distribution. You can subscribe to any region of the maize genome and receive BAC, gene prediction, and marker updates via your favorite feed reader. The XML generated can also be parsed by RSS enabled browsers such as Firefox 1.0/2.0, Safari 2.0, and Internet Explorer 7.0.</p>'
		       );  
    
    my $parameters = _initialize_parameters($object);

    for my $field (
        sort { $parameters->{$a}->{'order'} <=> $parameters->{$b}->{'order'} }
        keys %$parameters)
    {
        my %element_arguments = ();
        for my $key (keys %{ $parameters->{$field} }) {
            $element_arguments{$key} = $parameters->{$field}->{$key};
        }
        $form->add_element(%element_arguments);
    }

    $form->add_button('submit', 'RSS');
    $form->add_button('reset',  'Reset');

    $panel->print($form->render);
    return 1;
}

sub process {
    my $panel = shift;
    my ($object) = @_;
    
    $object->rss_xml;
    
    return 1;
}

sub is_required {
    my ($param) = @_;
    return 0 unless defined $FORM_ELEMENTS{$param};
    return ($FORM_ELEMENTS{$param}->{'required'} eq 'yes');
}

sub _initialize_parameters {
    my ($object) = @_;
    my %parameters = ();
    for my $field (keys %FORM_ELEMENTS) {
        $parameters{$field} = {};
        for my $key (keys %{ $FORM_ELEMENTS{$field} }) {
            $parameters{$field}->{$key} = $FORM_ELEMENTS{$field}->{$key};
        }
        $parameters{$field}->{'value'} ||= $object->param($field);
    }
    return \%parameters;
}

1;
