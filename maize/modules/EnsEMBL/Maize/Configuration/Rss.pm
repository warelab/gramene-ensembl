package EnsEMBL::Maize::Configuration::Rss;

use strict;
use EnsEMBL::Web::Configuration;
use Readonly;

use EnsEMBL::Maize::Component::Rss;

our @ISA = qw( EnsEMBL::Web::Configuration );

Readonly my $PAGE_TITLE => 'Maize Data Notification via RSS';

sub configure_rss {
    my $self       = shift;
    my $object     = $self->{object};
    my $print_form = 1;
    my $panel      = new EnsEMBL::Web::Document::Panel(
#						       'caption' => $PAGE_TITLE,
						       'object'  => $self->{object}
						       );
    if ($object->param('action') eq 'submit'){
	$panel->add_components(
			       qw(process_form EnsEMBL::Maize::Component::Rss::process)
			       );
	$print_form = 0;
    }
			       
    if ($print_form){
	$panel->add_components(
			       qw(show_form EnsEMBL::Maize::Component::Rss::show_form)
			       );
    }

    $self->{page}->content->add_panel($panel);
#    $self->{page}->title->set($PAGE_TITLE);
}
