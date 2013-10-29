package EnsEMBL::Web::Component::Search;

use EnsEMBL::Web::Component;
our @ISA = qw(EnsEMBL::Web::Component);

use strict;
use warnings;
no warnings "uninitialized";
use Data::Dumper;

use constant RESULTS_COUNT_IDX => 1;

sub no_results {
    my ($panel, $object) = @_;
    $panel->print(qq(<p>Your search returned no results</p>));
}

sub search_instructions {
    my ($panel, $object) = @_;
    $panel->print(
        qq(<p>You must enter a search term in the box at the top of the page</p>)
    );
}

sub results {
    my ($panel, $object) = @_;
    my $result_ref = $object->Obj->{'results'};
    
    # Sort by reverse counts
    for my $search_index (sort {
        $result_ref->{$b}->[RESULTS_COUNT_IDX] <=>
            $result_ref->{$a}->[RESULTS_COUNT_IDX]
    } keys %$result_ref) {
        my ($results, $count) = @{ $result_ref->{$search_index} };
        $panel->print(
            "<h3>Search results for $search_index</h3><p>$count entries matched your search strings.</p><ol>"
        );
        foreach my $result (@$results) {
            $panel->printf(qq(<li><strong>%s:</strong> <a href="%s">%s</a>),
                $result->{'subtype'}, $result->{'URL'}, $result->{'ID'});
            if ($result->{'URL_extra'}) {
                foreach my $E (@{ [ $result->{'URL_extra'} ] }) {
                    $panel->printf(qq( [<a href="%s" title="%s">%s</a>]),
                        $E->[2], $E->[1], $E->[0]);
                }
            }
            if ($result->{'desc'}) {
                $panel->printf(qq(<br />%s), $result->{'desc'});
            }
            $panel->print('</li>');
        }
        $panel->print('</ol>');
    }
    return 1;
}

1;
