package Maize::ClonepathTest;

=head1 NAME

Maize::ClonepathTest

=head1 SYNOPSIS

  use Maize::ClonepathTest

  
=head1 DESCRIPTION

Tests Maize::Clonpath

=cut

use strict;
use warnings;

use base qw(Test::Unit::TestCase);

use Maize::Clonepath;

=pod

=head2 new
    Constructor

=cut

sub new {
    my $self = shift()->SUPER::new(@_);
    return $self;
}

=pod

=head2 test_configuration
    Tests that configuration is loaded correctly

=cut

sub test_simple_configuration {
    my $self      = shift;
    my $clonepath = Maize::Clonepath->new;

    @ARGV = ("-h");

    my (%conf);

    $clonepath->load_configuration(\%conf);

    $self->assert(1, $conf{'help'}, "Help option should be set");
}

=pod

=head2 test_configuration_interpolates_constructor_variables
    Tests that a configuration option is evaluated based on constructor
    parameters.

=cut

sub test_configuration_interpolates_constructor_variables {
    my $self = shift;

    my $clonepath = Maize::Clonepath->new({ 'subst' => 'fusc' });

    @ARGV = ('-r=ob${subst}ate', "-b=Zea_mays");

    my (%conf);
    $clonepath->load_configuration(\%conf);

    $self->assert_equals('obfuscate', $conf{'registry_file'});
    $self->assert_equals('Zea_mays',  $conf{'bac_species'});
}

=pod

=head2 test_max_improved_sequence_length
    Tests that the maximum improved sequence length is as expected

=cut

sub test_max_improved_sequence_length {
    my $self = shift;

    my $clonepath = Maize::ClonepathStub->new;
    my $dbh       = $clonepath->get_clonepath_handle;
    $dbh->{mock_add_resultset} = [ [1000], [2000], [500] ];
    my $length = $clonepath->improved_sequence_max_length();
    $self->assert_equals(2000, $length);
}

1;

package Maize::ClonepathStub;

use base qw(Maize::Clonepath);

use DBI;

sub new {
    my $self = shift()->SUPER::new(@_);
    return $self;
}

=pod

=head2 _init_clonepath_handle
    Override subroutine to instead use DBD::Mock

=cut

sub _init_clonepath_handle {
    return DBI->connect('DBI:Mock:', '', '');
}

1;

=pod

=head1 AUTHOR

Shiran Pasternak E<lt>shiran@cshl.eduE<gt>

=cut

=head1 COPYRIGHT

Copyright (c) 2007 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
