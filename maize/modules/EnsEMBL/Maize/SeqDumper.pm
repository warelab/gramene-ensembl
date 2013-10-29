
=head1 NAME - EnsEMBL::Maize::SeqDumper

=head1 SYNOPSIS

  $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper;

  #don't dump snps or repeats
  $seq_dumper->disable_feature_type('repeat');
  $seq_dumper->disable_feature_type('variation');

  #dump EMBL format to STDOUT
  $seq_dumper->dump($slice, 'EMBL');

  #dump GENBANK format to a file
  $seq_dumper->dump($slice, 'GENBANK', 'out.genbank');

  #dump FASTA format to a file
  $seq_dumper->dump($slice, 'FASTA', 'out.fasta');

=head1 DESCRIPTION

  A relatively simple and lite-weight flat file dumper for Ensembl slices.
  The memory efficiency could be improved and this is currently
  not very good for dumping very large sequences such as whole chromosomes.

=head1 CONTACT

  Contact the Ensembl development list with questions: <ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;

package EnsEMBL::Maize::SeqDumper;

use vars qw(@ISA);

use POSIX qw(ceil);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use EnsEMBL::Web::SeqDumper;

our @ISA = qw(EnsEMBL::Web::SeqDumper);

=head2 new

  Arg [1]    : none
  Example    : $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper->new;
  Description: Creates a new SeqDumper 
  Returntype : Bio::EnsEMBL::Utils::SeqDumper
  Exceptions : none
  Caller     : general

=cut

sub new {
    my ($class, $slice) = @_;

    my $self = $class->SUPER::new(@_);
    for my $feature (qw/improved mdr/) {
        $self->{'feature_types'}->{$feature} = 1;

    }
    return $self;
}

sub _dump_feature_table {
    my $self   = shift;
    my $slice  = shift;
    my $FH     = shift;
    my $FORMAT = shift;

    $self->SUPER::_dump_feature_table($slice, $FH, $FORMAT);

    my @ff = ($FH, $FORMAT);

    if ($self->is_enabled('improved') && $slice->can('get_all_MiscFeatures'))
    {
        for my $improved_sequence (
            @{ $slice->get_all_MiscFeatures('improved_sequence') })
        {
            $self->write(@ff, 'misc_feature',
                $self->features2location([$improved_sequence]));
            $self->write(@ff, '', '/note="Improved sequence."');
        }
    }
    if ($self->is_enabled('mdr') && $slice->can('get_all_SimpleFeatures'))
    {
        for my $mdr (
            @{ $slice->get_all_SimpleFeatures() })
        {
            if ($mdr->analysis->logic_name =~ /^mdr_(.*)/) {
                my $rlevel = $1;
                my $repeats = ceil(10 ** $rlevel);
                $self->write(@ff, 'misc_feature',
                    $self->features2location([$mdr]));
                $self->write(@ff, '', qq(/note="Mathematically-Defined Repeat (>= $repeats copies)."));}
        }
    }
}

1;
