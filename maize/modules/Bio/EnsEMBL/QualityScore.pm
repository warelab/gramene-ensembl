#
# Ensembl module for Bio::EnsEMBL::QualityScore
#
# Cared for by Shiran Pasternak <shiran@cshl.edu>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::QualityScore - Stores conservation scores

=head1 SYNOPSIS

use Bio::EnsEMBL::QualityScore;
my $quality_score = new Bio::EnsEMBL::QualityScore(
    -window_size => $window_size, 
    -position    => $position, 
    -scores      => $scores,
);

SET VALUES
    $quality_score->window_size(10);
    $quality_score->position(1);
    $quality_score->scores($scores);
    $quality_score->y_axis_min(0);
    $quality_score->y_axis_max(100);


GET VALUES
    $window_size = $quality_score->window_size;
    $position = $quality_score->position;
    $scores = $quality_score->observed_scores;
    $y_axis_min = $quality_score->y_axis_min;
    $y_axis_max = $quality_score->y_axis_max;

=head1 DESCRIPTION

Object for storing quality scores. The scores are averaged over different
window sizes to speed up drawing over large regions. The scores are packed as 
floats and stored in a string. The scores can be stored and retrieved in 
either a packed or unpacked format. The unpacked format is as a space delimited
string eg ("75 90 22"). The packed format is a single precision float 
(4 bytes). It is recommended to use the unpacked format.

=head1 AUTHOR - Shiran Pasternak

This modules is part of the Ensembl project http://www.ensembl.org

Email shiran@cshl.edu

=head1 CONTACT

This modules is part of the EnsEMBL project (http://www.ensembl.org)

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::QualityScore;

use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning deprecate throw);

#store as 4 byte float
my $pack_size = 4;
my $pack_type = "f";

=head2 new (CONSTRUCTOR)
    Arg [-ADAPTOR] 
        : Bio::EnsEMBL::DBSQL::quality_score $adaptor
        (the adaptor for connecting to the database)
    Arg [-SLICE] (opt)
        : Bio::EnsEMBL::Slice $slice
        (the Slice on which the quality scores are called)
    Arg [-WINDOW_SIZE] (opt)
         : int $window_size
         (window size used to average the scores over)
    Arg [-POSITION] (opt)
        : int $position
        (position of the first score in genomic coordinates)
    Arg [-SCORES]
        : string $expected_score
        (packed or unpacked string of expected scores)
    Arg [-PACKED] (opt)
        : boolean $packed
        (whether the scores are packed (1) or unpacked (0))
    Arg [Y_AXIS_MIN] (opt)
	    : float $y_axis_min
	    (minimum score value used for display)
    Arg [Y_AXIS_MAX] (opt)
	    : float $y_axis_max
	    (maximum score value used for display)
    Example :
	my $quality_score = new Bio::EnsEMBL::QualityScore(
                    	    -slice       => $slice,
                            -window_size => $window_size, 
                            -position    => $position, 
                            -scores      => $scores, 
    );
       Description: Creates a new quality_score object
       Returntype : Bio::EnsEMBL::QualityScore
       Exceptions : none
       Caller     : general
       Status     : At risk

=cut

sub new {

    my ($class, @args) = @_;

    my $self = {};
    bless $self, $class;

    my ($adaptor, $slice, $window_size, $position, $scores, $packed,
        $y_axis_min, $y_axis_max)
        = rearrange(
        [   qw(ADAPTOR SLICE WINDOW_SIZE
                POSITION SCORES PACKED
                Y_AXIS_MIN Y_AXIS_MAX)
        ],
        @args
        );

    $self->adaptor($adaptor)         if (defined($adaptor));
    $self->slice($slice)             if (defined($slice));
    $self->window_size($window_size) if (defined($window_size));
    $self->position($position)       if (defined($position));
    $self->scores($scores)           if (defined($scores));
    $self->y_axis_min($y_axis_min)   if (defined($y_axis_min));
    $self->y_axis_max($y_axis_max)   if (defined($y_axis_max));

    if (defined($packed)) {
        $self->packed($packed);
    } else {
        $self->packed(0);
    }
    return $self;
}

=head2 new_fast

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: This is an ultra fast constructor which requires knowledge of
               the objects internals to be used.
  Returntype :
  Exceptions : none
  Caller     :

=cut

sub new_fast {
    my ($class, $hashref) = @_;

    return bless $hashref, $class;
}

=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::QualityScoreAdaptor $adaptor
  Example    : $quality_score->adaptor($adaptor);
  Description: Getter/Setter for the adaptor this object used for database
               interaction
  Returntype : Bio::EnsEMBL::DBSQL::QualityScoreAdaptor object
  Exceptions : thrown if the argument is not a
               Bio::EnsEMBL::DBSQL::QualityScoreAdaptor object
  Caller     : general
  Status     : At risk

=cut

sub adaptor {
    my ($self, $adaptor) = @_;

    if (defined($adaptor)) {
        throw(
            "$adaptor is not a Bio::EnsEMBL::DBSQL::QualityScoreAdaptor object"
            )
            unless (
            $adaptor->isa("Bio::EnsEMBL::DBSQL::QualityScoreAdaptor"));
        $self->{'adaptor'} = $adaptor;
    }

    return $self->{'adaptor'};
}

=head2 slice

  Arg [1]    : (optional) Bio::EnsEMBL::Slice $slice
  Example    : $seqname = $feature->slice()->name();
  Description: Getter/Setter for the Slice that is associated with this 
               feature.  The slice represents the underlying sequence that this
               feature is on.  Note that this method call is analagous to the
               old SeqFeature methods contig(), entire_seq(), attach_seq(),
               etc.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : thrown if an invalid argument is passed
  Caller     : general
  Status     : Stable

=cut

sub slice {
    my $self = shift;

    if (@_) {
        my $sl = shift;
        if (defined($sl) && (!ref($sl) || !$sl->isa('Bio::EnsEMBL::Slice'))) {
            throw('slice argument must be a Bio::EnsEMBL::Slice');
        }

        $self->{'slice'} = $sl;
    }

    return $self->{'slice'};
}


=head2 window_size

  Arg [1]    : (opt) integer window_size
  Example    : my $window_size = $quality_score->window_size();
  Example    : $quality_score->window_size(1);
  Description: Getter/Setter for the window_size of this conservation score
  Returntype : integer, Returns 1 if value not defined
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub window_size {
    my ($self, $window_size) = @_;

    if (defined $window_size) {
        $self->{'window_size'} = $window_size;
    }
    $self->{'window_size'} = '1' unless (defined($self->{'window_size'}));
    return $self->{'window_size'};
}

=head2 position

  Arg [1]    : (opt) integer
  Example    : $quality_score->position(1);
  Description: Getter/Setter for the genomic position of the first score
  Returntype : integer. Return 1 if value not defined
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub position {
    my ($self, $position) = @_;

    if (defined $position) {
        $self->{'position'} = $position;
    }

    $self->{'position'} ||= '1';
    return $self->{'position'};
}

=head2 start

  Example    : $quality_score->start();
  Description: Wrapper round position call 
  Returntype : integer. Return 1 if value not defined
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub start {
    my ($self, $start) = @_;
    if (defined $start) {
        $self->{'start'} = $start;
    }
    return $self->{'start'};
}

=head2 end

  Example    : $quality_score->end();
  Description: wrapper around position
  Returntype : integer. Return 1 if value not defined
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub end {
    my ($self, $end) = @_;
    if (defined $end) {
        $self->{'end'} = $end;
    }
    return $self->{'end'};
}

=head2 scores

  Arg [1]    : (opt) string of difference scores (expected - observed)
               (can be either packed or space delimited)
  Example    : $quality_score->scores("1.85 -2.54 1.56");
  Example    : my $scores = $quality_score->scores();
  Description: Getter/Setter for the difference score string
  Returntype : string (either packed or space delimited)
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub scores {
    my ($self, $scores) = @_;

    if (defined $scores) {
        $self->{'scores'} = $scores;
    }
    return $self->{'scores'};
}

=head2 score

  Arg [1]    : (opt) string of the aggregate score
               (can be either packed or space delimited)
  Example    : $quality_score->score("1.85 -2.54 1.56");
  Example    : my $score = $quality_score->score();
  Description: Getter/Setter for the score string (synonym of 'scores')
  Returntype : string (either packed or space delimited)
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub score {
    my ($self, $score) = @_;

    return $self->scores($score);
}

=head2 y_axis_min

  Arg [1]    : (opt) float
  Example    : $quality_score->y_axis_min(-0.5);
  Example    : $y_axis_min = $quality_score->y_axis_min;
  Description: Getter/Setter for the minimum score
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub y_axis_min {
    my ($self, $y_axis_min) = @_;

    if (defined $y_axis_min) {
        $self->{'y_axis_min'} = $y_axis_min;
    }
    return $self->{'y_axis_min'};
}

=head2 y_axis_max

  Arg [1]    : (opt) float
  Example    : $quality_score->y_axis_max(2.45);
  Example    : $y_axis_max = $quality_score->y_axis_min;
  Description: Getter/Setter for the maximum score
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub y_axis_max {
    my ($self, $y_axis_max) = @_;

    if (defined $y_axis_max) {
        $self->{'y_axis_max'} = $y_axis_max;
    }
    return $self->{'y_axis_max'};
}

=head2 packed

  Arg [1]    : (opt) boolean 
  Example    : $quality_score->packed(1);
  Example    : $packed = $quality_score->packed;
  Description: Getter/Setter for the whether the scores are packed or space
               delimited
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub packed {
    my ($self, $packed) = @_;

    if ($packed) {
        $self->{'packed'} = $packed;
    }
    return $self->{'packed'};
}

=head2 reverse

  Example    : $quality_score->reverse;
  Description: reverse scores and position in the quality_score object
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub reverse {
    my ($self) = @_;
    my $num_scores = 0;
    return if (!defined($self->scores));
    if ($self->packed) {
        $num_scores = length($self->scores) / $pack_size;
    } else {
        my @scores = split ' ', $self->scores;
        $num_scores = scalar(@scores);
    }

    #swap position orientation and reverse position in alignment
    my $end = $self->position + (($num_scores - 1) * $self->window_size);

    #10.10.06 +1 so position starts at 1 not 0
    #$self->position($self->genomic_align_block->length - $end);
    $self->position($self->slice->length - $end + 1);

    #swap position orientation and reverse position in alignment
    if (defined $self->position) {
        $end = $self->position
            + (($num_scores - 1) * $self->window_size);

        #10.10.06 +1 so position starts at 1 not 0
        #$self->position($self->genomic_align_block->length - $end);
        $self->position($self->slice->length - $end + 1);
    }

    #reverse score strings
    $self->scores(
        _reverse_score($self->scores, $num_scores, $self->packed));
}

#internal method used by reverse to reverse the score strings
#arguments:
#score_str : string, score string
#num_scores : integer, number of scores in the string
#packed : boolean, whether the scores are packed or not
sub _reverse_score {
    my ($score_str, $num_scores, $packed) = @_;

    my $rev_str;
    if ($packed) {
        for (my $i = $num_scores - 1; $i >= 0; $i--) {
            my $value = substr $score_str, $i * $pack_size, $pack_size;
            $rev_str .= $value;
        }
    } else {
        my @scores = split ' ', $score_str;
        my $rev_str;
        for (my $i = $num_scores - 1; $i >= 0; $i--) {
            $rev_str .= $scores[$i];
        }
    }
    return $rev_str;
}

#print the contents of the quality_score object
sub _print {
    my ($self, $FILEH) = @_;

    $FILEH ||= \*STDOUT;

    print $FILEH <<PRINT;
Bio::EnsEMBL::QualityScore object ($self)
  slice       = @{[ $self->slice ]}
  window_size = @{[ $self->window_size ]}
  position    = @{[ $self->position ]}
  scores      = @{[ $self->scores ]}
PRINT
}

1;
