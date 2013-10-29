#
# Ensembl module for Bio::EnsEMBL::DBSQL::QualityScoreAdaptor
#
# Cared for by Shiran Pasternak <shiran@cshl.edu>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::QualityScoreAdaptor - Object adaptor to access data in the quality_score table

=head1 SYNOPSIS

  Connecting to the database using the Registry

     use Bio::EnsEMBL::Registry;
 
     my $reg = "Bio::EnsEMBL::Registry";

      $reg->load_registry_from_db(-host=>"ensembldb.ensembl.org", -user=>"anonymous");

      my $quality_score_adaptor = $reg->get_adaptor(
         "core", "QualityScore");

  Store data in the database

     $quality_score_adaptor->store($quality_score);

  To retrieve score data from the database using the default display_size
     $quality_scores = $quality_score_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_link_species_set, $slice);

  To retrieve one score per base in the slice
     $quality_scores = $quality_score_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_link_species_set, $slice, $slice->end-$slice->start+1);
  Print the scores
   foreach my $score (@$quality_scores) {
      printf("position %d observed %.4f expected %.4f difference %.4f\n",  $score->position, $score->observed_score, $score->expected_score, $score->diff_score);
   }

  A simple example script for extracting scores from a slice can be found in ensembl-compara/scripts/examples/getConservationScores.pl

=head1 DESCRIPTION

This module is used to access data in the quality_score table.
Each score is represented by a Bio::EnsEMBL::QualityScore. The position and an observed, expected score and a difference score (expected-observed) is stored for each column in a multiple alignment. Not all bases in an alignment have a score (for example, if there is insufficient coverage) and termed here as 'uncalled'. 
In order to speed up processing of the scores over large regions, the scores are stored in the database averaged over window_sizes of 1 (no averaging), 10, 100 and 500. When retrieving the scores, the most appropriate window_size is estimated from the length of the alignment or slice and the number of scores requested, given by the display_size. There is no need to specify the window_size directly. If the number of scores requested (display_size) is smaller than the alignment length or slice length, the scores will be either averaged if display_type = "AVERAGE" or the maximum value taken if display_type = "MAX". Scores in uncalled regions are not returned. If a score for each column in an alignment is required, the display_size should be set to be the same size as the alignment length or slice length. 

=head1 AUTHOR - Kathryn Beal

This modules is part of the Ensembl project http://www.ensembl.org

Email kbeal@ebi.ac.uk

=head1 CONTACT

This modules is part of the EnsEMBL project (http://www.ensembl.org)

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::DBSQL::QualityScoreAdaptor;
use vars qw(@ISA);
use strict;

use POSIX qw(ceil floor);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::QualityScore;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info deprecate);

#global variables

#store as 4 byte float. If change here, must also change in
#ConservationScore.pm
my $_pack_size = 4;
my $_pack_type = "f";

my $_bucket;
my $_score_index = 0;

#my $_no_score_value = 0.0; #value if no score
my $_no_score_value = undef;    #value if no score

my $PACKED = 1;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_all_by_Slice

  Arg  1     : Bio::EnsEMBL::Slice $slice
  Arg  2     : (opt) integer $slice_length (default $slice->length)
  Arg  3     : (opt) integer $display_size (default 700)
  Arg  4     : (opt) string $display_type (one of "AVERAGE" or "MAX") (default "AVERAGE")
  Arg  5     : (opt) integer $window_size
  Example    : my $quality_scores =
                    $quality_score_adaptor->fetch_all_by_Slice(
                        $slice,
                        $slice_length,
                        $display_size,
                        $display_type
                    );
  Description: Retrieve the corresponding
               Bio::EnsEMBL::QualityScore objects. 
	       Each quality score object contains a position in genomic
               coordinates, and the raw scores.
               The $slice_length is the total length of the region to be 
               displayed.
               It is used to automatically calculate the window_size.
               Display_size is the number of scores that will be returned. If 
               the $slice_length is larger than the $display_size, the scores 
               will either be averaged if the display_type is "AVERAGE" or the 
               maximum taken if display_type is "MAXIMUM". 
	       Window_size defines which set of pre-averaged scores to use. 
	       Valid values are 1, 10, 100 or 500. There is no need to define 
               the window_size because the program will select the most 
               appropriate window_size to use based on the slice_length and the
               display_size. 
               The min and max y axis values for 
               the array of quality score objects are set in the first 
               quality score object (index 0). 
  Returntype : ref. to an array of Bio::EnsEMBL::QualityScore 
               objects. 
  Caller     : object::methodname
  Status     : At risk

=cut

sub fetch_all_by_Slice {
    my $self = shift;
    my ($slice, $slice_length, $display_size, $display_type, $window_size)
        = @_;

    my $scores = [];

    #default display_size is 700
    if (!defined $display_size) {
        $display_size = 700;
    }

    #default display_mode is AVERAGE
    if (!defined $display_type) {
        $display_type = "AVERAGE";
    }

    #default slice_length is the slice length
    if (!defined $slice_length) {
        $slice_length = $slice->length;
    }

    #set up bucket object for storing bucket_size number of scores
    my $bucket_size = ($slice_length) / $display_size;

    #default window size is the largest bucket that gives at least
    #display_size values ie get speed but reasonable resolution
    my @window_sizes = (1, 10, 100, 500);

    #check if valid window_size
    my $found = 0;
    if (defined $window_size) {
        foreach my $win_size (@window_sizes) {
            if ($win_size == $window_size) {
                $found = 1;
                last;
            }
        }
        if (!$found) {
            warning("Invalid window_size $window_size");
            return $scores;
        }
    }

    if (!defined $window_size) {

        #set window_size to be the largest for when for loop fails
        $window_size = $window_sizes[ scalar(@window_sizes) - 1 ];
        for (my $i = 1; $i < scalar(@window_sizes); $i++) {
            if ($bucket_size < $window_sizes[$i]) {
                $window_size = $window_sizes[ $i - 1 ];
                last;
            }
        }
    }

    my $quality_scores
        = $self->_fetch_all_by_Slice_WindowSize($slice, $window_size,
        $PACKED);

    if (scalar(@$quality_scores) == 0) {
        return $scores;
    }

    # reset _score_index for new quality_scores
    $_score_index = 0;

    $scores = $self->_get_display_scores(
        {   'quality_scores' => $quality_scores,
            'display_type'   => $display_type,
            'window_size'    => $window_size,
            'slice'          => $slice
        }
    );

    # $scores = $quality_scores;

    if (scalar(@$scores) == 0) {
        return $scores;
    }

    #Find the min and max scores for y axis scaling. Save in first
    #conservation score object
    my ($min_y_axis, $max_y_axis) = _find_min_max_score($scores);

    #add min and max scores to the first conservation score object
    if ((scalar @$scores) > 0) {
        $scores->[0]->y_axis_min($min_y_axis);
        $scores->[0]->y_axis_max($max_y_axis);
    }
    return ($scores);

}

=head2 fetch_raw_scores_for_Slice

  Arg  1     : Bio::EnsEMBL::Slice $slice
  Example    : my $quality_score =
                    $quality_score_adaptor->fetch_raw_scores_for_Slice(
                        $slice
                    );
  Description: Retrieve a singular Bio::EnsEMBL::QualityScore object
               corresponding to the quality scores of the input slice.
  Returntype : Bio::EnsEMBL::QualityScore 
  Caller     : object::methodname
  Status     : At risk


=cut

sub fetch_raw_scores_for_Slice {
    my $self = shift;
    my ($slice) = @_;

    my $scores = [];

    my $quality_scores
        = $self->_fetch_all_by_Slice_WindowSize($slice, 1, $PACKED);

    my $merged_quality_score = $self->_merge_scores($quality_scores);
    return $merged_quality_score;
}

=head2 store

  Arg [1]    : Bio::EnsEMBL::QualityScore $cs
  Example    : $csa->store($cs);
  Description: Stores a conservation score object in the compara database if
               it has not been stored already.  
  Returntype : none
  Exceptions : thrown if $slice is not a 
               Bio::EnsEMBL::Compara::GenomicAlignBlock object
  Exceptions : thrown if the argument is not a Bio::EnsEMBL::Compara:ConservationScore
  Caller     : general
  Status     : At risk
=cut

sub store {
    my ($self, $score) = @_;

    unless (defined $score
        && ref $score
        && $score->isa('Bio::EnsEMBL::QualityScore'))
    {
        $self->throw("Must have quality score argument [$score]");
    }

    my $slice = $score->slice;

    # Check to see if slice, window_size, and position have been defined
    # (should be unique)
    my @bad_fields = ();
    for my $field (qw/slice window_size position/) {
        push @bad_fields, $field unless defined($score->$field);
    }
    if (scalar @bad_fields) {
        $self->throw(
            "Quality score is missing the following required parameters: ",
            join(", ", @bad_fields));
    }

    #check if slice is valid
    if (!$slice->isa("Bio::EnsEMBL::Slice")) {
        throw("[$slice] is not a Bio::EnsEMBL::Slice");
    }

    # TODO: Check that slice is at seq level

    # Store the quality score
    my $sql = <<SQL;
replace into quality_score (seq_region_id, window_size, position, score)
     values (?, ?, ?, ?)
SQL
    my $sth       = $self->prepare($sql);
    my $sth_index = 1;
    $sth->bind_param($sth_index++, $slice->get_seq_region_id);
    $sth->bind_param($sth_index++, $score->window_size);
    $sth->bind_param($sth_index++, $score->position);
    $sth->bind_param($sth_index++, $score->score);

    $sth->execute();

    #update the quality_score object so that its adaptor is set
    $score->adaptor($self);
}

### Internal methods

#  Arg  1      : integer $slice
#  Arg  2      : integer $window_size
#  Arg  3      : (opt) boolean $packed (default 0)
#  Example     : my $quality_scores =
#                        $quality_score_adaptor->fetch_all_by_Slice($slice);
#  Description : Retrieve the corresponding
#                Bio::EnsEMBL::QualityScore objects.
#  Returntype  : ref. to an array of Bio::EnsEMBL::QualityScore objects. If
#                $packed is true, return the scores in a packed format given
#                by $_pack_size and $_pack_type.
#  Caller      : object::methodname

sub _fetch_all_by_Slice_WindowSize {
    my ($self, $slice, $window_size, $packed) = @_;
    my $quality_scores = [];
    my $scores;

    #whether to return the scores in packed or unpacked format
    #default to unpacked (space delimited string of floats)
    if (!defined $packed) {
        $packed = 0;
    }

    my $sql = <<SQL;
select seq_region_id, window_size, position, score
  from quality_score
 where seq_region_id = ?
   and window_size = ?
SQL

    my @slices = ();

    my $print_slice = sub {
        my $slice = shift;
        join(" ",
            $slice->get_seq_region_id, $slice->seq_region_name, $slice->start,
            $slice->end, $slice->coord_system_name),
            ;
    };

    my $seq_coord_system = $self->_get_sequence_level_coord_system();
    my $sth              = $self->prepare($sql);
    for my $segment (@{ $slice->project($seq_coord_system->name) }) {
        my $quality_slice = $segment->to_Slice();
        $sth->execute($quality_slice->get_seq_region_id, $window_size);

        while (my @values = $sth->fetchrow_array()) {
            if (!$packed) {
                $scores = _unpack_scores($values[3]);
            } else {
                $scores = $values[3];
            }

            # Lengths are 0-based, positions are 1-based
            my $position
                = ($segment->from_start - 1) 
                - $quality_slice->start
                + $values[2] - 1;
            my $quality_score = Bio::EnsEMBL::QualityScore->new_fast(
                {   'adaptor'     => $self,
                    'slice'       => $slice,
                    'window_size' => $values[1],
                    'position'    => $position,
                    'scores'      => $scores,
                    'packed'      => $packed
                }
            );
            push(@$quality_scores, $quality_score);
        }
    }

    #sort into numerical order based on position
    my @sorted_scores
        = sort { $a->{position} <=> $b->{position} } @$quality_scores;
    return \@sorted_scores;
}

=head2 _get_sequence_level_coord_system

    Finds the coordinate system on which the sequence is defined

=cut

sub _get_sequence_level_coord_system {
    my $self                 = shift;
    my $coord_system_adaptor = $self->db()->get_adaptor('CoordSystem');
    my ($seq_coord_system)
        = grep { $_->is_sequence_level }
        @{ $coord_system_adaptor->fetch_all() };
    return $seq_coord_system;
}

#find the min and max scores for y axis scaling
sub _find_min_max_score {
    my ($scores) = @_;
    my $min;
    my $max;

    foreach my $score (@$scores) {

        #find min and max of diff scores
        if (defined $score->scores) {

            #if min hasn't been defined yet, then define min and max
            unless (defined $min) {
                $min = $score->scores;
                $max = $score->scores;
            }
            if ($min > $score->scores) {
                $min = $score->scores;
            }
            if ($max < $score->scores) {
                $max = $score->scores;
            }
        }
    }

    return ($min, $max);
}

#reverse the conservation scores for complemented sequences
sub _reverse {
    my ($scores) = @_;

    #reverse each quality_score
    foreach my $s (@$scores) {
        $s->reverse;
    }

    #reverse array so position values go from small to large
    my @rev = reverse @$scores;

    return \@rev;
}

#unpack scores.
sub _unpack_scores {
    my ($scores) = @_;
    if (!defined $scores) {
        return "";
    }
    my $num_scores = length($scores) / $_pack_size;

    my $score = "";
    for (my $i = 0; $i < $num_scores * $_pack_size; $i += $_pack_size) {
        my $value = substr $scores, $i, $_pack_size;
        $score .= unpack($_pack_type, $value) . " ";
    }
    return $score;
}

#find the score index (row) that contains $pos in alignment coords
#use global variable $_score_index to keep track of where I am in the scores
#array
#$scores : array of conservation scores
#$num_scores : number of scores in the array
#$score_lengths : number of scores in each row of the array
#$pos : position to find
#$win_size : window size used from the database
sub _find_score_index {
    my ($scores, $num_scores, $score_lengths, $pos, $win_size) = @_;
    my $i;
    my $length;

    #special case for first window size
    if (   $pos < $scores->[0]->position
        && $pos > ($scores->[0]->position - $win_size))
    {
        return 0;
    }

    for ($i = $_score_index; $i < $num_scores; $i++) {
        $length = ($score_lengths->[$i] - 1) * $win_size;

        if (   $pos >= $scores->[$i]->position
            && $pos <= $scores->[$i]->position + $length)
        {
            $_score_index = $i;
            return ($i);
        }

        #smaller than end so there is no score for this position
        if ($pos < ($scores->[$i]->position + $length)) {
            $_score_index = $i;
            return -1;
        }
    }
    return -1;
}

#print scores (unpack first if necessary)
sub _print_scores {
    my ($scores, $packed) = @_;
    my $num_scores = scalar(@$scores);
    my $cnt;
    my ($start, $end);
    my $i;
    my @values;
    my $total_scores = 0;

    print "num scores $num_scores\n";
    for ($cnt = 0; $cnt < $num_scores; $cnt++) {
        if ($packed) {
            $end = (length($scores->[$cnt]->expected_score) / 4);
        } else {
            @values = split ' ', $scores->[$cnt]->diff_score;
            $end = scalar(@values);
        }
        print "row $cnt length $end\n";
        $total_scores += $end;
        for ($i = 0; $i < $end; $i++) {
            my $score;
            if ($packed) {
                my $value = substr $scores->[$cnt]->expected_score,
                    $i * $_pack_size, $_pack_size;
                $score = unpack($_pack_type, $value);
            } else {
                $score = $values[$i];
            }
            print "$i score $score \n";
        }
    }
    print "Total $total_scores\n";

}

use Bio::Seq::PrimaryQual;

sub _get_display_scores {
    my $self           = shift;
    my ($params)       = @_;
    my $quality_scores = $params->{'quality_scores'};
    my $display_type   = $params->{'display_type'};
    my $window_size    = $params->{'window_size'};
    my $slice          = $params->{'slice'};

    my $display_scores = [];

    for my $score (@$quality_scores) {
        my $slice = $score->slice;
        my $bioseq = Bio::Seq::PrimaryQual->new(-qual => $score->score());
        my $score_length = scalar @{ $bioseq->qual };
        my ($quality_position, $start_position, $end_position)
            = $self->_calculate_quality_positions($score, $score_length,
            $window_size);
        for (
            my $position = $start_position;
            $position < $end_position;
            $position += $window_size
            )
        {
            my $quality_score = $bioseq->qualat($quality_position);
            my $display_score = Bio::EnsEMBL::QualityScore->new_fast(
                {   'adaptor'     => $self,
                    'slice'       => $slice,
                    'window_size' => $window_size,
                    'position'    => $position,
                    'scores'      => $quality_score,
                    'packed'      => 0
                }
            );
            push @$display_scores, $display_score;
            $quality_position++;
        }
    }
    return $display_scores;
}

=head2 _calculate_quality_positions

    Calculates the start and end positions for a given quality score object in
    its containing slice

=cut

sub _calculate_quality_positions {
    my $self         = shift;
    my $score        = shift;
    my $score_length = shift;
    my $window_size  = shift || 1;

    my $slice = $score->slice();
    my ($quality_position, $start_position, $end_position);
    if ($score->position < 0) {
        $quality_position = ceil(-$score->position / $window_size);
        $start_position   = $slice->start;
        $end_position     = $start_position
            + ($score_length * $window_size + $score->position);
    } else {
        $quality_position = 1;
        $start_position   = $slice->start + $score->position - 1;
        $end_position     = $start_position + ($score_length * $window_size);
    }
    $end_position = $slice->end if $end_position > $slice->end;
    return ($quality_position, $start_position, $end_position);
}

sub _merge_scores {
    my $self           = shift;
    my $quality_scores = shift;

    my @raw_values = ();
    for my $score (@$quality_scores) {
        my $slice = $score->slice;
        my $bioseq = Bio::Seq::PrimaryQual->new(-qual => $score->score());
        my ($quality_position, $start_position, $end_position)
            = $self->_calculate_quality_positions($score,
            scalar @{ $bioseq->qual });
        for (
            my $position = $start_position;
            $position <= $end_position;
            $position++
            )
        {
            my $quality_score = 0;

            # Gracefully handle gaps: treat as phred 0
            if ($quality_position <= $bioseq->length()) {
                $quality_score = $bioseq->qualat($quality_position);
            }
            push @raw_values, $quality_score;
            $quality_position++;
        }
    }

    my $first_score = $quality_scores->[0];
    my $slice       = $first_score->slice();
    my $window_size = $first_score->window_size();
    my $merged      = Bio::EnsEMBL::QualityScore->new_fast(
        {   'adaptor'     => $self,
            'slice'       => $slice,
            'window_size' => $window_size,
            'position'    => $slice->start,
            'scores'      => join(' ', @raw_values),
        }
    );
    return $merged;
}

1;
