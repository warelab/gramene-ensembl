

=head1 NAME

    repeats_coding_exon_stats.pl -registry xxx -species species -logic_name repeat_logic_name or all by default

=cut

=head1 SYNOPSIS

     repeats_coding_exon_stats.pl
           -registry registry_file for connection to ensembl database 
           -species species_name_of_interest
           -logic_name [optional] repeat_logic_name or 
                     all repeats from all repeat annotation by default
    

=cut

=head1 DESCRIPTION

    This prgoram report the coverage of repeat analysis on the coding Exons

=head1 APPENDIX



=cut

use Getopt::Long;
use Pod::Usage;
use warnings;
use  Bio::EnsEMBL::Registry;

my ($help, $registry, $species, %sums, $logic_name, $debug);
   
GetOptions (
    "help"       => \$help,
    "registry=s" => \$registry,
    "species=s"  => \$species,
    "logic_name=s" => \$logic_name,
    "debug"     => \$debug,
    ) or pod2usage;

if( $help || !$registry || !$species){
    
    pod2usage(-verbose => 2);
}

Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( print( "\n[*DIE] No $species set in $registry\n\n" ) &&
	      pod2usage() );

my $rf_adaptor =  $ENS_DBA->get_RepeatFeatureAdaptor();
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor();

#Necessary to get unique bits of Y
my $slices = $slice_adaptor->fetch_all('toplevel');

warn("Fetched ", scalar @$slices, " toplevel seq regions");

my $totals;

foreach my $slice (@$slices) {
    
    print "\n\nprocessing slice ". $slice->seq_region_name. "\n" ;#if $debug;
#    next unless  $slice->seq_region_name eq '1';

    my $coding_exons = get_coding_exon_regions($slice);
    print "\ncoding_exons number " . @$coding_exons . "\n" ;#if($debug);
    # it will be 1 even if there is none, only the segment will be [0,-1]
    
    my $uncovered = 0;
    foreach my $coding_exon (@$coding_exons) {#all consolidated coding exons on this slice
	my ($start, $end) = @$coding_exon;	
	my $coding_exon_length = ($end - $start + 1);
	print "start=$start end=$end length=$coding_exon_length\n" if($debug);

	next unless ($coding_exon_length > 0);
	#my $slice_coding_exon = $slice_adaptor->fetch_by_region('toplevel', 8, 1, 1000000); #
	my $slice_coding_exon = $slice_adaptor->fetch_by_region('toplevel', $slice->seq_region_name, $start, $end);
	#Store the total coding_exon_length
	$totals->{'coding_exon_length'} += $coding_exon_length;
	print "coding_exon_length=".	$totals->{'coding_exon_length'} ."\n" if $debug;

	#restricted repeat_features 
	my $rfs = $logic_name ? 
	    $rf_adaptor->fetch_all_by_Slice( $slice_coding_exon, $logic_name):
	    $rf_adaptor->fetch_all_by_Slice( $slice_coding_exon);
	print "number of repeat feature overlapping this coding exon region is ". scalar @$rfs . "\n" if $debug;
	
	my $restricted_segs;	
	if (@$rfs >= 1) {# if there are overlapping repeat features for this coding exon
	    $restricted_segs = restrict_overlapping_segs($rfs, $slice_coding_exon);
	}else{
	    	$totals->{'uncovered'} += $coding_exon_length;
		next;
	}
	
	#skip if there is no repeat features overlapping this coding exon
	next unless ($restricted_segs && @$restricted_segs >= 1);
	print "num repeat segs " . @$restricted_segs . "\n" if($debug);
	
	my $rf_seg_covering_coding_exon = 0;
	foreach my $seg (@$restricted_segs) {
	    my $seg_length = $seg->[2]-$seg->[1]+1;
	    $rf_seg_covering_coding_exon += $seg_length;
	}#end of each coding exon overlapping repeat feature segments
	
	$totals->{'rf_seg_covering_coding_exon'} += $rf_seg_covering_coding_exon;
	
	#The uncovered portion of the coding_exon must be the total length of the coding_exon (in slice coords) minus the portion covered by the restricted reference genomic_aligns (there may be more than one)
	my $uncovered_coding_exon_len = $coding_exon_length-$rf_seg_covering_coding_exon;
	print "uncovered  ($coding_exon_length-$rf_seg_covering_coding_exon) " . ($coding_exon_length-$rf_seg_covering_coding_exon) . "\n\n" if($debug || $uncovered_coding_exon_len < 0);
	$totals->{'uncovered'} += ($coding_exon_length-$rf_seg_covering_coding_exon);
    }# end of foreach coding exon
   # last;

    for my $k (keys %$totals){
	
	print "$k = ", $totals->{$k}, "\n";
    }
    
    
}#end of each slide

print "\n\nSummary\n";
for my $k (keys %$totals){    
    print "$k = ", $totals->{$k}, "\n";
    printf ("%.2f percetage\n",  $totals->{$k}/$totals->{coding_exon_length}) unless $k eq 'coding_exon_length';
}





sub get_coding_exon_regions {
  my ($this_slice) = @_;
  my $regions = [];

  return undef if (!$this_slice);

  my $all_coding_exons = [];
  my $all_genes = $this_slice->get_all_Genes_by_type("protein_coding");
  foreach my $this_gene (@$all_genes) {
    my $all_transcripts = $this_gene->get_all_Transcripts();
    foreach my $this_transcript (@$all_transcripts) {
      push(@$all_coding_exons, @{$this_transcript->get_all_translateable_Exons()});
    }
  }

  #consolidate the exons
  my $last_start = 0;
  my $last_end = -1;
  foreach my $this_exon (sort {$a->seq_region_start <=> $b->seq_region_start} @$all_coding_exons) {
      print "exon start " . $this_exon->seq_region_start . " end " . $this_exon->seq_region_end . " last_start $last_start last_end $last_end\n" if $debug;

    if ($last_end < $this_exon->seq_region_start) {
      if ($last_end > 0) {
        push(@$regions, [$last_start, $last_end]);
      }
      $last_end = $this_exon->seq_region_end;
      $last_start = $this_exon->seq_region_start;
    } elsif ($this_exon->seq_region_end > $last_end) {
      $last_end = $this_exon->seq_region_end;
    }
  }

  #Add final region
  push (@$regions, [$last_start, $last_end]);
  return $regions;
}


#
#Need to deal with overlapping blocks:
#Sort gabs by dnafrag_start of the reference 
#Compare current gab ($gab) with previous gab (last_gab)
#If the there is no overlap, add last_gab to array to be returned (restricted_gabs)
#If the end of the current ref_ga is less than the end of previous ref_ga (last_end), current gab must be within prev gab so skip
#If the end of the current ref_ga is greater than the end of previous ref_ga (last_end) and the start is not the same, then restrict
#block 1 from last_start to current ref_ga start. If the start positions are the same, use the longer block 2.
#
sub restrict_overlapping_segs {
    my ($rfs, $slice_ce) = @_;

    my $start = $slice_ce->start;
    my $end = $slice_ce->end;
    my $restricted_segs;
    my $last_start = 0;
    my $last_end = -1;
    my $last_seg;

    #Sort on ref_ga->dnafrag_start
    foreach my $rf (sort {$a->seq_region_start <=> $b->seq_region_start} @$rfs) {
	warn"!!!!! not on the same sequence, repeat feature on ".$rf->seq_region_name . ", codinge exon on ".$slice_ce->seq_region_name
	&& exit    unless $rf->seq_region_name eq $slice_ce->seq_region_name;

        print "  rf_start " . $rf->seq_region_start . " rf_end " . $rf->seq_region_end . "\n" if ($debug);
       
	if ($last_end < 0) {
            #first time through
            $last_end = $rf->seq_region_end;
            $last_start = $rf->seq_region_start;
            $last_seg = [$rf->seq_region_name, 
			 $rf->seq_region_start < $start ? $start :$rf->seq_region_start, 
			 $rf->seq_region_end > $end ? $end : $rf->seq_region_end];
        } elsif ($rf->seq_region_start > $last_end) {
            #no overlap, no restriction necessary
            print "  OVER No overlap\n" if($debug);

            #Store the last_gab
            push @$restricted_segs, $last_seg;
	    print "pushin seg: ". (join ',', @$last_seg) . "\n" if $debug;

            #Set 'last' to new gab
            $last_end = $rf->seq_region_end;
            $last_start = $rf->seq_region_start;
            $last_seg = [$rf->seq_region_name, 
			 $rf->seq_region_start < $start ? $start :$rf->seq_region_start, 
			 $rf->seq_region_end > $end ? $end : $rf->seq_region_end];
	} elsif ($rf->seq_region_end <= $last_end) {
            #block 1 covers block 2, no restriction necessary
            print "  OVER block 2 covered by block 1\n" if($debug);

            #Don't set 'last'
        } elsif ($rf->seq_region_end > $last_end) {
            #block 2 extends beyond block 1
                print "  OVER block 2 extend block 1 end $last_end and block2 start " . $rf->seq_region_start . "\n" if($debug);
                $last_end = $rf->seq_region_end;
		$last_seg = [$rf->seq_region_name, 
			     $last_start < $start ? $start : $last_start, 
			     $last_end > $end ? $end : $last_end];
		#push @$restricted_segs, $last_seg;
		#print "pushin seg: ". (join ',', @$last_seg) . "\n" if $debug;
                #Set 'last' to new gab

        }else{
	    die "unexpected situation last seg [$last_start, $last_end], rf [" .$rf->seq_region_start, ", " . $rf->seq_region_end. "]\n" if ($debug);
	}
    }

    #Need to deal with last gab
    push @$restricted_segs, $last_seg;
    
    print "pushin seg: ". (join ',', @$last_seg) . "\n" if $debug;

    return $restricted_segs;
}

1;
