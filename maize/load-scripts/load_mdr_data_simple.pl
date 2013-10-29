#!/usr/local/bin/perl -w

=head1 NAME

load_corebins.pl - Populates a core Ensembl DB with copy number features

=head1 SYNOPSIS

perl load_corebins.pl [options] copy_number.file

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -n --no_insert
 -t --threshold
 -b --bridge

=head1 OPTIONS

Reads the B<corebins.file>, and uses its copy number data to 
load the simple_feature table. Chains all MDR blocks <= specified bridge apart.

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

B<-n --no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Populates a core Ensembl DB with copy number features given 
  a user-input rlevel log threshold  

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper); # For debug 

use DBI;
use FindBin qw($Bin) ;
use File::Basename qw(dirname);

use vars qw($BASEDIR);
BEGIN{
  # Set the perl libraries
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/bioperl-live';
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::SimpleFeature;
use Bio::Seq;
use Bio::SeqIO;

use vars qw($I $ENS_DBA);
my $help=0;
my $man=0;
my($threshold, $species, $file, $no_insert, $bridge);
GetOptions
    (
     "help|?"          => \$help,
     "man"             => \$man,
     "threshold=f"     => \$threshold,
     "species=s"       => \$species,
     "registry_file=s" => \$file,
     "no_insert"       => \$no_insert,
     "bridge=f"        => \$bridge,
     ) or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

$I = $no_insert ? 0 : 1; # Put stuff in the database?

# Validate file paths
$file    ||= $BASEDIR.'/conf/ensembl.registry';
my $copy_file = shift @ARGV;
$copy_file   || (warn("Need the path to a copy number\n") && pod2usage(1));

map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
} $file, $copy_file;

# Load the ensembl file
$species || ( warn( "Need a --species\n" ) && pod2usage(1) );
Bio::EnsEMBL::Registry->load_all( $file );
$ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( warn( "No core DB for $species set in $file\n" ) &&
	      pod2usage(1) );

# Process the copy number score file
my ($copy_coordz, $scorez) = copy_coords ($copy_file, $threshold);

# Manipulate the bins so that adjacent one (<=$bridge) are merged
my %copy_coords = merge_coords ($copy_coordz, $bridge);

# print table of contig coords for debugging
foreach my $contig (keys %copy_coords){
    my @ranges = @{$copy_coords{$contig}};
    foreach my $range (@ranges){
        my ($start, $end) = @{$range};
	my $length = ($end - $start) + 1;
	print "$contig\t$start\t$end\t$length\n";
    }
}

# Prepare some Ensembl adaptors and some static vars
my $slice_adaptor    = $ENS_DBA->get_adaptor('Slice');
my $simple_adaptor   = $ENS_DBA->get_adaptor('SimpleFeature');
my $analysis_adaptor = &fetch_analysis( "mdr_$threshold" );

my $display_label    = 'MDR_' . $threshold;

# Store the data in simple_feature
# omitting scores for now since chaining makes them weird
foreach my $contig (keys %copy_coords){
    warn "Working on $contig\n";
    
    my @ranges = @{$copy_coords{$contig}};
#    my @scorez = @{$scorez{$contig}};
    
    # get the contig slice
    my $slice = $slice_adaptor->fetch_by_region('contig', $contig);
    
#    my $scorez_count = 0;
    foreach my $range (@ranges){
	my ($start, $end) = @{$range};
#	my ($store_score) = @{$scorez[$scorez_count]};

#	print "$contig\t$start\t$end\n";

	# construct the data point
	my $mdr_feature = Bio::EnsEMBL::SimpleFeature->new
	    (-start         => $start,
	     -end           => $end,
	     -strand        => 1,
	     -slice         => $slice,
	     -analysis      => $analysis_adaptor,
	     -display_label => $display_label,
#	     -score         => $store_score,
	     -score         => 1,
	     );
	
	$simple_adaptor -> store ($mdr_feature) if $I;
#	$scorez_count++;
    }
}

####SUBS####

sub merge_coords {
    my $copy_coordz = shift;
    my $bridge      = shift;

    my %copy_coords;
    foreach my $contig (sort keys %{$copy_coordz}){
	warn "Merging bins for $contig\n";

	my @ranges = @{$copy_coordz->{$contig}};
	
	# prime the cache with the first bin                                    
	my $start_cache = $ranges[0]->[0];
	my $end_cache   = $ranges[0]->[1];

	# cycle through other bins and merge until >10 bp gap                        
	foreach my $range (@ranges){
	    my ($start, $end) = @{$range};
	    
	    # bail if we're looking at the first one                              
	    if ($start_cache == $start){
		next;
	    }
	    
	    # extend if <= bridge
	    elsif (($start - $end_cache) <= $bridge){
		$end_cache = $end;
	    }

	    # store if >=bridge
	    else {
		push @{$copy_coords{$contig}}, [ $start_cache, $end_cache ];
		$start_cache = $start;
		$end_cache   = $end;
	    }
	}
	
	# push the last cached values onto the stack                                      
	push @{$copy_coords{$contig}}, [ $start_cache, $end_cache ];
	
    }
    return (%copy_coords);
}
    
sub copy_coords {
    my ($file,$threshold) = @_;
    
    my %ranges;
    my %scores;
    my $copy = Bio::SeqIO -> new (-format => 'qual', -file => "$copy_file");
    while (my $copy_obj = $copy -> next_seq()){
	my $defline = $copy_obj -> display_id();
	my @def_components = split (/\:/, $defline);
	my $id = $def_components[2];

#	my $id = $copy_obj -> display_id();

	warn "Storing score data for $id\n";
	
	my @scores = @{$copy_obj -> qual()};

	my @copies;
	my $cpavg = 0;
	my @range;
	my $count = 1; #start coord                                                      
	foreach my $score (@scores){
	    $score++; #can't take log of zero, so inc by 1 everywhere                     
	    my $log_sc = log($score) / log(10);
	    if ($log_sc >= $threshold){
		push (@range, $count);
		push (@copies, $score);
	    } else {

		# calculate ranges
		push @{$ranges{$id}}, [ $range[0],$range[-1] ] if @range;
		@range = ();

		# calculate copy numbers for score field
		if (@copies){
		    my $elements = @copies;
		    my $sum = 0;
		    foreach my $cp (@copies){
			$sum+=$cp;
		    }
		    
		    $cpavg = $sum / $elements;

		    push @{$scores{$id}}, [ $cpavg ];
		}
		@copies =();
	    }

	    $count++;
	}
	
	# get the last bit if it goes to the edge
	push @{$ranges{$id}}, [ $range[0],$range[-1] ] if @range;
	push @{$scores{$id}}, [ $cpavg ] if @copies;
	
    }

    return (\%ranges, \%scores);
}

# Returns an Analysis object; either fetched from the DB if one exists,          
# or created fresh, in which case it is stored to the DB.                        
sub fetch_analysis{
    my $logic_name = shift || die("Need a logic_name" );
    my $db_file    = shift || '';
    my $adaptor = $ENS_DBA->get_adaptor('Analysis');

    my %args = ( -logic_name=>$logic_name,
                 $db_file ? (-db_file=>$db_file) : () );

    my $analysis;
    if( $analysis = $adaptor->fetch_by_logic_name($args{-logic_name}) ){
        # Analysis found in database already; use this.                          
        return $analysis;
    }

    # No analysis - create one from scratch                                      
    $analysis = Bio::EnsEMBL::Analysis->new(%args);
    $adaptor->store($analysis);
    return $analysis;
}
