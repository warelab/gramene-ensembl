#!/opt/local/bin/perl

# Purpose: This script interfaces with ensembl using dbi and assigns new analyses
#          to markers that are of the core bin variety using an input core bin marker
#          file provided my maizegdb 

# TODO: Use ensembl api

# Options:
#          -f is the core bin marker file
#          -d is the database
#          -s is the server
#          -u is the username
#          -p is the password

# use the following modules                                                 
use Getopt::Std;
use DBI;

# get inputs                                                                    
my %opts = ();
getopts ('f:d:s:u:p:', \%opts);
my $file   = $opts{'f'};
my $database = $opts{'d'};
my $server = $opts{'s'};
my $user = $opts{'u'};
my $pass = $opts{'p'};

# connect to database
my $dbh = DBI -> connect ('dbi:mysql:'.$database.':'.$server,
			  $user,
			  $pass,
			  {
			      RaiseError => 1,
			      AutoCommit => 1,
			  }
			  ) || die "Database connection no good: $DBI::errstr";

# store the core bin marker data
my %corebinmarker;
my %corebinprobe;
open (corebin, "$file");
while (my $line = <corebin>){
    chomp $line;
    my ( $bin,
	 $marker,
	 $probe,
	 $type,
	 $insert,
	 $enzyme,
	 $sequence ) = split (/\t/, $line);
    
    $corebinmarker{$marker} = 0;
    $corebinprobe{$probe}   = 0;
}
close (corebin);

# get the existing analyses
(my $analysis_id_max, $analyses_ref) = &analysis_cache();
my %analyses = %$analyses_ref;

# update the new core bin analysis type if it does not already exist
unless (exists ($analyses{"core_bin_marker"})){
    $analysis_id_max++;
    my $insert_logic_sql =
	"insert into analysis values (?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
    my $insert_analysis_sth = $dbh -> prepare_cached ($insert_logic_sql);
    $insert_analysis_sth -> execute($analysis_id_max, 
				    'NULL', 
				    'core_bin_marker', 
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL',
				    'NULL');
}

# update the analyses
(my $analysis_id_max, $analyses_ref) = &analysis_cache();
my %analyses_updated = %$analyses_ref;

# create marker_feature query and execute
my $marker_feature_sql = "select * from marker_feature";
my $marker_feature_sth = $dbh -> prepare ($marker_feature_sql);
$marker_feature_sth -> execute();

# cycle through the marker_feature data to transform the core bin markers
my %analysis_tracker;
my %markers_seen;
while (my $markerrow_ref = $marker_feature_sth -> fetchrow_hashref()){
    my $marker_id = $markerrow_ref -> {marker_id};
    
#    print stderr "$marker_id\t";
    
    # create marker_synonym sql and execute for the marker_id above
    my $marker_syn_sql = "select * from marker_synonym where marker_id='$marker_id'";
    my $marker_syn_sth = $dbh -> prepare ($marker_syn_sql);
    $marker_syn_sth -> execute();
 
    # each marker_id currently has only a single syn, it's original name,
    # if and when this changes, the following will have to be changed

    # check to see if there is more than one synonym for the marker_id
    my @syn_data = @{$marker_syn_sth -> fetchall_arrayref()};
    my $syn_data_scalar = @syn_data;
    if ($syn_data_scalar > 1){
	print stderr "more than one synonym for this marker -- check\n";
	next;
    }

    # if only one syn, continue processing
    $marker_syn_sth -> execute();
    my $syn_data = $marker_syn_sth -> fetchrow_hashref();
    my $name = $syn_data -> {name};
    my $source = $syn_data -> {source};
   
#    print stderr "$name\t";
    
    # see if this marker name corresponds to any of those in the %corebinmarker
    my $matches = 0;
    
    if (exists ($corebinmarker{$name})){

	$corebinmarker{$name}++;
#	print stderr "matched core bin marker\n";

	my $update_analysis_sql = 
	    "update marker_feature set analysis_id ='$analyses_updated{'core_bin_marker'}' where marker_id ='$marker_id';";
	my $update_analysis_sth = $dbh -> prepare_cached ($update_analysis_sql);
#	$update_analysis_sth -> execute();
		
    }
    else {
#	print stderr "unmatched\n";
	next;
    }
}

foreach my $corebinmarker (keys %corebinmarker){
    print "$corebinmarker\t$corebinmarker{$corebinmarker}\n";
}

$dbh -> disconnect();


##subs##

# create hash of analysis types already available in the analysis table                                   
sub analysis_cache {
    my $analysis_sql = "select * from analysis;";
    my $analysis_sth = $dbh -> prepare ($analysis_sql);
    $analysis_sth -> execute();
    
    my %analyses;
    my $analysis_id_max = 0;
    while (my $analysisrow_ref = $analysis_sth -> fetchrow_hashref()){
	my $analysis_id = $analysisrow_ref -> {analysis_id};
	my $logic_name = $analysisrow_ref -> {logic_name};
	
	($analysis_id_max = $analysis_id) if ($analysis_id > $analysis_id_max);
	$analyses{$logic_name} = $analysis_id;
    }
    return ($analysis_id_max, \%analyses);
}
    
