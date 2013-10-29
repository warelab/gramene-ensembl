#!/opt/local/bin/perl

# Purpose: This script interfaces with ensembl using dbi and assigns new analyses
#          to markers given a marker-source correspondence file

# TODO: use ensembl api

# Options:
#          -f is a file of new analysis types and their regexes
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

# store the new analysis types and their regexes
my %regexes;
open (regex, "$file");
while (my $line = <regex>){
    chomp $line;
    (my $analysis, my $regex) = split (/\t/, $line);
    $regexes{$analysis}{$regex} = 1;
}
close (regex);

# get the existing analyses
(my $analysis_id_max, $analyses_ref) = &analysis_cache();
my %analyses = %$analyses_ref;

# update the new analysis types if they do not already exist
foreach my $analysis (keys %regexes){
    unless (exists ($analyses{$analysis})){
	$analysis_id_max++;
	my $insert_logic_sql =
	    "insert into analysis values (?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	my $insert_analysis_sth = $dbh -> prepare_cached ($insert_logic_sql);
	$insert_analysis_sth -> execute($analysis_id_max, 'NULL', $analysis, 'NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL');
    }
}

# update the analyses
(my $analysis_id_max, $analyses_ref) = &analysis_cache();
my %analyses_updated = %$analyses_ref;

# create marker_feature query and execute
my $marker_feature_sql = "select * from marker_feature";
my $marker_feature_sth = $dbh -> prepare ($marker_feature_sql);
$marker_feature_sth -> execute();

# cycle through the marker_feature data
my %analysis_tracker;
my %markers_seen;
while (my $markerrow_ref = $marker_feature_sth -> fetchrow_hashref()){
    my $marker_id = $markerrow_ref -> {marker_id};
    
    print stderr "$marker_id\t";
    
    # create marker_synonym sql and execute for the marker_id above
    my $marker_syn_sql = "select * from marker_synonym where marker_id='$marker_id'";
    my $marker_syn_sth = $dbh -> prepare ($marker_syn_sql);
    $marker_syn_sth -> execute();
 
    # each marker_id currently has only a single syn, it's original name,
    # if and when this changes, the following will have to turn into a while loop
    my $syn_data = $marker_syn_sth -> fetchrow_hashref();
    my $name = $syn_data -> {name};
    my $source = $syn_data -> {source};
   
    print stderr "$name\t";
    
    # cycle through the regexes and see where the bastard falls
    my $matches = 0;
    foreach my $analysis (keys %regexes){
	foreach my $regex (keys %{$regexes{$analysis}}){
	    if ($name =~m/$regex/){
		my $analysis_id = $analyses_updated{$analysis};
		my $update_analysis_sql = 
		    "update marker_feature set analysis_id ='$analysis_id' where marker_id ='$marker_id';";
		my $update_analysis_sth = $dbh -> prepare_cached ($update_analysis_sql);
		$update_analysis_sth -> execute();
		
		# iterate for end stats; once for each marker
		($analysis_tracker{$analysis}++) unless (exists ($markers_seen{$marker_id}));

		# signify how many matches for this marker, and warn 
		$matches++;
	    }
	    else {
		next;
	    }
	}
    }
    
    # warn about how many matches tallied for any given marker
    if ($matches == 0){
	print stderr "unmatched\n";
    }
    else {
	print stderr "$matches matches\n";
    }
    
    # put into seen category
    $markers_seen{$marker_id} = 1;
    
}

# print out the stats
foreach my $analysis (keys %analysis_tracker){
    print stderr "$analysis\t$analysis_tracker{$analysis}\n";
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
    
