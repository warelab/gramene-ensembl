#!/opt/local/bin/perl

# Purpose: This script interfaces with ensembl using dbi and assigns new core bin positions
#          to misc features that are of the core bin variety using an input ranges file

# TODO: Use ensembl api

# Options:
#          -f is the ranges file
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

# store the ranges data
my %corebinstart;
my %corebinend;
open (VBINS, "$file");
while (my $line = <VBINS>){
    chomp $line;
    (my $bin, my $range) = split (/\t/, $line);
    (my $start, my $end) = split (/\-/, $range);
    (my $chr, my $binner) = split (/\./, $bin);

    $corebinstart{$bin} = $start;
    $corebinend{$bin}   = $end;
    
}
close (VBINS);

# create misc_feature query and execute
my $misc_feature_sql = "select mf.misc_feature_id, mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.seq_region_strand  from misc_feature mf, misc_feature_misc_set mm where mf.misc_feature_id=mm.misc_feature_id and misc_set_id=4";

my $misc_feature_sth = $dbh -> prepare ($misc_feature_sql);
$misc_feature_sth -> execute();

# cycle through the virtual core bins and apply the new ranges

while (my $markerrow_ref = $misc_feature_sth -> fetchrow_hashref()){
    my $mf_id = $markerrow_ref -> {misc_feature_id};
    my $seq_region_start = $markerrow_ref -> {seq_region_start};
    my $seq_region_end = $markerrow_ref -> {seq_region_end};
    
    # create misc_attrib sql and execute for the mf_id above
    my $ma_sql = "select * from misc_attrib where misc_feature_id='$mf_id'";
    my $ma_sth = $dbh -> prepare ($ma_sql);
    $ma_sth -> execute();
 
    my $ma_data = $ma_sth -> fetchrow_hashref();
    my $value = $ma_data -> {value};
   
    # see if any updates are necessary
    if (exists ($corebinstart{$value})){

	my $update_start_sql = 
	    "update misc_feature set seq_region_start ='$corebinstart{$value}' where misc_feature_id ='$mf_id';";
	my $update_end_sql =
	    "update misc_feature set seq_region_end ='$corebinend{$value}' where misc_feature_id ='$mf_id';";
	my $update_start_sth = $dbh -> prepare_cached ($update_start_sql);
	my $update_end_sth   = $dbh -> prepare_cached ($update_end_sql);

	$update_start_sth -> execute();
	$update_end_sth -> execute();
		
    }
    else {

	next;
    }
}

$dbh -> disconnect();


