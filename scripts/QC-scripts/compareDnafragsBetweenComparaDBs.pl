#!/usr/bin/env perl

use DBI;
use List::Compare;
use Data::Dump qw(dump);

my $user='weix';
my $pass='warelab';
my $host='cabot';

my %sqls;
$sqls{gdb} = 'select genome_db_id, name from genome_db where assembly_default=1';
$sqls{dnafrag} = "select dnafrag_id, name from dnafrag where genome_db_id = ?";


my $plant_compara_db = shift;
my $master_compara_db = shift;

my $dsn_plant = "dbi:mysql:${plant_compara_db}:$host";
my $dsn_master = "dbi:mysql:${master_compara_db}:$host";

my $dbh_plant = DBI->connect($dsn_plant, $user, $pass, 
			     { RaiseError => 1, AutoCommit => 0 }) or 
    die "cannot connect to $dsn_plant";

my $dbh_master = DBI->connect($dsn_master, $user, $pass) or 
    die "cannot connect to $dsn_master";

my $genome_db_ref = $dbh_plant->selectall_arrayref($sqls{gdb});

my $sth_df_plant = $dbh_plant->prepare($sqls{dnafrag}) or die "cannot prepare $sqls{dnafrag}";
my $sth_df_master = $dbh_master->prepare($sqls{dnafrag}) or die "cannot prepare $sqls{dnafrag}";;

my %gid2dnafrag;

for my $gdb( @{$genome_db_ref} ){

    print"$gdb->[0], $gdb->[1]\n";
    $sth_df_plant->execute($gdb->[0]) or die "cannot execute";

    #my @r;
    while(my @r = $sth_df_plant->fetchrow_array){

	#print dump(@r); #"$r[0], $r[1]";
	push @{$gid2dnafrag{$gdb->[0]}->{plant}}, [@r];
    }
    $sth_df_master->execute($gdb->[0]);

    while( my @r = $sth_df_master->fetchrow_array){
	push @{$gid2dnafrag{$gdb->[0]}->{master}}, [@r];
    }
    
    my @df_name_plant = map {$_->[1]} @{$gid2dnafrag{$gdb->[0]}->{plant}};
    my @df_name_master = map {$_->[1]} @{$gid2dnafrag{$gdb->[0]}->{master}};
    
    my $lc = List::Compare->new(\@df_name_plant, \@df_name_master);
    my @plant_only = $lc->get_unique;
    my @master_only = $lc->get_complement;
    my @intersection = $lc->get_intersection;

    print "There are ", scalar @intersection, " shared dnafrag names for $gdb->[1] (eg $intersection[0])\n";
    print "There are ", scalar @plant_only, " plant only dnafrag names for $gdb->[1] (eg $plant_only[0])\n";
    print "There are ", scalar @master_only, " master only dnafrag names for $gdb->[1] (eg $master_only[0])\n";
    #last;
}

$dbh_plant->disconnect;
$dbh_master->disconnect;
