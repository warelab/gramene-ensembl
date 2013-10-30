use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

my $database;
my $port = 3306;
my $user;
my $host;
my $pass;
my $upload = 0;
my $from_coord;
my $from_coord_version = undef;
my $to_coord;
my $to_coord_version = undef;

GetOptions('dbname=s'    => \$database,
           'port=s'        => \$port,
           'user=s'        => \$user,
           'host=s'        => \$host,
           'pass=s'        => \$pass,
           'from_coord=s'  => \$from_coord,
           'to_coord=s'    => \$to_coord, 
	   'from_coord_version=s' => \$from_coord_version,
	   'to_coord_version=s' => \$to_coord_version);

if(!defined($database) or !defined($user) or !defined($host) ){
  die "you must set database user host and pass\n";
}

#connect to the database

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (-host   => $host,
   -dbname  => $database,
   -group   => 'core',
   -user    => $user,
   -port    => $port,
   -pass    => $pass);

  my $sa = $db->get_SliceAdaptor();
  my $csa = $db->get_CoordSystemAdaptor();

# now project.
my $last_coord1_seq_region_id;
my $last_coord2_seq_region_id;
my $last_coord1_start;
my $last_coord2_start;
my $last_coord1_end;
my $last_coord2_end;
my $last_coord2_strand;

my $merge_count =0;

open (PROJ,">".$db->dbc->dbname()."_assembly_table.data") ||
  die "Could not open file ".$db->dbc->dbname()."_assembly_table.sql\n";

my $seperator = "|";

my $cs1_id = $csa->fetch_by_name($from_coord,$from_coord_version)->dbID;
print "coord id for $from_coord, $from_coord_version is $cs1_id\n";
my $sth = $sa->dbc->db_handle->prepare("select seq_region_id, name from seq_region where coord_system_id = ?");


$sth->execute($cs1_id);
my ($id, $name);
$sth->bind_columns(\$id,\$name);
my %coord1_2_seq_region_id;
while($sth->fetch){
  $coord1_2_seq_region_id{$name} = $id;
}

my $cs2_id = $csa->fetch_by_name($to_coord,$to_coord_version)->dbID;
print "coord id for $to_coord is $cs2_id\n";
$sth->execute($cs2_id);
$sth->bind_columns(\$id,\$name);
my %coord2_2_seq_region_id;
while($sth->fetch){
  $coord2_2_seq_region_id{$name} = $id;
}





foreach my $slice (@{$sa->fetch_all($from_coord,$from_coord_version)}){
  foreach my $seg (@{$slice->project($to_coord, $to_coord_version )}){
    my $clone = $seg->to_Slice();
    my $clone_start = $clone->start();
    my $clone_end = $clone->end();
    if($clone_start > $clone_end){
      my $temp = $clone_end;
      $clone_end = $clone_start;
      $clone_start = $temp;
    }
      if($last_coord1_seq_region_id eq $coord1_2_seq_region_id{$slice->seq_region_name()} and
         $last_coord2_seq_region_id eq $coord2_2_seq_region_id{$clone->seq_region_name()} and
         ($last_coord1_end+1) ==  $seg->from_start() and
         ($last_coord2_end+1) == $clone_start and
          ($last_coord2_strand ==  $clone->strand()) and
         (($seg->from_start()-$last_coord1_start) ==
           ($clone_end-$last_coord2_start)) and  # should be the same length
           $seperator ne "#"){
            $merge_count++;
            $last_coord1_end = $seg->from_end();
            $last_coord2_end = $clone_end;

      }
      else{
        if(defined($last_coord1_seq_region_id)){
          print PROJ  $last_coord1_seq_region_id,"\t",
                      $last_coord2_seq_region_id,"\t",
                      $last_coord1_start,"\t",$last_coord1_end,"\t",
                      $last_coord2_start,"\t",$last_coord2_end,"\t",
                      $last_coord2_strand,"\n";
        }

        $last_coord1_seq_region_id = $coord1_2_seq_region_id{$slice->seq_region_name()};
        $last_coord2_seq_region_id = $coord2_2_seq_region_id{$clone->seq_region_name()};
        $last_coord1_start = $seg->from_start();
        $last_coord1_end   =  $seg->from_end();
        $last_coord2_start = $clone_start;
        $last_coord2_end   = $clone_end;
        $last_coord2_strand= $clone->strand();
      }
    }
  }
print PROJ  $last_coord1_seq_region_id,"\t",
  $last_coord2_seq_region_id,"\t",
  $last_coord1_start,"\t",$last_coord1_end,"\t",
  $last_coord2_start,"\t",$last_coord2_end,"\t",
  $last_coord2_strand,"\n";

print "number of merges is $merge_count\n";


close PROJ;


#use the .sql to update the meta table
# it was assembly.default, I change it to assembly.mapping
print  "In mysql do :-";
print "LOAD DATA LOCAL INFILE '".$ENV{PWD}."/",$db->dbc->dbname()."_assembly_table.data' IGNORE INTO TABLE assembly;\n";
print "INSERT INTO meta (meta_key, meta_value) values('assembly.mapping','$from_coord:$from_coord_version#$to_coord:$to_coord_version')\n";
print "Possibaly, you need to do the following insertion too.\nINSERT INTO meta (meta_key, meta_value) values('assembly.mapping','chromosome:$to_coord_version|$to_coord:$to_coord_version|$from_coord:$from_coord_version')\n";
#use the .sql to update the meta table


##########################script end#####################################################

  
