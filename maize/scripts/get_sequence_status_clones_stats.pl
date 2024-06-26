#!/usr/local/bin/perl -w 

# (1) Run program, issue: perl ./get_sequence_status_clones_stats.pl -e ./ensembl_clonepath.registry -s Zea_mays
# (2) Run progra, issue: nice nohup perl ./get_sequence_status_clones_stats.pl -e ./ensembl_clonepath.registry -s Zea_mays >& get_sequence_status_clones_stats.out &
=pod

=head1 NAME

get_sequence_status_clones_stats.pl  - get statistics of sequence status clones Ensembl database

=head1 SYNOPSIS

perl get_sequence_status_clones_stats.pl [options]

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -f|--file             output file name

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s|--species>
  Use this species entry from the registry file [REQUIRED].

B<-f|--file>
  Use this file to name output filename

=head1 DESCRIPTION

  Add description here ...

=cut

use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

# Get today's date
my ($day, $month, $year) = (localtime)[3, 4, 5];
my $date = sprintf("%02d-%02d-%04d",$month+1, $day, $year+1900);

# Global variables
use vars qw( $ENS_DBA $file $date);

# Script arguments processing
BEGIN{
    my ($help, $man);
    my ($species, $reg);
   GetOptions( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "file=s"             => \$file
	       ) or pod2usage(2);
    pod2usage(-verbose => 2) if $man;
    pod2usage(1) if $help;

   # Load the ensembl file
    $species || ( print( "Need a --species\n" ) && pod2usage() );
    $reg    ||= './conf/ensembl.registry';
    -e $reg || ( print( "File $reg does not exist\n" )    && pod2usage(1) );
    -r $reg || ( print( "Cannot read $reg\n" )            && pod2usage(1) );
    -f $reg || ( print( "File $reg is not plain-text\n" ) && pod2usage(1) );
    -s $reg || ( print( "File $reg is empty\n" )          && pod2usage(1) );
    Bio::EnsEMBL::Registry->load_all( $reg ); # Connect to the database

   # Get DBAdaptor
    $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
    $ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) && pod2usage() );

    print( "Stats of Sequence Status Clones for $species\n" );
    my $pre_text = "  Target DB: ";
   # Use DBAdaptor's dbc method to get DBConnection, then use DBConnection's dbname, host and port methods 
   # to get database name, host and port
    foreach my $dba( $ENS_DBA ){
      print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
	     ":". $dba->dbc->port ."\n");
  }
} # end BEGIN

# Get stats for FULLTOP clones
my ($fulltop_clones, $fulltop_bases, $fulltop_min_bases, $fulltop_max_bases, $fulltop_ave_bases, $fulltop_median_bases, $fulltop_stddev_bases
   ) = &get_sequence_status_clones_stats('FULLTOP');

# Get stats for PREFIN clones                                                                                                                                                        
my ($prefin_clones, $prefin_bases, $prefin_min_bases, $prefin_max_bases, $prefin_ave_bases, $prefin_median_bases, $prefin_stddev_bases
    ) = &get_sequence_status_clones_stats('PREFIN');

# Get stats for ACTIVEFIN clones                                                                                                                                                     
my ($activefin_clones, $activefin_bases, $activefin_min_bases, $activefin_max_bases, $activefin_ave_bases, $activefin_median_bases, $activefin_stddev_bases
    ) = &get_sequence_status_clones_stats('ACTIVEFIN');

# Get stats for IMPROVED clones
my ($improved_clones, $improved_bases, $improved_min_bases, $improved_max_bases, $improved_ave_bases, $improved_median_bases, $improved_stddev_bases
    ) = &get_sequence_status_clones_stats('IMPROVED');

# Create stats table of sequence status clones
&create_stats_table($fulltop_clones, $fulltop_bases, $fulltop_min_bases, $fulltop_max_bases, $fulltop_ave_bases, $fulltop_median_bases, $fulltop_stddev_bases,
		    $prefin_clones, $prefin_bases, $prefin_min_bases, $prefin_max_bases, $prefin_ave_bases, $prefin_median_bases, $prefin_stddev_bases,
		    $activefin_clones, $activefin_bases, $activefin_min_bases, $activefin_max_bases, $activefin_ave_bases, $activefin_median_bases, $activefin_stddev_bases,
		    $improved_clones, $improved_bases, $improved_min_bases, $improved_max_bases, $improved_ave_bases, $improved_median_bases, $improved_stddev_bases
                   );

exit;

#================================================
#SUBROUTINES SECTION
#================================================
sub get_sequence_status_clones_stats {
    my ($clone_seq_status) = shift || die( "Need clone sequence status name" );
   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
   my $sth;

   # Get stats for the specified sequence status clone:
   my $clones_sql = "select count(*) from clone where status='".$clone_seq_status."'"; 
   my $bases_sql = "select sum(sequence_length) from clone where status='".$clone_seq_status."'";
   my $min_bases_sql = "select min(sequence_length) from clone where status='".$clone_seq_status."'";
   my $max_bases_sql = "select max(sequence_length) from clone where status='".$clone_seq_status."'";
   my $ave_bases_sql = "select avg(sequence_length) from clone where status='".$clone_seq_status."'";
   my $stddev_bases_sql = "select stddev(sequence_length) from clone where status='".$clone_seq_status."'";

   $sth = $ENS_DBA->dbc->prepare($clones_sql);
   $sth->execute();
   my $clones = $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($bases_sql);
   $sth->execute();
   my $bases = $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($min_bases_sql);
   $sth->execute();
   my $min_bases = $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($max_bases_sql);
   $sth->execute();
   my $max_bases = $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($ave_bases_sql);
   $sth->execute();
   my $ave_bases = sprintf "%.0f", $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($stddev_bases_sql);
   $sth->execute();
   my $stddev_bases = sprintf "%.0f", $sth->fetchrow_array();
   #================ 
   # Find the median:
   #================
   my ($median_bases_sql, $median_bases);
   my ($limit_value, @row);
   my @seq_length = ();
   my $seq_length_1 = 0;
   my $seq_length_2 = 0;
   # (1) If total clones is even:
   if ($clones % 2 == 0) {
      $limit_value = ($clones / 2) + 1;
      $median_bases_sql = "select sequence_length from clone where status='".$clone_seq_status."' ".
                          "order by sequence_length limit ".$limit_value;
      $sth = $ENS_DBA->dbc->prepare($median_bases_sql);
      $sth->execute();
      while ( @row = $sth->fetchrow_array() ) {
         push(@seq_length, $row[0]);
      }
      $seq_length_1 = pop(@seq_length);
      $seq_length_2 = pop(@seq_length);
      $median_bases = sprintf "%.0f", ($seq_length_1 + $seq_length_2) / 2;
   } # end if total clones is even
   # (2) if total clones is odd:
   else {
      $limit_value = ($clones + 1) / 2;
      $median_bases_sql = "select sequence_length from clone where status='".$clone_seq_status."' ".
                          "order by sequence_length limit ".$limit_value;
      $sth = $ENS_DBA->dbc->prepare($median_bases_sql);
      $sth->execute();
      while ( @row = $sth->fetchrow_array() ) {
         push(@seq_length, $row[0]);
      }
      $median_bases = sprintf "%.0f", pop(@seq_length);   
   } # end if total clones is odd

   $sth->finish();

   return ($clones, $bases, $min_bases, $max_bases, $ave_bases, $median_bases, $stddev_bases
          );     
}

sub create_stats_table {
   my ($fulltop_clones, $fulltop_bases, $fulltop_min_bases, $fulltop_max_bases, $fulltop_ave_bases, $fulltop_median_bases, $fulltop_stddev_bases,
       $prefin_clones, $prefin_bases, $prefin_min_bases, $prefin_max_bases, $prefin_ave_bases, $prefin_median_bases, $prefin_stddev_bases,
       $activefin_clones, $activefin_bases, $activefin_min_bases, $activefin_max_bases, $activefin_ave_bases, $activefin_median_bases, $activefin_stddev_bases,
       $improved_clones, $improved_bases, $improved_min_bases, $improved_max_bases, $improved_ave_bases, $improved_median_bases, $improved_stddev_bases) = @_;
      
   # get total clones and bases
   my $total_clones = $fulltop_clones + $prefin_clones + $activefin_clones + $improved_clones;
   my $total_bases = $fulltop_bases + $prefin_bases + $activefin_bases + $improved_bases;

   if (!$file) { $file = 'sequence_status_clones_stats'; }
   my $html_filename = $file."_".$date.".html";
   open (HTMLOUT, ">$html_filename") || die("Cannot Open HTML File");
   # Write to html file
   # print table header
   print(HTMLOUT "<center>"."\n");
   print(HTMLOUT "<table border=1 cellspacing=1 cellpadding=1>"."\n");
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<th colspan=8 align=center>Statistics of Clones by Sequence Status (".$date.")</th>"."\n");
   print(HTMLOUT "</tr>"."\n");
  
   # print table column headers
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<th valign=bottom>Status</th>"."\n");
   print(HTMLOUT "<th valign=bottom>Clones</th>"."\n");
   print(HTMLOUT "<th valign=bottom>Bases</th>"."\n");
   print(HTMLOUT "<th valign=bottom>Minimum (bases)</th>"."\n");
   print(HTMLOUT "<th valign=bottom>Maximum (bases)</th>"."\n");
   print(HTMLOUT "<th valign=bottom>Average (bases)</th>"."\n");
   print(HTMLOUT "<th valign=bottom>Median (bases)</th>"."\n");
   print(HTMLOUT "<th valign=bottom>Standard Deviation (bases)</th>"."\n");
   print(HTMLOUT "</tr>"."\n");

   # print FULLTOP row
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>FULLTOP</td>"."\n");
   print(HTMLOUT "<td>".$fulltop_clones."\n");
   print(HTMLOUT "<td>".$fulltop_bases."\n");
   print(HTMLOUT "<td>".$fulltop_min_bases."\n");
   print(HTMLOUT "<td>".$fulltop_max_bases."\n");
   print(HTMLOUT "<td>".$fulltop_ave_bases."\n");
   print(HTMLOUT "<td>".$fulltop_median_bases."\n");
   print(HTMLOUT "<td>".$fulltop_stddev_bases."\n");
   print(HTMLOUT "</tr>"."\n");
 
   # print PREFIN row
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>PREFIN</td>"."\n");
   print(HTMLOUT "<td>".$prefin_clones."\n");
   print(HTMLOUT "<td>".$prefin_bases."\n");
   print(HTMLOUT "<td>".$prefin_min_bases."\n");
   print(HTMLOUT "<td>".$prefin_max_bases."\n");
   print(HTMLOUT "<td>".$prefin_ave_bases."\n");
   print(HTMLOUT "<td>".$prefin_median_bases."\n");
   print(HTMLOUT "<td>".$prefin_stddev_bases."\n");
   print(HTMLOUT "</tr>"."\n");

   # print ACTIVEFIN row 
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>ACTIVEFIN</td>"."\n");
   print(HTMLOUT "<td>".$activefin_clones."\n");
   print(HTMLOUT "<td>".$activefin_bases."\n");
   print(HTMLOUT "<td>".$activefin_min_bases."\n");
   print(HTMLOUT "<td>".$activefin_max_bases."\n");
   print(HTMLOUT "<td>".$activefin_ave_bases."\n");
   print(HTMLOUT "<td>".$activefin_median_bases."\n");
   print(HTMLOUT "<td>".$activefin_stddev_bases."\n");
   print(HTMLOUT "</tr>"."\n");

   # print IMPROVED row
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>IMPROVED</td>"."\n");
   print(HTMLOUT "<td>".$improved_clones."\n");
   print(HTMLOUT "<td>".$improved_bases."\n");
   print(HTMLOUT "<td>".$improved_min_bases."\n");
   print(HTMLOUT "<td>".$improved_max_bases."\n");
   print(HTMLOUT "<td>".$improved_ave_bases."\n");
   print(HTMLOUT "<td>".$improved_median_bases."\n");
   print(HTMLOUT "<td>".$improved_stddev_bases."\n");
   print(HTMLOUT "</tr>"."\n");

   # print Total row
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<th>Total</th>"."\n");
   print(HTMLOUT "<td>".$total_clones."\n");
   print(HTMLOUT "<td>".$total_bases."\n");
   print(HTMLOUT "<td><br>"."\n");
   print(HTMLOUT "<td><br>"."\n");
   print(HTMLOUT "<td><br>"."\n");
   print(HTMLOUT "<td><br>"."\n");
   print(HTMLOUT "<td><br>"."\n");
   print(HTMLOUT "</tr>"."\n");   

   print(HTMLOUT "</table>"."\n");
   print(HTMLOUT "</center>"."\n");

   # print descriptions of FULLTOP, PREFIN, ACTIVEFIN, and IMPROVED
   print(HTMLOUT "<center>"."\n");
   print(HTMLOUT "<table>"."\n");
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>"."\n");
   print(HTMLOUT "<ul>"."\n");
   print(HTMLOUT "<li>"."\n");
   print(HTMLOUT "FULLTOP = 2 x 384 paired-end attempts (6X coverage); completed shotgun phase; initial assembly."."\n");
   print(HTMLOUT "</li>"."\n");
   print(HTMLOUT "<li>"."\n");
   print(HTMLOUT "PREFIN = Completed automated improvement phase (AutoFinish)."."\n");
   print(HTMLOUT "</li>"."\n");
   print(HTMLOUT "<li>"."\n");
   print(HTMLOUT "ACTIVEFIN = Active work being done by a finisher."."\n");
   print(HTMLOUT "</li>"."\n");
   print(HTMLOUT "<li>"."\n");
   print(HTMLOUT "IMPROVED = Finished sequence in gene regions; improved regions will be indicated; ".
         "once order and orientation of improved segments are confirmed, ".
         "a comment will be added to indicate this."."\n");
   print(HTMLOUT "</li>"."\n");
   print(HTMLOUT "</ul>"."\n");
   print(HTMLOUT "</td>"."\n");
   print(HTMLOUT "</tr>"."\n");
   print(HTMLOUT "</table>"."\n");
   print(HTMLOUT "</center>"."\n");

}

1;
