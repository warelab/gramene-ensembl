#!/usr/local/bin/perl -w 

=pod

=head1 NAME

generate_QC_report.pl  - generate QC report from Ensembl database

=head1 SYNOPSIS

perl generate_QC_report.pl [options]

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

    print( "QC report for $species\n" );
    my $pre_text = "  Target DB: ";
   # Use DBAdaptor's dbc method to get DBConnection, then use DBConnection's dbname, host and port methods 
   # to get database name, host and port
    foreach my $dba( $ENS_DBA ){
      print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
	     ":". $dba->dbc->port ."\n");
  }
} # end BEGIN

# Get clone statistics
my ($total_clones,
    $clone_min_size, $clone_max_size,
    $clone_ave_size, $clone_median_size,
    $clone_stddev_size) = &get_clones_stats('clone');

# Get contigs statistics
my ($total_contigs,
    $contig_min_size, $contig_max_size,
    $contig_ave_size, $contig_median_size,
    $contig_stddev_size) = &get_contigs_stats('contig');


# Get contigs per clone statistics
my ($contigsperclone_min_size, $contigsperclone_max_size,
    $contigsperclone_ave_size, $contigsperclone_median_size,
    $contigsperclone_stddev_size) = &get_contigs_perclone_stats;

#Get gene model statistics
my ($total_FGENESH_gene_models,
    $genemodel_min_size, $genemodel_max_size,
    $genemodel_ave_size, $genemodel_median_size,
    $genemodel_stddev_size) = &get_gene_models_stats('Fgenesh');

# Get gene model classes statistics
my ($total_FGENESH_TE_gene_models, $percent_TE_gene_models,
    $min_TE_gene_model_size, $max_TE_gene_model_size, $ave_TE_gene_model_size, $median_TE_gene_model_size, $std_TE_gene_model_size,
    $total_FGENESH_WH_gene_models, $percent_WH_gene_models,
    $min_WH_gene_model_size, $max_WH_gene_model_size, $ave_WH_gene_model_size, $median_WH_gene_model_size, $std_WH_gene_model_size,
    $total_FGENESH_NH_gene_models, $percent_NH_gene_models,
    $min_NH_gene_model_size, $max_NH_gene_model_size, $ave_NH_gene_model_size, $median_NH_gene_model_size, $std_NH_gene_model_size
#    $total_FGENESH_corrupted_translation_gene_models, $percent_corrupted_translation_gene_models,
#    $min_corrupted_translation_gene_model_size, $max_corrupted_translation_gene_model_size, $ave_corrupted_translation_gene_model_size, 
#    $median_corrupted_translation_gene_model_size, $std_corrupted_translation_gene_model_size 
    ) = &get_gene_model_classes_stats($total_FGENESH_gene_models);

# Get nucleotides statistics for gene model classes 
my ($total_gene_model_nucleotides, 
    $total_TE_gene_model_nucleotides, $percent_TE_gene_model_nucleotide,
    $total_WH_gene_model_nucleotides, $percent_WH_gene_model_nucleotide,
    $total_NH_gene_model_nucleotides, $percent_NH_gene_model_nucleotide,
    $total_corrupted_translation_gene_model_nucleotides, $percent_corrupted_translation_gene_model_nucleotide
    ) = &get_gene_models_nucleotides_stats;

# Get total nucleotides of the current BAC clones
my $total_nucleotides = &get_total_nucleotides;

# Get nucleotides statistics for MDRs
my ($total_mdr_025_nucleotides, $percent_mdr_025_nucleotides,
    $total_mdr_1_nucleotides, $percent__mdr_1_nucleotides,
    $total_mdr_2_nucleotides, $percent_mdr_2_nucleotides,
    $total_mdr_3_nucleotides, $percent_mdr_3_nucleotides) = &get_mdr_stats($total_nucleotides);

# Get Interpro statistics (NOTE: $interpro is reference to a hash)
my ($interpro) = &get_interpro_stats($total_FGENESH_gene_models);

# Create QC report table
&create_QC_report_table($total_clones, $clone_min_size, $clone_max_size, $clone_ave_size, $clone_median_size, $clone_stddev_size,
                        $total_contigs, $contig_min_size, $contig_max_size, $contig_ave_size, $contig_median_size, $contig_stddev_size,
                        $contigsperclone_min_size, $contigsperclone_max_size, $contigsperclone_ave_size, $contigsperclone_median_size, $contigsperclone_stddev_size,
			$total_FGENESH_gene_models, $genemodel_min_size, $genemodel_max_size, $genemodel_ave_size, $genemodel_median_size, $genemodel_stddev_size,
                        $total_FGENESH_TE_gene_models, $percent_TE_gene_models,
			$min_TE_gene_model_size, $max_TE_gene_model_size, $ave_TE_gene_model_size, $median_TE_gene_model_size, $std_TE_gene_model_size,
			$total_FGENESH_WH_gene_models, $percent_WH_gene_models,
			$min_WH_gene_model_size, $max_WH_gene_model_size, $ave_WH_gene_model_size, $median_WH_gene_model_size, $std_WH_gene_model_size,
			$total_FGENESH_NH_gene_models, $percent_NH_gene_models,
			$min_NH_gene_model_size, $max_NH_gene_model_size, $ave_NH_gene_model_size, $median_NH_gene_model_size, $std_NH_gene_model_size,
#			$total_FGENESH_corrupted_translation_gene_models, $percent_corrupted_translation_gene_models,
#			$min_corrupted_translation_gene_model_size, $max_corrupted_translation_gene_model_size, $ave_corrupted_translation_gene_model_size,
#			$median_corrupted_translation_gene_model_size, $std_corrupted_translation_gene_model_size,
			$total_gene_model_nucleotides,
			$total_TE_gene_model_nucleotides, $percent_TE_gene_model_nucleotide,
			$total_WH_gene_model_nucleotides, $percent_WH_gene_model_nucleotide,
			$total_NH_gene_model_nucleotides, $percent_NH_gene_model_nucleotide,
			$total_corrupted_translation_gene_model_nucleotides, $percent_corrupted_translation_gene_model_nucleotide,
			$total_nucleotides,
			$total_mdr_025_nucleotides, $percent_mdr_025_nucleotides,
			$total_mdr_1_nucleotides, $percent__mdr_1_nucleotides,
			$total_mdr_2_nucleotides, $percent_mdr_2_nucleotides,
			$total_mdr_3_nucleotides, $percent_mdr_3_nucleotides,
			$interpro
		       );

exit;

#================================================
#SUBROUTINES SECTION
#================================================
sub get_clones_stats {
    my ($coord_system_name) = shift || die( "Need coord system name" );
   
   # Get clone coord_system_id
    my $coord_system_id = &get_coord_system_id($coord_system_name);
 
   # Queries to find (1) total number, (2) min length, (3) max length, (4) ave length and (5) stddev 
   # of current BAC clones
   my $current_clones_sql = "(select s.seq_region_id, s.name, s.length AS LENGTH from ".
                            "seq_region s, seq_region_attrib sa ".
                            "where s.seq_region_id= sa.seq_region_id and s.coord_system_id=".$coord_system_id.
                            " and sa.value='current') tmpTable";
    my $total_sql = "select count(*) from ".$current_clones_sql;    
    my $min_sql = "select min(LENGTH) from ".$current_clones_sql;
    my $max_sql = "select max(LENGTH) from ".$current_clones_sql;
    my $ave_sql = "select avg(LENGTH) from ".$current_clones_sql;
    my $stddev_sql = "select stddev(LENGTH) from ".$current_clones_sql;
   
   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
    my $sth;

    $sth = $ENS_DBA->dbc->prepare($total_sql);
    $sth->execute();
    my $total = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($min_sql);
    $sth->execute();
    my $min = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($max_sql);
    $sth->execute();
    my $max = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($ave_sql);
    $sth->execute();
    my $ave = sprintf "%.0f", $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($stddev_sql);
    $sth->execute();
    my $stddev = sprintf "%.0f", $sth->fetchrow_array();
   
   #================ 
   # Find the median:
   #================

    my ($median_sql, $limit_value, @row, $median);
    my @length = ();
    my $lentgh_1 = 0; 
    my $lentgh_2 = 0; 
   # (2) If the total_entries is odd, :
    if ($total % 2 == 1) {
        $limit_value = ($total + 1) / 2; 
      $median_sql = "select s.seq_region_id, s.name, s.length AS LENGTH from seq_region s, seq_region_attrib sa ".
                    "where s.seq_region_id= sa.seq_region_id and s.coord_system_id=".$coord_system_id.
                    " and sa.value='current' order by s.length limit ".$limit_value;


        $sth = $ENS_DBA->dbc->prepare($median_sql);
        $sth->execute();
        while ( @row = $sth->fetchrow_array() ) {
            push(@length, $row[2]);
        }
        $median = sprintf "%.0f", pop(@length);
    }# end if total_entries is odd
         
   # (3) if the total entries is even
        else {
            $limit_value = ($total / 2) + 1;
      $median_sql = "select s.seq_region_id, s.name, s.length AS LENGTH from seq_region s, seq_region_attrib sa ".
                    "where s.seq_region_id= sa.seq_region_id and s.coord_system_id=".$coord_system_id.
                    " and sa.value='current' order by s.length limit ".$limit_value;

            $sth = $ENS_DBA->dbc->prepare($median_sql);
            $sth->execute();
            while ( @row = $sth->fetchrow_array() ) {
                push(@length, $row[2]);
            }
            $lentgh_1 = pop(@length);
            $lentgh_2 = pop(@length);
            $median = sprintf "%.0f", ($lentgh_1 + $lentgh_2) / 2;
        }
    $sth->finish();
 
   return ($total, $min, $max,
           $ave, $median, 
           $stddev);
}

sub get_contigs_stats {
    my ($coord_system_name) = shift || die( "Need coord system name" );
   
   # Get clone coord_system_id
    my $coord_system_id = &get_coord_system_id($coord_system_name);
 
   # Queries to find (1) total number, (2) min length, (3) max length, (4) ave length and (5) stddev of contigs
   my $contigs_sql = "seq_region a left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                     "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                     "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                     "where a.coord_system_id = ".
                     "(select coord_system_id from coord_system where name = 'clone') and ".
                     "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                     "v.value = 'current'";

    my $total_sql = "select count(*) from ".$contigs_sql;
    my $min_sql =  "select min(f.length) from ".$contigs_sql;
    my $max_sql = "select max(f.length) from ".$contigs_sql;
    my $ave_sql = "select avg(f.length) from ".$contigs_sql;
    my $stddev_sql = "select stddev(f.length) from ".$contigs_sql;

   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
    my $sth;
    $sth = $ENS_DBA->dbc->prepare($total_sql);
    $sth->execute();
    my $total = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($min_sql);
    $sth->execute();
    my $min = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($max_sql);
    $sth->execute();
    my $max = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($ave_sql);
    $sth->execute();
    my $ave = sprintf "%.0f", $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($stddev_sql);
    $sth->execute();
    my $stddev = sprintf "%.0f", $sth->fetchrow_array();

   
   #================ 
   # Find the median:
   #================
   my ($median_sql, $limit_value, @row, $median);
   my @length = ();
   my $lentgh_1 = 0; 
   my $lentgh_2 = 0;
   # (1) If total contigs is odd:
   if ($total % 2 == 1) {
      $limit_value = ($total + 1) / 2;
      $median_sql = "select f.length from ".$contigs_sql." order by f.length limit ".$limit_value;
      $sth = $ENS_DBA->dbc->prepare($median_sql);
      $sth->execute();
      while ( @row = $sth->fetchrow_array() ) {
	 push(@length, $row[0]);
         
      }
      $median = sprintf "%.0f", pop(@length);
   }# end if total contigs is odd    
   # (2) If total contigs is even:
   else {
      $limit_value = ($total / 2) + 1;
      $median_sql = "select f.length from ".$contigs_sql." order by f.length limit ".$limit_value;
      $sth = $ENS_DBA->dbc->prepare($median_sql);
      $sth->execute();
      while ( @row = $sth->fetchrow_array() ) {
	 push(@length, $row[0]);
      }
      $lentgh_1 = pop(@length);
      $lentgh_2 = pop(@length);
      $median = sprintf "%.0f", ($lentgh_1 + $lentgh_2) / 2;  
   }# end if total contigs is even 

   $sth->finish();
   return ($total, $min, $max, $ave, $median, $stddev);
}

sub get_contigs_perclone_stats {
   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle    
   my $contigs_perclone_sql = "(select count(*) AS COUNTING from seq_region a ".
                              "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                              "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                              "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                              "where a.coord_system_id = ".
                              "(select coord_system_id from coord_system where name = 'clone') and ".
                              "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                              "v.value = 'current' group by asm.asm_seq_region_id) tmpTable";

   my $min_sql = "select min(COUNTING) from ".$contigs_perclone_sql; 

   my $max_sql = "select max(COUNTING) from ".$contigs_perclone_sql;
    
   my $ave_sql = "select avg(COUNTING) from ".$contigs_perclone_sql;

   my $stddev_sql = "select stddev(COUNTING) from ".$contigs_perclone_sql; 

   my $sth;
   $sth = $ENS_DBA->dbc->prepare($min_sql);
   $sth->execute();
   my $min = $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($max_sql);
   $sth->execute();
   my $max = $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($ave_sql);
   $sth->execute();
   my $ave = sprintf "%.0f", $sth->fetchrow_array();

   $sth = $ENS_DBA->dbc->prepare($stddev_sql);
   $sth->execute();
   my $stddev = sprintf "%.0f", $sth->fetchrow_array();
   
   #================ 
   # Find the median:
   #================
   # (1) Find total entries
   my ($total_entries, $total_entries_sql);
   $total_entries_sql = "select count(COUNTING) from ".
                        "(select count(*) AS COUNTING from seq_region a ".
                        "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                        "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                        "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                        "where a.coord_system_id = ".
                        "(select coord_system_id from coord_system where name = 'clone') and ".
                        "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                        "v.value = 'current' group by asm.asm_seq_region_id) tmpTable";

   $sth = $ENS_DBA->dbc->prepare($total_entries_sql);
   $sth->execute();
   $total_entries = $sth->fetchrow_array();
   
   my ($median_sql, $limit_value, @row, $median);
   my @contigs_per_clone = ();
   my $contigs_per_clone_1 = 0;
   my $contigs_per_clone_2 = 0;
   # (2) If the total entries is even:
   if ($total_entries % 2 == 0) {
      $limit_value = ($total_entries / 2) + 1;
      $median_sql = "select count(*) AS COUNTING from seq_region a ".
                    "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                    "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                    "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                    "where a.coord_system_id = ".
                    "(select coord_system_id from coord_system where name = 'clone') and ".
                    "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                    "v.value = 'current' group by asm.asm_seq_region_id order by COUNTING limit ".$limit_value;

       $sth = $ENS_DBA->dbc->prepare($median_sql);
       $sth->execute();
       while ( @row = $sth->fetchrow_array() ) {
           push(@contigs_per_clone, $row[0]);
       }
       $contigs_per_clone_1 = pop(@contigs_per_clone);
       $contigs_per_clone_2 = pop(@contigs_per_clone);
       $median = sprintf "%.0f", ($contigs_per_clone_1 + $contigs_per_clone_2) / 2;
   }  # end if total_entries number is even
   # (3) if the total entries is odd:
   else {
      $limit_value = ($total_entries + 1) / 2;
      $median_sql = "select count(*) AS COUNTING from seq_region a ".
                    "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                    "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                    "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                    "where a.coord_system_id = ".
                    "(select coord_system_id from coord_system where name = 'clone') and ".
                    "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                    "v.value = 'current' group by asm.asm_seq_region_id order by COUNTING limit ".$limit_value;

       $sth = $ENS_DBA->dbc->prepare($median_sql);
       $sth->execute();
       while ( @row = $sth->fetchrow_array() ) {
           push(@contigs_per_clone, $row[0]);
       }
       $median = sprintf "%.0f", pop(@contigs_per_clone);
   } # end if total_entries number is odd     

   $sth->finish();
   return ($min, $max, $ave, $median, $stddev);
 
}

sub get_gene_models_stats {
    my $analysis_logic_name = shift || die( "Need analysis logic name" );
   # Get analysis_id for analysis
    my $analysis_id = &get_analysis_id($analysis_logic_name);

   # Queries to find (1) total gene models, (2) min, max, ave, median and stddev gene model size
   my $gene_models_sql = "gene where seq_region_id in (select f.seq_region_id from seq_region a ".
                         "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                         "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                         "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                         "where a.coord_system_id = ".
                         "(select coord_system_id from coord_system where name = 'clone') and ".
                         "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                         "v.value = 'current')";

    my $total_sql = "select count(*) from ".$gene_models_sql;
    my $min_sql = "select min(MODEL_LENGTH) from (select ((seq_region_end - seq_region_start) + 1) AS MODEL_LENGTH from ".$gene_models_sql.") tmpTable";
    my $max_sql = "select max(MODEL_LENGTH) from (select ((seq_region_end - seq_region_start) + 1) AS MODEL_LENGTH from ".$gene_models_sql.") tmpTable";
    my $ave_sql = "select avg(MODEL_LENGTH) from (select ((seq_region_end - seq_region_start) + 1) AS MODEL_LENGTH from ".$gene_models_sql.") tmpTable";
    my $stddev_sql = "select stddev(MODEL_LENGTH) from (select ((seq_region_end - seq_region_start) + 1) AS MODEL_LENGTH from ".$gene_models_sql.") tmpTable";
   
   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
    my $sth;
    $sth = $ENS_DBA->dbc->prepare($total_sql);
    $sth->execute();
    my $total = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($min_sql);
    $sth->execute();
    my $min = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($max_sql);
    $sth->execute();
    my $max = $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($ave_sql);
    $sth->execute();
    my $ave = sprintf "%.0f", $sth->fetchrow_array();

    $sth = $ENS_DBA->dbc->prepare($stddev_sql);
    $sth->execute();
    my $stddev = sprintf "%.0f", $sth->fetchrow_array();

   #================ 
   # Find the median:
   #================
    my ($median_sql, $limit_value, @row, $median);   
    my @gene_model_size = ();
    my $gene_model_size_1 = 0;
    my $gene_model_size_2 = 0;
   # (1) If total gene models is even:
    if ($total % 2 == 0) {
	$limit_value = ($total / 2) + 1;
	$median_sql = "select ((seq_region_end - seq_region_start) + 1) AS MODEL_LENGTH from ".$gene_models_sql." order by MODEL_LENGTH limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$gene_model_size_1 = pop(@gene_model_size);
	$gene_model_size_2 = pop(@gene_model_size);
	$median = sprintf "%.0f", ($gene_model_size_1 + $gene_model_size_2) / 2;
    } # end if total gene models is even

   # (2) If total gene models is odd:
    else {
	$limit_value = ($total + 1) / 2;
	$median_sql = "select ((seq_region_end - seq_region_start) + 1) AS MODEL_LENGTH from ".$gene_models_sql." order by MODEL_LENGTH limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$median = sprintf "%.0f", pop(@gene_model_size);       
    } # end if total gene models is odd 


    $sth->finish();
    return ($total, $min, $max, $ave, $median, $stddev);
 
}

sub get_gene_model_classes_stats {
    my ($total_FGENESH_gene_models) = shift || die( "Need total FGENESH gene models" );
  
   my $gene_models_query = "gene where seq_region_id in ".
                           "(select f.seq_region_id from seq_region a ".
                           "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                           "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                           "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                           "where a.coord_system_id = ".
                           "(select coord_system_id from coord_system where name = 'clone') and ".
                           "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                           "v.value = 'current')"; 

   # Queries to get total TE, WH, NH, Corrupted_translation gene models
    my $total_FGENESH_TE_gene_models_query = "select count(*) from ".$gene_models_query." and biotype='transposon_pseudogene'";
    my $total_FGENESH_WH_gene_models_query = "select count(*) from ".$gene_models_query." and biotype='protein_coding'";
    my $total_FGENESH_NH_gene_models_query = "select count(*) from ".$gene_models_query." and biotype='protein_coding_unsupported'";
#    my $total_FGENESH_corrupted_translation_gene_models_query = "select count(*) from ".$gene_models_query." and biotype='corrupted_translation'";

   # Queries to get min, max, ave, std gene model size
    my $min_gene_model_size_sql = "select min(NT) from (select ((seq_region_end - seq_region_start) + 1) AS NT from ".$gene_models_query;
    my $max_gene_model_size_sql = "select max(NT) from (select ((seq_region_end - seq_region_start) + 1) AS NT from ".$gene_models_query;
    my $ave_gene_model_size_sql = "select avg(NT) from (select ((seq_region_end - seq_region_start) + 1) AS NT from ".$gene_models_query;
    my $std_gene_model_size_sql = "select stddev(NT) from (select ((seq_region_end - seq_region_start) + 1) AS NT from ".$gene_models_query;
   
   # Queries to get min, max, ave, std TE gene model size
    my $min_TE_gene_model_size_sql = $min_gene_model_size_sql." and biotype='transposon_pseudogene') tmpTable"; 
    my $max_TE_gene_model_size_sql = $max_gene_model_size_sql." and biotype='transposon_pseudogene') tmpTable";
    my $ave_TE_gene_model_size_sql = $ave_gene_model_size_sql." and biotype='transposon_pseudogene') tmpTable";
    my $std_TE_gene_model_size_sql = $std_gene_model_size_sql." and biotype='transposon_pseudogene') tmpTable";

   # Queries to get min, max, ave, std WH gene model size
    my $min_WH_gene_model_size_sql = $min_gene_model_size_sql." and biotype='protein_coding') tmpTable"; 
    my $max_WH_gene_model_size_sql = $max_gene_model_size_sql." and biotype='protein_coding') tmpTable";
    my $ave_WH_gene_model_size_sql = $ave_gene_model_size_sql." and biotype='protein_coding') tmpTable";
    my $std_WH_gene_model_size_sql = $std_gene_model_size_sql." and biotype='protein_coding') tmpTable";

   # Queries to get min, max, ave, std NH gene model size
    my $min_NH_gene_model_size_sql = $min_gene_model_size_sql." and biotype='protein_coding_unsupported') tmpTable";
    my $max_NH_gene_model_size_sql = $max_gene_model_size_sql." and biotype='protein_coding_unsupported') tmpTable";
    my $ave_NH_gene_model_size_sql = $ave_gene_model_size_sql." and biotype='protein_coding_unsupported') tmpTable";
    my $std_NH_gene_model_size_sql = $std_gene_model_size_sql." and biotype='protein_coding_unsupported') tmpTable";

   # Queries to get min, max, ave, std corrupted_translation model size
 #   my $min_corrupted_translation_gene_model_size_sql = $min_gene_model_size_sql." and biotype='corrupted_translation') tmpTable";
 #   my $max_corrupted_translation_gene_model_size_sql = $max_gene_model_size_sql." and biotype='corrupted_translation') tmpTable";
 #   my $ave_corrupted_translation_gene_model_size_sql = $ave_gene_model_size_sql." and biotype='corrupted_translation') tmpTable";
 #   my $std_corrupted_translation_gene_model_size_sql = $std_gene_model_size_sql." and biotype='corrupted_translation') tmpTable";    

   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
    my $sth;
   
   # Get TE gene model statistics
    $sth = $ENS_DBA->dbc->prepare($total_FGENESH_TE_gene_models_query);
    $sth->execute();
    my $total_FGENESH_TE_gene_models = $sth->fetchrow_array();
    my $percent_TE_gene_models = sprintf "%.0f", ($total_FGENESH_TE_gene_models / $total_FGENESH_gene_models) * 100;
    $sth = $ENS_DBA->dbc->prepare($min_TE_gene_model_size_sql);
    $sth->execute();
    my $min_TE_gene_model_size = $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($max_TE_gene_model_size_sql);
    $sth->execute();
    my $max_TE_gene_model_size = $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($ave_TE_gene_model_size_sql);
    $sth->execute();
    my $ave_TE_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($std_TE_gene_model_size_sql);
    $sth->execute();
    my $std_TE_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();   

   # Get WH gene model statistics
    $sth = $ENS_DBA->dbc->prepare($total_FGENESH_WH_gene_models_query);
    $sth->execute();
    my $total_FGENESH_WH_gene_models = $sth->fetchrow_array();
    my $percent_WH_gene_models = sprintf "%.0f", ($total_FGENESH_WH_gene_models / $total_FGENESH_gene_models) * 100;
    $sth = $ENS_DBA->dbc->prepare($min_WH_gene_model_size_sql);
    $sth->execute();
    my $min_WH_gene_model_size = $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($max_WH_gene_model_size_sql);
    $sth->execute();
    my $max_WH_gene_model_size = $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($ave_WH_gene_model_size_sql);
    $sth->execute();
    my $ave_WH_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($std_WH_gene_model_size_sql);
    $sth->execute();
    my $std_WH_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();
   
   # Get NH gene model statistics
    $sth = $ENS_DBA->dbc->prepare($total_FGENESH_NH_gene_models_query);
    $sth->execute();
    my $total_FGENESH_NH_gene_models = $sth->fetchrow_array();
    my $percent_NH_gene_models = sprintf "%.0f", ($total_FGENESH_NH_gene_models / $total_FGENESH_gene_models) * 100;
    $sth = $ENS_DBA->dbc->prepare($min_NH_gene_model_size_sql);
    $sth->execute();
    my $min_NH_gene_model_size = $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($max_NH_gene_model_size_sql);
    $sth->execute();
    my $max_NH_gene_model_size = $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($ave_NH_gene_model_size_sql);
    $sth->execute();
    my $ave_NH_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();
    $sth = $ENS_DBA->dbc->prepare($std_NH_gene_model_size_sql);
    $sth->execute();
    my $std_NH_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();

   # Get corrupted_translation gene model statistics
 #   $sth = $ENS_DBA->dbc->prepare($total_FGENESH_corrupted_translation_gene_models_query);
 #   $sth->execute();
 #   my $total_FGENESH_corrupted_translation_gene_models = $sth->fetchrow_array();
 #   my $percent_corrupted_translation_gene_models = sprintf "%.0f", ($total_FGENESH_corrupted_translation_gene_models / $total_FGENESH_gene_models) * 100;
 #   $sth = $ENS_DBA->dbc->prepare($min_corrupted_translation_gene_model_size_sql);
 #   $sth->execute();
 #   my $min_corrupted_translation_gene_model_size = $sth->fetchrow_array();
 #   $sth = $ENS_DBA->dbc->prepare($max_corrupted_translation_gene_model_size_sql);
 #   $sth->execute();
 #   my $max_corrupted_translation_gene_model_size = $sth->fetchrow_array();
 #   $sth = $ENS_DBA->dbc->prepare($ave_corrupted_translation_gene_model_size_sql);
 #   $sth->execute();
 #   my $ave_corrupted_translation_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();
 #   $sth = $ENS_DBA->dbc->prepare($std_corrupted_translation_gene_model_size_sql);
 #   $sth->execute();
 #   my $std_corrupted_translation_gene_model_size = sprintf "%.0f", $sth->fetchrow_array();

   #===========================
   # Get median gene model size
   #===========================
    my $median_gene_model_size_sql = "select ((seq_region_end - seq_region_start) + 1) AS NT from ".$gene_models_query;
    my ($limit_value, @row);
    my @gene_model_size = ();
    my $gene_model_size_1 = 0;
    my $gene_model_size_2 = 0;

   # Get median TE gene model size
    my ($median_TE_gene_model_size_sql, $median_TE_gene_model_size);
   # (1) If total TE gene models is even:
    if ($total_FGENESH_TE_gene_models % 2 == 0) { 
	$limit_value = ($total_FGENESH_TE_gene_models / 2) + 1;
	$median_TE_gene_model_size_sql = $median_gene_model_size_sql." and biotype='transposon_pseudogene' order by NT limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_TE_gene_model_size_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$gene_model_size_1 = pop(@gene_model_size);
	$gene_model_size_2 = pop(@gene_model_size);
	$median_TE_gene_model_size = sprintf "%.0f", ($gene_model_size_1 + $gene_model_size_2) / 2;
    } # end if total TE gene models is even
   # (2) If total TE gene models is odd:
    else {
	$limit_value = ($total_FGENESH_TE_gene_models + 1) / 2;
	$median_TE_gene_model_size_sql = $median_gene_model_size_sql." and biotype='transposon_pseudogene' order by NT limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_TE_gene_model_size_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$median_TE_gene_model_size = sprintf "%.0f", pop(@gene_model_size);                
    } # end if total TE gene models is odd

   # Get median WH gene model size
    my ($median_WH_gene_model_size_sql, $median_WH_gene_model_size);
   # (1) If total WH gene models is even:
    if ($total_FGENESH_WH_gene_models % 2 == 0) {
	$limit_value = ($total_FGENESH_WH_gene_models / 2) + 1;
	$median_WH_gene_model_size_sql = $median_gene_model_size_sql." and biotype='protein_coding' order by NT limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_WH_gene_model_size_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$gene_model_size_1 = pop(@gene_model_size);
	$gene_model_size_2 = pop(@gene_model_size);
	$median_WH_gene_model_size = sprintf "%.0f", ($gene_model_size_1 + $gene_model_size_2) / 2;     
    } # end if total WH gene models is even
   # (2) If total WH gene models is odd:
    else {
	$limit_value = ($total_FGENESH_WH_gene_models + 1) / 2;
	$median_WH_gene_model_size_sql = $median_gene_model_size_sql." and biotype='protein_coding' order by NT limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_WH_gene_model_size_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$median_WH_gene_model_size = sprintf "%.0f", pop(@gene_model_size);   
    } # end if total WH gene models is odd

   # Get median NH gene model size
    my ($median_NH_gene_model_size_sql, $median_NH_gene_model_size);
   # (1) If total NH gene models is even:
    if ($total_FGENESH_NH_gene_models % 2 == 0) {
	$limit_value = ($total_FGENESH_NH_gene_models / 2) + 1;
	$median_NH_gene_model_size_sql = $median_gene_model_size_sql." and biotype='protein_coding_unsupported' order by NT limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_NH_gene_model_size_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$gene_model_size_1 = pop(@gene_model_size);
	$gene_model_size_2 = pop(@gene_model_size);
	$median_NH_gene_model_size = sprintf "%.0f", ($gene_model_size_1 + $gene_model_size_2) / 2;
    } # end if total NH gene models is even
   # (2) If total NH gene models is odd:
    else {
	$limit_value = ($total_FGENESH_NH_gene_models + 1) / 2;
	$median_NH_gene_model_size_sql = $median_gene_model_size_sql." and biotype='protein_coding_unsupported' order by NT limit ".$limit_value;
	$sth = $ENS_DBA->dbc->prepare($median_NH_gene_model_size_sql);
	$sth->execute();
	while ( @row = $sth->fetchrow_array() ) {
	    push(@gene_model_size, $row[0]);
	}
	$median_NH_gene_model_size = sprintf "%.0f", pop(@gene_model_size);      
    } # end if total NH gene models is odd

   # Get median corrupted_translation gene model size
#    my ($median_corrupted_translation_gene_model_size_sql, $median_corrupted_translation_gene_model_size);
   # (1) If total corrupted_translation gene models is even:
#    if ($total_FGENESH_corrupted_translation_gene_models % 2 == 0) {
#	$limit_value = ($total_FGENESH_corrupted_translation_gene_models / 2) + 1;
#	$median_corrupted_translation_gene_model_size_sql = $median_gene_model_size_sql." and biotype='corrupted_translation' order by NT limit ".$limit_value;
#	$sth = $ENS_DBA->dbc->prepare($median_corrupted_translation_gene_model_size_sql);
#	$sth->execute();
#	while ( @row = $sth->fetchrow_array() ) {
#	    push(@gene_model_size, $row[0]);
#	}
#	$gene_model_size_1 = pop(@gene_model_size);
#	$gene_model_size_2 = pop(@gene_model_size);
#	$median_corrupted_translation_gene_model_size = sprintf "%.0f", ($gene_model_size_1 + $gene_model_size_2) / 2;
#    } # end if corrupted_translation gene models is even
   # (2) If total corrupted_translation gene models is odd:
#    else {
#	$limit_value = ($total_FGENESH_corrupted_translation_gene_models + 1) / 2;
#	$median_corrupted_translation_gene_model_size_sql = $median_gene_model_size_sql." and biotype='corrupted_translation' order by NT limit ".$limit_value;
#	$sth = $ENS_DBA->dbc->prepare($median_corrupted_translation_gene_model_size_sql);
#	$sth->execute();
#	while ( @row = $sth->fetchrow_array() ) {
#	    push(@gene_model_size, $row[0]);
#	}
#	$median_corrupted_translation_gene_model_size = sprintf "%.0f", pop(@gene_model_size);  
#    } # end if total corrupted_translation gene models is odd

    $sth->finish();
   return ($total_FGENESH_TE_gene_models, $percent_TE_gene_models,
           $min_TE_gene_model_size, $max_TE_gene_model_size, $ave_TE_gene_model_size, $median_TE_gene_model_size, $std_TE_gene_model_size,
           $total_FGENESH_WH_gene_models, $percent_WH_gene_models,
           $min_WH_gene_model_size, $max_WH_gene_model_size, $ave_WH_gene_model_size, $median_WH_gene_model_size, $std_WH_gene_model_size,
           $total_FGENESH_NH_gene_models, $percent_NH_gene_models,
           $min_NH_gene_model_size, $max_NH_gene_model_size, $ave_NH_gene_model_size, $median_NH_gene_model_size, $std_NH_gene_model_size
#           $total_FGENESH_corrupted_translation_gene_models, $percent_corrupted_translation_gene_models, 
#           $min_corrupted_translation_gene_model_size, $max_corrupted_translation_gene_model_size, $ave_corrupted_translation_gene_model_size, 
#           $median_corrupted_translation_gene_model_size, $std_corrupted_translation_gene_model_size
          );
 
}

sub get_gene_models_nucleotides_stats {
   my $gene_model_nucleotides_query = "select ((seq_region_end - seq_region_start) + 1) AS NT from gene ".
                                      "where seq_region_id in ".
                                      "(select f.seq_region_id from seq_region a ".
                                      "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                                      "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                                      "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                                      "where a.coord_system_id = (select coord_system_id from coord_system where name = 'clone') and ".
                                      "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
                                      "v.value = 'current')";
   # Queries to get total gene model nucleotides and total nucleotides for TE, WH, NH and corrupted_translation gene models
   my $total_gene_model_nucleotides_query = "select sum(NT) from (".$gene_model_nucleotides_query.") tmpTable";
   my $total_TE_gene_model_nucleotides_query = "select sum(NT) from (".$gene_model_nucleotides_query." and biotype='transposon_pseudogene') tmpTable";
   my $total_WH_gene_model_nucleotides_query = "select sum(NT) from (".$gene_model_nucleotides_query." and biotype='protein_coding') tmpTable";
   my $total_NH_gene_model_nucleotides_query = "select sum(NT) from (".$gene_model_nucleotides_query." and biotype='protein_coding_unsupported') tmpTable";
   my $total_corrupted_translation_gene_model_nucleotides_query = "select sum(NT) from (".$gene_model_nucleotides_query." and biotype='corrupted_translation') tmpTable";
   
   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
   my $sth;
   
   # Get total gene model nucleotides
   $sth = $ENS_DBA->dbc->prepare($total_gene_model_nucleotides_query);
   $sth->execute();
   my $total_gene_model_nucleotides = $sth->fetchrow_array();
     
   # Get total TE gene model nucleotides
   $sth = $ENS_DBA->dbc->prepare($total_TE_gene_model_nucleotides_query);
   $sth->execute();
   my $total_TE_gene_model_nucleotides = $sth->fetchrow_array();
   my $percent_TE_gene_model_nucleotide = sprintf "%.0f", ($total_TE_gene_model_nucleotides / $total_gene_model_nucleotides) * 100;

   # Get total WH gene model nucleotides
   $sth = $ENS_DBA->dbc->prepare($total_WH_gene_model_nucleotides_query);
   $sth->execute();
   my $total_WH_gene_model_nucleotides = $sth->fetchrow_array();
   my $percent_WH_gene_model_nucleotide = sprintf "%.0f", ($total_WH_gene_model_nucleotides / $total_gene_model_nucleotides) * 100;
   
   # Get total NH gene model nucleotides
   $sth = $ENS_DBA->dbc->prepare($total_NH_gene_model_nucleotides_query);
   $sth->execute();
   my $total_NH_gene_model_nucleotides = $sth->fetchrow_array();
   my $percent_NH_gene_model_nucleotide = sprintf "%.0f", ($total_NH_gene_model_nucleotides / $total_gene_model_nucleotides) * 100;

   # Get total corrupted_translation gene model nucleotides
   $sth = $ENS_DBA->dbc->prepare($total_corrupted_translation_gene_model_nucleotides_query);
   $sth->execute();
   my $total_corrupted_translation_gene_model_nucleotides = $sth->fetchrow_array();
   my $percent_corrupted_translation_gene_model_nucleotide = sprintf "%.0f", ($total_corrupted_translation_gene_model_nucleotides / $total_gene_model_nucleotides) * 100;

   return ($total_gene_model_nucleotides,
           $total_TE_gene_model_nucleotides, $percent_TE_gene_model_nucleotide,
           $total_WH_gene_model_nucleotides, $percent_WH_gene_model_nucleotide,
           $total_NH_gene_model_nucleotides, $percent_NH_gene_model_nucleotide,
           $total_corrupted_translation_gene_model_nucleotides, $percent_corrupted_translation_gene_model_nucleotide);
}

sub get_total_nucleotides {
   # Get total number of nucleotides of the current BAC clones
   my $total_nucleotides_sql = "select sum(sr.length) from seq_region sr, seq_region_attrib sra ".
                               "where sr.seq_region_id=sra.seq_region_id and ".
                               "sr.coord_system_id=(select coord_system_id from coord_system where name = 'clone') and ".
                               "sra.value='current'";
   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
   my $sth = $ENS_DBA->dbc->prepare($total_nucleotides_sql);
   $sth->execute();
   my $total_nucleotides = $sth->fetchrow_array();
   $sth->finish();
   return ($total_nucleotides);   
}

sub get_mdr_stats {
    my ($total_nucleotides) = shift || die( "Need total nucleotides" );
   
   # Queries to get total nucleotides for various MDRs
   my $mdr_nucleotides_sql = "(select ((sf.seq_region_end - sf.seq_region_start) + 1) AS NT from simple_feature sf ".
                             "where sf.seq_region_id in ".
                             "(select f.seq_region_id from seq_region a ".
                             "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
                             "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
                             "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
                             "where a.coord_system_id = (select coord_system_id from coord_system where name = 'clone') and ".
                             "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') ".
                             "and v.value = 'current')";
   my $total_mdr_025_nucleotides_sql = "select sum(NT) from ".$mdr_nucleotides_sql." and analysis_id=".
       "(select analysis_id from analysis where logic_name='mdr_0.25')) tmpTable";
   my $total_mdr_1_nucleotides_sql = "select sum(NT) from ".$mdr_nucleotides_sql." and analysis_id=".
       "(select analysis_id from analysis where logic_name='mdr_1')) tmpTable";
   my $total_mdr_2_nucleotides_sql = "select sum(NT) from ".$mdr_nucleotides_sql." and analysis_id=".
       "(select analysis_id from analysis where logic_name='mdr_2')) tmpTable";
   my $total_mdr_3_nucleotides_sql = "select sum(NT) from ".$mdr_nucleotides_sql." and analysis_id=".
       "(select analysis_id from analysis where logic_name='mdr_3')) tmpTable";
   
   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
    my $sth;

   # Get total nucleotides for mdr_025
    $sth = $ENS_DBA->dbc->prepare($total_mdr_025_nucleotides_sql);
    $sth->execute();
    my $total_mdr_025_nucleotides = $sth->fetchrow_array();
    my $percent_mdr_025_nucleotides = sprintf "%.0f", ($total_mdr_025_nucleotides / $total_nucleotides) * 100;

   # Get total nucleotides for mdr_1
    $sth = $ENS_DBA->dbc->prepare($total_mdr_1_nucleotides_sql);
    $sth->execute();
    my $total_mdr_1_nucleotides = $sth->fetchrow_array();
    my $percent__mdr_1_nucleotides = sprintf "%.0f", ($total_mdr_1_nucleotides / $total_nucleotides) * 100;

   # Get total nucleotides for mdr_2
    $sth = $ENS_DBA->dbc->prepare($total_mdr_2_nucleotides_sql);
    $sth->execute();
    my $total_mdr_2_nucleotides = $sth->fetchrow_array();
    my $percent_mdr_2_nucleotides = sprintf "%.0f", ($total_mdr_2_nucleotides / $total_nucleotides) * 100;
    
   # Get total nucleotides for mdr_3
    $sth = $ENS_DBA->dbc->prepare($total_mdr_3_nucleotides_sql);
    $sth->execute();
    my $total_mdr_3_nucleotides = $sth->fetchrow_array(); 
    my $percent_mdr_3_nucleotides = sprintf "%.0f", ($total_mdr_3_nucleotides / $total_nucleotides) * 100;

    $sth->finish();
   return ($total_mdr_025_nucleotides, $percent_mdr_025_nucleotides,
           $total_mdr_1_nucleotides, $percent__mdr_1_nucleotides,
           $total_mdr_2_nucleotides, $percent_mdr_2_nucleotides,
           $total_mdr_3_nucleotides, $percent_mdr_3_nucleotides);
}

sub get_interpro_stats {
    my ($total_FGENESH_gene_models) = shift || die( "Need total Fgenesh gene models" );
   
   # Query to get number of gene models having a hit to each database in interpro
   my $sql = "select p.analysis_id, a.logic_name, count(distinct p.translation_id) AS COUNT from protein_feature p, analysis a ".
             "where p.translation_id in ".
             "(select gene_id from gene where seq_region_id in ".
             "(select f.seq_region_id from seq_region a ".
             "left join assembly asm on asm.asm_seq_region_id = a.seq_region_id ".
             "left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id ".
             "left join seq_region_attrib v on v.seq_region_id = a.seq_region_id ".
             "where a.coord_system_id = (select coord_system_id from coord_system where name = 'clone') and ".
             "f.coord_system_id = (select coord_system_id from coord_system where name = 'contig') and ".
             "v.value = 'current')) ".
             "and p.analysis_id=a.analysis_id group by p.analysis_id order by COUNT desc";

   # Get DBConnection. Then use DBConnection's prepare method to prepare sql statement and get DBI statement handle
    my $sth;
    $sth = $ENS_DBA->dbc->prepare($sql);
    $sth->execute();

    my %interpro = ();
    my @row;
    my $logic_name = '';
    my $hits = 0;
    my $percent_hit = 0;
    while ( @row = $sth->fetchrow_array() ) {
	$logic_name = $row[1];
	$hits = $row[2];
	$percent_hit = sprintf "%.2f", ($hits / $total_FGENESH_gene_models) * 100;
	$interpro{$logic_name}->{ 'HITS' } = $hits;
	$interpro{$logic_name}->{ 'PERCENT_HITS' } = $percent_hit;
    }#end while

	$sth->finish();
    return (\%interpro);
}

sub create_QC_report_table {
   my ($total_clones, $clone_min_size, $clone_max_size, $clone_ave_size, $clone_median_size, $clone_stddev_size, 
       $total_contigs, $contig_min_size, $contig_max_size, $contig_ave_size, $contig_median_size, $contig_stddev_size, 
       $contigsperclone_min_size, $contigsperclone_max_size, $contigsperclone_ave_size, $contigsperclone_median_size, $contigsperclone_stddev_size,
       $total_FGENESH_gene_models, $genemodel_min_size, $genemodel_max_size, $genemodel_ave_size, $genemodel_median_size, $genemodel_stddev_size, 
       $total_FGENESH_TE_gene_models, $percent_TE_gene_models,
       $min_TE_gene_model_size, $max_TE_gene_model_size, $ave_TE_gene_model_size, $median_TE_gene_model_size, $std_TE_gene_model_size,
       $total_FGENESH_WH_gene_models, $percent_WH_gene_models,
       $min_WH_gene_model_size, $max_WH_gene_model_size, $ave_WH_gene_model_size, $median_WH_gene_model_size, $std_WH_gene_model_size,
       $total_FGENESH_NH_gene_models, $percent_NH_gene_models,
       $min_NH_gene_model_size, $max_NH_gene_model_size, $ave_NH_gene_model_size, $median_NH_gene_model_size, $std_NH_gene_model_size,
#       $total_FGENESH_corrupted_translation_gene_models, $percent_corrupted_translation_gene_models,
#       $min_corrupted_translation_gene_model_size, $max_corrupted_translation_gene_model_size, $ave_corrupted_translation_gene_model_size, 
#       $median_corrupted_translation_gene_model_size, $std_corrupted_translation_gene_model_size,
       $total_gene_model_nucleotides,
       $total_TE_gene_model_nucleotides, $percent_TE_gene_model_nucleotide,
       $total_WH_gene_model_nucleotides, $percent_WH_gene_model_nucleotide,
       $total_NH_gene_model_nucleotides, $percent_NH_gene_model_nucleotide,
       $total_corrupted_translation_gene_model_nucleotides, $percent_corrupted_translation_gene_model_nucleotide,
       $total_nucleotides,
       $total_mdr_025_nucleotides, $percent_mdr_025_nucleotides,
       $total_mdr_1_nucleotides, $percent__mdr_1_nucleotides,
       $total_mdr_2_nucleotides, $percent_mdr_2_nucleotides,
       $total_mdr_3_nucleotides, $percent_mdr_3_nucleotides,
       $interpro) = @_; 
   
   my $db_host = $ENS_DBA->dbc->host;
   my $db_name = $ENS_DBA->dbc->dbname;
   my $db_port = $ENS_DBA->dbc->port;

   if (!$file) { $file = 'QC_report'; }
   my $report_text_filename = $file."_".$date.".txt"; 
   my $report_html_filename = $file."_".$date.".html";

   open (TXTOUT, ">$report_text_filename") || die("Cannot Open Text File");
   open (HTMLOUT, ">$report_html_filename") || die("Cannot Open HTML File");

   # Write to the plain text file
   print(TXTOUT "QC Report Filename: ".$report_text_filename."\n\n");
   print(TXTOUT $db_host.": ".$db_name."(".$date.")\n");
   print(TXTOUT "=========================================\n");

   print(TXTOUT "BAC Clones: ".$total_clones."\n");
   print(TXTOUT "BAC Clone Minimum Size (bp): ".$clone_min_size."\n");
   print(TXTOUT "BAC Clone Maximum Size (bp): ".$clone_max_size."\n");
   print(TXTOUT "BAC Clone Average Size (bp): ".$clone_ave_size."\n");
   print(TXTOUT "BAC Clone Median Size (bp): ".$clone_median_size."\n");
   print(TXTOUT "BAC BAC Clone Standard deviation (bp): ".$clone_stddev_size."\n\n");

   print(TXTOUT "Contigs: ".$total_contigs."\n");
   print(TXTOUT "Contig Minimum Size (bp): ".$contig_min_size."\n");
   print(TXTOUT "Contig Maximum Size (bp): ".$contig_max_size."\n");
   print(TXTOUT "Contig Average Size (bp): ".$contig_ave_size."\n");
   print(TXTOUT "Contig Median Size (bp): ".$contig_median_size."\n");
   print(TXTOUT "Contig Standard deviation (bp): ".$contig_stddev_size."\n");
   print(TXTOUT "Minimum Contigs Per Clone: ".$contigsperclone_min_size."\n");
   print(TXTOUT "Maximum Contigs Per Clone: ".$contigsperclone_max_size."\n");
   print(TXTOUT "Average Contigs Per Clone: ".$contigsperclone_ave_size."\n");
   print(TXTOUT "Median Contigs Per Clone: ".$contigsperclone_median_size."\n");
   print(TXTOUT "Standard deviation Contigs Per Clone: ".$contigsperclone_stddev_size."\n\n");

   print(TXTOUT "Fgenesh Gene Models: ".$total_FGENESH_gene_models."\n");
   print(TXTOUT "Fgenesh Gene Model Minimum Size (bp): ".$genemodel_min_size."\n");
   print(TXTOUT "Fgenesh Gene Model Maximum Size (bp): ".$genemodel_max_size."\n");
   print(TXTOUT "Fgenesh Gene Model Average Size (bp): ".$genemodel_ave_size."\n");
   print(TXTOUT "Fgenesh Gene Model Median Size (bp): ".$genemodel_median_size."\n");
   print(TXTOUT "Fgenesh Gene Model Standard deviation Size (bp): ".$genemodel_stddev_size."\n\n");

   print(TXTOUT "TE in Gene Models: ".$total_FGENESH_TE_gene_models."\n");
   print(TXTOUT "% Gene Models in TE: ".$percent_TE_gene_models."\n");
   print(TXTOUT "Minimum TE Gene Model Size (bp): ".$min_TE_gene_model_size."\n");
   print(TXTOUT "Maximum TE Gene Model Size (bp): ".$max_TE_gene_model_size."\n");
   print(TXTOUT "Average TE Gene Model Size (bp): ".$ave_TE_gene_model_size."\n");
   print(TXTOUT "Median TE Gene Model Size (bp): ".$median_TE_gene_model_size."\n");
   print(TXTOUT "Standard Deviation TE Gene Model Size (bp): ".$std_TE_gene_model_size."\n\n");
   
   print(TXTOUT "WH in Gene Models: ".$total_FGENESH_WH_gene_models."\n");
   print(TXTOUT "% Gene Models in WH: ".$percent_WH_gene_models."\n");
   print(TXTOUT "Minimum WH Gene Model Size (bp): ".$min_WH_gene_model_size."\n");
   print(TXTOUT "Maximum WH Gene Model Size (bp): ".$max_WH_gene_model_size."\n");
   print(TXTOUT "Average WH Gene Model Size (bp): ".$ave_WH_gene_model_size."\n");
   print(TXTOUT "Median WH Gene Model Size (bp): ".$median_WH_gene_model_size."\n");
   print(TXTOUT "Standard Deviation WH Gene Model Size (bp): ".$std_WH_gene_model_size."\n\n");

   print(TXTOUT "NH in Gene Models: ".$total_FGENESH_NH_gene_models."\n");
   print(TXTOUT "% Gene Models in NH: ".$percent_NH_gene_models."\n");
   print(TXTOUT "Minimum NH Gene Model Size (bp): ".$min_NH_gene_model_size."\n");
   print(TXTOUT "Maximum NH Gene Model Size (bp): ".$max_NH_gene_model_size."\n");
   print(TXTOUT "Average NH Gene Model Size (bp): ".$ave_NH_gene_model_size."\n");
   print(TXTOUT "Median NH Gene Model Size (bp): ".$median_NH_gene_model_size."\n");
   print(TXTOUT "Standard Deviation NH Gene Model Size (bp): ".$std_NH_gene_model_size."\n\n");

 #  print(TXTOUT "Corrupted_translation in Gene Models: ".$total_FGENESH_corrupted_translation_gene_models."\n");
 #  print(TXTOUT "% Gene Models in Corrupted_translation: ".$percent_corrupted_translation_gene_models."\n");
 #  print(TXTOUT "Minimum Corrupted_translation Gene Model Size (bp): ".$min_corrupted_translation_gene_model_size."\n");
 #  print(TXTOUT "Maximum Corrupted_translation Gene Model Size (bp): ".$max_corrupted_translation_gene_model_size."\n");
 #  print(TXTOUT "Average Corrupted_translation Gene Model Size (bp): ".$ave_corrupted_translation_gene_model_size."\n");
 #  print(TXTOUT "Median Corrupted_translation Gene Model Size (bp): ".$median_corrupted_translation_gene_model_size."\n");
 #  print(TXTOUT "Standard Deviation Corrupted_translation Gene Model Size (bp): ".$std_corrupted_translation_gene_model_size."\n\n");

   print (TXTOUT "Nucleotides (bp) in Gene Models: ".$total_gene_model_nucleotides."\n");
   print (TXTOUT "Nucleotides (bp) in TE Gene Models: ".$total_TE_gene_model_nucleotides."\n");
   print (TXTOUT "% Gene Model Nucleotides in TE Gene Models: ".$percent_TE_gene_model_nucleotide."\n");
   print (TXTOUT "Nucleotides (bp) in WH Gene Models: ".$total_WH_gene_model_nucleotides."\n");
   print (TXTOUT "% Gene Model Nucleotides in WH Gene Models: ".$percent_WH_gene_model_nucleotide."\n");
   print (TXTOUT "Nucleotides (bp) in NH Gene Models: ".$total_NH_gene_model_nucleotides."\n");
   print (TXTOUT "% Gene Model Nucleotides in NH Gene Models: ".$percent_NH_gene_model_nucleotide."\n");
   print (TXTOUT "Nucleotides (bp) in corrupted_translation Gene Models: ".$total_corrupted_translation_gene_model_nucleotides."\n");
   print (TXTOUT "% Gene Model Nucleotides in corrupted_translation Gene Models: ".$percent_corrupted_translation_gene_model_nucleotide."\n\n");
    
   print (TXTOUT "Total Nucleotides (bp): ".$total_nucleotides."\n\n");

   print (TXTOUT "Nucleotides (bp) in MDR >= 2 Copies: ".$total_mdr_025_nucleotides."\n");
   print (TXTOUT "% Nucleotides in MDR >= 2 Copies: ".$percent_mdr_025_nucleotides."\n");
   print (TXTOUT "Nucleotides (bp) in MDR >= 10 Copies: ".$total_mdr_1_nucleotides."\n");
   print (TXTOUT "% Nucleotides in MDR >= 10 Copies: ".$percent__mdr_1_nucleotides."\n");
   print (TXTOUT "Nucleotides (bp) in MDR >= 100 Copies: ".$total_mdr_2_nucleotides."\n");
   print (TXTOUT "% Nucleotides in MDR >= 100 Copies: ".$percent_mdr_2_nucleotides."\n\n");
   print (TXTOUT "Nucleotides (bp) in MDR >= 1000 Copies: ".$total_mdr_3_nucleotides."\n");
   print (TXTOUT "% Nucleotides in MDR >= 1000 Copies: ".$percent_mdr_3_nucleotides."\n\n");

   print (TXTOUT "InterPro:\n");
   foreach my $key (sort (keys(%$interpro))) { #dereference to hash reference
       print (TXTOUT $key." Hits in Gene Models: ".$interpro->{$key}->{ 'HITS' }."\n");  
       print (TXTOUT "% Gene Models with ".$key." Hits: ".$interpro->{$key}->{ 'PERCENT_HITS'  }."\n");
   }

   close(TXTOUT);

   # Write to html file
   print(HTMLOUT "<table border=1 cellspacing=5 cellpadding=5>"."\n");
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<th colspan=2>".$db_host.": ".$db_name."(".$date.")</th>"."\n");
   print(HTMLOUT "</tr>"."\n");
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<th align=left>Status</th>"."\n");
   print(HTMLOUT "<th align=left>Count</th>"."\n");
   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>BAC Clones<br>"."\n");
   print(HTMLOUT "Minimum Size (bp)<br>"."\n");
   print(HTMLOUT "Maximum Size (bp)<br>"."\n");
   print(HTMLOUT "Average Size (bp)<br>"."\n");
   print(HTMLOUT "Median Size (bp)<br>"."\n");
   print(HTMLOUT "Standard deviation (bp)</td>"."\n");
   print(HTMLOUT "<td>".$total_clones."<br>"."\n");
   print(HTMLOUT $clone_min_size."<br>"."\n");
   print(HTMLOUT $clone_max_size."<br>"."\n");
   print(HTMLOUT $clone_ave_size."<br>"."\n");
   print(HTMLOUT $clone_median_size."<br>"."\n");
   print(HTMLOUT $clone_stddev_size."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");
  
   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>Contigs<br>"."\n");
   print(HTMLOUT "Minimum Size (bp)<br>"."\n");
   print(HTMLOUT "Maximum Size (bp)<br>"."\n");
   print(HTMLOUT "Average Size (bp)<br>"."\n");
   print(HTMLOUT "Median Size (bp)<br>"."\n");
   print(HTMLOUT "Standard deviation (bp)<br><br>"."\n");
   print(HTMLOUT "Minimum Contigs Per Clone<br>"."\n");
   print(HTMLOUT "Maximum Contigs Per Clone<br>"."\n");
   print(HTMLOUT "Average Contigs Per Clone<br>"."\n");
   print(HTMLOUT "Median Contigs Per Clone<br>"."\n");
   print(HTMLOUT "Standard deviation Contigs Per Clone</td>"."\n");
   print(HTMLOUT "<td>".$total_contigs."<br>"."\n");
   print(HTMLOUT $contig_min_size."<br>"."\n");
   print(HTMLOUT $contig_max_size."<br>"."\n");
   print(HTMLOUT $contig_ave_size."<br>"."\n");
   print(HTMLOUT $contig_median_size."<br>"."\n");
   print(HTMLOUT $contig_stddev_size. "<br><br>"."\n");
   print(HTMLOUT $contigsperclone_min_size."<br>"."\n");
   print(HTMLOUT $contigsperclone_max_size."<br>"."\n");
   print(HTMLOUT $contigsperclone_ave_size."<br>"."\n");
   print(HTMLOUT $contigsperclone_median_size."<br>"."\n");
   print(HTMLOUT $contigsperclone_stddev_size."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>Fgenesh Gene Models<br>"."\n");
   print(HTMLOUT "Minimum Size (bp)<br>"."\n");
   print(HTMLOUT "Maximum Size (bp)<br>"."\n");
   print(HTMLOUT "Average Size (bp)<br>"."\n");
   print(HTMLOUT "Median Size (bp)<br>"."\n");
   print(HTMLOUT "Standard deviation Size (bp)</td>"."\n");
   print(HTMLOUT "<td>".$total_FGENESH_gene_models."<br>"."\n");
   print(HTMLOUT $genemodel_min_size."<br>"."\n");
   print(HTMLOUT $genemodel_max_size."<br>"."\n");
   print(HTMLOUT $genemodel_ave_size."<br>"."\n");
   print(HTMLOUT $genemodel_median_size."<br>"."\n");
   print(HTMLOUT $genemodel_stddev_size."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>TE in Gene Models<br>"."\n");
   print(HTMLOUT "% Gene Models in TE<br>"."\n");
   print(HTMLOUT "Minimum TE Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Maximum TE Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Average TE Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Median TE Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Standard Deviation TE Gene Model Size (bp)</td>"."\n");
   print(HTMLOUT "<td>".$total_FGENESH_TE_gene_models."<br>"."\n");
   print(HTMLOUT $percent_TE_gene_models."<br>"."\n");
   print(HTMLOUT $min_TE_gene_model_size."<br>"."\n");
   print(HTMLOUT $max_TE_gene_model_size."<br>"."\n");
   print(HTMLOUT $ave_TE_gene_model_size."<br>"."\n");
   print(HTMLOUT $median_TE_gene_model_size."<br>"."\n");
   print(HTMLOUT $std_TE_gene_model_size."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>WH in Gene Models<br>"."\n");
   print(HTMLOUT "% Gene Models in WH<br>"."\n");
   print(HTMLOUT "Minimum WH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Maximum WH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Average WH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Median WH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Standard Deviation WH Gene Model Size (bp)</td>"."\n");
   print(HTMLOUT "<td>".$total_FGENESH_WH_gene_models."<br>"."\n");
   print(HTMLOUT $percent_WH_gene_models."<br>"."\n");
   print(HTMLOUT $min_WH_gene_model_size."<br>"."\n");
   print(HTMLOUT $max_WH_gene_model_size."<br>"."\n"); 
   print(HTMLOUT $ave_WH_gene_model_size."<br>"."\n");
   print(HTMLOUT $median_WH_gene_model_size."<br>"."\n");
   print(HTMLOUT $std_WH_gene_model_size."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");   

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>NH in Gene Models<br>"."\n");
   print(HTMLOUT "% Gene Models in NH<br>"."\n");
   print(HTMLOUT "Minimum NH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Maximum NH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Average NH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Median NH Gene Model Size (bp)<br>"."\n");
   print(HTMLOUT "Standard Deviation NH Gene Model Size (bp)</td>"."\n");
   print(HTMLOUT "<td>".$total_FGENESH_NH_gene_models."<br>"."\n");
   print(HTMLOUT $percent_NH_gene_models."<br>"."\n");
   print(HTMLOUT $min_NH_gene_model_size."<br>"."\n");
   print(HTMLOUT $max_NH_gene_model_size."<br>"."\n");
   print(HTMLOUT $ave_NH_gene_model_size."<br>"."\n");
   print(HTMLOUT $median_NH_gene_model_size."<br>"."\n");
   print(HTMLOUT $std_NH_gene_model_size."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");

#   print(HTMLOUT "<tr>"."\n");
#   print(HTMLOUT "<td>Corrupted_translation in Gene Models<br>"."\n");
#   print(HTMLOUT "% Gene Models in Corrupted_translation<br>"."\n");
#   print(HTMLOUT "Minimum Corrupted_translation Gene Model Size (bp)<br>"."\n");
#   print(HTMLOUT "Maximum Corrupted_translation Gene Model Size (bp)<br>"."\n");
#   print(HTMLOUT "Average Corrupted_translation Gene Model Size (bp)<br>"."\n");
#   print(HTMLOUT "Median Corrupted_translation Gene Model Size (bp)<br>"."\n");
#   print(HTMLOUT "Standard Deviation Corrupted_translation Gene Model Size (bp)</td>"."\n");
#   print(HTMLOUT "<td>".$total_FGENESH_corrupted_translation_gene_models."<br>"."\n");
#   print(HTMLOUT $percent_corrupted_translation_gene_models."<br>"."\n");
#   print(HTMLOUT $min_corrupted_translation_gene_model_size."<br>"."\n");
#   print(HTMLOUT $max_corrupted_translation_gene_model_size."<br>"."\n");
#   print(HTMLOUT $ave_corrupted_translation_gene_model_size."<br>"."\n");
#   print(HTMLOUT $median_corrupted_translation_gene_model_size."<br>"."\n");
#   print(HTMLOUT $std_corrupted_translation_gene_model_size."</td>"."\n");
#   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>Nucleotides (bp) in Gene Models<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in TE Gene Models<br>"."\n");
   print(HTMLOUT "% Gene Model Nucleotides in TE Gene Models<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in WH Gene Models<br>"."\n");
   print(HTMLOUT "% Gene Model Nucleotides in WH Gene Models<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in NH Gene Models<br>"."\n");
   print(HTMLOUT "% Gene Model Nucleotides in NH Gene Models<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in corrupted_translation Gene Models<br>"."\n");
   print(HTMLOUT "% Gene Model Nucleotides in corrupted_translation Gene Models</td>"."\n");
   print(HTMLOUT "<td>".$total_gene_model_nucleotides."<br>"."\n");
   print(HTMLOUT $total_TE_gene_model_nucleotides."<br>"."\n");
   print(HTMLOUT $percent_TE_gene_model_nucleotide."<br>"."\n");
   print(HTMLOUT $total_WH_gene_model_nucleotides."<br>"."\n");
   print(HTMLOUT $percent_WH_gene_model_nucleotide."<br>"."\n");
   print(HTMLOUT $total_NH_gene_model_nucleotides."<br>"."\n");
   print(HTMLOUT $percent_NH_gene_model_nucleotide."<br>"."\n");
   print(HTMLOUT $total_corrupted_translation_gene_model_nucleotides."<br>"."\n");
   print(HTMLOUT $percent_corrupted_translation_gene_model_nucleotide."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>Total Nucleotides (bp)</td>"."\n");
   print(HTMLOUT "<td>".$total_nucleotides."</td>"."\n");
   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>Mathematically-Defined Repeats<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in MDR >= 2 Copies<br>"."\n");
   print(HTMLOUT "% Nucleotides in MDR >= 2 Copies<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in MDR >= 10 Copies<br>"."\n");
   print(HTMLOUT "% Nucleotides in MDR >= 10 Copies<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in MDR >= 100 Copies<br>"."\n");
   print(HTMLOUT "% Nucleotides in MDR >= 100 Copies<br>"."\n");
   print(HTMLOUT "Nucleotides (bp) in MDR >= 1000 Copies<br>"."\n");
   print(HTMLOUT "% Nucleotides in MDR >= 1000 Copies</td>"."\n");
   print(HTMLOUT "<td><br>"."\n");
   print(HTMLOUT $total_mdr_025_nucleotides."<br>"."\n");
   print(HTMLOUT $percent_mdr_025_nucleotides."<br>"."\n");
   print(HTMLOUT $total_mdr_1_nucleotides."<br>"."\n");
   print(HTMLOUT $percent__mdr_1_nucleotides."<br>"."\n");
   print(HTMLOUT $total_mdr_2_nucleotides."<br>"."\n");
   print(HTMLOUT $percent_mdr_2_nucleotides."<br>"."\n");
   print(HTMLOUT $total_mdr_3_nucleotides."<br>"."\n");
   print(HTMLOUT $percent_mdr_3_nucleotides."</td>"."\n"); 

   print(HTMLOUT "<tr>"."\n");
   print(HTMLOUT "<td>"."\n");
   print(HTMLOUT "Interpro<br>"."\n");
   foreach my $key (sort (keys(%$interpro))) { #dereference to hash reference
       print(HTMLOUT $key." Hits in Gene Models<br>"."\n");
       print(HTMLOUT "% Gene Models with ".$key." Hits<br>"."\n");
   }
   print(HTMLOUT "</td>"."\n");
   print(HTMLOUT "<td>"."\n");
   print(HTMLOUT "<br>"."\n");
   foreach my $key (sort (keys(%$interpro))) { #dereference to hash reference
       print(HTMLOUT $interpro->{$key}->{ 'HITS' }."<br>"."\n");
       print(HTMLOUT $interpro->{$key}->{ 'PERCENT_HITS'  }."<br>"."\n");    
   }
   print(HTMLOUT "</td>"."\n");
   print(HTMLOUT "</tr>"."\n");

   print(HTMLOUT "</table>"."\n");
   close(HTMLOUT);
}

sub get_coord_system_id {
    my $coord_system_name = shift || die( "Need a coord_system_name" );
   # Get CoordSystemAdaptor
    my $csa = $ENS_DBA->get_CoordSystemAdaptor();
   # Use CoordSystemAdaptor's fetch_by_name method to retrieve a CoordSystem object by its name
    my $cs = $csa->fetch_by_name($coord_system_name);
   # Use the CoordSystem object to get the object's dbID variable value
    return my $coord_system_id = $cs->dbID();
}

sub get_analysis_id {
    my $analysis_logic_name = shift || die( "Need a logic_name" );
   # Get AnalysisAdaptor
    my $aa = $ENS_DBA->get_AnalysisAdaptor;
   #Use AnalysisAdaptor's fetch_by_logic_name method to retrieve an Analysis object by its logic name
    my $a  = $aa->fetch_by_logic_name($analysis_logic_name);
   # Use the Analysis object to get the object's dbID variable value
    return my $analysis_id = $a->dbID;
}

1;
