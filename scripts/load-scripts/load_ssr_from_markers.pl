#!/usr/local/bin/perl


BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/';
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/';
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" }
        qw ( 
             ensembl/modules);

#use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;


use DBI;

use Bio::EnsEMBL::DBLoader;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::Registry;

#use lib "/home/weix/scripts/markers/lib/gramene/lib/perl";
use Gramene::Marker::DB;
use List::MoreUtils qw( uniq );

$| = 1;
my $marker_source = 'markersdb';
my $registry_file;
my $logic_name;
my $ensembl_species=$ENV{ENSEMBL_SPECIES};
my $synonyms_file;
my $feature_file;
my $ENS_DBA;
my $mdb_analysis;
my $MDB;
my $conf_file;

{  #Argument Processing
  my $help=0;
  my $man=0;
  GetOptions
      (
       "help|?"          => \$help,
       "man"             => \$man,
       "species=s"       => \$ensembl_species,
       "registry_file=s" => \$registry_file,
       "logic_name=s"    => \$logic_name,
       # "synonyms_file=s"    => \$synonyms_file,
       # "feature_file=s"    => \$feature_file
       "mdb_analysis=i"  => \$mdb_analysis,
       "conf_file=s"     => \$conf_file,
      )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;


  $registry_file    ||= $ENV{GrameneEnsemblDir}.'/conf/ensembl.registry';
  $conf_file  && do{$ENV{GrameneConfPath} = $conf_file};

  # Load the ensembl file
  $ensembl_species || ( warn( "Need a --species\n" ) && pod2usage(1) );
  Bio::EnsEMBL::Registry->load_all( $registry_file )
            or die "load_all($registry_file) failed";

  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $ensembl_species, 'core' );
  $ENS_DBA 
      || ( warn( "No core DB for $ensembl_species set in $registry_file\n" ) &&
                pod2usage(1) );

  $MDB     = Gramene::Marker::DB->new;
print "created MDB\n";

  $MDB
    || ( warn( "No core Markers DB object created\n" ) &&
	 pod2usage(1) );
}

my $slice_adaptor    = $ENS_DBA->get_adaptor('Slice');
my $analysis_adaptor = $ENS_DBA->get_adaptor('Analysis');

my $db = $ENS_DBA->dbc();
my $analysis = &fetch_analysis($logic_name);
my $analysis_id=$analysis->dbID();

#retrieve mapmappingspings for this analysis
my $mapping_objs = $MDB->search_Mapping(analysis_id => $mdb_analysis) 
  || die "no records found for analysis_id $mdb_analysis in mapping table\n";

my %marker_ids;
my ($markersdb_mapping_cnt, $ensembl_marker_cnt, $ensembl_marker_feature_cnt);
while( my $mapping_obj = $mapping_objs->next){
  #print "mapping_obj id = " .  $mapping_obj->mapping_id . "\n";
    
  ++$markersdb_mapping_cnt;
  
  my $feat_chr_start = $mapping_obj->start;
  my $feat_chr_end   = $mapping_obj->end;
  my $map_name       = $MDB->retrieve_Map( map_id => $mapping_obj->map_id )->map_name;
  my ($chr, $cs);
  if($map_name =~ m/Chr[.-_]?\s*(\S+)/xmsi){
    $cs = 'Chromosome';
    $chr = $1;
  }elsif( $map_name =~ m/(\S+)\s*:\s*(\S+)/xmsi ){
    $cs = $1;
    $chr = $2;
  }

  my $marker_id  = $mapping_obj->marker_id;  
  my $marker_obj = $MDB->retrieve_Marker(marker_id => $marker_id);
  my $marker_display_syn_id = $marker_obj->display_synonym_id;
  my $marker_display_name   = $MDB->retrieve_MarkerSynonym(marker_synonym_id => $marker_display_syn_id)->marker_name;
  my $marker_syns_obj          = $MDB->search_MarkerSynonym(marker_id => $marker_id);
  my @syn_names = ($marker_display_name); #put display_name the 1st in the list,
                                          # uniq preserve the list order

  while( my $marker_syn = $marker_syns_obj->next ){
    push @syn_names, $marker_syn->marker_name;
  }

  my @uniq_syn_names = uniq @syn_names;

  unless( $marker_ids{$marker_id} ){
    my $ensembl_marker_id = add_ensembl_marker(\@uniq_syn_names);
    $marker_ids{$marker_id} = $ensembl_marker_id;
  }

  my $sth=$db->prepare("select seq_region_id from seq_region sr join coord_system cs using (coord_system_id) where cs.name = ? and sr.name = ?" );
  warn("get seq_region_id for $cs, $chr\n");
  $sth->execute( $cs, $chr ) or die "error fetching seqregion_id:".$db->db_handle()->errstr;
  my ($seq_region_id)=$sth->fetchrow_array or die "error fetching seq_region_id:".$db->db_handle->errstr;

  my $marker_feature_id = add_ensembl_marker_feature($marker_ids{$marker_id}, $seq_region_id, 
			       $feat_chr_start, $feat_chr_end, $analysis_id);
  if($marker_feature_id){
    ++$ensembl_marker_feature_cnt;
  }else{
    print "Cannot add_ensembl_marker_feature for: $marker_ids{$marker_id}, $seq_region_id, $feat_chr_start, $feat_chr_end, $analysis_id\n";
  }
}

$ensembl_marker_cnt = scalar keys %marker_ids;
print "
SSR mappings from markers db: $markersdb_mapping_cnt, 
Uniq SSRs from markers db   : $ensembl_marker_cnt, 
SSRs in ensembl marker feature :$ensembl_marker_feature_cnt
";

sub add_ensembl_marker{

  my $syn_names_ref = shift;

  #2 SSRs can share the same name
  #my $ensembl_marker_id = get_marker_id($syn_names_ref->[0]);
  #return $ensembl_marker_id if $ensembl_marker_id;

  $ensembl_marker_id = insert_marker_id();
  insert_marker_synonym($ensembl_marker_id, $syn_names_ref);
  return $ensembl_marker_id;

}

sub insert_marker_synonym {
  my($marker_id, $syn_names_ref) = @_;
  my $sth;
  
  $sth=$db->prepare("select max(marker_synonym_id) from marker_synonym") ;
  $sth->execute or die "max marker_synonym_id:".$db->db_handle->errstr;
  my ($marker_synonym_id)=$sth->fetchrow_array or die "max marker_synonym_id:".$db->db_handle->errstr;
  $marker_synonym_id ||=0;
  $marker_synonym_id++;
  
  my @syn_names = @{$syn_names_ref};
  
  my $display_name = shift @syn_names;
  $sth=$db->prepare("insert into marker_synonym values (?,?,?,?)");
  $sth->execute($marker_synonym_id,$marker_id,$marker_source,$display_name) 
    or die "cannot insert marker_synonym:($marker_synonym_id,$marker_id,$marker_source,$display_name)".$db->db_handle->errstr;
  $sth=$db->prepare("update marker set display_marker_synonym_id=$marker_synonym_id where marker_id=?");
  $sth->execute($marker_id) or 
    die "cannot update display marker synonym to $marker_synonym_id in marker table".$db->db_handle->errstr;
  
  for my $marker_name (@syn_names) {
    $marker_synonym_id++;
    $sth=$db->prepare("insert into marker_synonym values (?,?,?,?)");
    $sth->execute($marker_synonym_id,$marker_id,$marker_source,$marker_name) or 
      die "cannot insert marker_synonym:$marker_name".$db->db_handle->errstr;
  } 

}

sub get_marker_id {
   my($marker_name)=@_;
  my $sql_statement = "select marker_id
                       from marker_synonym
                       where marker_synonym.name = ? and source = ?";
  #print "$sql_statement\n";
  my $sth_contigs=$db->prepare($sql_statement) or die "marker query:$DBI::errstr";
  $sth_contigs->execute($marker_name, $marker_source);
  my($marker_id)=$sth_contigs->fetchrow_array or return undef;
  return $marker_id;
}

sub get_marker_feature_id {
  my($marker_id, $seq_region_id, $seq_region_start, $seq_region_end, $analysis_id)=@_;
  my $sql_statement = "select marker_feature_id
                       from marker_feature
                       where marker_id = ? and 
                             seq_region_id = ? and
                             seq_region_start = ? and
                             seq_region_end = ? and
                             analysis_id    = ?";
#marker_feature_id | marker_id | seq_region_id | seq_region_start | seq_region_end | analysis_id | map_weight |
  #print "$sql_statement\n";
  my $sth_contigs=$db->prepare($sql_statement) or die "marker query:$DBI::errstr";
  $sth_contigs->execute($marker_id, $seq_region_id, $seq_region_start, $seq_region_end, $analysis_id);
  my $marker_feature_id=$sth_contigs->fetchrow_array or return undef;
  return $marker_feature_id;
}



sub add_ensembl_marker_feature {
  my($marker_id, $contig_id, $ctgstart, $ctgstop, $analysis_id)=@_;
  
  my $marker_feature_id;
  if( $marker_feature_id = get_marker_feature_id(
						     $marker_id, 
						     $contig_id, 
						     $ctgstart, 
						     $ctgstop, 
						     $analysis_id
						    )
    ){
    print "Already in marker feature: $marker_id, $contig_id, $ctgstart, $ctgstop, $analysis_id\n";
    return $marker_feature_id;
  }

  # now insert marker_feature
  my $sth=$db->prepare("select max(marker_feature_id) from marker_feature");
  $sth->execute or die "max marker_feature_id:".$db->db_handle->errstr;

  $marker_feature_id = $sth->fetchrow_array ;
  $marker_feature_id ||=0;
  $marker_feature_id++;
  
  $sth = $db->prepare("insert into marker_feature
           (marker_feature_id,marker_id,seq_region_id,seq_region_start,
                                      seq_region_end,analysis_id,map_weight)
          values (?, ?, ?,?,?,?,?)"
      ) and $sth->execute(marker_feature_id,
			  $marker_id,
			  $contig_id,
			  $ctgstart,
			  $ctgstop,
			  $analysis_id,
			  1
			 )
      or die "cannot insert marker_feature".$db->db_handle->errstr;
   print "Marker_feature inserted : $marker_feature_id \n";

  return $marker_feature_id;
}

sub insert_marker_id {
  my $sth=$db->prepare("select max(marker_id) from marker") ;
      $sth->execute or die "max marker_id:".$db->db_handle->errstr;
  my ($marker_id)=$sth->fetchrow_array ;
    $marker_id ||=0;
  $marker_id++;
  $sth=$db->prepare("insert into marker 
           (marker_id,display_marker_synonym_id,priority,type)
        values ($marker_id,'0','100','microsatellite')") 
        and $sth->execute
      or die "cannot insert marker_id:".$db->db_handle->errstr;
  return $marker_id;
}




#======================================================================
# Returns an Analysis object; either fetched from the DB if one exists,
# or created fresh, in which case it is stored to the DB.
sub fetch_analysis{
  my $logic_name = shift || die("Need a logic_name" );
  my $db_file    = shift || '';
  my $adaptor = $ENS_DBA->get_adaptor('Analysis');

  my %args = ( -logic_name=>$logic_name,
               $db_file ? (-db_file=>$db_file) : () );

  # Hard-coded nastyness to make analyses correspond to Ensembl
  if( $logic_name eq 'MyAnalysis' ){
    $args{-logic_name} = 'MyLogic';
    $args{-db}         = 'MyDB';
    $args{-program}    = 'MyProgram';
  }

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

#======================================================================

=head1 SYNOPSIS
  
  "help|?"          => \$help,
  "man"             => \$man,
  "species=s"       => \$ensembl_species,
  "registry_file=s" => \$registry_file,
  "logic_name=s"    => \$logic_name,
  "mdb_analysis=i"  => \$mdb_analysis,
  "conf_file=s"     => \$conf_file,

=cut
