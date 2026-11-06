#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 VCF info format

 INFO=<ID=RS_VALIDATED,Number=0,Type=Flag,Description="RS validated flag, present when the RS was validated
 validation status">
 INFO=<ID=SID,Number=.,Type=String,Description="Identifiers of studies that report a variant">
 INFO=<ID=SS_VALIDATED,Number=1,Type=Integer,Description="Number of submitted variants clustered in an RS t
cated by the dbSNP validation status">

=cut

## collect variant synonyms & load to db
##  - initially from PharmGKB database

use lib "/usr/local/ensembl-108/ensembl-variation/modules/";
use lib "/usr/local/ensembl-108/ensembl/modules/";
use lib "/usr/local/ensembl-108/ensembl-compara/modules/";

use strict;
use warnings;
use HTTP::Tiny;
use Getopt::Long;
use DBI qw(:sql_types);
use Time::Piece;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Source;

our $DEBUG = 0;

my ($data_file, $registry_file, $species, $source_name, $clean, $source_version, $source_url, $source_description, $host,
   $port, $user, $pass, $db_name, $no_allele);

GetOptions ("data_file=s"          => \$data_file,
            "species=s"            => \$species,
            "source_name=s"        => \$source_name,
            "source_url=s"         => \$source_url,
            "source_version=s"     => \$source_version,
            "source_description=s" => \$source_description,
            "registry=s"           => \$registry_file,            
            "host=s"               => \$host,
            "port=s"               => \$port,
            "user=s"               => \$user,
            "pass=s"               => \$pass,
            "db_name=s"            => \$db_name,
            "clean"                => \$clean,
	    "no_allele"		   => \$no_allele 
            );

usage() unless defined $registry_file && defined $source_name;

usage() if($source_name eq 'EVA' && !$data_file);

$species ||= 'homo_sapiens';


my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($registry_file);
my $dba = $reg->get_DBAdaptor($species, 'variation');
my $core_dba = $reg->get_DBAdaptor($species, 'core');
my $slice_adaptor = $core_dba->get_SliceAdaptor;

## include failed variants to avoid missing any links
$dba->include_failed_variations(1);


## collect synonyms by source
if($source_name eq 'PharmGKB') { 
  my $source = get_source($species, $dba, $source_name, $source_version, $source_url, $source_description, 'germline');

  my $synonyms;

  ## collect data files if not already available
  $data_file = download_file('PharmGKB', 'rsid.zip') unless $data_file;

  ## extract synonyms from file
  $synonyms  = extract_PharmGKB($data_file);
  
  ## add synonyms to the database
  import_synonyms($synonyms, $source, $dba, $species); 
}

elsif($source_name eq 'dbSNP') {
  ### import synonyms from production database ### Access directly database because this database doesn't have adaptors 
  my $prod_dbc = Bio::EnsEMBL::DBSQL::DBConnection->new
    ( '-host'    => $host,
      '-port'    => $port, 
      '-user'    => $user,
      '-pass'    => $pass,
      '-dbname'  => $db_name, 
      '-driver'  =>  'mysql' 
    ); 

  my $sql_query = "SELECT rs_name,synonym_name,source_id FROM human_synonym"; 
  my $synonyms_output = query_synonyms_data($prod_dbc, $sql_query); 

  my %data; # synonym_name => {rs_name} & {source_id} 

  # get data from human_synonym table (in production db) 
  # some synonyms names are duplicated - from different sources 
  foreach my $result (@{$synonyms_output}){  
    $data{$result->[0].' '.$result->[1]}{rs_name} = $result->[0]; 
    $data{$result->[0].' '.$result->[1]}{source_id} = $result->[2]; 
  }

  # get all available sources from production database  
  my $sources = get_all_sources($dba, $prod_dbc);

  import_synonyms_all_sources($sources, $dba, $species, \%data);  
} 

elsif($source_name eq 'UniProt') {
  # download UniProt file if it's not provided  
  $data_file = download_file('UniProt', 'humsavar') unless $data_file;
 
  my $result = parse_uniprot($data_file); 

  # Source UniProt - variation_set table 
  my %source_uniprot = (
  "UniProt"     => {
    description => "Variants with annotations provided by UniProt",
    url         => "http://www.uniprot.org/",
    set         => "Uniprot variants",
    type        => "Variation",
    status      => 'mixed',
  }, 
  ); 

  # get source info 
  $source_description = $source_uniprot{$source_name}->{description}; 
  $source_url = $source_uniprot{$source_name}->{url}; 
  my $set = $source_uniprot{$source_name}->{set}; 
  my $object_type = $source_uniprot{$source_name}->{type}; 
  my $source_status = $source_uniprot{$source_name}->{status}; 

  # get source and returns source_id 
  my $source = get_source($species, $dba, $source_name, $source_version, $source_url, $source_description, $source_status); 
  my $source_id = $source->dbID();
  my $source_name = $source->name(); 

  my @ids; 
  my %synonym;
  
  if (exists($result->{'synonyms'})) {
    %synonym = %{$result->{'synonyms'}};
    # To get all the ids of the source (Uniprot)
    @ids = keys(%synonym);
    warn "There are no UniProt synonyms to import: check if data file has the correct format.\n" if !@ids;
  }

  # Add the synonyms
  my $variation_set_id = add_set_uniprot($dba,$set,$source_description);
  add_synonyms_uniprot(\%synonym, $source_name, $variation_set_id, $dba, $species); 
}

elsif($source_name eq 'EVA') {
  my $result = parse_eva_rsid_file($data_file);

  # Source EVA
  my %source_eva = (
  "EVA" => {
    description => "Short variant data imported from EVA",
    url         => "https://www.ebi.ac.uk/eva/",
    status      => "germline",
  }, 
  );

  # get source info 
  $source_description ||= $source_eva{$source_name}->{description};
  $source_url ||= $source_eva{$source_name}->{url};
  my $source_status = $source_eva{$source_name}->{status};

  my $source = get_source($species, $dba, $source_name, $source_version, $source_url, $source_description, $source_status);

  import_rsid($result, $source, $dba, $species);
}

elsif($source_name eq 'pig_chip') {
  my $int_dba = $reg->get_DBAdaptor('multi', 'intvar') || warn "Error getting intvar db adaptor\n";;

  # Source Pig SNP Consortium
  my %source_data = (
    "name"           => "Pig SNP Consortium",
    "version"        => undef,
    "description"    => "PorcineSNP60 BeadChip",
    "url"            => "https://www.animalgenome.org/repository/pig/",
    "somatic_status" => "germline"
  );

  my $source = get_source($species, $dba, $source_data{name}, $source_data{version}, $source_data{url}, $source_data{description}, $source_data{somatic_status});

  add_synonym_pig_chip($species, $dba, $int_dba, $source->dbID());
}

else{
    die "data source : $source_name not supported\n";
} 

=head2 download_file

collect current export from PharmGKB and UniProt site

=cut

sub download_file {
  my $source = shift; 
  my $data_file = shift;
  my $output_file; 

  my $http = HTTP::Tiny->new();
  my $url; 

  if($source eq "PharmGKB") {
    $url = "https://api.pharmgkb.org/v1/download/file/data/" . $data_file; 
  }
  elsif($source eq "UniProt") {
    $url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/" . $data_file;
  }
  else {
    die "Can't download file for data source: $data_file\n"; 
  }
  
  my $response = $http->get($url); 
  die "Failed to pick up file: $data_file\n" unless $response->{success};

  open (my $out, ">" . $data_file) or die ("Failed to write data locally : $!\n");
  print $out  $response->{content};

  if($source eq "PharmGKB") {
    eval{
    `unzip rsid.zip`
    };
    die "Failed to unzip data :$@\n" unless $@ eq '';
    $output_file = "rsid.tsv"; 
  }
  else {
    $output_file = $data_file; 
  }

  return $output_file; 
} 

=head2 parse_uniprot

parse UniProt file 

=cut 

sub parse_uniprot {
  my $infile = shift;
  
  my %synonym;
    
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open (IN, "zcat $infile |") or die ("Could not open $infile for reading");
  }
  else {
    open (IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    # A regexp to catch the meta information in the header. Just echo this to the stdout for logging purposes
    if ($_ =~ m/^\s*(Description|Name|Release)\:/) {
      print STDOUT $_ . "\n";
    }
    
    # Main regexp to extract relevant variation information
    if ($_ =~ m/^(\S+)\s+\S+\s+(VAR\_\d+)\s+\w\.\S+\s+(LP\/P|LB\/B|US)\s+(\-|rs\d*)\s*(.*)$/) {
      
      # Get the data that was caught by the regexp
      my $gene = $1;
      my $uniprot_id = $2;
      my $rs_id = $4;
      my $phenotype = $5;
      
      # If no rsId was given, will attempt to get one by looking up the Uniprot id in the synonym table
      if ($rs_id ne '-') {
        push(@{$synonym{$rs_id}},$uniprot_id);
      }
    }
  }
  close(IN);
  
  my %result = ('synonyms' => \%synonym);
  return \%result;
}

=head2 add_set_uniprot

add UniProt variation set 

=cut

sub add_set_uniprot {  
  my $db_adaptor      = shift;
  my $set_name        = shift; 
  my $set_description = shift;

  my $attrib_name = 'ph_uniprot';

  # check if UniProt already exists 
  my $stmt = qq{
    SELECT
      variation_set_id
    FROM
      variation_set
    WHERE
      name = '$set_name'
    LIMIT 1
  };

  my $sth = $db_adaptor->dbc->prepare($stmt);
  $sth->execute();
  my $variation_set_id;
  $sth->bind_columns(\$variation_set_id);
  $sth->fetch();

  if (!defined($variation_set_id)) {  
    my $stmt_attrib = qq{
      SELECT
        attrib_id
      FROM
        attrib
      WHERE
        value = '$attrib_name'
      LIMIT 1
    };
    my $sth_attrib = $db_adaptor->dbc->prepare($stmt_attrib);
    $sth_attrib->execute();
    my $attrib_id;
    $sth_attrib->bind_columns(\$attrib_id);
    $sth_attrib->fetch();

    if(!defined($attrib_id)){
      die "Error: attrib '$attrib_name' is not defined in attrib table.\n";
    }

    my $ins_stmt = qq{
      INSERT IGNORE INTO
        variation_set (
        name,
        description,
        short_name_attrib_id
        )
      VALUES (
        '$set_name',
        '$set_description',
        '$attrib_id' 
      )
    };

    $db_adaptor->dbc->do($ins_stmt); 
    $variation_set_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
 
  } 
  
  return $variation_set_id; 
}

=head2 add_synonyms_uniprot

import UniProt synonyms;

=cut

sub add_synonyms_uniprot {
  my $synonyms           = shift;
  my $source_name        = shift;
  my $variation_set_id   = shift; 
  my $db_adaptor         = shift;
  my $species            = shift; 

  # If we actually didn't get any synonyms, just return
  return if (!defined($synonyms) || !scalar(keys(%{$synonyms})));
 
  my $variation_adaptor = $db_adaptor->get_VariationAdaptor($species, 'variation', );

  my $stmt_check_var_set = qq{
    SELECT 
      variation_set_id
    FROM 
      variation_feature 
    WHERE 
      variation_id=?
  }; 
 
  my $stmt_update_after_check = qq{
    UPDATE
      variation_feature 
    SET 
      variation_set_id = ('$variation_set_id')
    WHERE
      variation_id=?
  }; 

  my $stmt_update = qq{
    UPDATE
      variation_feature 
    SET 
      variation_set_id = CONCAT_WS(',',variation_set_id,'$variation_set_id')  
    WHERE
      variation_id=?
  }; 
  
  my $stmt_insert_var_set = qq{
    INSERT IGNORE INTO
      variation_set_variation (variation_id, variation_set_id)
    VALUES (?, ?)
  };

  my $check_variation_set_sth = $db_adaptor->dbc->prepare($stmt_check_var_set); 
  my $update_after_check_sth = $db_adaptor->dbc->prepare($stmt_update_after_check);
  my $update_variation_set_sth = $db_adaptor->dbc->prepare($stmt_update);
  my $insert_var_set = $db_adaptor->dbc->prepare($stmt_insert_var_set);
  
  foreach my $rs_id (keys %{$synonyms}) {
    my $variation = $variation_adaptor->fetch_by_name($rs_id);
    unless($variation){
      warn "Variant $rs_id not found\n"; 
      next; 
    }

    my $var_dbid = $variation->dbID();

    $variation->add_synonym($source_name, @{$synonyms->{$rs_id}}); 
    $variation_adaptor->store_synonyms($variation); 
    
    $check_variation_set_sth->bind_param(1, $var_dbid, SQL_VARCHAR);
    $check_variation_set_sth->execute();
    my $var_set_id;

    $check_variation_set_sth->bind_columns(\$var_set_id);
    $check_variation_set_sth->fetch();

    if(!defined($var_set_id)){
      $update_after_check_sth->execute($var_dbid) || die ("Failed to update variation_set_id: variation_id $var_dbid\n");
    }
    else{
      $update_variation_set_sth->execute($var_dbid) || die ("Failed to update variation_set_id: variation_id $var_dbid\n");
    }

    # Insert into variation_set_variation
    $insert_var_set->execute($var_dbid, $variation_set_id) || die ("Failed to insert $rs_id into variation_set_variation\n");
  }

}

=head2 extract_PharmGKB

extract data from PharmGKB file

=cut

sub extract_PharmGKB {

  my $data_file = shift;

  my %synonyms;

  open my $rslist, $data_file ||die "Failed to open synonym list to load: $!\n";
  while(<$rslist>) {
    next if/RSID/;

    my @a = split/\t/;
    $synonyms{$a[0]} = $a[3];
  }
  return \%synonyms;
}

=head2 parse_eva_file

parse EVA data file

=cut

sub parse_eva_rsid_file {
  my $infile = shift;
  
  my %pos_rsid;
  my $rsid;

  # Open the input file for reading
  my $fh;
  if($infile =~ /gz$/) {
    open ($fh, "zcat $infile |") or die ("Could not open $infile for reading\n");
  }
  else {
    open ($fh, '<', $infile) or die ("Could not open $infile for reading\n");
  }

  # Read through the file and parse out the desired fields: 
  # current rs id and synonym
  while (my $row = <$fh>) {
    chomp $row;
    next if($row =~ /^#/);
   #1       74      rs5413863902    G       C       .       .       SID=PRJEB51985;SS_VALIDATED=0;VC=SO:0001483
   #1       84      rs5413863040    TA      T       .       .       SID=PRJEB51985;SS_VALIDATED=0;VC=SO:0000159
   #6       46567780        rs5437819034    AAT     A       .       .       SID=PRJEB51985;SS_VALIDATED=0;VC=SO:0000159
   #6       46567783        rs5437818934    T       TCGAAA  .       .       SID=PRJEB51985;SS_VALIDATED=0;VC=SO:0000667
   #6       46567787        rs5437819023    C       CGCAA   .       .       SID=PRJEB51985;SS_VALIDATED=0;VC=SO:0000667
   # the context base is used in case of Indels, in the case of ref alle and alt alle not same in length, need to add 1 to the pos
 
    my @row_data = split /\t/, $row;

    my $chr = $row_data[0];
    my $coord = $row_data[1];
    
    $rsid = $row_data[2];

    my $ref = uc $row_data[3];
    my $alt = uc $row_data[4];
    my $attr = uc $row_data[7];

    if( length $ref != length $alt ){
	$coord += 1;
	if( length $ref > length $alt){ #deletion
		$ref =~ s/^$alt//i ;
		$alt = '-';
	}else{ #insertion
		$alt =~ s/^$ref//i;
		$ref = '-';
	}
    }

    my $pos_key = join (':', ($chr, $coord));
    $pos_rsid{$pos_key} = [$rsid, $ref, $alt, $attr];

warn("DEBUG pos_key=$pos_key, $rsid, $ref, $alt, $attr\n");
  }
  close($fh);

  return \%pos_rsid;
}

=head2 import_rsid

import rsid from refhash and source object;

=cut

sub import_rsid {

  my $pos_rsid = shift;
  my $source   = shift;
  my $dba      = shift;
  my $species  = shift;

  my $variation_adaptor = $dba->get_VariationAdaptor($species, 'variation' );
  my $variationfeature_adaptor = $dba->get_VariationFeatureAdaptor( $species, 'variation_feature' );

  my $update_cnt = 0;
  my $already_in = 0;
  my $new = 0;
  my %newv; 
  my $dup = 0;
  my $inconsist = 0;

SNPL:  foreach my $pos (keys %{$pos_rsid}) {

    my ($rsid, $ref_allele, $alt_allele, $info) = @{$pos_rsid->{$pos}};
    my $attrib;

    map{my @t = split '='; $attrib->{$t[0]} = $t[1]}(split ';', $info);

    my ($chr, $chr_start) = split /:/, $pos;	
    my $loc_id = "$pos:${ref_allele}_${alt_allele}";
    my $var_features = $variationfeature_adaptor->fetch_all_by_location_identifier($pos);

#use Data::Dumper;
#warn( "warn var_features  is ", Dumper($var_features ) ,"\n" );

    unless (@{$var_features}){
    	warn("No snps on $pos found, variant $rsid is new, create new variation");
	create_variation_record($pos, $rsid, $ref_allele, $alt_allele, $info, $source, $variation_adaptor, $variationfeature_adaptor, \%newv) && $new++; 
	next SNPL;	
    }

 #   next SNPL;
	
    my $var_feature = undef;  
    unless( scalar @{$var_features} == 1 ){

	#let's find out if one of the feature does match our rsID pos, ref and alt
	my $match = 0;
	VF: foreach my $vf( @{$var_features} ){
		my $start = $vf->seq_region_start();
		unless ($start == $chr_start){
			warn ("DB start $start does not match vcf pos $loc_id, keep looking for next db vf\n");
			next VF; 
		}

		my $db_ref_alle =  $vf->ref_allele_string;
		my $db_alt_alleles = $vf->alt_alleles;
		unless ($db_ref_alle and $ref_allele =~ /^$db_ref_alle$/i ){
		   warn("WARN: vcf ref $ref_allele (variant:$rsid:$pos) does not match db ref allele $db_ref_alle, keep looking for next db vf\n");
		   next VF;
		}
		my $match_altallele = 0;
		map{
                        ${match_altallele}++ && last if ($_ and $_ =~ /^$alt_allele$/i );
                    }@${db_alt_alleles};

		if ($match_altallele){
		   $var_feature = $vf;
		   last VF;
		}
	}
	unless (defined $var_feature){
       	   warn("More than one snp on $loc_id found, none of them match db record, create new variant $rsid");
	   #$dup++;
           create_variation_record($pos, $rsid, $ref_allele, $alt_allele, $info, $source, $variation_adaptor, $variationfeature_adaptor, \%newv ) && $new++;
	   next SNPL; 
	}
    }

    $var_feature = shift @{$var_features} unless defined $var_feature;
    warn "var_feature is undefined" && return undef unless defined $var_feature;
    my $var = $var_feature->variation;

    if (lc($var->name) eq lc($rsid)){
       warn( "$rsid already in the database, no need to update, skip\n" );
       $already_in++;
       next SNPL;
    }
    # santiy check
    unless ( $no_allele){
    	my @all_alleles = map{ $_->allele } @{$var->get_all_Alleles}; 
    	my ($found_ref, $found_alt);
    	map{ $found_ref = 1 if (uc $_ eq uc $ref_allele) } @all_alleles; 
    	map{ $found_alt = 1 if (uc $_ eq uc $alt_allele) } @all_alleles; 
    	unless( $found_ref && $found_alt){
		warn("WARN: cannot find ref_allele($ref_allele) or alt_allele($alt_allele) for variant:$rsid:$pos in db, create new record\n");
		create_variation_record($pos, $rsid, $ref_allele, $alt_allele, $info, $source, $variation_adaptor, $variationfeature_adaptor, \%newv) && $new++;
		next SNPL;
    	}
    }
    #update variation_feature    

    $var->add_synonym($var->source_name(), $var->name);
    $variation_adaptor->store_synonyms($var);

    my $vid = $var->dbID;
    my $srid = $source->dbID;

    eval{
	my $sth_vf = $variationfeature_adaptor->prepare("update variation_feature set variation_name=?, source_id=? where variation_id=?");
    	$sth_vf->execute($rsid, $srid, $vid);
	# or die "Cannot update variation_feature $vid with name $rsid and source_id $srid";
    	$sth_vf->finish;

    	my $sth_v = $variation_adaptor->prepare("update variation set name=?, source_id=? where variation_id=?");
    	$sth_v->execute($rsid, $srid, $vid);
	# or die "Cannot update variation $vid with name $rsid and source_id $srid";
    	$sth_v->finish; 
    };
    if( $@ ) { 
	warn("$@");
	next SNPL;
    };
    $update_cnt++;
  } #end of SNPL

    print "Updated $update_cnt rsIDs\nSkipped $already_in already_in, $new new, $newv{skipped} new skipped, $newv{stored} new inserted, $dup duplicated, $inconsist inconsist SNPs\n";

}


sub create_variation_record {
	
	my ( $pos, $rsid, $ref_allele, $alt_allele, $info, $source, $variation_adaptor, $variationfeature_adaptor, $newcnt_hash) = @_;

	my $v = Bio::EnsEMBL::Variation::Variation->new_fast({
						 name   => $rsid,
                                                 _source_id => $source->dbID,
                                                 is_somatic       => 0
						});

        my ($chr, $coord) = split ':', $pos;
        my $chr_slice = $slice_adaptor->fetch_by_region( 'toplevel', $chr );
        my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
               (-start   => $coord,
                 -end     => $coord,
                 -strand  => 1,
                 -slice   => $chr_slice,
                 -allele_string => "${ref_allele}/${alt_allele}",
                 -variation_name => $rsid,
                 -map_weight  => 1,
                 -variation => $v,
                 );

        eval{ 
		$variation_adaptor->store($v);
		#$vf->variation( $v );
		$vf->source( $v->source );
		$vf->is_somatic($v->is_somatic);
		$variationfeature_adaptor->store($vf); 
	 };
        if( $@ ){
                warn("Error store variation and feature: $rsid, $ref_allele, $alt_allele, $info\n$@");
                $newcnt_hash->{skipped}++;
                return undef;
        }

        $newcnt_hash->{stored}++;
	return $v;

}


=head2 import_synonyms

import synonyms from refhash and source object;

=cut

sub import_synonyms {

  my $synonyms = shift;
  my $source   = shift;
  my $dba      = shift;
  my $species  = shift;

  my $variation_adaptor = $dba->get_VariationAdaptor($species, 'variation', );

  foreach my $var_name (keys %{$synonyms}) {

    my $var = $variation_adaptor->fetch_by_name($var_name);
    unless($var){
      warn "variant $var_name not found\n";
      next;
    }
    $var->add_synonym($source->name(), $synonyms->{$var_name});
    $variation_adaptor->store_synonyms($var);
  }
}

=head2 import_synonyms_all_sources

import all synonyms from production database  

=cut

sub import_synonyms_all_sources {

  my $sources           = shift;  
  my $dest_db           = shift; 
  my $species           = shift; 
  my $data              = shift;

  my $variation_adaptor = $dest_db->get_VariationAdaptor($species, 'variation');  

  foreach my $synonym_name (keys %{$data}) { 
    my $var_name = $data->{$synonym_name}->{rs_name}; 
    my $var = $variation_adaptor->fetch_by_name($var_name); 

    unless($var) {
      warn "variant $var_name not found\n";
      next; 
    }
  
    my $source_id = $data->{$synonym_name}->{source_id};
    my $source = $sources->{$source_id}; 
    # Clean synonym name before insert in db - remove rsid from name  
    $synonym_name =~ s/.*\s//g; 
    $var->add_synonym($source->name(), $synonym_name);
    $variation_adaptor->store_synonyms($var); 

  } 
} 

=head2 add_synonym_pig_chip

add synonyms for Pig SNP Consortium

=cut

sub add_synonym_pig_chip {
  my $species = shift;
  my $dba = shift;
  my $int_dba = shift;
  my $source_id = shift;

  my $variation_adaptor = $dba->get_VariationAdaptor($species, 'variation', );

  my $synonym_ext_sth = $int_dba->dbc->prepare(qq[ SELECT rs_name, synonym_name FROM pig_synonym ]);
  $synonym_ext_sth->execute();
  my $all_syn = $synonym_ext_sth->fetchall_arrayref();

  my $synonym_ins_sth   = $variation_adaptor->dbc->prepare(qq[ INSERT IGNORE INTO variation_synonym (variation_id, name, source_id)
                                                     VALUES (?,?,?) ]);

  foreach my $synonym(@{$all_syn}){
    my $var = $variation_adaptor->fetch_by_name($synonym->[0]);
    if(defined $var){
      $synonym_ins_sth->execute( $var->dbID, $synonym->[1], $source_id );
    }
    else{
      warn "Variant $synonym->[0] not found\n";
    }
  }
}

=head2 get_source

get or add source object

=cut

sub get_source {

  my $species       = shift;
  my $dba           = shift;
  my $source_name   = shift;
  my $version       = shift;
  my $url           = shift;
  my $description   = shift;
  my $source_status = shift;

  my $source_adaptor = $dba->get_SourceAdaptor($species, 'variation', );
  my $source = $source_adaptor->fetch_by_name($source_name);

  if (defined $source) {
    warn("DEBUG source_name=$source_name, defined source object\n");
    ## do we need to update the version of an existing source?
    if(defined $version) { 
      $source->version($version);
      $source_adaptor->update_version($source);   
    }
    # Update version with current date 
    else{ 
      $version = localtime->strftime('%Y%m%d');
      $source->version($version);
      $source_adaptor->update_version($source);
    }
  }
  else{
    ## update enter new source
    print "Source information not held for $source_name - adding supplied info\n" unless defined $source ;
    $source = Bio::EnsEMBL::Variation::Source->new
       (-name           => $source_name,
        -url            => $url           || undef,
        -version        => $version       || undef,
        -description    => $description   || undef,
        -somatic_status => $source_status || undef,
        -data_types     => ['variation', 'variation_synonym']
      );
    eval{
      $source_adaptor->store($source);
    };
    die "ERROR storing source: $@\n" unless $@ eq ''; 
  }
  return $source;
}

=head2 get_all_sources

get sources object 

=cut

sub get_all_sources {

  my $dest_db  = shift;
  my $prod_dbc = shift; 

  # my %source_names; # hash with source id and corresponding source name 
  my %source_list; 

  my $source_adaptor = $dest_db->get_SourceAdaptor('human', 'variation'); 

  # Get sources from production that are linked to human synonym data 
  my $sql_query = "SELECT source_id,name,version,description,url FROM source WHERE source_id IN
                   (SELECT DISTINCT(source_id) FROM human_synonym)";
  my $source_output = query_synonyms_data($prod_dbc, $sql_query);  

  foreach my $result_from_source (@{$source_output}) {  
    # my $source_id = $result_from_source->[0]; 
    my $source_name = $result_from_source->[1]; 
    my $source = $source_adaptor->fetch_by_name($source_name);

    if(!defined($source)) { 
      # get info for source  
      $source = Bio::EnsEMBL::Variation::Source->new
         (-name        => $source_name, 
          -url         => $result_from_source->[4],
          -version     => $result_from_source->[2],
          -description => $result_from_source->[3], 
          -data_types  => ['variation_synonym']
        );

      eval{
        $source_adaptor->store($source);
      };
      die "ERROR storing source: $@\n" unless $@ eq ''; 
    } 
    $source_list{$result_from_source->[0]} = $source; 
  } 

  return \%source_list; 
} 

=head2 query_synonyms_data 

fetch data from production database   

=cut

sub query_synonyms_data {

  my $prod_dbc   = shift;
  my $sql_query  = shift; 

  my $statement = $prod_dbc->prepare($sql_query);  
  $statement->execute(); 
  my $source_output = $statement->fetchall_arrayref();  
  $statement->finish(); 

  return $source_output; 
} 

sub usage{

  die "\nUsage : import_variant_synonyms -registry [registry file] -source_name [name]
\n\tOptions:
\t         -data_file          [name of file to load]  - mandatory for EVA
\t         -source_name        [source]                - 'PharmGKB', 'dbSNP', 'UniProt', 'EVA' or 'pig_chip'
\t         -source_version     [version]
\t         -source_url         [url]
\t         -source_description [longer description]
\t         -species            [species]               - defaults to human
\t         -registry           [registry file]
\t         -db_name            [production db name]
\t         -host               [production db host]
\t         -port               [production db port]
\t         -user               [user]
\t         -pass               [pass]
\t         -clean                                      - remove old data
\t	   -no_allele				       - there is no data in allele table (this is the case when the vcf loaded does not have genotypes\n\n";
}
