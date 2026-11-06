#!/usr/bin/env perl

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
            "registry=s"           => \$registry_file,            
            );

$species ||= 'homo_sapiens';

open my $fh, $data_file or die "cannot open $data_file to read";
my %vfid2allelestr = map {chomp; split /\t/;} <$fh>;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($registry_file);
my $dba = $reg->get_DBAdaptor($species, 'variation');

my $variationfeature_adaptor = $dba->get_VariationFeatureAdaptor( $species, 'variation_feature' );
my $sth_vf = $variationfeature_adaptor->prepare("update variation_feature set allele_string=? where variation_feature_id=?");


foreach my $vfid (keys %vfid2allelestr) {
#	print "key=$vfid, $vfid2allelestr{$vfid}\n";

	my $allele_string = $vfid2allelestr{$vfid};
	$allele_string =~ s=_=/=g;
	$sth_vf->execute($allele_string, $vfid) or die "Cannot update variation_feature $vfid with allele_string $allele_string";
} 

$sth_vf->finish;


__END__

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all($registry_file);
my $dba = $reg->get_DBAdaptor($species, 'variation');

my $variation_adaptor = $dba->get_VariationAdaptor($species, 'variation' );
my $variationfeature_adaptor = $dba->get_VariationFeatureAdaptor( $species, 'variation_feature' );

open my $fh, $data_file or die "cannot open $data_file to read";
my @vfs =  
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
    	warn("No snps on $pos found, variant $rsid is new, skip");
	$new++;
	next SNPL;
    }
    
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

		my $db_ref_alle = uc $vf->ref_allele_string;
		my $db_alt_alleles = $vf->alt_alleles;
		unless ($db_ref_alle and $ref_allele eq $db_ref_alle){
		   warn("WARN: vcf ref $ref_allele (variant:$rsid:$pos) does not match db ref allele $db_ref_alle, keep looking for next db vf\n");
		   next VF;
		}
		my $match_altallele = 0;
		map{
                        ${match_altallele}++ && last if ($_ and uc $alt_allele eq uc $_);
                    }@${db_alt_alleles};
		if ($match_altallele){
		   $var_feature = $vf;
		   last VF;
		}
		#my $db_allele_str = uc $vf->allele_string();
		#my @alleles = split '/', $db_allele_str;
		#my $match_allele = 0;
		#map{   
		#    my $a = $_;
		#    $match_allele = 0;
		#    map{
		#	$match_allele++ && last if ($a eq $_);
		#    }@allele;
		#	warn("WARN: cannot find vcf $a (variant:$rsid:$pos) in one of db($db_allele_str), keep looking for next db vf\n") && (next VF)
		#	unless $match_allele; 
		#}( ${ref_allele}, ${alt_allele});
		#$var_feature = $vf && last if $match_allele;
	}
	unless (defined $var_feature){
       	   warn("More than one snp on $loc_id found, none of them match db record, variant $rsid is skipped");
	   $dup++;
	   next SNPL; 
	}
    }

    $var_feature = shift @{$var_features} unless defined $var_feature;
    die "var_feature is undefined" unless defined $var_feature;
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
    	unless( $found_ref ){
		warn("WARN: cannot find ref_allele $ref_allele for variant:$rsid:$pos in db, skip\n");
		$inconsist++;
		next SNPL;
    	}
    	unless( $found_alt ){
		warn("WARN: cannot find alt_allele $alt_allele for variant:$rsid:$pos in db, skip\n");
		$inconsist++;
		next SNPL;
    	}
    }
    #update variation_feature    

    $var->add_synonym($var->source_name(), $var->name);
    $variation_adaptor->store_synonyms($var);

    my $vid = $var->dbID;
    my $srid = $source->dbID;

    #$var_feature->name($rsid); 
    #$var_feature->source( $source );
    #$variationfeature_adaptor->update($var_feature);
    #$variationfeature_adaptor->store;
    my $sth_vf = $variationfeature_adaptor->prepare("update variation_feature set variation_name=?, source_id=? where variation_id=?");
    $sth_vf->execute($rsid, $srid, $vid) or die "Cannot update variation_feature $vid with name $rsid and source_id $srid";
    $sth_vf->finish;

    #$var->name($rsid);
    #$var->source( $source );
    #$var->update_attributes( $attrib );
    my $sth_v = $variation_adaptor->prepare("update variation set name=?, source_id=? where variation_id=?");
    eval {$sth_v->execute($rsid, $srid, $vid) or die "Cannot update variation $vid with name $rsid and source_id $srid";
    $sth_v->finish; };
    if( $@ ) { warn("$@")};
	
    #$variation_adaptor->update($var);
    #$variation_adaptor->store;
    
    $update_cnt++;
  }
    print "Updated $update_cnt rsIDs\nSkipped $already_in already_in, $new new, $dup duplicated, $inconsist inconsist SNPs\n";
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
