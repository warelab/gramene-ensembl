=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

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

package EnsEMBL::Web::Document::HTML::FTPtable;

### This module outputs a table of links to the FTP site

use strict;
use warnings;

use HTML::Entities qw(encode_entities);

use EnsEMBL::Web::Document::Table;

use base qw(EnsEMBL::Web::Document::HTML);

sub render {
  my $self = shift;

  my $hub             = $self->hub;
  my $species_defs    = $hub->species_defs;

  ## Links go to the EnsemblPlants FTP archive of the release this build is based on (set in gramene/conf/SiteDefs.pm).
  ## Local database names carry the Gramene release (e.g. _core_70_116_11); the archive uses the EG release (_core_63_116_11).
  my $rel = $species_defs->EG_FTP_RELEASE_VERSION;
  (my $ftp_pub = $species_defs->EG_FTP_URL) =~ s{/[^/]+/?$}{};   # .../pub/plants -> .../pub, for pan_ensembl
  my $eg_db = sub { (my $db = shift // '') =~ s/^(.+?_(?:core|otherfeatures|funcgen|variation))_\d+_/${1}_${rel}_/; return $db; };

  my @species = $species_defs->valid_species;
  my %title = (
    dna       => 'Masked and unmasked genome sequences associated with the assembly (contigs, chromosomes etc.)',
    cdna      => 'cDNA sequences for protein-coding genes',
    prot      => 'Protein sequences for protein-coding genes',
    rna       => 'Non-coding RNA gene predictions',
    embl      => 'Ensembl Genomes database dumps in EMBL nucleotide sequence database format',
    genbank   => 'Ensembl Genomes database dumps in GenBank nucleotide sequence database format',
    gtf       => 'Gene sets for each species in GTF format. These files include annotations of both coding and non-coding genes',
    gff3      => 'Gene sets for each species in GFF3 format. These files include annotations of both coding and non-coding genes',
    emf       => 'Alignments of resequencing data from the Compara database',
    gvf       => 'Variation data in GVF format',
    vcf       => 'Variation data in VCF format',
    vep       => 'Cache files for use with the VEP script',
    coll      => 'Additional regulation data (not in database)',
    bed       => 'Constrained elements calculated using GERP',
    files     => 'Additional release data stored as flat files rather than MySQL for performance reasons',
    ancestral => 'Ancestral Allele data in FASTA format',
    bam       => 'Alignments against the genome',
    core      => '%s core data export',
    otherfeatures     => '%s other features data export',
    variation      => '%s variation data export',
    funcgen   => '%s funcgen data export',
    pan       => 'Pan-taxomic Compara data export',
    compara   => '%s Compara data export',
    mart      => '%s BioMart data export',
    tsv       => 'Tab separated files containing selected data for individual species and from comparative genomics',

  );
  $title{$_} = encode_entities($title{$_}) for keys %title;

  my @rows;
  foreach my $spp (sort @species) {
    (my $sp_name = $spp) =~ s/_/ /;
     # FTP folders are named by production name (e.g. avena_sativa_gca947311555v1cm), not the URL name (Avena_sativa_aslak)
     my $sp_dir = $species_defs->get_config($spp, 'SPECIES_PRODUCTION_NAME') || lc($spp);
     my $sp_var = lc($spp).'_variation';
     my $common = $species_defs->get_config($spp, 'SPECIES_COMMON_NAME');

    my $genomic_unit = $species_defs->get_config($spp, 'GENOMIC_UNIT');
    my $collection;
    if($genomic_unit =~ /bacteria/i){
      my $group = lc $species_defs->get_config($spp, 'SPECIES_GROUP');
      $collection = $group . '_collection/'  if $group;
    }
    my $ftp_base_path_stub = "$ftp_pub/$genomic_unit/release-$rel";
    my @mysql;
    foreach my $db( qw/core otherfeatures funcgen variation/){
      my $db_config =  $species_defs->get_config($spp, 'databases')->{'DATABASE_' . uc($db)};
      if($db_config){
        my $title = sprintf($title{$db}, $sp_name);
        my $db_name = $eg_db->($db_config->{NAME});
        push(@mysql, qq{<a rel="external" title="$title" href="$ftp_base_path_stub/mysql/$db_name/">MySQL($db)</a>});
      }
    }
        
    my $data = {
species    => qq{<strong><i>$sp_name</i></strong>},
dna        => qq{<a rel="external"  title="$title{'dna'}" href="$ftp_base_path_stub/fasta/$sp_dir/dna/">FASTA</a> (DNA)},
cdna       => qq{<a rel="external"  title="$title{'cdna'}" href="$ftp_base_path_stub/fasta/$sp_dir/cdna/">FASTA</a> (cDNA)},
prot       => qq{<a rel="external"  title="$title{'prot'}" href="$ftp_base_path_stub/fasta/$sp_dir/pep/">FASTA</a> (protein)},
embl       => qq{<a rel="external"  title="$title{'embl'}" href="$ftp_base_path_stub/embl/} . $sp_dir . qq{/">EMBL</a>},
genbank    => qq{<a rel="external"  title="$title{'genbank'}" href="$ftp_base_path_stub/genbank/} . $sp_dir . qq{/">GenBank</a>},
gtf        => qq{<a rel="external"  title="$title{'gtf'}" href="$ftp_base_path_stub/gtf/$sp_dir/">GTF</a>},
gff3       => qq{<a rel="external"  title="$title{'gff3'}" href="$ftp_base_path_stub/gff3/$sp_dir/">GFF3</a>},
mysql      => join('<br/>',@mysql),
tsv        => qq{<a rel="external"  title="$title{'tsv'}" href="$ftp_base_path_stub/tsv/$sp_dir/">TSV</a>},
vep        => qq{<a rel="external"  title="$title{'vep'}" href="$ftp_base_path_stub/variation/indexed_vep_cache/">VEP</a>},
    };
    # VCF-only variation species (variation_source.database = 0) have no GVF/VCF dumps on the FTP site
    my $var_db   = $species_defs->get_config($spp, 'databases')->{'DATABASE_VARIATION'};
    my $var_meta = $var_db && $var_db->{'meta_info'} && $var_db->{'meta_info'}{1};
    my $vcf_only = $var_meta && $var_meta->{'variation_source.database'} && $var_meta->{'variation_source.database'}[0] eq '0';
    if ($var_db && !$vcf_only) {
      $data->{'gvf'} = qq{<a rel="external" title="$title{'gvf'}" href="$ftp_base_path_stub/variation/gvf/$sp_dir/">GVF</a>};
      $data->{'vcf'} = qq{<a rel="external" title="$title{'vcf'}" href="$ftp_base_path_stub/variation/vcf/$sp_dir/">VCF</a>};
    }
    push(@rows, $data);
  }

  my $genomic_unit = $species_defs->GENOMIC_UNIT;

  my $g_units = {
   'bacteria' => 'Bacterial',
   'fungi'    => 'Fungal',
   'metazoa'  => 'Metazoa',
   'plants'   => 'Plants',
   'protists' => 'Protists',
  };


  my $table    = EnsEMBL::Web::Document::Table->new(
    [
      {key=>'species',    sort=>'html', title=>'Species'},
      {key => 'dna',      sort=>'none', title => 'DNA'},    
      {key => 'cdna',     sort=>'none', title => 'CDNA'},   
      {key => 'prot',     sort=>'none', title => 'Protein'},    
      {key => 'embl',     sort=>'none', title => 'EMBL'},   
      {key => 'genbank',  sort=>'none', title => 'GENBANK'},
      {key => 'mysql',    sort=>'none', title => 'MySQL'},  
      {key => 'tsv',      sort=>'none', title => 'TSV'},    
      {key => 'gtf',      sort=>'none', title => 'GTF'},    
      {key => 'gff3',     sort=>'none', title => 'GFF3'},   
      {key => 'gvf',      sort=>'none', title => 'GVF'},    
      {key => 'vcf',      sort=>'none', title => 'VCF'},    
    ],
    \@rows,
    { data_table=>1, exportable=>0 }
  );
  $table->code = 'FTPtable::'.scalar(@rows);
  $table->{'options'}{'data_table_config'} = {iDisplayLength => 10};

  my $eg_compara  = sub { (my $db = shift // '') =~ s/_\d+_(\d+)$/_${rel}_$1/; return $db; };
  my $pan_compara = $eg_compara->($species_defs->get_config('MULTI', 'databases')->{DATABASE_COMPARA_PAN_ENSEMBL}->{NAME});
  my $compara     = $eg_compara->($species_defs->get_config('MULTI', 'databases')->{DATABASE_COMPARA}->{NAME});
  my $multi_sp = $g_units->{$genomic_unit};
  my $multi_table    = EnsEMBL::Web::Document::Table->new(
    [
      {key=>'database',    sort=>'html',title=>'Database'},
      {key => 'mysql',    sort=>'none',title => 'MySQL' },  
      {key => 'tsv',      sort=>'none',title => 'TSV'   },    
      {key => 'emf',      sort=>'none',title => 'EMF'   },    
    ],
    [
      {
      database => qq{<strong>Pan_compara Multi-species</strong>},
      mysql   => qq{<a rel="external" title="$title{pan}" href="$ftp_pub/pan_ensembl/release-$rel/mysql/$pan_compara/">MySQL</a>},
      emf     => qq{<a rel="external" title="$title{emf}" href="$ftp_pub/pan_ensembl/release-$rel/emf/ensembl-compara/homologies/">EMF</a>},
      tsv     => qq{<a rel="external" title="$title{tsv}" href="$ftp_pub/pan_ensembl/release-$rel/tsv/ensembl-compara/">TSV</a>},
      },
      {
      database => qq{<strong>$multi_sp Multi-species</strong>},
      mysql   => sprintf(qq{<a rel="external" title="%s" href="$ftp_pub/$genomic_unit/release-$rel/mysql/$compara/">MySQL</a>},sprintf($title{compara},ucfirst $genomic_unit)),
      emf     => qq{<a rel="external" title="$title{emf}" href="$ftp_pub/pan_ensembl/release-$rel/emf/ensembl-compara/homologies/">EMF</a>},
      tsv     => qq{<a rel="external" title="$title{tsv}" href="$ftp_pub/$genomic_unit/release-$rel/tsv/ensembl-compara/">TSV</a>},
      },
      {
      database => qq{<strong>Ensembl Mart</strong>},
      mysql   => sprintf(qq{<a rel="external" title="%s" href="$ftp_pub/$genomic_unit/release-$rel/mysql/$genomic_unit\_mart_$rel/">MySQL</a>},sprintf($title{mart},ucfirst $genomic_unit)),
      }
    ],
    { data_table=>0 }
  );
  $table->add_option('class','no_col_toggle');

  return sprintf(qq{<h3>Multi-species data</h3>%s<h3>Single species data</h3><div id="species_ftp_dl" class="js_panel"><input type="hidden" class="panel_type" value="Content"/>%s</div>},
    $multi_table->render,$table->render);
      
}



1;
