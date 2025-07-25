=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package EnsEMBL::Web::Document::HTML::ToolsTable;

### Allows easy removal of items from template

use strict;

use EnsEMBL::Web::Document::Table;

use base qw(EnsEMBL::Web::Document::HTML);

sub render { 
  my $self        = shift;
  my $hub         = $self->hub;
  my $sd          = $hub->species_defs;
  my $sp          = $sd->ENSEMBL_PRIMARY_SPECIES;
  my $img_url     = $sd->img_url;
  my $is_bacteria = $sd->GENOMIC_UNIT eq 'bacteria'; 
  my $url;

  my $sitename = $sd->ENSEMBL_SITETYPE;
  my $html = '<h2>Processing your data</h2>';

  ## Table for online tools
  my $table = EnsEMBL::Web::Document::Table->new([
      { key => 'name',  title => 'Name',            width => '20%', align => 'left' },
      { key => 'desc',  title => 'Description',     width => '40%', align => 'left' },
      { key => 'tool',  title => 'Online tool',     width => '10%', align => 'center' },
      { key => 'limit', title => 'Upload limit',    width => '10%', align => 'center' },
      { key => 'code',  title => 'Download script', width => '10%', align => 'center' },
      { key => 'docs',  title => 'Documentation',   width => '10%', align => 'center' },
    ], [], { cellpadding => 4 }
  );

  my $tools_limit = '50MB';

  ## VEP
  if ($sd->ENSEMBL_VEP_ENABLED) {
    my $vep_link = $hub->url({'species' => $sp, qw(type Tools action VEP)});
    $table->add_row({
      'name'  => sprintf('<a href="%s" class="nodeco"><b>Variant Effect Predictor</b><br /><img src="%svep_logo_sm.png" alt="[logo]" /></a>', $vep_link, $img_url),
      'desc'  => 'Analyse your own variants and predict the functional consequences of known and unknown variants via our Variant Effect Predictor (VEP) tool.',
      'limit' => $tools_limit.'*',
      'tool'  => sprintf('<a href="%s" class="nodeco"><img src="%s16/tool.png" alt="Tool" title="Go to online tool" /></a>', $vep_link, $img_url),
      'code'  => sprintf('<a href="https://github.com/Ensembl/ensembl-tools/archive/release/%s.zip" rel="external" class="nodeco"><img src="%s16/download.png" alt="Download" title="Download Perl script" /></a>', $sd->ENSEMBL_VERSION, $img_url),
      'docs'  => sprintf('<a href="https://ensembl.org/info/docs/tools/vep/index.html"><img src="%s16/info.png" alt="Documentation" /></a>', $img_url)
    });
  }

  ## HMMER
  if ($sd->ENSEMBL_HMMER_ENABLED) {
    my $link = 'https://www.ebi.ac.uk/Tools/hmmer/search/phmmer';
    $table->add_row({
      'name'  => sprintf('<b><a class="nodeco" href="%s">HMMER</a></b>', $link),
      'desc'  => 'Quickly search our genomes for your protein sequence.',
      'limit' => '',
      'tool'  => sprintf('<a href="%s" class="nodeco"><img src="%s16/tool.png" alt="Tool" title="Go to online tool" /></a>', $link, $img_url),
      'code'  => '',
      'docs'  => ''
    });
  }

  ## BLAST
  if ($sd->ENSEMBL_BLAST_ENABLED) {
    my $link = $hub->url({'species' => $sp, qw(type Tools action Blast)});
    $table->add_row({
      'name' => sprintf('<b><a class="nodeco" href="%s">BLAST/BLAT</a></b>', $link),
      'desc' => 'Search our genomes for your DNA or protein sequence.',
      'tool' => sprintf('<a href="%s" class="nodeco"><img src="%s16/tool.png" alt="Tool" title="Go to online tool" /></a>', $link, $img_url),
      'limit' => $tools_limit,
      'code' => '',
      'docs' => sprintf('<a href="https://ensembl.org/Help/View?db=core;id=451" class="popup"><img src="%s16/info.png" alt="Documentation" /></a>', $img_url)
    });
  }


  ## ASSEMBLY CONVERTER
  if ($sd->ENSEMBL_AC_ENABLED) {
    my $link = $hub->url({'species' => $sp, qw(type Tools action AssemblyConverter)});
    $table->add_row({
      'name'  => sprintf('<b><a class="nodeco" href="%s">Assembly Converter</a></b>', $link),
      'desc'  => "Map (liftover) your data's coordinates to the current assembly.",
      'tool'  => sprintf('<a href="%s" class="nodeco"><img src="%s16/tool.png" alt="Tool" title="Go to online tool" /></a>', $link, $img_url),
      'limit' => $tools_limit,
      'code'  => '',
      'docs'  => '',
    });
  }

  ## ID HISTORY CONVERTER
  if ($sd->ENSEMBL_IDM_ENABLED) {
    my $link = $hub->url({'species' => $sp, qw(type Tools action IDMapper)});
    $table->add_row({
      'name'  => sprintf('<b><a class="nodeco" href="%s">ID History Converter</a></b>', $link),
      'desc'  => 'Convert a set of Ensembl IDs from a previous release into their current equivalents.',
      'tool'  => sprintf('<a href="%s" class="nodeco"><img src="%s16/tool.png" alt="Tool" title="Go to online tool" /></a>', $link, $img_url),
      'limit' => $tools_limit,
      'code'  => sprintf('<a href="https://github.com/Ensembl/ensembl-tools/tree/release/%s/scripts/id_history_converter" rel="external" class="nodeco"><img src="%s16/download.png" alt="Download" title="Download Perl script" /></a>', $sd->ENSEMBL_VERSION, $img_url),
      'docs'  => '',
    });
  }

  ## Allele frequency
  if ($sd->ENSEMBL_AF_ENABLED) {
    my $link = $hub->url({'species' => $sp, qw(type Tools action AlleleFrequency)});
    $table->add_row({
      'name'  => sprintf('<b><a class="nodeco" href="%s">Allele frequency calculator</a></b>', $link),
      'desc'  => "This tool calculates population-wide allele frequency for sites within the chromosomal region defined from a VCF file and populations defined in a sample panel file.",
      'tool'  => sprintf('<a href="%s" class="nodeco"><img src="%s16/tool.png" alt="Tool" title="Go to online tool" /></a>', $link, $img_url),
      'limit' => '',
      'code'  => '',
      'docs'  => '',
    });
  }

  ## VCF to PED
  if ($sd->ENSEMBL_VP_ENABLED) {
    my $link = $hub->url({'species' => $sp, qw(type Tools action VcftoPed)});
    $table->add_row({
      'name'  => sprintf('<b><a class="nodeco" href="%s">VCF to PED converter</a></b>', $link),
      'desc'  => "Parse a vcf file to create a linkage pedigree file (ped) and a marker information file, which together may be loaded into ld visualization tools like Haploview.",
      'tool'  => sprintf('<a href="%s" class="nodeco"><img src="%s16/tool.png" alt="Tool" title="Go to online tool" /></a>', $link, $img_url),
      'limit' => '',
      'code'  => sprintf('<a href="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/browser/vcf_to_ped_converter/version_1.1/vcf_to_ped_convert.pl" rel="external" class="nodeco"><img src="%s16/download.png" alt="Download" title="Download Perl script" /></a>', $img_url),
      'docs'  => '',
    });
  }


  $html .= $table->render;

  $html .= '* For larger datasets we provide an API script that can be downloaded (you will also need to install our Perl API, below, to run the script).';

  ## Table of other tools

  $html .= qq(<h2 class="top-margin">Accessing $sitename data</h2>);

  $table = EnsEMBL::Web::Document::Table->new([
      { key => 'name', title => 'Name', width => '20%', align => 'left' },
      { key => 'desc', title => 'Description', width => '30%', align => 'left' },
      { key => 'from', title => 'Get it from:', width => '30%', align => 'center' },
      { key => 'docs', title => 'Documentation', width => '10%', align => 'center' },
    ], [], { cellpadding => 4 }
  );

  ## BIOMART
  if ($sd->ENSEMBL_MART_ENABLED) {
    $table->add_row({
      'name' => '<b><a href="/biomart/martview">BioMart</a></b>',
      'desc' => "Use this data-mining tool to export custom datasets from $sitename.",
      'from' => qq(<a href="/biomart/martview">$sitename BioMart</a>),
      'docs' => sprintf('<a href="/info/data/biomart/index.html" class="popup"><img src="%s16/info.png" alt="Documentation" /></a>', $img_url)
    });
  }

  ## PERL API 
  my $ftp = $sd->ENSEMBL_FTP_URL;
  $table->add_row({
    'name' => '<b>Ensembl Perl API</b>',
    'desc' => 'Programmatic access to all Ensembl data using simple Perl scripts',
    'from' => qq(<a href="https://github.com/Ensembl">GitHub</a> or <a href="$ftp/ensembl-api.tar.gz" rel="external">FTP download</a> (current release only)),
    'docs' => sprintf('<a href="/info/docs/api/"><img src="%s16/info.png" alt="Documentation" /></a>', $img_url)
  });

  ## REST
  if (my $rest_url = $sd->ENSEMBL_REST_URL) {
    $table->add_row({
      "name" => sprintf("<b><a href=%s>Ensembl Genomes REST server</a></b>", $rest_url),
      'desc' => 'Access Ensembl data using your favourite programming language',
      "tool" => sprintf("<a href='%s' class='nodeco'><img src='%s16/tool.png' alt='Tool' title='Go to online tool' /></a>", $rest_url, $img_url),
      'code' => sprintf('<a href="https://github.com/EnsemblGenomes/eg-rest" rel="external" class="nodeco"><img src="%s16/download.png" alt="Download" title="Download source code" /></a>', $img_url),
      'docs' => sprintf('<a href="%s"><img src="%s16/info.png" alt="Documentation" /></a>', $sd->ENSEMBL_REST_DOC_URL || $rest_url, $img_url)
    });
  }
  $html .= $table->render;

  return $html;
}

1;
