# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::SyntenyTable
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::SyntenyTable

=head1 SYNOPSIS

  my $syntable = Bio::EnsEMBL::Analysis::RunnableDB::SyntenyTable->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -analysis => $analysis,
     );
  $syntable->fetch_input;
  $syntable->run;
  $syntable->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable SyntenyTable which wraps the program synteny_table.pl

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the marker_feature and marker tables

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::SyntenyTable;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning); 
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::SyntenyTable;
use Data::Dumper;

use vars qw(@ISA $BAND_TO_BASEPAIR);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);
$BAND_TO_BASEPAIR = 4096;

=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Arg [3]   : Bio::EnsEMBL::Analysis
  Function  : create a Bio::EnsEMBL::Analysis::RunnableDB object
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB
  Exceptions: throws if not passed either a dbadaptor, input id or
  an analysis object
  Example   : $rdb = $perl_path->new( -analysis => $self->analysis,
                                      -input_id => $self->input_id,
                                      -db => $self->adaptor->db );

=cut

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($target_db, $compara_db, $ort, $dist_thresh) = rearrange([ 'TARGET_DB', 'COMPARA_DB', 'ORT', 'DIST_THRESH' ], @args);
  if (!$target_db) {
	  throw("[*DIE] You must supply target db");
  }
  if (!$compara_db) {
	  throw("[*DIE] You must supply compara db");
  }

  $self->compara_db($compara_db);
  $self->target_db($target_db);
  $self->dist_thresh($dist_thresh);
  $self->{ort}=$ort;

  my ($query_species, $target_species, $query_taxon_id, $target_taxon_id);
  $query_species=$self->db->get_MetaContainer->get_scientific_name;
  $query_species=~s/\s/_/g;
  $target_species=$self->target_db->get_MetaContainer->get_scientific_name;
  $target_species=~s/\s/_/g;
  $query_taxon_id=$self->db->get_MetaContainer->get_taxonomy_id;
  $target_taxon_id=$self->target_db->get_MetaContainer->get_taxonomy_id;

  $self->{species_meta}={
	  query => {
		  name  => $query_species, #'Zea_mays',
		  core_adaptor => $self->db,
		  taxon_id => $query_taxon_id,
		  abbr => join('', map { (substr $_, 0, 1).(substr $_, -1)} split /_/, $query_species),
	  },
	  target => {
		  name  => $target_species, #'Oryza_sativa',
		  core_adaptor => $self->target_db,
		  taxon_id => $target_taxon_id,
		  abbr => join('', map {(substr $_, 0, 1).(substr $_, -1)} split /_/, $target_species),
	  },
  };
  return $self;
}

=head2 compara_db

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Compara::DBAdaptor
  Function  : container for dbadaptor
  Returntype: Bio::EnsEMBL::Compara::DBSQL::DBAdaptor
  Exceptions: throws if not passed a Bio::EnsEMBL::Compara::DBSQL::DBAdaptor
  Example   : 

=cut

sub compara_db {
  my $self = shift;
  my $compara_db = shift;
  if($compara_db){
    throw("[*DIE] Must pass RunnableDB:compara_db a Bio::EnsEMBL::Compara::DBSQL::DBAdaptor ".
          "not a ".$compara_db) 
      unless($compara_db->isa('Bio::EnsEMBL::Compara::DBSQL::DBAdaptor'));
    $self->{'compara_db'} = $compara_db;
  }
  return $self->{'compara_db'};
}

=head2 target_db

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : container for dbadaptor
  Returntype: Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions: throws if not passed a Bio::EnsEMBL::DBSQL::DBAdaptor
  Example   : 

=cut

sub target_db {
  my $self = shift;
  my $target_db = shift;
  if($target_db){
    throw("[*DIE] Must pass RunnableDB:target_db a Bio::EnsEMBL::DBSQL::DBAdaptor ".
          "not a ".$target_db) 
      unless($target_db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{'target_db'} = $target_db;
  }
  return $self->{'target_db'};
}

=head2 containers

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string
  Function  : container for specified variable. This pod refers to the
  methods below dist_thresh. These are simple
  containers which dont do more than hold and return an given value
  Returntype: string
  Exceptions: none
  Example   : my $ort = $self->ort;

=cut

sub dist_thresh {
  my $self = shift;
  $self->{'dist_thresh'} = shift if (@_);
  return $self->{'dist_thresh'};
}

=head2 fetch_input
    
  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::SyntenyTable
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut

sub fetch_input {
    my ($self) = @_;

	my $db=$self->db;
	my $compara_db=$self->compara_db;

	my $query_species=$self->{species_meta}{query}{query_species};
	my $target_species=$self->{species_meta}{target}{target_species};

	my @ort;
	if ( $self->{ort} ) {
		@ort=@{ $self->{ort} };
	} else {
	push @ort, join( "\t",
		 	'gene_stable_id',
		 	'sequence_name',
		 	'gene1_start',
		 	'gene1_end',
		 	'specie_id_1',
		 	'homolog2_member_id',
		 	'sequence_name2',
		 	'gene2_start',
		 	'gene2_end',
		 	'homology_relationship',
		 	'specie_id_2',
		); 

	##########################################################################
	# The following fetch data logic is adapted from Liya's script
	# brie:/home/liya/ensembl_script/dump_all_ortholog_genes.pl
	##########################################################################
	
	# Create the ADAPTORS that we will be using to fetch data from the database
	# A member can be a gene or a protein. Most of them are in the core databases of their correspondient specie.

	my $gene_adaptor = $db->get_adaptor('Gene') || throw( "[*DIE] Cannot db->get_adaptor('Gene')" );
	my $member_adaptor = $compara_db->get_adaptor('Member') || throw( "[*DIE] Cannot compara_db->get_adaptor('Member')" );
	my $HomologyAdaptor = $compara_db->get_adaptor('Homology') || throw("[*DIE] Cannot compara_db->get_adaptor('Homology')");

	my $dbname = $db->dbc->dbname;

	warn( "[iNFO] Collecting all genes from $dbname\n" );
	my $gene_id_arrayref = $gene_adaptor->list_stable_ids;
	my $total_genes = scalar @$gene_id_arrayref;
	warn( "[INFO] Collected genes $total_genes; processing...\n" );

	my( $num_genes, $num_coding, $num_chr_genes, $chr_specie_genes, $ignored_species) = (0,0,0,0,0);
	my $ignored=0;
	
	while( my $gene_id = shift @$gene_id_arrayref ){
		$num_genes ++; #COUNT
		
		# Get the gene object, and skip unless protein coding
		my $gene = $gene_adaptor->fetch_by_stable_id( $gene_id )
		|| warn( "[*DIE] Cannot fetch_by_stable_id gene $gene_id" );
		next unless $gene->biotype eq 'protein_coding';

		$num_coding ++; # COUNT
    	my ($homologue1,$homologue2);#added
    
    	######### TEST CODE
    	#liya : comment this line to get all genes instead of chromosome genes only
    	#unless($gene->seq_region_name =~ /^\d/){$ignored++; next;} # MANUAL MODIFICATION -turn this two "ifs" off when analizing maize
		$num_chr_genes++; #COUNT
    	
		#if( my $member =  $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$gene_id) ){		
		if( my $member =  $member_adaptor->fetch_by_stable_id($gene_id) ){		
 			#Obtain the homologues of the gene	
 			my $homologies = $HomologyAdaptor->fetch_all_by_Member($member);
 			foreach my $homology (@{$homologies}) {
				# if the relationship involves orthology
 				if ($homology->description =~ /ortholog/ ) {
 					($homologue1,$homologue2) = @{$homology->gene_list()};
					
					# To make sure that the $homologue2 gene is different from
					# the gene that it is been analyzed
  					unless ($homologue1->stable_id eq $gene_id) {
    					($homologue1,$homologue2) = ($homologue2,$homologue1);
  					}
  				
					if( $homologue2->taxon_id == $self->{species_meta}{target}{taxon_id} ){
					 	# Make sure the gene is on the longest assembled sequence
					 	$gene->project('toplevel');
					 
					 	$chr_specie_genes++;
						my $homologue_gene2=$homologue2->get_Gene();
						push @ort, join( "\t",
							 	$gene_id,
							 	$gene->seq_region_name,
							 	$gene->seq_region_start,
							 	$gene->seq_region_end,
							 	$homologue1->taxon_id,
							 	$homologue2->stable_id,
                                $homologue_gene2->seq_region_name,
                                $homologue_gene2->seq_region_start,
                                $homologue_gene2->seq_region_end,
							 	$homology->description,
							 	$homologue2->taxon_id, 
               			);
					}
					else{$ignored_species++;}
 				}
 			}
	  	}
  		# Update progress each 1000 genes
  		if( $num_genes % 1000 == 0 ){
      		warn("[INFO] processed $num_genes of $total_genes\n");
  		}
  	}

	##########################################################################
	}
	
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::SyntenyTable->new(
		-analysis => $self->analysis,
		#-program => "/data/ware/gramene-pipeline/bin/synteny_table.pl",
		-program => $self->analysis->program,
		-ort => \@ort,
		-query_species => $self->{species_meta}{query}{name},
		-target_species => $self->{species_meta}{target}{name},
		-dist_thresh => $self->dist_thresh,
    );
    $self->runnable($runnable);
    return 1;
}

=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : set analysis and slice on each feature
  Returntype: 1
  Exceptions: none
  Example   : 

=cut


sub write_output{
  my ($self) = @_;
  my $query_db = $self->db;
  my $target_db = $self->target_db;

  my $query_db_name=$query_db->dbc->dbname;
  my $target_db_name=$target_db->dbc->dbname;
#warn("DEBUG query_sp=$query_db_name, target_db=$target_db_name\n");

  my $query_gene_adaptor=$query_db->get_GeneAdaptor(); 
  my $target_gene_adaptor=$target_db->get_GeneAdaptor();

  my $switched = test_switched($query_gene_adaptor, $target_gene_adaptor, $self->output);
  
  my (%query_syn, %target_syn); 
  while (my $line= shift @{ $self->output }) {
	  next if $line=~/syn:type/i;
	  my @parts = (split( /\s+/, $line))[0, 3, 6, 8];
  	  my @empty = grep{ ! length($_) } @parts[0 .. 3];

      	  if( @empty ){
		  # Sanity check
		  warn("[WARN] Incomplete line in synteny file: $line\n");
       	  }

	  $query_syn{$parts[0]}||=[];
	  push @{ $query_syn{$parts[0]} }, [@parts[1 .. 3]];
	  $target_syn{$parts[1]}||=[];
	  push @{ $target_syn{$parts[1]} }, [@parts[0, 2, 3]];
  }

  my $query_attribute_adaptor = $query_db->get_AttributeAdaptor();
  my $target_attribute_adaptor = $target_db->get_AttributeAdaptor();
  ## check if query/target is switched
  if ($switched) {
  	  my %tmp=%query_syn;
  	  %query_syn=%target_syn;
  	  %target_syn=%tmp;
  }

  #warn ("DEBUG query_gene_adaptor is ", Dumper($query_gene_adaptor), "\n");
  #warn ("DEBUG target_gene_adaptor is ", Dumper($target_gene_adaptor), "\n");
  #exit;

  local $query_db->dbc->db_handle->{'AutoCommit'};
  local $target_db->dbc->db_handle->{'AutoCommit'};
  eval {
	  while (my ($key, $value) = each %query_syn) {
		  my $attr_value="$target_db_name;" . join("&", map { join "+", @$_ } @$value);
		  my $gene=$query_gene_adaptor->fetch_by_stable_id( $key );
		  warn ("DEBUG: query key=$key, attr_value=$attr_value\n");
		   #warn ("DEBUG query gene = $key", Dumper($gene), "\n");
		  #exit;
		  unless ($gene) {
		  	warn("[WARN] Can't find gene $key\n");
			next;
		  }
		  my $attribute=Bio::EnsEMBL::Attribute->new(
			  -CODE   => 'syn-gene-pairs',
			  -NAME   => 'syntenic gene pairs',
			  -DESCRIPTION    => 'syntenic gene relationship from gramene DAGchainer pipeline',
			  -VALUE  => $attr_value,
		  );
		  print "$key\t$gene\t$attr_value\n";
		  $query_attribute_adaptor->store_on_Gene($gene, [$attribute]);
	  }
	  print "Done loading query genes attrib to cores\n";
	  #exit;
	  while (my ($key, $value) = each %target_syn) {
		  my $attr_value="$query_db_name;" . join("&", map { join "+", @$_ } @$value);;
		  my $gene=$target_gene_adaptor->fetch_by_stable_id( $key );
		  unless ($gene) {
		  	warn("[WARN] Can't find gene $key\n");
			next;
		  }
		  my $attribute=Bio::EnsEMBL::Attribute->new(
			  -CODE   => 'syn-gene-pairs',
			  -NAME   => 'syntenic gene pairs',
			  -DESCRIPTION    => 'syntenic gene relationship from gramene DAGchainer pipeline',
			  -VALUE  => $attr_value,
		  );
		  print "$key\t$gene\t$attr_value\n";
		  $target_attribute_adaptor->store_on_Gene($gene, [$attribute]);
	  }
	  print "Done loading target genes attrib to cores\n";
	  $query_db->dbc->db_handle->commit;
	  $target_db->dbc->db_handle->commit;
  };
  if ($@) {
	  throw("[*DIE} Transaction aborted because $@");
	  eval {
		  $query_db->dbc->db_handle->rollback;
		  $target_db->dbc->db_handle->rollback;
	  };
  }
  return 1;
}

sub test_switched{

  my($query_gene_adaptor, $target_gene_adaptor, $output) = @_;

  my $switched = 0;

  my $first_line =  $output->[1];
  my ($g0, $s1, $s3) = (split( /\s+/, $first_line))[0, 1, 4];

  #we can only tell query from target apart by gene
  my $test_gene=$query_gene_adaptor->fetch_by_stable_id( $g0 );
  unless( $test_gene ){
        $test_gene=$target_gene_adaptor->fetch_by_stable_id( $g0 );
        if( $test_gene ){
                $switched = 1;
        }else{  
                throw("Cannot find gene $g0 in either query or target \n");
        }
  }
	
   return $switched;
}

1;
