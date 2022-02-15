#!/usr/local/bin/perl 

=head1 NAME

set_canonical_transcripts_from_TRaCE.pl

TRaCE output format is a tab delimited file with columns
gene.stable_id
transcript.stable_id
rank
transcript length
CDS length
domain coverage (aa)
modified AED scores for each sample

This script only needs the first 3 columns. When rank is 1, set the canonical transcript id of the gene
To speed things up, fetch stable id to primary id lookup tables
=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use autodie;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use JSON;

=head1 SYNOPSIS

load_interpro_json.pl  [options] interproscan.1.json interproscan.2.json
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species
    --debug
    --nowrite don't actually update the database

=head1 OPTIONS

=over 4

=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name in the registry file

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=item B<--nowrite> 

   don't update the database

=back

=head1 ARGUMENTS

   interproscan json output files

=cut

my ($species, $registry);
my ($debug, $nowrite);
my $margin=undef;
{  							#Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
	      ,"nowrite"=>\$nowrite
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  							#pod2usage(2) if $margin<0;
}
# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );
my $dbh = $ENS_DBA->dbc->db_handle;

# get the interpro table
print STDERR "reading interpro table\n" if $debug;
my %iprids;
my $sql = "select interpro_ac,id from interpro";
my $sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($ipr,$id) = @$row;
  $iprids{$ipr}{$id} = 1;
}
$sth->finish;
# get analysis table
print STDERR "reading analysis table\n" if $debug;
my %analysisIdLUT;
$sql = "select analysis_id, logic_name from analysis";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($id,$logic_name) = @$row;
  print STDERR "analysis: $id, $logic_name\n" if $debug;
  $analysisIdLUT{uc($logic_name)} = $id;
}
$sth->finish;

my %library2logic_name = (
	CDD => 'CDD',
	COILS => 'NCOILS',
	PANTHER => 'HMMPANTHER',
	PFAM => 'PFAM',
	PRINTS => 'PRINTS',
	PROSITE_PROFILES => 'PFSCAN',
	SUPERFAMILY => 'SUPERFAMILY'
);
my @needed_analyses = qw(interpro2go pfam cdd ncoils hmmpanther superfamily prints pfscan);
my @missing_analyses;
for my $ln (@needed_analyses) {
	push @missing_analyses, $ln unless exists $analysisIdLUT{uc($ln)};
}
my %analysis_defaults = (
	program => 'InterProScan',
	program_version => '5.52-86.0'
);

my %analysis_fields = (
	interpro2go	=> {
		db => 'InterPro2Go'
	},
	pfam => {
		db => 'Pfam',
		db_version => '33.1'
	},
	cdd => {
		db => 'CDD',
		db_version => '3.18'
	},
	ncoils => {
		db => 'ncoils',
		db_version => '2.2.1'
	},
	hmmpanther => {
		db => 'PANTHER',
		db_version => '15.0'
	},
	superfamily => {
		db => 'SuperFamily',
		db_version => '1.75'
	},
	prints => {
		db => 'PRINTS',
		db_version => '42.0'
	},
	pfscan => {
		db => 'Prosite_profiles',
		db_version => '2021_01'
	}	
);
my %analysis_description_fields = (
	interpro2go => {
		description => 'InterPro2GO mapping, defined by InterPro.',
		display_label => 'InterPro2GO mapping',
		displayable => 0
	},
	pfam => {
		description => 'Protein domains and motifs from the <a rel="external" href="http://pfam.xfam.org">Pfam</a> database.',
		display_label => 'Pfam',
		web_data => '{"type":"domain"}'
	},
	cdd => {
		description => '<a href="https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml">Conserved Domain Database</a> models.',
		display_label => 'CDD',
		web_data => '{"type":"domain"}'
	},
	ncoils => {
		description => 'Coiled-coil regions predicted by <a rel="external" href="http://www.sciencemag.org/cgi/reprint/252/5009/1162">Ncoils</a>.',
		display_label => 'Coiled-coils (Ncoils)',
		web_data => '{"type":"feature"}'
	},
	hmmpanther => {
		description => '<a href="http://www.pantherdb.org/">PANTHER</a> families.',
		display_label => 'PANTHER',
		web_data => '{"type":"domain"}'
	},
	superfamily => {
		description => 'Protein domains and motifs from the <a rel="external" href="http://supfam.org/SUPERFAMILY">SUPERFAMILY</a> database.',
		display_label => 'Superfamily',
		web_data => '{"type":"domain"}'
	},
	prints => {
		description => 'Protein fingerprints (groups of conserved motifs) from the <a href="http://www.bioinf.manchester.ac.uk/dbbrowser/PRINTS/index.php">PRINTS</a> database.',
		display_label => 'Prints',
		web_data => '{"type":"domain"}'
	},
	pfscan => {
		description => 'Protein domains and motifs from the <a rel="external" href="http://prosite.expasy.org">PROSITE</a> profiles database.',
		display_label => 'PROSITE profiles',
		web_data => '{"type":"domain"}'
	}
);
if (@missing_analyses) {
	$nowrite and die "Error: missing analyses [" . join(',',@missing_analyses) . "]\n";
	for my $ln (@missing_analyses) {
		# first insert into the analysis table
		my @fields = ('logic_name');
		my @values = ($ln);
		for my $k (keys %analysis_defaults) {
			push @fields, $k;
			push @values, $analysis_defaults{$k};
		}
		for my $k (keys %{$analysis_fields{$ln}}) {
			push @fields, $k;
			push @values, $analysis_fields{$ln}{$k};
		}
		my $qmarks = join(',', map {'?'} @fields);
		my $insert_sql = "insert into analysis (".join(',',@fields).") VALUES ($qmarks)";
		my $insert_sth = $dbh->prepare($insert_sql) or die "cannot prepare $insert_sql\n";
		$insert_sth->execute(@values);
		my $analysis_id = $insert_sth->{mysql_insertid};
		$analysisIdLUT{uc($ln)} = $analysis_id;
		$insert_sth->finish;
		# now that we have an analysis_id, insert into analysis_description
		@fields = ('analysis_id');
		@values = ($analysis_id);
		for my $k (keys %{$analysis_description_fields{$ln}}) {
			push @fields, $k;
			push @values, $analysis_description_fields{$ln}{$k};
		}
		$qmarks = join(',', map {'?'} @fields);
		$insert_sql = "insert into analysis_description (".join(',',@fields).") VALUES ($qmarks)";
		$insert_sth = $dbh->prepare($insert_sql) or die "cannot prepare $insert_sql\n";
		$insert_sth->execute(@values);
		$insert_sth->finish;
	}
}
# get external db ids for GO and Interpro
my %external_db_id;
$sql = "select external_db_id, db_name from external_db where db_name IN ('GO','Interpro')";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
	my ($id,$db) = @$row;
	$external_db_id{$db} = $id;
}
$sth->finish;
$external_db_id{GO} and $external_db_id{Interpro} or die "missing external db for GO or Interpro\n";
# get translation lookup table
print STDERR "reading translation table\n" if $debug;
my %translationIdLUT;
$sql = "select translation_id, stable_id from translation";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($tid,$stable_id) = @$row;
  $translationIdLUT{$stable_id} = $tid;
}
$sth->finish;
# get xref lookup table for GO and Interpro
my %xref;
$sql = "select x.* from external_db e, xref x where x.external_db_id = e.external_db_id and e.db_name = 'GO' and x.info_type = 'DIRECT' and x.info_text = 'InterPro'";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
	my ($id, $external_db_id, $dbprimary_acc, $display_label, $version, $description, $info_type, $info_text) = @$row;
	$xref{GO}{$dbprimary_acc} = $row;
}
$sth->finish;
$sql = "select x.* from external_db e, xref x where x.external_db_id = e.external_db_id and e.db_name = 'Interpro' and x.info_type = 'DIRECT'";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
	my ($id, $external_db_id, $dbprimary_acc, $display_label, $version, $description, $info_type, $info_text) = @$row;
	$xref{Interpro}{$dbprimary_acc} = $row;
}
$sth->finish;
# also download xrefs for objects and ontologies
$sql = "select object_xref_id, ensembl_id, xref_id from object_xref where ensembl_object_type = 'Translation' and analysis_id = ". $analysisIdLUT{INTERPRO2GO};
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
	my ($id, $ensembl_id, $xref_id) = @$row;
	$xref{object}{$ensembl_id}{$xref_id} = $id;
}
$sth->finish;


my $insert_xref_sql = "insert into xref (external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text) VALUES (?,?,?,?,?,?,?)";
my $insert_interpro_sql = "insert into interpro (interpro_ac, id) VALUES (?,?)";
my @pf_cols = qw(translation_id seq_start seq_end hit_start hit_end hit_name analysis_id score evalue);
my $pf_cols_str = join(',',@pf_cols);
my $question_marks = join(',', map {'?'} @pf_cols);
my $insert_protein_feature_sql = "insert into protein_feature ($pf_cols_str) VALUES ($question_marks)";
my $insert_object_xref_sql = "insert into object_xref (ensembl_id, ensembl_object_type, xref_id, analysis_id) VALUES (?,?,?,?)";
my $insert_ontology_xref_sql = "insert into ontology_xref (object_xref_id, source_xref_id, linkage_type) VALUES (?,?,?)";

my $insert_interpro_sth;
my $insert_protein_feature_sth;
my $insert_xref_sth;
my $insert_object_xref_sth;
my $insert_ontology_xref_sth;
 
unless ($nowrite){
  $insert_interpro_sth = $dbh->prepare($insert_interpro_sql) or die "cannot prepare $insert_interpro_sql\n";
  $insert_protein_feature_sth = $dbh->prepare($insert_protein_feature_sql) or die "cannot prepare $insert_protein_feature_sql\n";
	$insert_xref_sth = $dbh->prepare($insert_xref_sql) or die "cannot prepare $insert_xref_sql\n";
	$insert_object_xref_sth = $dbh->prepare($insert_object_xref_sql) or die "cannot prepare $insert_object_xref_sql\n";
	$insert_ontology_xref_sth = $dbh->prepare($insert_ontology_xref_sql) or die "cannot prepare $insert_ontology_xref_sql\n";
}

for my $jsonfile (@ARGV) {
  print STDERR "Process $jsonfile\n";
  open(my $fh, "<", $jsonfile);
  my $file_content = do {local $/; <$fh> };
  $file_content =~ s/}{/},{/g;
  my $ipr = decode_json $file_content;
  for my $res (@{$ipr->{results}}) {
    for my $xref (@{$res->{xref}}) {
      my $tid = $translationIdLUT{$xref->{id}};
      $tid or die "failed to get translation id for xref " . Dumper($xref);
      for my $match (@{$res->{matches}}) {
				# skip matches that have no interpro entry
				my $entry = $match->{signature}{entry};
				next unless $entry->{accession};
				# also skip if no library
				next unless $match->{signature}{signatureLibraryRelease}{library};
        my $hit_name = $match->{"model-ac"};
        my $ipr = $entry->{accession};
        if (not exists $iprids{$ipr}{$hit_name}) {
          print STDERR "insert interpro $ipr $hit_name\n" if $debug;
          $insert_interpro_sth->execute($ipr, $hit_name) if $insert_interpro_sth;
          $iprids{$ipr}{$hit_name}=1;
        }
				if (not exists $xref{Interpro}{$ipr}) {
					if ($insert_xref_sth) {
						my @values = ($external_db_id{Interpro}, $ipr, $entry->{name}, 0, $entry->{description}, "DIRECT","");
	          print STDERR "insert xref @values\n" if $debug;
						$insert_xref_sth->execute(@values);
						my $xref_id = $insert_xref_sth->{mysql_insertid};
						$xref{Interpro}{$ipr} = [$xref_id,@values];
					}
				}
				my $ipr_xref = $xref{Interpro}{$ipr}[0];
        my $score = $match->{score};
        my $analysis_id = $analysisIdLUT{$library2logic_name{$match->{signature}{signatureLibraryRelease}{library}}};
        $analysis_id or die "failed to get analysis_id for signature " . Dumper($match->{signature}{signatureLibraryRelease}{library});
				# add goXrefs
				for my $goXRef (@{$entry->{goXRefs}}) {
					my $go = $goXRef->{id};
					if (not exists $xref{GO}{$go} and $goXRef->{databaseName} eq 'GO') {
						if ($insert_xref_sth) {
							my @values = ($external_db_id{GO}, $go, $go, 0, $goXRef->{name}, "DIRECT", "InterPro");
							$insert_xref_sth->execute(@values);
		          print STDERR "insert xref @values\n" if $debug;
							my $xref_id = $insert_xref_sth->{mysql_insertid};
							$xref{GO}{$go} = [$xref_id,@values];
						}
					}
					my $go_xref = $xref{GO}{$go}[0];
					# create object xref
					if (not exists $xref{object}{$tid}{$go_xref}) {
						if ($insert_object_xref_sth) {
							my @values = ($tid, "Translation", $go_xref, $analysisIdLUT{INTERPRO2GO});
							$insert_object_xref_sth->execute(@values);
		          print STDERR "insert object_xref @values\n" if $debug;
							my $object_xref_id = $insert_object_xref_sth->{mysql_insertid};
							$xref{object}{$tid}{$go_xref}= $object_xref_id;
						}
					}
					# create ontology xref
					my $object_xref_id = $xref{object}{$tid}{$go_xref};
					if (not exists $xref{ontology}{$object_xref_id}{$ipr_xref}) {
						$insert_ontology_xref_sth->execute($object_xref_id, $ipr_xref, 'IEA');
	          print STDERR "insert ontology_xref $object_xref_id, $ipr_xref, 'IEA'\n" if $debug;
						$xref{ontology}{$object_xref_id}{$ipr_xref} = 'IEA';
					}
				}
        for my $location (@{$match->{locations}}) {
          my @pf = (
            $tid,
            $location->{start},
            $location->{end},
            $location->{hmmStart} || 0,
            $location->{hmmEnd} || 0,
            $hit_name,
            $analysis_id,
            $score || undef,
            $location->{evalue} || undef
          );
          # print STDERR "insert protein_feature @pf\n" if $debug;
          $insert_protein_feature_sth->execute(@pf) if $insert_protein_feature_sth;
        }
      }
    }
  }
}


$insert_interpro_sth->finish if $insert_interpro_sth;
$insert_protein_feature_sth->finish if $insert_protein_feature_sth;
$insert_xref_sth->finish if $insert_xref_sth;
$insert_object_xref_sth->finish if $insert_object_xref_sth;
$insert_ontology_xref_sth->finish if $insert_ontology_xref_sth;

$dbh->disconnect;

  
########################## subroutines ######################################

__END__


=head1 OUTPUT


=head1 AUTHOR

   Andrew Olson <olson@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

