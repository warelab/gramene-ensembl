# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::SyntenyTable
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME
    
Bio::EnsEMBL::Analysis::Runnable::SyntenyTable
    
=head1 SYNOPSIS

  my $synteny_table = Bio::EnsEMBL::Analysis::Runnable::SyntenyTable->
  new(
    -program => 'synteny_table.pl',
    -ort => $ort, #ort data in arrayref
	-query_species => "rice", 
	-target_species => "sorghum", 
    );
  $synteny_table->run;
  my @blocks = @{$synteny_table->output};

=head1 DESCRIPTION

Runs the program SyntenyTable.pl and produces synteny blocks 

=head1 CONTACT
    
Zhenyuan Lu &lt;luj@cshl.edu&gt;
Shiran Pasternak &lt;shiran@cshl.edu&gt;

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk
    
=cut

package Bio::EnsEMBL::Analysis::Runnable::SyntenyTable;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use File::Basename;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	
	#SETTING DEFAULTS#
	$self->program('synteny_table.pl') if (!$self->program);

	my ($ort, $query_species, $target_species, $dist_thresh) = rearrange([ 'ORT', 'QUERY_SPECIES', 'TARGET_SPECIES', 'DIST_THRESH' ], @args);

	if(!$ort){
		throw("[*DIE] You must supply ortholog data in an arrayref");
	}
	if(!$query_species){
		throw("[*DIE] You must supply query_species");
	}
	if(!$target_species){
		throw("[*DIE] You must supply target_species");
	}

	$dist_thresh ||= 5;
	$self->ort($ort);
	$self->query_species($query_species);
	$self->target_species($target_species);
	$self->dist_thresh($dist_thresh);
    return $self;
}

=head2 containers

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string
  Function  : container for specified variable. This pod refers to the
  four methods below ort, query_species, target_species. These are simple
  containers which dont do more than hold and return an given value
  Returntype: string
  Exceptions: none
  Example   : my $ort = $self->ort;

=cut

sub ort{
	my $self = shift;
	my $ort=shift;
	if ($ort) {
		throw("[*DIE] Must pass Runnable::SyntenyTable::ort an arrayref not a ".
			$ort) unless (ref($ort) eq 'ARRAY');
		$self->{'ort'} = $ort;
	}
	return $self->{'ort'};
}

sub query_species {
  my $self = shift;
  $self->{'query_species'} = shift if (@_);
  return $self->{'query_species'};
}

sub target_species {
  my $self = shift;
  $self->{'target_species'} = shift if (@_);
  return $self->{'target_species'};
}

sub dist_thresh {
  my $self = shift;
  $self->{'dist_thresh'} = shift if (@_);
  return $self->{'dist_thresh'};
}

=head2 ortfile

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::SyntenyTable
  Arg [2]   : string, filename
  Function  : will hold a given filename or if one is requested but none
  defined it will use the create_filename method to create a filename
  Returntype: string, filename
  Exceptions: none
  Example   : 

=cut


sub ortfile{
  my ($self, $filename) = @_;
  if($filename){
    $self->{'ortfile'} = $filename;
  }
  if (!$self->{'query_species'}) {
		throw("[*DIE] You must supply query species");
  }
  if (!$self->{'target_species'}) {
		throw("[*DIE] You must supply target species");
  }
  if(!$self->{'ortfile'}){
	my $name=join "_", map { join '', map {uc(substr ($_, 0, 1)).lc(substr ($_, 1))} split /_|\s/, $self->{$_} } ('query_species', 'target_species');
    $self->{'ortfile'} = $self->create_filename($name, 'ort');
  }
  if(!$self->resultsfile){
	my $resultsfile = $self->{'ortfile'};
	my $dist_thresh=$self->dist_thresh;
	$resultsfile=~s/\.ort$/.dist$dist_thresh.syn_table/;
	$self->resultsfile($resultsfile);
  }
  return $self->{'ortfile'};
}

=head2 write_ort_file

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::SyntenyTable
  Arg [2]   : arraryref
  Arg [3]   : filename
  Function  : This dumps an arrayref of ortholog data to a ort file
  Returntype: string, filename
  Exceptions: throw if failed to write sequence
  Example   : 

=cut


sub write_ort_file{
  my ($self, $ort, $filename) = @_;
 
  if(!$ort){
    $ort = $self->ort;
  }
  if(!$filename){
    $filename = $self->ortfile;
  }
  $filename = write_ortfile($ort, $filename);
  return $filename;
}

=head2 write_ort_file

  Arg [1]   : arraryref of ortholog data
  Arg [3]   : filename
  Function  : This dumps an arrayref of ortholog data to a ort file
  Returntype: string, filename
  Exceptions: throw if failed to write sequence
  Example   : 

=cut

sub write_ortfile{
	my ($ort, $filename)=@_;
	if (!$ort || !$filename) {
		throw("[*DIE] Must pass Runnable::SyntenyTable::ort an arrayref of ortholog data and a filename");
	} else {
		throw("[*DIE] Must pass Runnable::SyntenyTable::ort an arrayref not a ".
			$ort) unless (ref($ort) eq 'ARRAY');
		open FH, ">$filename" or throw("[*DIE] FAILED write to $filename");
		print FH "$_\n" foreach @$ort;
		close FH;
	}
	$filename
}

=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::SyntenyTable
  Arg [2]   : string, directory
  Function  : a generic run method. This checks the directory specifed
  to run it, write the query sequence to file, marks the query sequence
  file and results file for deletion, runs the analysis parses the 
  results and deletes any files
  Returntype: 1
  Exceptions: throws if no query sequence is specified
  Example   : 

=cut

sub run{
  my ($self, $dir) = @_;
  $self->workdir($dir) if($dir);
  throw("[*DIE] Can't run ".$self." without a ortholog ort") unless($self->ort);
  $self->checkdir();
  my $filename = $self->write_ort_file();
  my ($name, $workdir)=fileparse($filename, ".ort");
  $self->run_analysis();
  $self->parse_results;
  foreach my $file (glob("$workdir$name*")) {
      $self->files_to_delete($file);
  }
  $self->delete_files;
  return 1;
}

=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Essr
  Arg [2]   : string, program name
  Function  : constructs a commandline and runs the program passed
  in, the generic method in Runnable isnt used as essr.pl doesnt
  fit this module
  Returntype: none
  Exceptions: throws if run failed because system doesnt
  return 0 or the output file doesnt exist
  Example   : 

=cut

sub run_analysis {
    my ($self) = @_;

    my $cmd = join(" ", $self->program, "--distance", $self->dist_thresh, $self->ortfile);

    print "Running analysis " . $cmd . "\n";
	system "cp", $self->ortfile, "$ENV{HOME}";
    system($cmd) == 0
        or throw("[*DIE] FAILED to run " . $cmd . " SyntenyTable.pm:run_analysis");
	system "cp", $self->resultsfile, "$ENV{HOME}";

    if (!-e $self->resultsfile) {
        throw(    "[*DIE] FAILED to run " . $self->program . " on "
                . $self->ortfile . " "
                . $self->resultsfile
                . " has not been produced "
                . "SyntenyTable.pm:run_analysis");
    }
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::SyntenyTable
  Arg [2]   : string, filename
  Function  : open and parse the results file into synteny block data
  features
  Returntype: none 
  Exceptions: throws on failure to open or close output file
  Example   : 

=cut

sub parse_results {
    my ($self) = @_;

    # get our raw materials
	# my $query           = $self->queryfile;
    my $results         = $self->resultsfile;

    open(OUT, "$results")
        or throw("[*DIE] FAILED to open $results");

    while (my $line = <OUT>) {
        chomp $line;
        $self->output([$line]);
    }
    close(OUT) or throw("[*DIE] FAILED to close $results");
}

1;
