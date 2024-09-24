#!/usr/bin/env perl

package Script;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Bio::Seq;
use File::Find;
use File::Spec;
use File::Path qw/mkpath/;
use Getopt::Long qw/:config no_ignore_case auto_version bundling_override/;
use Pod::Usage;

my $rcsid = '$Revision$';
our ($VERSION) = $rcsid =~ /(\d+\.\d+)/;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->check();
  $self->setup();
  $self->process();
  return;
}

sub args {
  my ($self) = @_;
  my $opts = {};
  GetOptions(
    $opts, qw/
      dry
			overwrite
			registry|reg|r=s
      out_dir=s
      verbose 
      help
      man
      /
  ) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub opts {
  my ($self) = @_;
  return $self->{'opts'};
}

sub check {
  my ($self) = @_;
  my $o = $self->opts();

  my @required_params = qw/registry out_dir/;
	
  foreach my $r (@required_params) {
    if (!$o->{$r}) {
      pod2usage(
        -message => "-${r} has not been given at the command line but is a required parameter",
        -verbose => 1,
        -exitval => 1
      );
    }
  }
  
  foreach my $key (qw/out_dir/) {
    my $dir = $o->{$key};
    if(! -d $dir) {
      pod2usage(
        -message => "-${key} given location '${dir}' does not exist",
        -verbose => 1,
        -exitval => 1
      );
    }
  }
  
  return;
}

sub setup {
  my ($self) = @_;
  my $o = $self->opts();

	my $registry = $o->{'registry'};
	if (! -e $registry) {
		pod2usage(
			-message => "registry file '$registry' not found",
			-verbose => 1,
			-exitval => 1
		);
	}
	my $loaded = Bio::EnsEMBL::Registry->load_all( $registry );
	$self->v('Loaded %d DBAdaptor(s)', $loaded);
  
  return;
}


sub process {
  my ($self) = @_;
  my $dbas = $self->_get_dbs('core');
 
  foreach my $schema_type(keys %{$dbas}){
  
    while (my $dba = shift @{$dbas->{$schema_type}}) {
      $self->_process_dba($dba, $schema_type);
    }
  }
  
  return;
}

sub _process_dba {
  my ($self, $dba, $schema_type) = @_;
  my $species = $dba->get_MetaContainer()->get_production_name();
	my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
	my $dbName = $dba->get_MetaContainer()->single_value_by_key('species.display_name');
	my $o = $self->opts();
	my $species_path = join('/',$o->{'out_dir'},'fasta',$species);
  $self->v('Working with prod_name: %s assembly: %s', $species, $assembly);
	$self->v('path for fasta files: %s', $species_path);
	my @file_types = qw/cdna cds dna ncrna pep/;
	for my $type (@file_types) {
		my $target_dir = join("/",$species_path,$type);
		my $suffix = $type eq 'dna' ? 'toplevel' : 'all';
		my $file_name = ucfirst(join(".",$species,$assembly,$type,$suffix));
		if ($type eq 'ncrna') {
			$file_name = ucfirst(join(".",$species,$assembly,$type));
		}
		$self->_ensure_fasta($dba, $file_name, $target_dir, $type);
		$self->_ensure_blast($file_name, $species_path, $type, $dbName);
	}
  $dba->dbc()->disconnect_if_idle();
  return;
}

sub _ensure_fasta {
  my ($self, $dba, $datafile, $target_dir, $type) = @_;

  if(! -d $target_dir) {
    if($self->opts->{dry}) {
      $self->v("\tWould have created directory '%s'", $target_dir);
    }
    else {
      $self->v("\tCreating directory '%s'", $target_dir);
      mkpath($target_dir) or die "Cannot create the directory $target_dir: $!";
    }
  }

	my $file_path = join("/",$target_dir, $datafile . ".fa");
	if (not $self->opts->{overwrite} and (-e $file_path or -e "$file_path.gz")) {
		$self->v("skipping $file_path");
		return;
	}

	$self->v("dump $type into $file_path");
	$self->_dump_fasta($dba, $file_path, $type) unless ($self->opts->{dry});
	return;
}

sub _dump_fasta {
	my ($self, $dba, $outfile, $type) = @_;
	my $seqio_fh = new Bio::SeqIO(-format => 'fasta', -file => ">$outfile");

	my $slice_adaptor = $dba->get_SliceAdaptor;
	my $slices = $slice_adaptor->fetch_all('toplevel');
	if ($type eq 'dna') {
		
		my $outfile_sm = $outfile;
		$outfile_sm =~ s/dna\.toplevel/dna_sm.toplevel/;
		my $outfile_rm = $outfile;
		$outfile_rm =~ s/dna\.toplevel/dna_rm.toplevel/;

		my $seqio_unmasked = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile");
		my $seqio_sm = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile_sm");
		my $seqio_rm = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile_rm");

		foreach my $slice (@$slices) {
			my $seqstr = $slice->get_repeatmasked_seq(['RepeatMask'],1)->seq;
			$seqstr ||= $slice->seq;
			my $seq = Bio::Seq->new(
				-display_id => $slice->seq_region_name,
				-seq => $seqstr
			);
			$seqio_sm->write_seq($seq);
			$seqstr = $slice->get_repeatmasked_seq(['RepeatMask'],0)->seq;
			$seqstr ||= $slice->seq;
			$seq = Bio::Seq->new(
				-display_id => $slice->seq_region_name,
				-seq => $seqstr
			);
			$seqio_rm->write_seq($seq);
			$seqstr = $slice->seq;
			$seq = Bio::Seq->new(
				-display_id => $slice->seq_region_name,
				-seq => $seqstr
			);
			$seqio_unmasked->write_seq($seq);
		}
	}
	else {
		foreach my $slice (@$slices) {
            my $gs = $slice->get_all_Genes();
            for my $gene (@$gs) {
			    next if ($type eq 'ncrna' and $gene->biotype eq 'protein_coding');
			    next if ($type ne 'ncrna' and $gene->biotype ne 'protein_coding');
			    my @transcripts;
	            eval { @transcripts = @{ $gene->get_all_Transcripts } };
	            print STDERR "$@" && next if $@;

	            foreach my $trans (@transcripts) {
				    my $id = $trans->stable_id;
				    next if ($type eq 'ncrna' and $trans->biotype eq 'protein_coding');
				    next if ($type ne 'ncrna' and $trans->biotype ne 'protein_coding');

                    my $cdna_seq = $trans->spliced_seq;
				
			        if ($type eq 'cdna' or $type eq 'ncrna') {
				        $seqio_fh->write_seq(Bio::Seq->new(
						    -display_id => $id,
						    -seq => $cdna_seq
					        ));
				    }
				    if ($type eq 'cds') {
                        my $cdna_coding_start = $trans->cdna_coding_start;
                        my $cdna_coding_end   = $trans->cdna_coding_end;
                        my $seq_obj_cds       = Bio::Seq->new(
                            -display_id => $id,
                            -seq        => substr(
                                $cdna_seq,
                                $cdna_coding_start - 1,
                                $cdna_coding_end - $cdna_coding_start + 1
                                )
                                );

                        $seqio_fh->write_seq($seq_obj_cds);
				    }
				    if ($type eq 'pep') {
					    my $aa_obj = $trans->translate;
					    $seqio_fh->write_seq($aa_obj);
				    }
			    }
            }
		}
	}
}

sub _ensure_blast {
  my ($self, $datafile, $target_dir, $type, $dbName, $masked) = @_;

	# use the softmasked version if dna
	if ($type eq 'dna' and not $masked) {
		my $sm_file = $datafile;
		$sm_file =~ s/\.dna\./.dna_sm./;
		my $rm_file = $datafile;
		$rm_file =~ s/\.dna\./.dna_rm./;
		$self->_ensure_blast($sm_file, $target_dir, $type, $dbName . ' softmasked', 1);
		$self->_ensure_blast($rm_file, $target_dir, $type, $dbName . ' hardmasked', 1);
	}
	
	my $file_path = join("/",$target_dir, "$datafile");
	if (not $self->opts->{overwrite} and (-e "$file_path.nhr" or -e "$file_path.phr")) {
		$self->v("skipping BLAST $file_path");
		return;
	}

	$self->v("build blastdb $file_path");
	if (not $self->opts->{dry}) {
		my $fa_file = join("/",$target_dir,$type, $datafile . ".fa");
		if (not (-e $fa_file or -e "$fa_file.gz")) {
			$self->v("fasta file for blast db missing $fa_file");
			return
		}
        if (-z $fa_file) {
            $self->v("skipping BLAST : empty fasta file $fa_file");
            return
        }
		my $dbtype = $type eq 'pep' ? 'prot' : 'nucl';
		my $title = "$dbName $type";
		my $cmd = -e $fa_file
		? "makeblastdb -out $file_path -title \"$title\" -in $fa_file -dbtype $dbtype"
		: "gzip -cd $fa_file.gz | makeblastdb -out $file_path -title \"$title\" -in - -dbtype $dbtype";
		$self->v($cmd);
		system($cmd);
	}
	return;
}

sub _get_dbs {
  my ($self, $schema_type) = @_;
  my $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors();
  my %final_dbas;
  my %required_types = (core => 'core_like', funcgen => 'funcgen');  # Merge with _dba_to_ftp_type & make package var?

  if($schema_type){
    if(! exists $required_types{$schema_type}){
      die("Cannot _get_dbs schema_type is not supported:\t".$schema_type);
    }
    %required_types = ($schema_type => $required_types{$schema_type});
  }

  while(my $dba = shift @{$dbas}) {
    next if $dba->species() eq 'multi';
    next if lc($dba->species()) eq 'ancestral sequences';
    next if $dba->dbc()->dbname() =~ /^.+_userdata$/xms;
    
    my $type = $dba->get_MetaContainer()->single_value_by_key('schema_type');
    $dba->dbc()->disconnect_if_idle();
    next unless $type;

    if(exists $required_types{$type}){
      $final_dbas{$required_types{$type}} ||= [];
      push @{$final_dbas{$required_types{$type}}}, $dba;
    }
  }

  foreach my $rtype(values %required_types){
    my $num_dbs = (defined $final_dbas{$rtype}) ? scalar(@{$final_dbas{$rtype}}) : 0;   
    $self->v('Found %d '.$rtype.' like database(s)', $num_dbs);
  }

  return \%final_dbas;
}

sub v {
  my ($self, $msg, @params) = @_;
  return unless $self->opts()->{verbose};
  printf(STDERR $msg."\n", @params);
  return;
}

Script->run();

1;
__END__

=pod

=head1 NAME

ensure_blast.pl

=head1 SYNOPSIS

  #BASIC
  ./ensure_blast.pl -registry REG.PM -out_dir DIR [-dry] [-verbose] [-help | -man]
  
=head1 DESCRIPTION

A script which will export fasta files and build blast databases if necessary

=head1 OPTIONS

=over 8

=item B<--registry | --reg | -r>

REQUIRED. Ensembl registry file to process

=item B<--out_dir>

REQUIRED. Target directory

=item B<--verbose>

Makes the program give more information about what is going on. Otherwise
the program is silent.

=item B<--dry>

If specified the script will inform of the types of commands and actions it 
would have performed.

=item B<--overwrite>

If specified the script will overwrite files that are already present in the output directory.

=item B<--help>

Help message

=item B<--man>

Man page

=back

=head1 REQUIREMENTS

=over 8

=item Perl 5.8+

=item Bio::EnsEMBL

=item Post 66 databases

=back

=end

