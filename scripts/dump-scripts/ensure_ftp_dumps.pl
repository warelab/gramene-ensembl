#!/usr/bin/env perl

package Script;

use strict;
use warnings;
use lib '.';
use ExportView::GFF3ExporterMod;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::SequenceAdaptor;
use Bio::DB::HTS::Faidx;
use Bio::SeqIO;
use Bio::Seq;
use File::Find;
use File::Spec;
use File::Path qw/mkpath/;
use Getopt::Long qw/:config no_ignore_case auto_version bundling_override/;
use Pod::Usage;

# Serve genomic sequence from the bgzipped dna.toplevel FASTA this script has just written,
# instead of from the core's dna table. The API reads sequence 256 kb at a time with
# SUBSTR() on one LONGTEXT row per seq_region, and MyISAM reads the WHOLE row to do that:
# 84 s per chunk on a 640 Mb wheat chromosome, so cdna/cds/pep for a large genome took days
# per file and loaded the server for every other client. See _use_local_dna.
{
  no warnings 'redefine';
  my $db_fetch = \&Bio::EnsEMBL::DBSQL::SequenceAdaptor::_fetch_raw_seq;
  *Bio::EnsEMBL::DBSQL::SequenceAdaptor::_fetch_raw_seq = sub {
    my ($sa, $id, $start, $length) = @_;
    my $loc = $sa->{_gramene_local_dna};
    my $map = $loc ? $loc->{map}{$id} : undef;
    return $db_fetch->(@_) unless $map;
    # $map: [toplevel name, region length, toplevel start of region pos 1, ori]
    my ($name, $len, $off, $ori) = @$map;
    my $end = $start + $length - 1;
    $end = $len if $end > $len;
    my $seq = '';
    if ($start <= $end) {    # SUBSTR() past the end returns '', so do the same
      my ($s, $e) = $ori == 1
        ? ($off + $start - 1, $off + $end - 1)
        : ($off - $end + 1, $off - $start + 1);    # $off is then the toplevel END of pos 1
      $seq = $loc->{fai}->get_sequence2_no_length($name, $s - 1, $e - 1);    # htslib: 0-based, inclusive
      defined $seq and length($seq) == $e - $s + 1
        or die "local dna fetch failed for $name:$s-$e from $loc->{file}\n";
      $seq = uc $seq;
      if ($ori == -1) { $seq = reverse $seq; $seq =~ tr/ACGTRYKMBVDH/TGCAYRMKVBHD/; }
    }
    return \$seq;
  };
}

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
			repeat_logic=s
			bgzip_threads=i
      out_dir=s
      verbose 
      help
      man
      /
  ) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $opts->{repeat_logic} //= '^(repeatmask|dust|trf)';
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
      # One bad core must not strand the rest of the list: log it, carry on, fail at the end.
      next if eval { $self->_process_dba($dba, $schema_type); 1 };
      my $err = $@ || 'unknown error';
      my $db = eval { $dba->dbc->dbname } || '?';
      $self->v("GENOME FAILED %s: %s", $db, $err);
      push @{ $self->{failed} }, $db;
      eval { $dba->dbc->disconnect_if_idle };
    }
  }
  if (my @f = @{ $self->{failed} || [] }) {
    die sprintf("%d genome(s) failed: %s\n", scalar(@f), join(' ', @f));
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
	# Every sequence-level region needs a dna row, or the dna dump dies part-way (after
	# writing what it could) and no cdna/cds/pep can be made for genes on the missing regions.
	my $h = $dba->dbc->sql_helper;
	my $n_seq = $h->execute_single_result(-SQL => q{SELECT COUNT(*) FROM seq_region sr
		JOIN coord_system cs USING (coord_system_id)
		WHERE cs.species_id = ? AND FIND_IN_SET('sequence_level', cs.attrib)}, -PARAMS => [ $dba->species_id ]);
	my $n_dna = $h->execute_single_result(-SQL => q{SELECT COUNT(*) FROM dna});
	die "dna table has $n_dna rows for $n_seq sequence-level regions; this core is missing sequence\n"
		if $n_dna < $n_seq;
    # my @file_types = qw/cdna cds dna ncrna pep/;
	# dna first: once it is on disk, cdna/cds/pep read their sequence from it (_use_local_dna).
	my @file_types = qw/dna cdna cds pep/;
	for my $type (@file_types) {
		my $target_dir = join("/",$species_path,$type);
		my $suffix = $type eq 'dna' ? 'toplevel' : 'all';
		my $file_name = ucfirst(join(".",$species,$assembly,$type,$suffix));
		if ($type eq 'ncrna') {
			$file_name = ucfirst(join(".",$species,$assembly,$type));
		}
		$self->_ensure_fasta($dba, $file_name, $target_dir, $type);
		$self->_ensure_blast($file_name, $species_path, $type, $dbName);
		$self->_use_local_dna($dba, "$target_dir/$file_name.fa.gz")
			if ($type eq 'dna' and not $self->opts->{dry});
	}
    $species_path = join('/',$o->{'out_dir'},'gff3',$species);
    my $gff_name = ucfirst(join(".",$species,$assembly,"gff3"));
    $self->_ensure_gff3($dba, $gff_name, $species_path);
    delete $dba->get_SequenceAdaptor->{_gramene_local_dna};
    $dba->dbc()->disconnect_if_idle();
    return;
}

# Point this core's SequenceAdaptor at the finished dna.toplevel FASTA. A sequence-level
# region is served from it when it either is a toplevel region itself, or is placed whole
# (one assembly row covering it end to end) into a toplevel region; e.g. a contig that makes
# up a supercontig. The toplevel region must have the same name and length in the FASTA.
# Anything else (multi-level assemblies, partial placements) still comes from the database.
sub _use_local_dna {
    my ($self, $dba, $fa) = @_;
    my $sa = $dba->get_SequenceAdaptor;
    delete $sa->{_gramene_local_dna};
    unless (-s $fa and -s "$fa.fai" and -s "$fa.gzi") {
        $self->v("\tno indexed $fa; cdna/cds/pep sequence will come from the database");
        return;
    }
    my $fai = Bio::DB::HTS::Faidx->new($fa);
    my %fa_len = map { $_ => $fai->length($_) } $fai->get_all_sequence_ids;
    my $h = $dba->dbc->sql_helper;
    my $seqlevel = q{
        FROM seq_region c
        JOIN coord_system ccs ON ccs.coord_system_id = c.coord_system_id
             AND ccs.species_id = ? AND FIND_IN_SET('sequence_level', ccs.attrib)};
    my $is_top = q{
        JOIN seq_region_attrib ta ON ta.seq_region_id = t.seq_region_id
        JOIN attrib_type tt ON tt.attrib_type_id = ta.attrib_type_id AND tt.code = 'toplevel'};
    my %map;
    my $all = $h->execute_single_result(-SQL => "SELECT COUNT(*), IFNULL(SUM(c.length),0) $seqlevel",
        -PARAMS => [ $dba->species_id ], -CALLBACK => sub { [ @{$_[0]} ] });
    # 1. sequence-level regions that are toplevel themselves
    for my $r (@{ $h->execute(-SQL => "SELECT c.seq_region_id, c.name, c.length $seqlevel
            JOIN seq_region t ON t.seq_region_id = c.seq_region_id $is_top", -PARAMS => [ $dba->species_id ]) }) {
        my ($id, $name, $l) = @$r;
        $map{$id} = [ $name, $l, 1, 1 ] if defined $fa_len{$name} and $fa_len{$name} == $l;
    }
    # 2. sequence-level regions placed whole into a toplevel region
    for my $r (@{ $h->execute(-SQL => "SELECT c.seq_region_id, c.length, t.name, t.length,
                a.asm_start, a.asm_end, a.cmp_start, a.cmp_end, a.ori $seqlevel
            JOIN assembly a ON a.cmp_seq_region_id = c.seq_region_id
            JOIN seq_region t ON t.seq_region_id = a.asm_seq_region_id $is_top", -PARAMS => [ $dba->species_id ]) }) {
        my ($id, $clen, $tname, $tlen, $as, $ae, $cs, $ce, $ori) = @$r;
        next if $map{$id};
        next unless $cs == 1 and $ce == $clen and $ae - $as == $ce - $cs and abs($ori) == 1;
        next unless defined $fa_len{$tname} and $fa_len{$tname} == $tlen and $ae <= $tlen;
        $map{$id} = [ $tname, $clen, ($ori == 1 ? $as : $ae), $ori ];
    }
    # Self-check: 1 kb from the FASTA must equal 1 kb from the dna table, for one region of
    # each orientation. On any mismatch use the database for the whole genome.
    my %probe;
    for my $id (sort { $a <=> $b } keys %map) {
        $probe{ $map{$id}[3] } //= $id;    # one forward, one reverse (if any)
    }
    $sa->{_gramene_local_dna} = { fai => $fai, map => \%map, file => $fa };
    for my $id (values %probe) {
        my $len = $map{$id}[1];
        my $st = $len > 2000 ? int($len / 2) : 1;
        my $n = $len - $st + 1 < 1000 ? $len - $st + 1 : 1000;
        my $local = ${ $sa->_fetch_raw_seq($id, $st, $n) };
        delete local $sa->{_gramene_local_dna};
        my $db = ${ $sa->_fetch_raw_seq($id, $st, $n) };
        next if defined $db and $local eq uc $db;
        $self->v("\tWARNING local dna differs from the dna table for seq_region_id %d at %d; using the database", $id, $st);
        delete $sa->{_gramene_local_dna};
        return;
    }
    delete $sa->{_gramene_local_dna};
    my $bp = 0; $bp += $_->[1] for values %map;
    $self->v("\tlocal dna: %d of %d sequence-level regions (%.1f of %.1f Mb) served from %s",
        scalar(keys %map), $all->[0], $bp / 1e6, $all->[1] / 1e6, $fa);
    return unless %map;
    $sa->clear_cache;
    $sa->{_gramene_local_dna} = { fai => $fai, map => \%map, file => $fa };
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
    my @to_index = ($file_path);
    if ($type eq 'dna') {
        (my $sm = $file_path) =~ s/dna\.toplevel/dna_sm.toplevel/;
        (my $rm = $file_path) =~ s/dna\.toplevel/dna_rm.toplevel/;
        push @to_index, $sm, $rm;
    }
	# Finish indexing that was interrupted: a renamed .fa is complete, but if the run died during
	# bgzip/faidx it is still uncompressed (perhaps beside a partial .fa.gz), or the .fai is missing.
	if (not $self->opts->{overwrite} and not $self->opts->{dry}
		and grep { -s $_ or (-s "$_.gz" and not -s "$_.gz.fai") } @to_index) {
		for my $f (@to_index) {
			if (-s $f) {
				unlink "$f.gz", "$f.gz.fai", "$f.gz.gzi";
				$self->v("finishing interrupted index of $f");
				$self->_index_fasta($f);
			}
			elsif (-s "$f.gz" and not -s "$f.gz.fai") {
				$self->v("finishing interrupted faidx of $f.gz");
				system("samtools faidx $f.gz") == 0 or die "samtools faidx $f.gz failed\n";
			}
		}
	}
	if (not $self->opts->{overwrite} and ((-e $file_path and not -z $file_path) or (-e "$file_path.gz" and not -z "$file_path.gz"))) {
		$self->v("skipping $file_path");
		return;
	}

	$self->v("dump $type into $file_path");
    # Dump to <name>.part and rename only once complete: the skip test above treats any
    # non-empty file as done, so a dump killed part-way (e.g. the DB server restarting)
    # must never leave a file under its final name.
    if (not $self->opts->{dry}) {
        $self->_dump_fasta($dba, "$file_path.part", $type);
        for my $f (@to_index) {
            rename("$f.part", $f) or die "cannot rename $f.part -> $f: $!";
        }
    }
    for my $f (@to_index) {
        $self->v("index $f");
        $self->_index_fasta($f) unless ($self->opts->{dry});
    }
	return;
}

sub _dump_fasta {
	my ($self, $dba, $outfile, $type) = @_;
	my $seqio_fh = new Bio::SeqIO(-format => 'fasta', -file => ">$outfile");

	my $slice_adaptor = $dba->get_SliceAdaptor;
	my $slices = $slice_adaptor->fetch_all('toplevel');
	if ($type eq 'dna') {
		# Mask with every repeat analysis whose logic_name matches -repeat_logic. The old
		# hardcoded ['RepeatMask'] matches no analysis in current cores (they have
		# repeatmask_redat, repeatmask_repbase, dust, trf, ...), so the masked files came out
		# as byte-identical unmasked copies without any warning.
		my $logic_re = $self->opts->{repeat_logic};
		my @logic_names = grep { /$logic_re/i }
			map { $_->logic_name } @{ $dba->get_AnalysisAdaptor->fetch_all_by_feature_class('RepeatFeature') };
		if (@logic_names) {
			$self->v("\trepeat masking with: %s", join(', ', @logic_names));
		}
		else {
			$self->v("\tWARNING no repeat analyses match /%s/: dna_sm and dna_rm will be UNMASKED copies", $logic_re);
		}

		my $outfile_sm = $outfile;
		$outfile_sm =~ s/dna\.toplevel/dna_sm.toplevel/;
		my $outfile_rm = $outfile;
		$outfile_rm =~ s/dna\.toplevel/dna_rm.toplevel/;

		my $seqio_unmasked = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile");
		my $seqio_sm = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile_sm");
		my $seqio_rm = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile_rm");

		# Repeats come straight from SQL (no Perl object per feature) as toplevel intervals,
		# also kept as <name>.repeats.toplevel.bed.gz. Each region's sequence is read once and
		# both masked copies are made from it in memory, the same way as
		# Slice::_mask_features (clip to the region; soft = lowercase, hard = N).
		(my $bed = $outfile) =~ s/dna\.toplevel\.fa\.part$/repeats.toplevel.bed/;
		my $iv = @logic_names ? $self->_repeat_intervals($dba, \@logic_names, $bed) : {};
		foreach my $slice (@$slices) {
			my $plain = $slice->seq;
			my $name = $slice->seq_region_name;
			my @iv;
			if ($iv) {
				@iv = unpack('N*', $iv->{$name} // '');
			}
			else {    # repeats that SQL cannot place: let the API project them
				for my $f (map { @{ $slice->get_all_RepeatFeatures($_) } } @logic_names) {
					push @iv, $f->start, $f->end;
				}
			}
			my $len = length $plain;
			for my $m ([$seqio_sm, 1], [$seqio_rm, 0]) {
				my $masked = $plain;
				for (my $k = 0; $k < @iv; $k += 2) {
					my ($fs, $fe) = @iv[$k, $k + 1];
					next if $fe < 1 or $fs > $len;
					$fs = 1 if $fs < 1;
					$fe = $len if $fe > $len;
					my $l = $fe - $fs + 1;
					next if $l < 1;
					if ($m->[1]) { substr($masked, $fs - 1, $l) = lc substr($masked, $fs - 1, $l) }
					else         { substr($masked, $fs - 1, $l, 'N' x $l) }
				}
				$m->[0]->write_seq(Bio::Seq->new(-display_id => $name, -seq => $masked));
			}
			$seqio_unmasked->write_seq(Bio::Seq->new(-display_id => $name, -seq => $plain));
		}
		if ($iv and -e "$bed.gz.part") {
			rename("$bed.gz.part", "$bed.gz") or die "cannot rename $bed.gz.part: $!";
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
                        if ($aa_obj) {
                            $seqio_fh->write_seq($aa_obj);
                        } else {
                            $self->v("no translation for transcript $id");
                        }
				    }
			    }
            }
		}
	}
}

# Repeat features of the given analyses as 1-based toplevel intervals, packed per toplevel
# region name (pack 'N2' per interval), streamed from SQL. Same rows as
# RepeatFeatureAdaptor::fetch_all_by_Slice (inner join to repeat_consensus). A feature is placed
# if its seq_region is toplevel, or is placed whole (one assembly row, end to end) into a
# toplevel region. If any feature cannot be placed that way, returns undef and the caller
# masks with the API instead. Also writes a sorted, bgzipped BED to "$bed.gz.part".
sub _repeat_intervals {
    my ($self, $dba, $logic_names, $bed) = @_;
    my $h = $dba->dbc->sql_helper;
    my $top = q{
        JOIN seq_region_attrib ta ON ta.seq_region_id = t.seq_region_id
        JOIN attrib_type tt ON tt.attrib_type_id = ta.attrib_type_id AND tt.code = 'toplevel'
        JOIN coord_system tcs ON tcs.coord_system_id = t.coord_system_id AND tcs.species_id = ?};
    my %map;    # seq_region_id => [toplevel name, region length, offset, ori, toplevel length, is toplevel]
    for my $r (@{ $h->execute(-SQL => "SELECT t.seq_region_id, t.name, t.length FROM seq_region t $top",
            -PARAMS => [ $dba->species_id ]) }) {
        $map{ $r->[0] } = [ $r->[1], $r->[2], 1, 1, $r->[2], 1 ];
    }
    for my $r (@{ $h->execute(-SQL => "SELECT c.seq_region_id, c.length, t.name, t.length,
                a.asm_start, a.asm_end, a.cmp_start, a.cmp_end, a.ori
            FROM assembly a
            JOIN seq_region c ON c.seq_region_id = a.cmp_seq_region_id
            JOIN seq_region t ON t.seq_region_id = a.asm_seq_region_id $top",
            -PARAMS => [ $dba->species_id ]) }) {
        my ($id, $clen, $tname, $tlen, $as, $ae, $cs, $ce, $ori) = @$r;
        next if $map{$id};
        next unless $cs == 1 and $ce == $clen and $ae - $as == $ce - $cs and abs($ori) == 1 and $ae <= $tlen;
        $map{$id} = [ $tname, $clen, ($ori == 1 ? $as : $ae), $ori, $tlen, 0 ];
    }

    my $dbh = $dba->dbc->db_handle;
    my $ph = join(',', ('?') x @$logic_names);
    my $sth = $dbh->prepare(qq{SELECT r.seq_region_id, r.seq_region_start, r.seq_region_end
        FROM repeat_feature r
        JOIN repeat_consensus rc ON rc.repeat_consensus_id = r.repeat_consensus_id
        JOIN analysis a ON a.analysis_id = r.analysis_id
        WHERE a.logic_name IN ($ph)}, { mysql_use_result => 1 });
    $sth->execute(@$logic_names);
    open(my $out, '>', "$bed.unsorted") or die "cannot write $bed.unsorted: $!";
    my (%iv, $n, $bad);
    while (my ($id, $fs, $fe) = $sth->fetchrow_array) {
        my $m = $map{$id};
        if (not $m or $fs > $fe or (not $m->[5] and ($fs < 1 or $fe > $m->[1]))) { $bad++; next; }
        my ($tname, undef, $off, $ori, $tlen) = @$m;
        my ($ts, $te) = $ori == 1 ? ($off + $fs - 1, $off + $fe - 1) : ($off - $fe + 1, $off - $fs + 1);
        next if $te < 1 or $ts > $tlen;    # outside the region: the API skips these
        $ts = 1 if $ts < 1;
        $te = $tlen if $te > $tlen;
        $iv{$tname} .= pack('N2', $ts, $te);
        print $out "$tname\t", $ts - 1, "\t$te\n";
        $n++;
    }
    $sth->finish;
    close($out) or die "cannot write $bed.unsorted: $!";
    if ($bad) {
        $self->v("\tWARNING %d repeat features are not on a toplevel region or a whole placement; masking via the API", $bad);
        unlink "$bed.unsorted";
        return;
    }
    system("set -o pipefail; LC_ALL=C sort -S 2G -k1,1 -k2,2n -k3,3n $bed.unsorted | bgzip -c > $bed.gz.part") == 0
        or die "cannot sort/bgzip $bed.unsorted\n";
    unlink "$bed.unsorted";
    $self->v("\t%d repeat intervals on %d toplevel regions from SQL", $n // 0, scalar(keys %iv));
    return \%iv;
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
	# A DB over makeblastdb's 3 GB volume limit is written as .00.nhr, .01.nhr, ... plus a .nal
	# (.pal) alias and has no plain .nhr, so test the alias too; otherwise every large genome is
	# rebuilt on every run. A build that was interrupted leaves its .inprogress marker behind.
	my $built = grep { -e "$file_path.$_" } qw/nhr phr nal pal/;
	if (not $self->opts->{overwrite} and $built and not -e "$file_path.inprogress") {
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
		if (open(my $mark, '>', "$file_path.inprogress")) { close($mark); }
		if (system("set -o pipefail; $cmd") == 0) {
			unlink "$file_path.inprogress";
		}
		else {
			$self->v("WARNING makeblastdb failed for $file_path; left $file_path.inprogress so it is rebuilt next run");
		}
	}
	return;
}

sub _index_fasta {
    my ($self, $fa_file) = @_;
    # -@: single-threaded bgzip takes ~55 min per 14.7 GB wheat file, three files per genome.
    my $cmd = "bgzip -f -@ " . ($self->opts->{bgzip_threads} || 8) . " $fa_file";
    $self->v($cmd);
    system($cmd) == 0 or die "$cmd failed\n";
    $cmd = "samtools faidx $fa_file.gz";
    $self->v($cmd);
    system($cmd) == 0 or die "$cmd failed\n";
}

sub _ensure_gff3 {
  my ($self, $dba, $out_file, $target_dir) = @_;

  if(! -d $target_dir) {
    if($self->opts->{dry}) {
      $self->v("\tWould have created directory '%s'", $target_dir);
    }
    else {
      $self->v("\tCreating directory '%s'", $target_dir);
      mkpath($target_dir) or die "Cannot create the directory $target_dir: $!";
    }
  }

	my $file_path = join("/",$target_dir, $out_file);
	if (not $self->opts->{overwrite} and ((-e $file_path and not -z $file_path) or (-e "$file_path.gz" and not -z "$file_path.gz"))) {
		$self->v("skipping $file_path");
		return;
	}

	$self->v("dump gff3 into $file_path");
	if (not $self->opts->{dry}) {
    	my $slice_adaptor = $dba->get_SliceAdaptor;
        $slice_adaptor or die "failed to get SlickAdaptor\n";
    	my $slices = $slice_adaptor->fetch_all('toplevel');
	    open (my $gff, "| sort -k1,1 -k4,4n | bgzip -c > $file_path.gz.part");
        eval {
            my $exporter = ExportView::GFF3ExporterMod->new('debug' => 0);
            $exporter->header($gff, $slice_adaptor->db->get_MetaContainer->get_common_name(), $slice_adaptor->db->get_MetaContainer->get_genebuild());
            $exporter->export_genes_from_slices($gff, @$slices);
        };
        my $err = $@;
        close $gff;
        if ($err) {
            unlink "$file_path.gz.part";
            $self->v("BLAMMO! DUMP FAILED: $err\n");
            die $err;    # let the driver retry this genome rather than silently moving on
        }
        rename("$file_path.gz.part", "$file_path.gz") or die "cannot rename $file_path.gz.part: $!";
        # TBI cannot hold sequences > 512 Mb (most wheat chromosomes); fall back to CSI.
        if (system("tabix -p gff $file_path.gz") != 0) {
            $self->v("tbi index failed, trying csi for $file_path.gz");
            system("tabix -C -p gff $file_path.gz") == 0
                or $self->v("WARNING could not index $file_path.gz (e.g. a feature spanning a circular origin)");
        }
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

=item B<--repeat_logic>

Case-insensitive regex selecting the repeat analyses (by logic_name) used to build
dna_sm and dna_rm. Default C<^(repeatmask|dust|trf)>, i.e. RepeatMasker libraries plus
low-complexity (dust) and tandem repeats (trf); excludes e.g. tallymer and repeatdetector.

=item B<--bgzip_threads>

Threads for each bgzip of a FASTA file (bgzip -@). Default 8.

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

