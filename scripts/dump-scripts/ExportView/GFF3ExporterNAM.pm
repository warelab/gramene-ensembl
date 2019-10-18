package ExportView::GFF3ExporterNAM;

=pod

=head1 GFF3 Dumper

Library to export ensembl databases in GFF3 format

=head1 Example usage

my $db_adaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new( ...various connection settings... );

my $slice_adaptor = $db_adaptor->get_SliceAdaptor;

ExportView::GFF3Exporter->header(
	\*STDOUT,
	$slice_adaptor->db->get_MetaContainer->get_Species()->common_name(),
	$slice_adaptor->db->get_MetaContainer->get_genebuild()
);

my $slices = $adaptor->fetch_all('toplevel');

my $exporter = ExportView::GFF3Exporter->new('debug => 1);

$exporter->export_genes_from_slices(\*STDOUT, @$slices, $logicname);

=cut

use strict;
use warnings;

our %analysis2logicname;

#convert ensembl style strands -> GFF3 style strands
my $strand_map = {
	'-1'	=> '-',
	'0'		=> '.',
	'1'		=> '+',
};

my $feature_method_map = {
	'dna_align'				=> 'get_all_DnaAlignFeatures',
	'marker'				=> 'get_all_MarkerFeatures',
	'repeat'				=> 'get_all_RepeatFeatures',
	'affy'					=> 'get_all_AffyFeatures',
	'assembly_exception'	=> 'get_all_AssemblyExceptionFeatures',
	'ditag'					=> 'get_all_DitagFeatures',
	'external'				=> 'get_all_ExternalFeatures',
	'oligo'					=> 'get_all_OligoFeatures',
	'qtl'					=> 'get_all_QtlFeatures',
	'simple'				=> 'get_all_SimpleFeatures',
	'protein_align'			=> 'get_all_ProteinAlignFeatures',
};

my $feature_type_map = {
	'dna_align'				=> 'nucleotide_match',
	'marker'				=> 'repeat_region',
	'repeat'				=> 'repeat_region',
	'affy'					=> 'repeat_region',
	'assembly_exception'	=> 'repeat_region',
	'ditag'					=> 'repeat_region',
	'external'				=> 'repeat_region',
	'oligo'					=> 'repeat_region',
	'qtl'					=> 'repeat_region',
	'simple'				=> 'repeat_region',
	'protein_align'			=> 'repeat_region',
};

=pod

=head1 METHODS

=item new

=over

Standard constructor. Can pass in a boolean value to the 'debug' attribute to print out debugging information as you go along.

Otherwise, just sets up the internal counters and data structures necessary to dump the database.

=cut

sub new {

	my $class = shift;
	my %init = @_;

	my $self = bless {}, $class;
	
	return $self->initialize(%init);
}

=pod

=item initialize

Sets all the internal counters back to their starting values. Called automatically by the constructor, must be
explicitly called if you dump more than once from the same exporter object

=cut

sub initialize {
	my $self = shift;
	
	my %init = @_;
	
	my $default = {
		'id_counter'				=> 0,
		'overlapping_id'			=> {},
		'duplicate_transcript_id'	=> {},
		'duplicate_gene_id'			=> {},
		'duplicate_sequence_id'		=> {},
		'debug'						=> $init{'debug'} || 0,
		'source'			=> $init{'source'} || '',
		'export_transcripts'		=>  1,
		'export_introns'			=> defined $init{'export_introns'} ? $init{'export_introns'} : 0,
		'export_exons'				=> defined $init{'export_exons'} ? $init{'export_exons'} : 1,
		'export_cdses'				=> defined $init{'export_cdses'} ? $init{'export_cdses'} : 1,
		'export_utrs'                          => defined $init{'export_utr'} ? $init{'export_utrs'} : 1,
	};
	
	@{$self}{keys %$default} = values %$default;
	
	return $self;
}

=pod

=item header(filehandle, species, genebuild)

prints out the GFF header. Expects to be given the file handle to print to, the species being dumped, and the genebuild used

=cut

sub header {
	my $self = shift;
	my $fh = shift;
	my $species = shift;
	my $genebuild = shift;
	
	#print the header
print $fh qq{##gff-version 3
# generated on @{[scalar localtime]} by $0
# for species $species
# genebuild $genebuild
};
}

=pod

=item id_counter

Mainly used internally. Increments and returns the internal ID counter

=cut

sub id_counter {
	my $self = shift;
	return ++$self->{'id_counter'};
}

sub reset_id_counter{
	my $self = shift;
	$self->{'id_counter'} = 0;
}

sub source{
        my $self = shift;
        return $self->{'source'};
}

=pod

=item export_genes_from_slices(FILEHANDLE, \@slices, $logicname);

given a filehandle and an array of slices, exports all genes/transcripts/additional (as specified) from those slices

=cut

sub export_genes_from_slices {
	my $self = shift;
	my $fh = shift;
	my $slices_ref = shift;
	my $logicname = shift;
	
	my @slices = @$slices_ref;
	
	 select((select($fh), $| = 1)[0]);

	
	my %overlapping_id = ();
	my %duplicate_gene_id = ();
	my %duplicate_transcript_id = ();
	my %duplicate_sequence_id = ();
	
    print STDERR "retrieving ", scalar(@slices), " slices\n" if $self->{'debug'};
    my $slice_counter = 0;
	while (@slices) {

		my $slice = shift @slices;
        print STDERR "@ slice $slice_counter (", $slice->name, ") @ ", scalar(localtime), "\n" if ! $slice_counter++ % 500 && $self->{'debug'};
	
        my $gene_counter = 0;
        my $num_genes = scalar(@{$slice->get_all_Genes($logicname)});

        print STDERR "retrieving $num_genes $logicname genes @ ", scalar(localtime), "\n" if $num_genes > 0 && $self->{'debug'};
       
#print "debug before export_genes_from_slice"; 
        $self->export_genes_from_slice($fh, $slice, $logicname);
        
	}	#end while @slices
	
}

=pod

=item export_genes_from_slice(FILEHANDLE, $slice, $logicname);

given a filehandle and a single slice, exports all genes/transcripts/additional (as specified) from that slice

=cut

sub export_genes_from_slice {

	my $self	= shift;
	my $fh		= shift;
	my $slice	= shift;
	my $logicname   = shift;

	#foreach my $gene ( @{ $slice->get_all_Genes_by_type('avec_rnaseq_avec_sbi') } ) {
	foreach my $gene ( @{ $slice->get_all_Genes($logicname)} ) {
print "Debug inside export_genes_from_slice\n";
		my $analysis_obj = $gene->analysis;
		my $analysis_id = $analysis_obj->dbID;
		
		next if ($logicname && (lc $analysis_obj->logic_name ne lc $logicname) );
		unless($analysis2logicname{$analysis_id} ){	
			$analysis2logicname{$analysis_id} = $analysis_obj->logic_name;
		}
		my $extra_attributes = '';

		my $gene_id = $gene->stable_id;
		my $gid = $gene->dbID;
		my $gene_name = $gene_id || '';
die "Not found geneID for $gid" unless $gene_id;		
		if ($self->{'overlapping_id'}->{$gene_name}++) {
			$gene_id .= '.gene';
		}
		
		if ($self->{'duplicate_gene_id'}->{$gene_name}++) {
			$gene_id .= '.' . $self->{'duplicate_gene_id'}->{$gene_name};
			$extra_attributes =
				";Note=This feature has the same ID as another feature named "
				. $self->escape_attribute($gene_name);
		}

		my $gene_seqid = (split(':', $gene->seqname))[2];
		
		unless ($self->{'duplicate_sequence_id'}->{$gene->seqname}++) {
			print $fh join("\t",
				(map { $self->escape($_) }
				(
					$gene_seqid,
					'assembly', #$gene->source,
					$slice->coord_system->name,
					1,
					$gene->seq_region_length,
					'.',
					'.',
					'.',
				)),
				"ID=" . $self->escape_attribute($gene_seqid)
				. ";Name=" . $self->escape_attribute($gene->seqname),
			), "\n";
		}
		
		my $src = $self->source ||        $analysis2logicname{$analysis_id};

print "DEBUG source=$src\n";

		print $fh join("\t",
			(map { $self->escape($_) }
			(
				$gene_seqid,					#seqid
				$src,         #$gene->source,					#source
				'gene', 						#type		
				$gene->seq_region_start,		#start
				$gene->seq_region_end,			#end
				'.',							#score
			)),
			$strand_map->{$gene->strand},	#strand
			'.',							#phase
			"ID=gene:" . $self->escape_attribute($gene_id)		#attributes
			. ";biotype=" . $gene->biotype
			. ";logic_name=" . $analysis2logicname{$analysis_id}
			. $extra_attributes,
			
		), "\n";
		
		my $transcript_counter = 0;
		my $num_transcripts = scalar(@{$gene->get_all_Transcripts});
		print STDERR "Now retrieving $num_transcripts transcripts @ ", scalar(localtime), "\n" if $num_transcripts > 1 && $self->{'debug'};
		
		$self->export_transcripts(
			$fh,
			$gene,
			$gene_id,
			$gene->get_all_Transcripts,
			$src,
		) ;
	}# end foreach genes

}

=pod

=item export_transcripts(FILEHANDLE, $gene, $gene_id, $transcripts);

Mainly used internally by export_genes_from_slice.

Given a file handle, a gene, the gene's id (whatever it should be, and an arrayref of transcripts), exports those transcripst and specified
sub info to the file handle

=cut

sub export_transcripts {

	my $self		= shift;
	my $fh			= shift;
	my $gene		= shift;
	my $gene_id 	= shift;
	my $transcripts	= shift;
	my $source = shift;

print "DEBUG source=$source\n";
	foreach my $transcript (@$transcripts) {
		my $extra_attributes = '';
	
		my $transcript_id = $transcript->stable_id;
		
		my $transcript_name = $transcript_id || '';
		
		if ($self->{'overlapping_id'}->{$transcript_name}++) {
			$transcript_id .= '.transcript';
		}
		
		if ($self->{'duplicate_transcript_id'}->{$transcript_name}++) {
			$transcript_id .= '.' . $self->{'duplicate_transcript_id'}->{$transcript_name};
			$extra_attributes =
				";Note=This feature has the same ID as another feature named "
					. $self->escape_attribute($transcript_name);
		}
		
		
		my $transcript_seqid = (split(':', $gene->seqname))[2];
		
		print $fh join("\t",
			(map { $self->escape($_) }
			(
				$transcript_seqid,					#seqid
				$source,						#source
				'mRNA',								#type
				$transcript->seq_region_start,		#start
				$transcript->seq_region_end,		#end
				'.',								#score
			)),
			$strand_map->{$transcript->strand}, #strand
			'.',								#phase
			"ID=transcript:" . $self->escape_attribute($transcript_id)				#attributes
			. ";Parent=gene:" . $self->escape_attribute($gene_id)
			. ";biotype=" . $transcript->biotype
			. ";transcript_id=" . $self->escape_attribute($transcript_id)
			. $extra_attributes,
			
		), "\n";

		$self->export_introns(
			$fh,
			$gene,
			$transcript,
			$transcript_id,
			$transcript->get_all_Introns
		) if $self->{'export_introns'};

#$transcript->get_all_five_prime_UTRs

	       $self->export_Futrs(
                        $fh,
                        $gene,
                        $transcript,
                        $transcript_id,
                        $transcript->get_all_five_prime_UTRs,
                        $source,
                ) if $self->{'export_utrs'};

		$self->export_exons(
			$fh,
			$gene,
			$transcript,
			$transcript_id,
			$transcript->get_all_Exons,
			$source,
		) if $self->{'export_exons'};

		
		$self->export_cdses(
			$fh, 
			$gene,
			$transcript,
			$transcript_id,
			$transcript->get_all_Exons,
			$source,
		) if $self->{'export_cdses'};

		$self->export_Tutrs(
                        $fh,
                        $gene,
                        $transcript,
                        $transcript_id,
                        $transcript->get_all_three_prime_UTRs,
			$source,
                ) if $self->{'export_utrs'};
	}

}

=pod

=item export_introns(FILEHANDLE, $gene, $transcript, $transcript_id, $introns);

Mainly used internally by export_transcripts.

Given a file handle, a gene, a transcript, the transcript's id (whatever it should be), and an arrayref of introns, exports those introns to the file handle

=cut

sub export_introns {

	my $self			= shift;
	my $fh				= shift;
	my $gene			= shift;
	my $transcript		= shift;
	my $transcript_id	= shift;
	my $introns			= shift;

	$self->reset_id_counter;	
	foreach my $intron (sort {$a->seq_region_start <=> $b->seq_region_start } @$introns) {

		my $intron_id = 'intron.' . $self->id_counter;
	
		my $intron_seqid = (split(':', $gene->seqname))[2];
	
	#if ($intron->seq_region_end > $intron->seq_region_start) {
		print $fh join("\t",
			(map { $self->escape($_) }
			(
				$intron_seqid,					#seqid
				$gene->source,					#source
				'intron',							#type
				$intron->seq_region_start,		#start
				$intron->seq_region_end,			#end
				'.',							#score
			)),
			$strand_map->{$intron->strand},	#strand
			'.',							#phase
			"Parent=" . $self->escape_attribute($transcript_id)		#attributes
			. ";Name=" . $self->escape_attribute($intron_id),
			
		), "\n";
	}

}

=pod

=item export_exons(FILEHANDLE, $gene, $transcript, $transcript_id, $exons);

Mainly used internally by export_transcripts.

Given a file handle, a gene, a transcript, the transcript's id (whatever it should be), and an arrayref of exons, exports those exons to the file handle

=cut

sub export_exons {
	my $self			= shift;
	my $fh				= shift;
	my $gene			= shift;
	my $transcript		= shift;
	my $transcript_id	= shift;
	my $exons			= shift;
	my $source = shift;

	$self->reset_id_counter;
	foreach my $exon ( sort {$a->seq_region_start <=> $b->seq_region_start } @$exons) {


		my $exon_id = 'exon.' . $self->id_counter;

		my $exon_seqid = (split(':', $gene->seqname))[2];
		
		print $fh join("\t",
		#	(map { $self->escape($_) }
		#	(
				$exon_seqid,					#seqid
				$source,					#source
				'exon',							#type
				$exon->seq_region_start,		#start
				$exon->seq_region_end,			#end
				'.',							#score
		#	)),
			$strand_map->{$exon->strand},	#strand
			'.',							#phase
			"Parent=transcript:$transcript_id"		#attributes
			. ";Name=${transcript_id}.$exon_id"
			. ";ensembl_end_phase=". $exon->end_phase
			. ";ensembl_phase=".$exon->end_phase
			. ";exon_id=${transcript_id}.$exon_id"
			. ";rank=" . $exon->rank($transcript),
			
		), "\n";
	}

}

=pod

=item export_cdses(FILEHANDLE, $gene, $transcript, $transcript_id, $cdses);

Mainly used internally by export_transcripts.

Given a file handle, a gene, a transcript, the transcript's id (whatever it should be), and an arrayref of cdses, exports those cdses to the file handle

=cut

sub export_cdses {
	my $self			= shift;
	my $fh				= shift;
	my $gene			= shift;
	my $transcript		= shift;
	my $transcript_id	= shift;
	my $cdses			= shift;
	my $source = shift;

	$self->reset_id_counter;
	my $translation_id = $transcript->translation->stable_id;

	foreach my $cds ( sort {$a->seq_region_start <=> $b->seq_region_start } @$cdses) {

		my $cds_id = "CDS:$translation_id";
	
		my $cds_seqid = (split(':', $gene->seqname))[2];
		
		if (defined $cds->coding_region_start($transcript) && defined $cds->coding_region_end($transcript)) {
	
			print $fh join("\t",
				(map { $self->escape($_) }
				(
					$cds_seqid, 							#seqid
					$source,							#source
					'CDS',									#type
					$cds->coding_region_start($transcript),	#start
					$cds->coding_region_end($transcript),	#end
					'.',									#score
				)),
				$strand_map->{$cds->strand},			#strand
				$cds->phase >= 0 ? $cds->phase : '.',	#phase
				
				"ID=".$self->escape_attribute($cds_id)
				.";Parent=transcript:" . $self->escape_attribute($transcript_id)				#attributes
				. ";protein_id=" . $self->escape_attribute($translation_id),
				
			), "\n";
		}
	}

}

sub export_Futrs {
        my $self                        = shift;
        my $fh                          = shift;
        my $gene                        = shift;
        my $transcript          = shift;
        my $transcript_id       = shift;
        my $Futrs                       = shift;
	my $source = shift;

	$self->reset_id_counter;	
        foreach my $Futr ( sort {$a->seq_region_start <=> $b->seq_region_start } @$Futrs) {

                my $Futr_id = '5UTR.' . $self->id_counter;

                my $Futr_seqid = (split(':', $gene->seqname))[2];

                if (defined $Futr->seq_region_start && defined $Futr->seq_region_end) {

                        print $fh join("\t",
                                (map { $self->escape($_) }
                                (
                                        $Futr_seqid,                                                     #seqid
                                        $source,                                                        #source
                                        'five_prime_UTR',                                                                  #type
                                        $Futr->seq_region_start, #start
                                        $Futr->seq_region_end,   #end
                                        '.',                                                                    #score
                                )),
                                $strand_map->{$Futr->strand},                    #strand
				'.',                                                    #phase
                                "Parent=transcript:" . $self->escape_attribute($transcript_id),                             #attributes
                               

                        ), "\n";
                }
        }

}

sub export_Tutrs {
        my $self                        = shift;
        my $fh                          = shift;
        my $gene                        = shift;
        my $transcript          = shift;
        my $transcript_id       = shift;
        my $Tutrs                       = shift;
        my $source = shift;
	
	$self->reset_id_counter;
        foreach my $Tutr ( sort {$a->seq_region_start <=> $b->seq_region_start } @$Tutrs) {

                my $Tutr_id = '3UTR.' . $self->id_counter;

                my $Tutr_seqid = (split(':', $gene->seqname))[2];

                if (defined $Tutr->seq_region_start && defined $Tutr->seq_region_end) {

                        print $fh join("\t",
                                (map { $self->escape($_) }
                                (
                                        $Tutr_seqid,                                                     #seqid
                                        $source,                                                        #source
                                        'three_prime_UTR',                                                                  #type
                                        $Tutr->seq_region_start, #start
                                        $Tutr->seq_region_end,   #end
                                        '.',                                                                    #score
                                )),
                                $strand_map->{$Tutr->strand},                    #strand
				'.',                                                    #phase
                                "Parent=transcript:" . $self->escape_attribute($transcript_id),                             #attributes
                              

                        ), "\n";
                }
        }

}


sub export_features_from_slices {

	my $self		= shift;
	my $fh			= shift;
	my $slices		= shift;
	my $features	= shift;

	my $id_counter = 0;
	my %duplicate_sequence_id = ();

    my $slice_counter = 0;
	while (@$slices) {

        print STDERR "@ slice $slice_counter @ ", scalar(localtime), "\n" if ! $slice_counter++ % 500 && $self->{'debug'};
		my $slice = shift @$slices;

		next unless $self->should_process_slice($slice);
print "debug before feature\n";
		foreach my $feature (@$features) {

			my $feature_method = $feature_method_map->{$feature};
		
		print STDERR "# fetches $feature_method\n" if $self->{'debug'};
			#my $fc_counter = 0;
	
			foreach my $dna_feature (@{ $slice->$feature_method() } ) {
			
				#data bug. Skip all features with a seq_region_start of 0.
				next if $dna_feature->seq_region_start == 0;
			
				my $seqname_method = $dna_feature->can('hseqname') ? 'hseqname' : 'seqname';
			
				my $dna_feature_id = $dna_feature->$seqname_method() . '.' . ++$id_counter;
				
				#my $gene_seqid = join(':', (split(':', $gene->seqname))[0,1,2]);
				#my $dna_feature_seqid = join('_', (split(':', $dna_feature->seqname))[0,2]);
				my $dna_feature_seqid = (split(':', $dna_feature->seqname() ))[2];
	
				#print "\tgene: $gene ", $gene->display_id, "(", $gene->strand, ")\n";
				
				#XXX - HACK FOR SCOTT. THIS SHOULD BE REVISITED.
#				next unless defined $dna_feature->analysis;
				
				my @output = ((map { $self->escape($_) }
					(
						$dna_feature_seqid, # . '(' . $dna_feature->dbID . ')', #seqid
						$dna_feature->analysis->logic_name, #source
						$feature_type_map->{$feature},
						$dna_feature->seq_region_start, #start
						$dna_feature->seq_region_end, #end
						#'.', #score
                        $dna_feature->can('score') ? $dna_feature->score : '.', #score
					)),
					$strand_map->{$dna_feature->strand}, #strand
					'.', #phase
					#XXX - commented out for scott...for now
					"ID=" . $self->escape($dna_feature->dbID)
					. ";Name=" . $self->escape($dna_feature_id)
					#'' #XXX placeholder because of scott
					. $self->extra_feature_attributes($dna_feature, $feature),#attributes
					
				);
				
				print join ("\t", @output), "\n";
				
				#last if $fc_counter++ > 10;
			}
		}
	}

}

=pod

=item escape($string)

Mainly used internally. Escapes strings to allow for GFF3 encoding

=cut

sub escape {
	my $self = shift;
	my $string = shift;
	
	return '' unless defined $string;

	#$string =~ s/([\t\r\n])/sprintf("%%%02x",ord($1))/eg;
	$string =~ s/([^a-zA-Z0-9.:^*$@!+_?-|])/sprintf("%%%02x",ord($1))/eg;
	

	return $string;
}

=pod

=item escape_attribute($string)

Mainly used internally. Escapes strings to allow for GFF3 encoding as part of attributes. Different characters are escaped

=cut

sub escape_attribute {
	my $self = shift;
	my $string = shift;
	
	return '' unless defined $string;

	#$string =~ s/([\t\r\n])/sprintf("%%%02x",ord($1))/eg;
	$string =~ s/([,=;\t])/sprintf("%%%02x",ord($1))/eg;

	return $string;
}

=pod

=item extra_features_attribute($feature)

Used for subclasses to display extra attributes for a given feature. This class returns an empty string

=cut

sub extra_feature_attributes {
	my $self = shift;
	my $feature = shift;
	my $feature_type = shift;
	
	if ($feature_type eq 'dna_align') {
		#return sprintf( "Parent=%s;Target=%s %s %s %s",
		return sprintf( "Target=%s %s %s %s",
		       #$feature->dbID,
		       $feature->hseqname,
		       $feature->hstart,
		       $feature->hend,
		       ($feature->hstrand > 0 ? '+' : '-' ))
	}
	
	return '';
}

=pod

=item should_process_slice($slice)

Used for subclasses to determine if a given slice should be marched through. This class always returns true.

=cut

sub should_process_slice {
	return 1;
}

=pod

=back

=cut

1;
