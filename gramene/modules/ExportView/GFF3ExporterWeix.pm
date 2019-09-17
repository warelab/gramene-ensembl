package ExportView::GFF3ExporterWeix;

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

$exporter->export_genes_from_slices(\*STDOUT, @$slices);

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
		'export_transcripts'		=> defined $init{'export_transcripts'} ? $init{'export_transcripts'} : 1,
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

=pod

=item export_genes_from_slices(FILEHANDLE, @slices);

given a filehandle and an array of slices, exports all genes/transcripts/additional (as specified) from those slices

=cut

sub export_genes_from_slices {
	my $self = shift;
	my $fh = shift;
	my @slices = @_;
	
	
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
        my $num_genes = scalar(@{$slice->get_all_Genes});

        print STDERR "retrieving $num_genes genes @ ", scalar(localtime), "\n" if $num_genes > 0 && $self->{'debug'};
        
        $self->export_genes_from_slice($fh, $slice);
        
	}	#end while @slices
	
}

=pod

=item export_genes_from_slice(FILEHANDLE, $slice);

given a filehandle and a single slice, exports all genes/transcripts/additional (as specified) from that slice

=cut

sub export_genes_from_slice {

	my $self	= shift;
	my $fh		= shift;
	my $slice	= shift;

	#foreach my $gene ( @{ $slice->get_all_Genes_by_type('avec_rnaseq_avec_sbi') } ) {
	foreach my $gene ( @{ $slice->get_all_Genes} ) {

		my $analysis_obj = $gene->analysis;
		my $analysis_id = $analysis_obj->dbID;
		
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
					$gene->source,
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
		
		print $fh join("\t",
			(map { $self->escape($_) }
			(
				$gene_seqid,					#seqid
				$analysis2logicname{$analysis_id},         #$gene->source,					#source
				'gene', 						#type		
				$gene->seq_region_start,		#start
				$gene->seq_region_end,			#end
				'.',							#score
			)),
			$strand_map->{$gene->strand},	#strand
			'.',							#phase
			"ID=" . $self->escape_attribute($gene_id)		#attributes
			. ";Name=" . $self->escape_attribute($gene_name)
			. ";biotype=" . $gene->biotype
			. $extra_attributes,
			
		), "\n";
		
		my $transcript_counter = 0;
		my $num_transcripts = scalar(@{$gene->get_all_Transcripts});
		print STDERR "retrieving $num_transcripts transcripts @ ", scalar(localtime), "\n" if $num_transcripts > 1 && $self->{'debug'};
		
		$self->export_transcripts(
			$fh,
			$gene,
			$gene_id,
			$gene->get_all_Transcripts,
			$analysis2logicname{$analysis_id},
		) if $self->{'export_transcripts'};
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
			"ID=" . $self->escape_attribute($transcript_id)				#attributes
			. ";Parent=" . $self->escape_attribute($gene_id)
			. ";Name=" . $self->escape_attribute($transcript_name)
			. ";biotype=" . $transcript->biotype
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
	
	foreach my $intron (@$introns) {

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

	foreach my $exon (@$exons) {


		my $exon_id = 'exon.' . $self->id_counter;

		my $exon_seqid = (split(':', $gene->seqname))[2];
		
		print $fh join("\t",
			(map { $self->escape($_) }
			(
				$exon_seqid,					#seqid
				$source,					#source
				'exon',							#type
				$exon->seq_region_start,		#start
				$exon->seq_region_end,			#end
				'.',							#score
			)),
			$strand_map->{$exon->strand},	#strand
			'.',							#phase
			"Parent=" . $self->escape_attribute($transcript_id)		#attributes
			. ";Name=" . $self->escape_attribute($exon_id),
			
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

	foreach my $cds (@$cdses) {

		my $cds_id = 'CDS.' . $self->id_counter;
	
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
				
				"Parent=" . $self->escape_attribute($transcript_id)				#attributes
				. ";Name=" . $self->escape_attribute($cds_id),
				
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

        foreach my $Futr (@$Futrs) {

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
                                "Parent=" . $self->escape_attribute($transcript_id)                             #attributes
                                . ";Name=" . $self->escape_attribute($Futr_id),

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

        foreach my $Tutr (@$Tutrs) {

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
                                "Parent=" . $self->escape_attribute($transcript_id)                             #attributes
                                . ";Name=" . $self->escape_attribute($Tutr_id),

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

=item export_pairwise_wga_from_slice() 
 
export whole genome alignment into gff3 format

=cut

sub export_pairwise_wga_from_slice {

        my $self        = shift;
        my $fh          = shift;
        my $slice       = shift;
        my $gabAdaptor  = shift;
        my $mlss        = shift;


        my $seqid = $slice->seq_region_name;
        my $source = $mlss->method->type;
        my $type = $feature_type_map->{'dna_align'};
        my( $score, $phase) = qw(. .);

        foreach my $gab ( @{ $gabAdaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice) } ) {


            my($start, $end, $strand, $attributes);

            $attributes->{GroupID} = $gab->group_id;

            for my $ga( @{$gab->genomic_align_array} ){

                my $ga_species = lc $ga->dnafrag->genome_db->name;
                $ga_species =~ s/\s+/_/g;
                #my $gffref_species = lc $slice->adaptor->db->get_MetaContainer->get_short_name;
                my $gffref_species = lc $slice->adaptor->db->get_MetaContainer->get_production_name;
                $gffref_species =~ s/\s+/_/g;

                #print "In ga $ga_species ? $gffref_species\n";
                if( $ga_species eq $gffref_species ){
                    #this is the reference species in GFF3


                    $start = $ga->dnafrag_start;
                    $end = $ga->dnafrag_end;
                    $strand = $strand_map->{$ga->dnafrag_strand};
                    $attributes->{ID} = $ga->dbID;
                    $attributes->{Name} = $mlss->name . "-". $attributes->{ID};

                }else{
                    #this is the target genome in GFF3

                    $attributes->{Target} = sprintf( "%s %s %s %s",
                       $ga->dnafrag->name,
                       $ga->dnafrag_start,
                       $ga->dnafrag_end,
                       $strand_map->{$ga->dnafrag_strand});
                    #print "target =  $attributes->{Target}\n";

                }
            }
	    
	    my $attrib_string = join ';', (map{ "$_=$attributes->{$_}"} keys %{$attributes} );
	    print $fh join("\t",
			   $seqid,					#seqid
			   $source,         #$gene->source,					#source
			   $type, 						#type		
			   $start,		#start
			   $end,			#end
			   '.',							#score
			   $strand,	#strand
			   '.',							#phase
			   $attrib_string,
			   "\n");
	    
	}
}

sub wga_header {
        my $self = shift;
        my $fh = shift;
        my $species = shift;
        my $method = shift;
        my $genome_info = shift;

        unless( $genome_info ){ 

        #print the header
            print $fh qq{##gff-version 3
## generated on @{[scalar localtime]} by $0
## for species $species
## method $method
                         };

        }else{

            my $ref_species = $genome_info->{ref};
            my $nonref_species = $genome_info->{nonref};
            my $dbname = $genome_info->{dbname};

            my $ref_dnafrag = join "\n## ",
            (map {join ' ',
                  ('sequence-region-reference', $_->name, '1', $_->length ) }                             @{$genome_info->{genomes}->{$ref_species}->{dnafrag}}
             );

            my $nonref_dnafrag = join "\n## ",
            (map {join ' ',
                  ('sequence-region-target', $_->name, '1', $_->length ) }                                @{$genome_info->{genomes}->{$nonref_species}->{dnafrag}}
             );

            print $fh qq{##gff-version 3
## generated on @{[scalar localtime]} by $0
## method $method
## $dbname
## reference species = $ref_species
## target species = $nonref_species
## $ref_dnafrag
## $nonref_dnafrag
                         };

        }

	print $fh "\n";
}

1;

__END__

