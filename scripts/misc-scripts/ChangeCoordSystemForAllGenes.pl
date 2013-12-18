use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use String::Similarity;
use Getopt::Long;
use Pod::Usage;

my $n = qq{\n};


=pod

=head1 Name

    ChangeCoordSystemForAllGenes - project gene features to a new coordinate system

=head1 Synopsis

    ChangeCoordSystemForAllGenes.pl [OPTIONS]

=head1 Options

    --registry registry file
    --species
    --from_cs  the coordinate system the feature is on
    --to_cs    the coordinate system to project to
    --version         the coordinate system version
    --help 



=cut


    my ($registry, $species, $transform_cs, $from_cs, $version, $help);

GetOptions(

	   "--registry=s" =>  \$registry,
	   "--species=s"  =>  \$species,
	   "--to_cs=s"    =>  \$transform_cs,
	   "--from_cs=s"  =>  \$from_cs,
	   "--version=s"  =>  \$version,
	   "--help"       =>  \$help,
	   );

pod2usage( -verbose=>2 ) if $help;
unless ( $registry && $species && $transform_cs && $from_cs){
    pod2usage ( -verbose => 2);
}

$version ||= 1;

Bio::EnsEMBL::Registry->load_all($registry);

my $dbAd = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

### gene, transcript, exon must ALL be on same seqid - thus you must collect the lists to 
### modify first as once a gene is changed the API won't pick up the transcripts as being
### from that gene

### don't take risks with slices - need to make sure we get ALL genes

my $csa = $dbAd->get_CoordSystemAdaptor();
my $asmad = $dbAd->get_AssemblyMapperAdaptor() or die;
my $sad = $dbAd->get_SliceAdaptor() or die;
my $geneAd = $dbAd->get_GeneAdaptor;
my $transAd = $dbAd->get_TranscriptAdaptor;
my $exonAd = $dbAd->get_ExonAdaptor;
my $DBI = $dbAd->dbc;



### where to move them to
my $new_cs_obj = $csa->fetch_by_name($transform_cs,$version) or die;
#print qq{\nnew cs: }, $new_cs_obj->name();

### grab the gene

#/ fetch_all for exon is wrong - i.e. has the constraint for genes/transcripts that includes biotype!?!
my $genes_ref = $geneAd->fetch_all() or die;
my $trans_ref = $transAd->fetch_all() or die;
my $exons_ref = $exonAd->generic_fetch() or die;
#my $gene = $exonAd->fetch_by_dbID(265) or die;
#my $gene = $geneAd->fetch_by_dbID(2050);



&wrapper($genes_ref,'gene');
&wrapper($trans_ref,'transcript');
&wrapper($exons_ref,'exon');
#&wrapper([$gene],'exon');
#&wrapper([$gene],'transcript');
#&wrapper([$gene],'gene');




sub wrapper {

    my ($list, $type) = @_;

    my $c = 0;
    for my $feat (@{$list}) {
	
	my $f_cs = $feat->slice->coord_system->name;
	next if ( lc $f_cs ne lc $from_cs);

        if (lc $f_cs eq  lc $transform_cs) {
            print qq{\n}.$type.' is already in this cs (dbID: '.$feat->dbID.')';
        }
        else {
            print qq{\n[}.++$c.'] converting '.$type.' cs: (dbID: '.$feat->dbID.').';
            &transform_feature_in_situ($feat,$type);
	
        }
     }

	print "transformed $c $type features\n";
    return;
}

sub transform_feature_in_situ {

    # the adaptors are file scoped
    my ($gene, $type) = @_;

    #y the api does know what transcripts a gene has even though the mapping is from transcript (fK) to gene (pk)
    #y as it creates a transcript array entry
    #y for seq_region name you need to map via slice

    #print qq{\nGENE DETAILS:\n};
    #print qq{\nfeature stable id: }, $gene->stable_id();
    #print qq{\nfeature stable id: }, $gene->dbID();
    #print qq{\nfeature seq region: }, $gene->slice()->seq_region_name();
    #print qq{\nfeature start: }, $gene->start();
    #print qq{\nfeature end: }, $gene->end();
    #print qq{\nfeature strand: }, $gene->strand();
    #print qq{\nold cs type: },$gene->slice()->coord_system()->name();
    #print qq{\n\nold seq region id: }, $sad->get_seq_region_id($gene->slice);

    my $slice_gene_ad = $gene->slice()->adaptor();
    #### precisa comecar a usar versao - pode ter mais de um sistema de assembly no mesmo tempo!?!
    #b nao tem muito sentido - i.e. genes devem ser differentes.

    #print qq{\nold cs },$gene->slice()->coord_system()->name(),qq{\n\n};
    my $old_cs_obj = $gene->slice()->coord_system();
    my $asm_mapper = $asmad->fetch_by_CoordSystems($new_cs_obj, $old_cs_obj);

    #The ordering of the coodinate systems is arbitrary.The
    #The following two statements are equivalent:
    #$mapper = $asma->fetch_by_CoordSystems($cs1,$cs2);
    #$mapper = $asma->fetch_by_CoordSystems($cs2,$cs1);
    #my $asm_mapper1 = $asma1->fetch_by_CoordSystems($new_cs_obj, $old_cs_obj);

    #print ($gene->slice->seq_region_name, $gene->start, $gene->end, $gene->strand, $old_cs_obj);
    #use Data::TreeDraw; draw($asm_mapper->map($gene->slice->seq_region_name, $gene->start, $gene->end,
    #$gene->strand, $old_cs_obj));

    # @ctg_ids = $asm_mapper->list_ids( '13', 1_000_000, 1, $chr_cs );

    ### returns in list context - the strangest thing?!?
    my ($mapper_obj) = $asm_mapper->map(
        $gene->slice->seq_region_name, $gene->start, $gene->end,
        $gene->strand, $old_cs_obj
    );

    #use Data::TreeDraw; draw(\@mapper_obj); 

    #print qq{\n\n}, @mapper_obj;

    #print qq{\n\nnew id: }, $mapper_obj->id;
    #print qq{\nnew start: }, $mapper_obj->start;
    #print qq{\nnew end: }, $mapper_obj->end;
    #print qq{\nnew strand: }, $mapper_obj->strand;
    #print qq{\nnew length: }, $mapper_obj->length;

    my $geneid = $gene->dbID;
    my $new_seq = $mapper_obj->id;
    my $new_start = $mapper_obj->start;
    my $new_end = $mapper_obj->end;
    my $new_strand = $mapper_obj->strand;

    
    print qq{\n\n > UPDATE $type\n > SET seq_region_id = '$new_seq', seq_region_start = '$new_start',
    seq_region_end = '$new_end', seq_region_strand = '$new_strand'\n > WHERE $type.id = '$geneid'; };

    #return;
    #exit;

    my $update = $DBI->prepare(qq{
        UPDATE 
            $type
        SET
            seq_region_id = '$new_seq',
            seq_region_start = '$new_start',
            seq_region_end = '$new_end',
            seq_region_strand = '$new_strand'
        WHERE 
            ${type}_id = '$geneid';
    });

    $update->execute or die;

    return
}



### the below is silly!?! just send in complete lists of things - who cares about relations - they will be
### maintained later

__END__
exit;


for my $gene (@genes) {

    #y we separate these as we don't want to skip transcripts or exons in case they aren't on the same level

    #r/ must get these first as once the gene is changed the api won't find the transcripts as its offspring anymore
    my @trans = @{$gene->get_all_Transcripts};

    if ($gene->slice->coord_system->name eq $transform_cs) {
        print qq{\nGene is already in this cs (dbID: }.$gene->dbID.')';
    }
    else {
        print qq{\n[}.++$g.'] converting gene cs: (dbID: '.$gene->dbID.').';
        &transform_feature_in_situ($gene,'gene');
    }

    print qq{\nGene has }, scalar @trans, qq{ transcripts.};

    for my $tran (@trans) {

        #r/ again must get all exons before modifying the transcript!?!
        my @exons = @{$tran->get_all_Exons};

        if ($tran->slice->coord_system->name eq $transform_cs) {
            print qq{\nTranscript is already in this cs (dbID: }.$tran->dbID.')';
        }
        else {
            print qq{\n[}.++$t.'] converting transcript cs: (dbID: '.$tran->dbID.').';
            &transform_feature_in_situ($tran,'transcript');
        }

        print qq{\nTranscript has }, scalar @exons, qq{ exons.};

        for my $exon (@exons) {

            if ($exon->slice->coord_system->name eq $transform_cs) {
                print qq{\nExon is already in this cs (dbID: }.$exon->dbID.')';
            }
            else {
                print qq{\n[}.++$e.'] converting exon cs: (dbID: '.$exon->dbID.').';
                &transform_feature_in_situ($tran,'exon');
            }


    }

    }



}







#my   my $update_gene_sql = qq(
#       UPDATE gene
#          SET biotype = ?,
#              analysis_id = ?,
#              display_xref_id = ?,
#              status = ?,
#              description = ?,
#              is_current = ?,
#              canonical_transcript_id = ?,
#              canonical_annotation = ?
#        WHERE gene_id = ?
#  );
#
#  my $sth = $self->prepare( $update_gene_sql );
#  $sth->execute();
#  $sth->bind_param( 1, $gene->biotype(),        SQL_VARCHAR );
#  $sth->bind_param( 2, $gene->analysis->dbID(), SQL_INTEGER );
#  $sth->bind_param( 3, $display_xref_id,        SQL_INTEGER );
#  $sth->bind_param( 4, $gene->status(),         SQL_VARCHAR );
#  $sth->bind_param( 5, $gene->description(),    SQL_VARCHAR );
#  $sth->bind_param( 6, $gene->is_current(),     SQL_TINYINT );


### The bind_param_array method is used to bind an array of values to a placeholder embedded in the prepared
### statement which is to be executed with "execute_array

# we will just update it straight off

__END__

print qq{\n\n\nnew seq region id: }, $mapper_obj->id;
print qq{\nstart: }, $mapper_obj->start;
print qq{\nend: }, $mapper_obj->end;
print qq{\nstrand: }, $mapper_obj->strand;


__END__
#print qq{\nfirst $asm_mapper1};
#print qq{\nsecond $asm_mapper2};

#use Data::TreeDraw; draw($asma1);
#use Data::TreeDraw; draw($asm_mapper2);




  # decompose this slice into its symlinked components.
  # this allows us to handle haplotypes and PARs
  my $normal_slice_proj =
    $slice_adaptor->fetch_normalized_slice_projection($self);
  foreach my $segment (@$normal_slice_proj) {
    my $normal_slice = $segment->[2];

    $slice_cs = $normal_slice->coord_system();

    my $asma = $db->get_AssemblyMapperAdaptor();
    my $asm_mapper = $asma->fetch_by_CoordSystems($slice_cs, $cs);

    # perform the mapping between this slice and the requested system
    my @coords;

    if( defined $asm_mapper ) {
     @coords = $asm_mapper->map($normal_slice->seq_region_name(),
				 $normal_slice->start(),
				 $normal_slice->end(),
				 $normal_slice->strand(),
				 $slice_cs);
    } else {
      $coords[0] = Bio::EnsEMBL::Mapper::Gap->new( $normal_slice->start(),
						   $normal_slice->end());
    }


    my $last_rank = 0;
    #construct a projection from the mapping results and return it
    foreach my $coord (@coords) {
      my $coord_start  = $coord->start();
      my $coord_end    = $coord->end();
      my $length       = $coord_end - $coord_start + 1;
#  my $db = $slice_adaptor->db();
#  my $csa = $db->get_CoordSystemAdaptor();
#  my $cs = $csa->fetch_by_name($cs_name, $cs_version);
#  my $slice_cs = $self->coord_system();

__END__

#use Data::TreeDraw; draw($gene); 
my $gene_replacement = $gene->transform('toplevel');
#my $gene_replacement = $gene->transform('supercontig');
#print qq{\n\n$stable\n\n};
use Data::TreeDraw; draw($gene_replacement); 
#$gene_replacement->stable_id($stable);
print qq{\ngene start: }, $gene_replacement->stable_id();
print qq{\ngene seq region: }, $gene_replacement->slice()->seq_region_name();
print qq{\ngene start: }, $gene_replacement->start();
print qq{\ngene end: }, $gene_replacement->end();
print qq{\ngene strand: }, $gene_replacement->strand();


sub project {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;

  throw('Coord_system name argument is required') if(!$cs_name);

  my $slice_adaptor = $self->adaptor();

  if(!$slice_adaptor) {
    warning("Cannot project without attached adaptor.");
    return [];
  }

  if(!$self->coord_system()) {
    warning("Cannot project without attached coord system.");
    return [];
  }


  my $db = $slice_adaptor->db();
  my $csa = $db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_by_name($cs_name, $cs_version);
  my $slice_cs = $self->coord_system();

  if(!$cs) {
    throw("Cannot project to unknown coordinate system " .
          "[$cs_name $cs_version]");
  }

  # no mapping is needed if the requested coord system is the one we are in
  # but we do need to check if some of the slice is outside of defined regions
  if($slice_cs->equals($cs)) {
    return $self->_constrain_to_region();
  }

  my @projection;
  my $current_start = 1;

  # decompose this slice into its symlinked components.
  # this allows us to handle haplotypes and PARs
  my $normal_slice_proj =
    $slice_adaptor->fetch_normalized_slice_projection($self);
  foreach my $segment (@$normal_slice_proj) {
    my $normal_slice = $segment->[2];

    $slice_cs = $normal_slice->coord_system();

    my $asma = $db->get_AssemblyMapperAdaptor();
    my $asm_mapper = $asma->fetch_by_CoordSystems($slice_cs, $cs);

    # perform the mapping between this slice and the requested system
    my @coords;

    if( defined $asm_mapper ) {
     @coords = $asm_mapper->map($normal_slice->seq_region_name(),
				 $normal_slice->start(),
				 $normal_slice->end(),
				 $normal_slice->strand(),
				 $slice_cs);
    } else {
      $coords[0] = Bio::EnsEMBL::Mapper::Gap->new( $normal_slice->start(),
						   $normal_slice->end());
    }


    my $last_rank = 0;
    #construct a projection from the mapping results and return it
    foreach my $coord (@coords) {
      my $coord_start  = $coord->start();
      my $coord_end    = $coord->end();
      my $length       = $coord_end - $coord_start + 1;
      
#      if( $last_rank != $coord->rank){
#	$current_start = 1;
#	print "LAST rank has changed to ".$coord->rank."from $last_rank \n";
#     }
#      $last_rank = $coord->rank;

      #skip gaps
      if($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {

        my $coord_cs     = $coord->coord_system();

        # If the normalised projection just ended up mapping to the
        # same coordinate system we were already in then we should just
        # return the original region.  This can happen for example, if we
        # were on a PAR region on Y which refered to X and a projection to
        # 'toplevel' was requested.
        if($coord_cs->equals($slice_cs)) {
          # trim off regions which are not defined
          return $self->_constrain_to_region();
        }
	#create slices for the mapped-to coord system
        my $slice = $slice_adaptor->fetch_by_seq_region_id(
                                                    $coord->id(),
                                                    $coord_start,
                                                    $coord_end,
                                                    $coord->strand());

	my $current_end = $current_start + $length - 1;

        push @projection, bless([$current_start, $current_end, $slice],
                                "Bio::EnsEMBL::ProjectionSegment");
      }

      $current_start += $length;
    }
  }

  return \@projection;
}


__END__

sub project {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;

  my $slice = $self->{'slice'};

  if(!$slice) {
    warning("Feature cannot be projected without attached slice.");
    return [];
  }


  #get an adaptor from the attached slice because this feature may not yet
  #be stored and may not have its own adaptor
  my $slice_adaptor = $slice->adaptor();

  if(!$slice_adaptor) {
    throw("Cannot project feature because associated slice does not have an " .
          " adaptor");
  }

  my $strand = $self->strand() * $slice->strand();
  #fetch by feature always gives back forward strand slice:
  $slice = $slice_adaptor->fetch_by_Feature($self);
  $slice = $slice->invert if($strand == -1);
  return $slice->project($cs_name, $cs_version);
}
  my $db = $slice_adaptor->db();
  my $csa = $db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_by_name($cs_name, $cs_version);
  my $slice_cs = $self->coord_system();




sub _transform_between_Slices {
  my ($self, $to_slice) = @_;

  my $from_slice = $self->contig;

  $self->throw("New contig [$to_slice] is not a Bio::EnsEMBL::Slice")
   unless ($to_slice->isa("Bio::EnsEMBL::Slice") or $to_slice->isa("Bio::EnsEMBL::LRGSlice") );

  if ((my $c1 = $from_slice->chr_name) ne (my $c2 = $to_slice->chr_name)) {
    $self->warn("Can't transform between chromosomes: $c1 and $c2");
    return;
  }

  my($start, $end, $strand);

  #first convert to assembly coords
  if($from_slice->strand == 1) {
    $start  = $from_slice->chr_start + $self->start - 1;
    $end    = $from_slice->chr_start + $self->end   - 1;
    $strand = $self->strand;
  } else {
    $start  = $from_slice->chr_end - $self->end   + 1;
    $end    = $from_slice->chr_end - $self->start + 1;
    $strand = $self->strand;
  }

  #now convert to the other slice's coords 
  if($to_slice->strand == 1) {
    $self->start ($start - $to_slice->chr_start + 1); 
    $self->end   ($end   - $to_slice->chr_start + 1); 
    $self->strand($strand);
  } else {
    $self->start ($to_slice->chr_end - $end   + 1);
    $self->end   ($to_slice->chr_end - $start + 1);
    $self->strand($strand * -1);
  }

  $self->contig($to_slice);

  return $self;
}

###### the actual mechanics of a CS transfer!?!
#
#  my($start, $end, $strand);
#
#  #first convert to assembly coords
#  if($from_slice->strand == 1) {
#    $start  = $from_slice->chr_start + $self->start - 1;
#    $end    = $from_slice->chr_start + $self->end   - 1;
#    $strand = $self->strand;
#  } else {
#    $start  = $from_slice->chr_end - $self->end   + 1;
#    $end    = $from_slice->chr_end - $self->start + 1;
#    $strand = $self->strand;
#  }
#
#  #now convert to the other slice's coords 
#  if($to_slice->strand == 1) {
#    $self->start ($start - $to_slice->chr_start + 1); 
#    $self->end   ($end   - $to_slice->chr_start + 1); 
#    $self->strand($strand);
#  } else {
#    $self->start ($to_slice->chr_end - $end   + 1);
#    $self->end   ($to_slice->chr_end - $start + 1);
#    $self->strand($strand * -1);
#  }
#

sub update {
  my ($self, $gene) = @_;
  my $update = 0;

  if ( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
    throw("Must update a gene object, not a $gene");
  }

  my $update_gene_sql = qq(
       UPDATE gene
          SET biotype = ?,
              analysis_id = ?,
              display_xref_id = ?,
              status = ?,
              description = ?,
              is_current = ?,
              canonical_transcript_id = ?,
              canonical_annotation = ?
        WHERE gene_id = ?
  );

  my $display_xref = $gene->display_xref();
  my $display_xref_id;

  if ( $display_xref && $display_xref->dbID() ) {
    $display_xref_id = $display_xref->dbID();
  } else {
    $display_xref_id = undef;
  }

  my $sth = $self->prepare( $update_gene_sql );
  $sth->execute();

  $sth->bind_param( 1, $gene->biotype(),        SQL_VARCHAR );
  $sth->bind_param( 2, $gene->analysis->dbID(), SQL_INTEGER );
  $sth->bind_param( 3, $display_xref_id,        SQL_INTEGER );
  $sth->bind_param( 4, $gene->status(),         SQL_VARCHAR );
  $sth->bind_param( 5, $gene->description(),    SQL_VARCHAR );
  $sth->bind_param( 6, $gene->is_current(),     SQL_TINYINT );

  if ( defined( $gene->canonical_transcript() ) ) {
    $sth->bind_param( 7, $gene->canonical_transcript()->dbID(),
      SQL_INTEGER );
  } else {
    $sth->bind_param( 7, 0, SQL_INTEGER );
  }

  $sth->bind_param( 8, $gene->canonical_annotation(), SQL_VARCHAR );
  $sth->bind_param( 9, $gene->dbID(), SQL_INTEGER );

  $sth->execute();

  # maybe should update stable id ???
}

sub transform_isSitu {
  my $self = shift;

  # @_ is what is left of @_ after initial shift
  my $new_gene = $self->SUPER::transform(@_);

  if( exists $self->{'_transcript_array'} ) {
    my @new_transcripts;
    my ( $strand, $slice );
    my $low_start = POSIX::INT_MAX;
    my $hi_end = POSIX::INT_MIN;
    for my $old_transcript ( @{$self->{'_transcript_array'}} ) {
      my $new_transcript = $old_transcript->transform(@_);
      # this can fail if gene transform failed  
      
      return undef unless $new_transcript;

      if( ! defined $new_gene ) {
	if( $new_transcript->start() < $low_start ) {
	  $low_start = $new_transcript->start();
	}
	if( $new_transcript->end() > $hi_end ) {
	  $hi_end = $new_transcript->end();
	}
	$slice = $new_transcript->slice();
	$strand = $new_transcript->strand();
      }
      push( @new_transcripts, $new_transcript );
    }

    if( ! defined $new_gene ) {
      %$new_gene = %$self;
      bless $new_gene, ref( $self );

      $new_gene->start( $low_start );
      $new_gene->end( $hi_end );
      $new_gene->strand( $strand );
      $new_gene->slice( $slice );
    }

    $new_gene->{'_transcript_array'} = \@new_transcripts;
  }
  return $new_gene;
}




use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use String::Similarity;
my $n = qq{\n};

$ARGV[0] or die $n.'give me a gff3 file';
open my $FH, '<', $ARGV[0] or die $n.'could not open file!?!';

$ARGV[1] or die $n.'give me a database name';

my $dbAd = Bio::EnsEMBL::DBSQL::DBAdaptor->new( 
    -host   => 'localhost',
    -user   => 'dsth',
    -dbname => $ARGV[1],
);

my $geneAd = $dbAd->get_GeneAdaptor;

#y dodgy approach for a big file but...
my @lines = <$FH>;

my $dbAdtrans = $dbAd->get_TranscriptAdaptor();
my $bored = $dbAd->get_TranslationAdaptor() or die;

my @ids;
for my $line (@lines) {
    #print $line;
    if ($line =~ /^\S+?\t\S+?\tpolypeptide\t.+ID=(\S+?);.*translation=(\S+?);/) { 

        #print $1.$n; # the ids
        #print $2.$n; # the translation
        push @ids, [$1,$2] 
    }
}

#my @bastards;
#my @fuckwits;
my @different;
my @multiple_transcripts;
my $same = 0;
my @no_gene;

for my $i (0..$#ids) {

    #y temp
    #last if $i == 4;
    #exit if $i == 4;

    my $gene = $ids[$i]->[0]; 
    $gene =~ s/:pep//;
    #print $gene.$n;
    
    #$gene =~ s/.\d+$//;

    my $dbGene = $geneAd->fetch_by_stable_id($gene);

    if (!$dbGene) {
        warn qq{\ncould not get gene $gene};
        push @no_gene, $gene;
        next;
    }

    my @trans = @{$dbGene->get_all_Transcripts()};

    if (scalar @trans == 1) {
   


        my $seq = $trans[0]->spliced_seq();

#print qq{\nseq $seq};

        #y/ let's clarify things - mixing Bio::Perl and API is a headache
        # my $BioSeq = $trans[0]->seq;
        # translated bioseq is still just a bioseq object
        #print $n.'raw BioSeq: '.$BioSeq.$n;
        #print $'.'BioSeq translated: '.$BioSeq->translate.$n;
        #print $'.'BioSeq translated seq: '.$BioSeq->translate->seq.$n;
        #y for other frames!?!
        # $trans->seq->trunc(1,$trans->length)->translate->seq;

        #y let's get the peptide seq directly (translateable seq -> dna, translate -> directly to peptide string

        #/ to get bioseq for translation
        #my $ApiPeptideSeq = $trans[0]->translate();
        #/ if we want the translation without going through translateable_seq need to grab the translation object
        my $TrlObj = $trans[0]->translation();
        my $ApiPeptideSeq = $TrlObj->seq;
        chomp ($ApiPeptideSeq);

        #$seqViaBP =~ s/\*$//;

        #print $n.'seq via api: '.$ApiPeptideSeq;
        
        #y if you want the transcript via the translation id - a re-map to check integrity/non-canonical=-translations
        #my $APItrans = $bored->fetch_by_stable_id($trans[0]->stable_id);

        #/ must apply uc here!?!
        my $GffPeptideSeq = uc($ids[$i]->[1]);
        chomp ($GffPeptideSeq);
        my $GffPepLen = length $GffPeptideSeq;
        #my $blah2 = length $seqViaBP;
        my $ApiPepLen = length $ApiPeptideSeq;
#        push @fuckwits, $gene if ($GffPepLen != $ApiPepLen);

        #y to retive the translation directly
        #print $trans[0]->translation().$n;
       
        #/ twat the uc was applied in situ and not as an  lvalue?!?
        print '================================================================'.$n;
        print 'Checking gene '.$gene.$n.$n;
        print 'seq from gff file (has length: '.$GffPepLen.'): '.$n.$GffPeptideSeq.$n.$n;
        #print 'seq from gff file (has length: '.$GffPepLen.'): '.$n.uc($GffPeptideSeq).$n.$n;
        print 'translation via api (has length: '.$ApiPepLen.'): '.$n.$ApiPeptideSeq.$n;

        my $sim = similarity $GffPeptideSeq, $ApiPeptideSeq;
        print $n.'Algorithms for Approximate String Matching gives score (0 -> totally diff, 1 -> the same): ';

        if ($GffPeptideSeq eq $ApiPeptideSeq) {
            print $n.'They are identical!'.$n;
            $same++;
        }
        else {
            push @different, [$gene,$sim];
        }
    }
    else { 
        print 'has multiple transcripts'.$n; 
        push @multiple_transcripts, $gene;
    }
}

print $n.$n.'#######################################################################';
print $n.'there were '.scalar @multiple_transcripts.' genes with more than one transcript - if they exist'
  .' then we need to do this on a per transcript basis?!?';
print $n.'identical: '.$same;
print $n.'there were '.scalar @different.' api OTF translations different to Gff polypeptides';

print $n.$n.'could not fetch: '.scalar @no_gene.' genes: '.@no_gene;
print qq{\n\n};


use Data::TreeDraw; draw(\@different); 

