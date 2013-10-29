package EnsEMBL::Web::Object::Gene;

use strict;
use warnings;
no warnings "uninitialized";

use EnsEMBL::Web::Object;
use EnsEMBL::Web::Proxy::Object;
use EnsEMBL::Web::Proxy::Factory;


@EnsEMBL::Web::Object::Gene::ISA = qw(EnsEMBL::Web::Object);

sub get_slice_object {
  my $self = shift;
  my $slice = $self->Obj->feature_Slice->expand( $self->param('flank5_display'), $self->param('flank3_display') );
  my $T = new EnsEMBL::Web::Proxy::Object( 'Slice', $slice, $self->__data );
  $T->highlight_display( $self->Obj->get_all_Exons );
  return $T;
}

sub get_Slice {
  my( $self, $context, $ori ) = @_;
  my $db  = $self->get_db ;
  my $dba = $self->DBConnection->get_DBAdaptor($db);
  my $slice = $self->Obj->feature_Slice;
  if( $context =~ /(\d+)%/ ) {
    $context = $slice->length * $1 / 100;
  }
  if( $ori && $slice->strand != $ori ) {
    $slice = $slice->invert();
  }
  return $slice->expand( $context, $context );
}

sub type_name         { my $self = shift; return $self->species_defs->translate('Gene'); }
sub gene              { my $self = shift; return $self->Obj;             }
sub stable_id         { my $self = shift; return $self->Obj->stable_id;  }
sub feature_type      { my $self = shift; return $self->Obj->type;       }
sub source            { my $self = shift; return $self->Obj->source;     }
sub version           { my $self = shift; return $self->Obj->version;    }
sub logic_name        { my $self = shift; return $self->Obj->analysis->logic_name; }
sub coord_system      { my $self = shift; return $self->Obj->slice->coord_system->name; }
sub seq_region_type   { my $self = shift; return $self->coord_system;    }
sub seq_region_name   { my $self = shift; return $self->Obj->slice->seq_region_name; }
sub seq_region_start  { my $self = shift; return $self->Obj->start;      }
sub seq_region_end    { my $self = shift; return $self->Obj->end;        }
sub seq_region_strand { my $self = shift; return $self->Obj->strand;     }
sub feature_length    { my $self = shift; return $self->Obj->feature_Slice->length; }

sub get_external_id {
  my( $self, $type ) = @_; 
  my $links = $self->get_database_matches($self->gene);
  my $ext_id;
  foreach my $link (@$links) {
    $ext_id = $link->primary_id if ($link->database eq $type);
  }
  return $ext_id;
}

sub get_database_matches {
  my $self = shift;
  my @DBLINKS;
  eval { @DBLINKS = @{$self->Obj->get_all_DBLinks};};
  return \@DBLINKS  || [];
}

sub get_all_transcripts{
  my $self = shift;
  unless ($self->{'data'}{'_transcripts'}){
    foreach my $transcript (@{$self->gene()->get_all_Transcripts}){
      my $transcriptObj = EnsEMBL::Web::Proxy::Object->new(
        'Transcript', $transcript, $self->__data
      );
      $transcriptObj->gene($self->gene);
      push @{$self->{'data'}{'_transcripts'}} , $transcriptObj;
    }
  }
  return $self->{'data'}{'_transcripts'};
}

sub chromosome {
  my $self = shift;
  return undef if lc($self->coord_system) ne 'chromosome';
  return $self->Obj->slice->seq_region_name;
}

sub display_xref {
  my $self = shift; 
  my $trans_xref = $self->Obj->display_xref();
  return ($trans_xref->display_id, $trans_xref->dbname, $trans_xref->primary_id, $trans_xref->db_display_name, $trans_xref->info_text ) if $trans_xref;
}

sub gene_description {
  my $self = shift;
  my %description_by_type = ( 'bacterial_contaminant' => "Probable bacterial contaminant" );
  return $self->Obj->description() || $description_by_type{ $self->Obj->biotype } || 'No description';
}

sub mod_date {
  my $self = shift;
  my $time = $self->gene()->modified_date;
  return $self->date_format( $time,'%d/%m/%y' ), $self->date_format( $time, '%y/%m/%d' );
}

sub created_date {
  my $self = shift;
  my $time = $self->gene()->created_date;
  return $self->date_format( $time,'%d/%m/%y' ), $self->date_format( $time, '%y/%m/%d' );
}

sub get_db {
  my $self = shift;
  my $db = $self->param('db') || 'core';
  return $db eq 'est' ? 'otherfeatures' : $db;
}


sub gene_type {
  my $self = shift;
  my $db = $self->get_db;
  my $type = '';
  if( $db eq 'core' ){
    $type = $self->logic_name;
    $type ||= $self->db_type;
  } else {
    $type = $self->db_type;
    $type ||= $self->logic_name;
  }
  $type ||= $db;
  if( $type !~ /[A-Z]/ ){ $type = ucfirst($type) } #All lc, so format
  return $type;
}


sub date_format { 
  my( $self, $time, $format ) = @_;
  my( $d,$m,$y) = (localtime($time))[3,4,5];
  my %S = ('d'=>sprintf('%02d',$d),'m'=>sprintf('%02d',$m+1),'y'=>$y+1900);
  (my $res = $format ) =~s/%(\w)/$S{$1}/ge;
  return $res;
}


sub location_string {
  my $self = shift;
  return sprintf( "%s:%s-%s", $self->seq_region_name, $self->seq_region_start, $self->seq_region_end );
}

sub readable_location {
  my $self = shift;
  my $gene = $self->gene;
  return sprintf "%s: %s",
     $self->neat_sr_name( $self->coord_system, $gene->slice->seq_region_name ),
     $self->round_bp( $gene->start );
}

sub get_contig_location {
  my $self    = shift;
  my ($pr_seg) = @{$self->Obj->project('seqlevel')};
  return undef unless $pr_seg;
  return (
    $self->neat_sr_name( $pr_seg->[2]->coord_system->name, $pr_seg->[2]->seq_region_name ),
    $pr_seg->[2]->seq_region_name,
    $pr_seg->[2]->start
  );
}

sub get_alternative_locations {
  my $self = shift;
  my @alt_locs = map { [ $_->slice->seq_region_name, $_->start, $_->end, $_->slice->coord_system->name ] }
     @{$self->Obj->get_all_alt_locations};
  return \@alt_locs;
}

sub get_homology_matches{
  my $self = shift;
  my $homology_source = shift;
  my $homology_description = shift;
  $homology_source = "ENSEMBL_HOMOLOGUES" unless (defined $homology_source);
  $homology_description= "ortholog" unless (defined $homology_description);

# warn("Getting homology matches");
  
  my %homologues = %{$self->fetch_homology_species_hash($homology_source, $homology_description)};
# warn("Got @{[scalar keys %homologues]} for (src=$homology_source desc=$homology_description)");
  return unless keys %homologues;
  my $gene = $self->Obj;
  my $geneid = $gene->stable_id;
  my %homology_list;
  my $adaptor_call = $self->param('gene_adaptor') || 'get_GeneAdaptor';

  # hash to convert descriptions into more readable form

  my %desc_mapping_sr7 = ('ortholog_one2one' => '1 to 1', 'apparent_ortholog_one2one' => '1 to 1 (apparent)', 'ortholog_one2many' => '1 to many', 'between_species_paralog' => 'paralogue (between species)', 'ortholog_many2many' => 'many to many', 'within_species_paralog' => 'paralog (within species)');


  foreach my $displayspp (keys (%homologues)){
# warn("Homologue for: $displayspp");
    ( my $spp = $displayspp ) =~ tr/ /_/;
    my $order_sr7=0;
    foreach my $homology (@{$homologues{$displayspp}}){


      my ($homologue, $homology_desc, $homology_subtype) = @{$homology};
      next unless ($homology_desc =~ /$homology_description/);
      my $homologue_id = $homologue->stable_id;
      my $homology_desc_sr7= $desc_mapping_sr7{$homology_desc};   # mapping to more readable form
      $homology_desc_sr7= "no description" unless (defined $homology_desc_sr7);
      $homology_list{$displayspp}{$homologue_id}{'homology_desc'} = $homology_desc_sr7 ;
      $homology_list{$displayspp}{$homologue_id}{'homology_subtype'} = $homology_subtype ;
      $homology_list{$displayspp}{$homologue_id}{'spp'} = $displayspp ;
      $homology_list{$displayspp}{$homologue_id}{'description'} = $homologue->description ;
      
      $homology_list{$displayspp}{$homologue_id}{'order'} = $order_sr7 ;
      
      if ($self->species_defs->valid_species($spp)){
        my $database_spp = $self->DBConnection->get_databases_species( $spp, 'core') ;
        unless( $database_spp->{'core'} ) {
          warn "NO CORE DB CONNECTION ($spp)";
          next;
        }
        my $geneadaptor_spp = $database_spp->{'core'}->$adaptor_call;
                my $gene_spp = $geneadaptor_spp->fetch_by_stable_id( $homologue->stable_id, 1 );
        unless ( $gene_spp ) {
          warn "Gene @{[$homologue->stable_id]} not in core database for $spp";
          next;
        }
        my $display_xref = $gene_spp->display_xref;
        my $display_id = $display_xref ?  $display_xref->display_id() : 'Novel Ensembl prediction';
        $homology_list{$displayspp}{$homologue_id}{'display_id'} = $display_id;
        $homology_list{$displayspp}{$homologue_id}{'description'} = $gene_spp->description || 'No description';
        $homology_list{$displayspp}{$homologue_id}{'location'}= $gene_spp->feature_Slice->name;
        $database_spp->{'core'}->dbc->disconnect_if_idle();
      }
      $order_sr7++;
      
    }
  }
  return \%homology_list;
}


sub fetch_homology_species_hash {
  my $self = shift;
  my $homology_source = shift;
  my $homology_description = shift;
  
  
  $homology_source = "ENSEMBL_HOMOLOGUES" unless (defined $homology_source);
  $homology_description= "ortholog" unless (defined $homology_description);
  
  my $geneid = $self->stable_id;
  my $databases = $self->database('compara') ;
  my %homologues;

  return {} unless $databases;

  my $member_adaptor = $databases->get_MemberAdaptor;
  my $query_member = $member_adaptor->fetch_by_source_stable_id("ENSEMBLGENE",$geneid);

  return {} unless defined $query_member ;
  my $homology_adaptor = $databases->get_HomologyAdaptor;
  my $homologies_array = $homology_adaptor->fetch_all_by_Member_method_link_type($query_member,$homology_source);
# print STDERR "HOMOLOGIES: ", Data::Dumper::Dumper($homologies_array), "\n";
  
  my $query_taxon = $query_member->taxon;
  my %classification;
  my $idx = 1;
  my $node = $query_taxon;
  while ($node){
    $classification{$node->get_tagvalue('scientific name')} = $idx;
    $node = $node->parent;
    $idx++;


  }
 
 
 foreach my $homology (@{$homologies_array}){
# warn("Looking at homology $homology");
    next unless ($homology->description =~ /$homology_description/);
    foreach my $member_attribute (@{$homology->get_all_Member_Attribute}) {
      my ($member, $attribute) = @{$member_attribute};
      next if ($member->stable_id eq $query_member->stable_id);
     # push (@{$homologues{$member->genome_db->name}}, [ $member, $homology->description, $homology->dnds_ratio ]);

     # warn("Added homologue @{[$member->genome_db->name]}");
      push (@{$homologues{$member->genome_db->name}}, [ $member, $homology->description, $homology->subtype ]);
      
    }
  }

  foreach my $species_name (keys %homologues){
    @{$homologues{$species_name}} = sort {$classification{$a->[2]} <=> $classification{$b->[2]}} @{$homologues{$species_name}};

  }
  
  return \%homologues;
}


sub get_disease_matches{
  my $self = shift;
  my %disease_list;
  my $disease_adaptor;
  return undef unless ($disease_adaptor = $self->database('disease'));
  my %omim_disease = ();
  my @diseases = $disease_adaptor->disease_name_by_ensembl_gene($self->gene());
  foreach my $disease (@diseases){
    next unless $disease;
    my $desc = $disease->name;
    foreach my $loc ($disease->each_Location){
      my $omim_id = $loc->db_id;
      push @{$omim_disease{$desc}}, $omim_id;
    }
  }
  return \%omim_disease ;
}

#----------------------------------------------------------------------

sub get_das_factories {
   my $self = shift;
   return [ $self->__data->{_object}->adaptor()->db()->_each_DASFeatureFactory ];
}

sub get_das_features_by_name {
  my $self = shift;
  my $name  = shift || die( "Need a source name" );
  my $scope = shift || '';
  my $data = $self->__data;     
  my $cache = $self->Obj;
  $cache->{_das_features} ||= {}; # Cache
  my %das_features;
  foreach my $dasfact( @{$self->get_das_factories} ){
    my $type = $dasfact->adaptor->type;
    next if $dasfact->adaptor->type =~ /^ensembl_location/;
    my $name = $dasfact->adaptor->name;
    next unless $name;
    my $dsn = $dasfact->adaptor->dsn;
    my $url = $dasfact->adaptor->url;

# Construct a cache key : SOURCE_URL/TYPE
# Need the type to handle sources that serve multiple types of features

    my $key = $url || $dasfact->adaptor->protocol .'://'.$dasfact->adaptor->domain;
    if ($key =~ m!/das$!) {
	$key .= "/$dsn";
    }
    $key .= "/$type";
    unless( $cache->{_das_features}->{$key} ) { ## No cached values - so grab and store them!!
      my $featref = ($dasfact->fetch_all_by_ID($data->{_object}, $data ))[1];
      $cache->{_das_features}->{$key} = $featref;
    }
    $das_features{$name} = $cache->{_das_features}->{$key};
  }
  return @{ $das_features{$name} || [] };
}

sub get_das_features_by_slice {
  my $self = shift;
  my $name  = shift || die( "Need a source name" );
  my $slice = shift || die( "Need a slice" );
  
  my $cache = $self->Obj;     

  $cache->{_das_features} ||= {}; # Cache
  my %das_features;
    
  foreach my $dasfact( @{$self->get_das_factories} ){
    my $type = $dasfact->adaptor->type;
    next unless $dasfact->adaptor->type =~ /^ensembl_location/;
    my $name = $dasfact->adaptor->name;
    next unless $name;
    my $dsn = $dasfact->adaptor->dsn;
    my $url = $dasfact->adaptor->url;

# Construct a cache key : SOURCE_URL/TYPE
# Need the type to handle sources that serve multiple types of features

    my $key = $url || $dasfact->adaptor->protocol .'://'.$dasfact->adaptor->domain;
    $key .= "/$dsn/$type";

    unless( $cache->{_das_features}->{$key} ) { ## No cached values - so grab and store them!!
      my $featref = ($dasfact->fetch_all_Features( $slice, $type ))[0];
      $cache->{_das_features}->{$key} = $featref;
    }
    $das_features{$name} = $cache->{_das_features}->{$key};
  }

  return @{ $das_features{$name} || [] };
}

sub get_gene_slices {
  my( $self, $master_config, @slice_configs ) = @_;
  foreach my $array ( @slice_configs ) {
    if($array->[1] eq 'normal') {
      my $slice= $self->get_Slice( $array->[2], 1 );
      $self->__data->{'slices'}{ $array->[0] } = [ 'normal', $slice, [], $slice->length ];
    } else {
      $self->__data->{'slices'}{ $array->[0] } = $self->get_munged_slice( $master_config, $array->[2], 1 );
    }
  }
}


# Calls for GeneSNPView

# Valid user selections
sub valids {
  my $self = shift;
  my %valids = ();    ## Now we have to create the snp filter....
  foreach( $self->param() ) {
    $valids{$_} = 1 if $_=~/opt_/ && $self->param( $_ ) eq 'on';
  }
  return \%valids;
}

sub getVariationsOnSlice {
  my( $self, $slice, $subslices, $gene ) = @_;
  my $sliceObj = EnsEMBL::Web::Proxy::Object->new(
        'Slice', $slice, $self->__data
       );

  my ($count_snps, $filtered_snps) = $sliceObj->getFakeMungedVariationFeatures($subslices,$gene);
  $self->__data->{'sample'}{"snp_counts"} = [$count_snps, scalar @$filtered_snps];
  $self->__data->{'SNPS'} = $filtered_snps;
  return ($count_snps, $filtered_snps);
}


sub get_source {
  my $self = shift;
  my $default = shift;

  my $vari_adaptor = $self->Obj->adaptor->db->get_db_adaptor('variation');
  unless ($vari_adaptor) {
    warn "ERROR: Can't get variation adaptor";
    return ();
  }

  if ($default) {
    return  $vari_adaptor->get_VariationAdaptor->get_default_source();
  }
  else {
    return $vari_adaptor->get_VariationAdaptor->get_all_sources();
  }

}


sub store_TransformedTranscripts {
  my( $self ) = @_;

  my $offset = $self->__data->{'slices'}{'transcripts'}->[1]->start -1;
  foreach my $trans_obj ( @{$self->get_all_transcripts} ) {
    my $transcript = $trans_obj->Obj;
    my $raw_coding_start = defined( $transcript->coding_region_start ) ? $transcript->coding_region_start : $transcript->start;
       $raw_coding_start -= $offset;
    my $raw_coding_end   = defined( $transcript->coding_region_end )   ? $transcript->coding_region_end : $transcript->end;
       $raw_coding_end -= $offset;

    my $coding_start = $raw_coding_start + $self->munge_gaps( 'transcripts', $raw_coding_start );
    my $coding_end   = $raw_coding_end   + $self->munge_gaps( 'transcripts', $raw_coding_end );

    my $raw_start = $transcript->start;
    my $raw_end   = $transcript->end  ;
    my @exons = ();
    foreach my $exon (@{$transcript->get_all_Exons()}) {
      my $es = $exon->start - $offset; 
      my $ee = $exon->end   - $offset;
      my $O = $self->munge_gaps( 'transcripts', $es );
      push @exons, [ $es + $O, $ee + $O, $exon ];
    }
    $trans_obj->__data->{'transformed'}{'exons'}        = \@exons;
    $trans_obj->__data->{'transformed'}{'coding_start'} = $coding_start;
    $trans_obj->__data->{'transformed'}{'coding_end'}   = $coding_end;
    $trans_obj->__data->{'transformed'}{'start'}        = $raw_start;
    $trans_obj->__data->{'transformed'}{'end'}          = $raw_end;
  }
}

sub store_TransformedSNPS {
  my $self = shift;
  my $valids = $self->valids;
  foreach my $trans_obj ( @{$self->get_all_transcripts} ) {
    my $T = $trans_obj->stable_id;
    my $snps = {};
    foreach my $S ( @{$self->__data->{'SNPS'}} ) {
      foreach( @{$S->[2]->get_all_TranscriptVariations||[]} ) {
        if( $T eq $_->transcript->stable_id ) {
          $snps->{ $S->[2]->dbID } = $_ if $valids->{'opt_'.lc($_->consequence_type)};
        }
      }
    }
    $trans_obj->__data->{'transformed'}{'snps'} = $snps;
  }
}

sub store_TransformedDomains {
  my( $self, $key ) = @_;
  my %domains;
  my $offset = $self->__data->{'slices'}{'transcripts'}->[1]->start -1;
  foreach my $trans_obj ( @{$self->get_all_transcripts} ) {
    my $transcript = $trans_obj->Obj;
    next unless $transcript->translation;
    foreach my $pf ( @{$transcript->translation->get_all_ProteinFeatures($key)} ) {
## rach entry is an arry containing the actual pfam hit, and mapped start and end co-ordinates
      my @A = ($pf);
      foreach( $transcript->pep2genomic( $pf->start, $pf->end ) ) {
        my $O = $self->munge_gaps( 'transcripts', $_->start - $offset, $_->end - $offset) - $offset;
        push @A, $_->start + $O, $_->end + $O;
      }
      push @{$trans_obj->__data->{'transformed'}{lc($key).'_hits'}}, \@A;
    }
  }
}

sub munge_gaps {
  my( $self, $slice_code, $bp, $bp2  ) = @_;
  my $subslices = $self->__data->{'slices'}{ $slice_code }[2];
  foreach( @$subslices ) {

    if( $bp >= $_->[0] && $bp <= $_->[1] ) {
      return defined($bp2) && ($bp2 < $_->[0] || $bp2 > $_->[1] ) ? undef : $_->[2] ;
    }
  }
  return undef;
}

sub get_munged_slice {
  my $self = shift;
  my $master_config = shift;
  my $slice = $self->get_Slice( @_ );
  my $gene_stable_id = $self->stable_id;

  my $length = $slice->length();
  my $munged  = '0' x $length;
  my $CONTEXT = $self->param( 'context' );
  my $EXTENT  = $CONTEXT eq 'FULL' ? 1000 : $CONTEXT;
  ## first get all the transcripts for a given gene...
  my @ANALYSIS = ( $self->get_db() eq 'core' ? (lc($self->species_defs->AUTHORITY)||'ensembl') : 'otter' );
  @ANALYSIS = qw(ensembl havana ensembl_havana_gene) if $ENV{'ENSEMBL_SPECIES'} eq 'Homo_sapiens';
warn ">>>>> @ANALYSIS <<<<<<";
  my $features = [map { @{ $slice->get_all_Genes($_)||[]} } @ANALYSIS ];
  my @lengths;
  if( $CONTEXT eq 'FULL' ) {
    @lengths = ( $length );
  } else {
    foreach my $gene ( grep { $_->stable_id eq $gene_stable_id } @$features ) {
      foreach my $X ('1') { ## 1 <- exon+flanking, 2 <- exon only
        my $extent = $X == 1 ? $EXTENT : 0;
        foreach my $transcript (@{$gene->get_all_Transcripts()}) {
          foreach my $exon (@{$transcript->get_all_Exons()}) {
            my $START    = $exon->start            - $extent;
            my $EXON_LEN = $exon->end-$exon->start + 1 + 2 * $extent;
            substr( $munged, $START-1, $EXON_LEN ) = $X x $EXON_LEN;
          }
        }
      }
    }
    @lengths = map { length($_) } split /(0+)/, $munged;
  }
  ## @lengths contains the sizes of gaps and exons(+- context)

  $munged = undef;

  my $collapsed_length = 0;
  my $flag = 0;
  my $subslices = [];
  my $pos = 0;
  foreach(@lengths,0) {
    if ($flag=1-$flag) {
      push @$subslices, [ $pos+1, 0, 0 ] ;
      $collapsed_length += $_;
    } else {
      $subslices->[-1][1] = $pos;
    }
    $pos+=$_;
  }

## compute the width of the slice image within the display
  my $PIXEL_WIDTH =
    $self->param('image_width') -
        ( $master_config->get( '_settings', 'label_width' ) || 100 ) -
    3 * ( $master_config->get( '_settings', 'margin' )      ||   5 );

## Work out the best size for the gaps between the "exons"
  my $fake_intron_gap_size = 11;
  my $intron_gaps  = ((@lengths-1)/2);
  if( $intron_gaps * $fake_intron_gap_size > $PIXEL_WIDTH * 0.75 ) {
     $fake_intron_gap_size = int( $PIXEL_WIDTH * 0.75 / $intron_gaps );
  }
## Compute how big this is in base-pairs
  my $exon_pixels  = $PIXEL_WIDTH - $intron_gaps * $fake_intron_gap_size;
  my $scale_factor = $collapsed_length / $exon_pixels;
  my $padding      = int($scale_factor * $fake_intron_gap_size) + 1;
  $collapsed_length += $padding * $intron_gaps;

## Compute offset for each subslice
  my $start = 0;
  foreach(@$subslices) {
    $_->[2] = $start - $_->[0];
    $start += $_->[1]-$_->[0]-1 + $padding;
  }

  return [ 'munged', $slice, $subslices, $collapsed_length+2*$EXTENT ];

}

sub generate_query_hash {
  my $self = shift;
  return {
    'gene' => $self->stable_id,
    'db'   => $self->get_db,
  };

}


# Calls for GeneRegulationView 

sub features {
  my $self = shift;
  return $self->gene->get_all_regulatory_features(1) || [];
}


=head2 vega_projection

 Arg[1]	     : EnsEMBL::Web::Proxy::Object
 Arg[2]	     : Alternative assembly name
 Example     : my $v_slices = $object->ensembl_projection($alt_assembly)
 Description : map an object to an alternative (vega) assembly
 Return type : arrayref

=cut

sub vega_projection {
	my $self = shift;
	my $alt_assembly = shift;
	my $alt_projection = $self->Obj->feature_Slice->project('chromosome', $alt_assembly);
	my @alt_slices = ();
	foreach my $seg (@{ $alt_projection }) {
		my $alt_slice = $seg->to_Slice;
		push @alt_slices, $alt_slice;
	}
	return \@alt_slices;
}


1;

__END__


sub features
  Input:       EnsEMBL::Web::Gene object
  Description:  Returns all the features that regulate this gene
  Output:      Array ref

