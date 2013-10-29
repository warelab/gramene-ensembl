package Bio::EnsEMBL::GlyphSet::gramene_dnafeature;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_feature;
@ISA = qw(Bio::EnsEMBL::GlyphSet_feature);

sub my_label { 
    my($self) = @_; 
    return $self->{'config'}->get($self->{'extras'}{'row'},'text') ; 
}

sub features {
    my ($self) = @_;

#    print STDERR "gramene_dnafeature: extras is a ",ref($self->{'extras'}),"\n";
#    print STDERR "gramene_dnafeature: extras = {"
#	,join(", " ,
#	     map { "$_=>".$self->{'extras'}->{$_} } keys %{$self->{'extras'}}
#	     )
#	,"}\n";
#    print STDERR "gramene_dnafeature: row=",$self->{'extras'}{'row'},"\n";
    my $feature=$self->{'config'}->get($self->{'extras'}{'row'},'feature') ; 

    my $ff=
      $self->{'container'}->get_all_DnaAlignFeatures($feature,0);
				   #analysis.logic_name  ^
    print STDERR "gramene_dnafeature: ",scalar(@$ff)," $feature hits\n";

    return $ff;
}

#method in GlyphSet_feature.pm is ok (actually better, since includes label)
#sub colour {
#    my ($self,$id)=@_; #don't care about id
#    return $self->{'feature_colour'}; # i.e.  col part of config
#}


sub href {
    my ($self, $id ) = @_;
    return ''; # Kiran -- Removing the link when some one clicks on the feature. it somehow seems to go to the show_alignment.pl 
    ( my $estid = $id ) =~ s/(.*?)\.\d+/$1/ ;
    my $row=$self->{'extras'}{'row'};
    my $Config=$self->{'config'};
    my $linkout=$Config->get($row,'linkout');
    return $self->{'config'}{'ext_url'}->get_url( 'EST', $estid ) 
         if exists $linkout->{'GenBank'} ;
    my $z=  $self->zmenu;
    my @zk = grep { $_ ne 'caption' } keys %$z;
    my @zkng = grep { $_ !~ 'gramene\.org' } @zk;
    return $z->{$zkng[0]} if @zkng;
    return $z->{$zk[0]} if @zk;
    return '' ; 
}

sub zmenu {
    my ($self, $id, $start, $end ) = @_;
    my $row=$self->{'extras'}{'row'};
    my $Config=$self->{'config'};
	  my $has_origin = undef;
	  my %current_linkout=();
	  #$istem deletes the last period and everything after it from the 
	  #id.  So "AU108727.1" becomes "AU108727".  This is important in 
	  #linking out to other sites such as GrainGenes.
	  (my $istem=$id)=~s/\.[^.]*$//;

	  my $linkout=$Config->get($row,'linkout');
	  
	  foreach my $key (keys %$linkout){
	    my $link_url= $linkout->{$key};
	    $link_url=~s/FEATUREIDSTEM/$istem/;
	    $link_url=~s/FEATUREID/$id/;
	    $current_linkout{$key}=$link_url;
	  }

	  # Blast link:
	  my $blastdb=$Config->get($row,'blastdb');
	  my $slice_start= $self->{'container'}->chr_start;
	  my $slice_end= $self->{'container'}->chr_end;
	  my ($chr_start,$chr_end);	#start and end of feature on chrom
	  if( $self->{'container'}->strand >=0) {
	      $chr_start=$slice_start+($start-1);
	      $chr_end=$slice_start+($end-1);
	  } else {
	      $chr_start=$slice_end-($start-1);
	      $chr_end=$slice_end-($end-1);
	  }
#########################
# Somehow, the above Chromosome coordinates of the feature are dependent on the zoom level that the user is in the browser
# this means, at different zoom levels, smaller or larger chromosome segements are selected for blast. May not be a problem for large zoom
# levels but could be a problem for small zoom levels (zooming in compeltely) when the feature spans more than the zoom level.
# Padding added to these chromosome coordinates in the script showalignment.pl may alleviate the situation.. but this needs to be looked into

###########################
	  #$current_linkout{alignment}="http://www.gramene.org/perl/blast/show_alignment.pl?chr="
	  $current_linkout{alignment}="http://blast.wormbase.org/db/searches/gramene_blast/show_alignment.pl?chr="
	      .$self->{'container'}->chr_name
	      ."&chr_start=$chr_start&chr_end=$chr_end"
	      ."&feature_id=$id&blastdb=$blastdb" 
	      ."&species=".$ENV{ENSEMBL_SPECIES}
	     if $blastdb;




#	  print STDERR "$istem  chr:"
#	   ,join(", "
#		  , $self->{'container'}->chr_name
#		  , $self->{'container'}->chr_end
#		  , $self->{'container'}->strand
#		)
#	  ,"  feature: $start, $end","\n";

    return { 'caption' => "$row $id", 
# 		   "$id" => $self->href( $id )
#	    ,'GenBank' => qq(http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=n&form=1&field=Sequence+ID&term=$id) -should be in linkout if valid
	    ,%current_linkout
    };
}
1;

__END__
### Blast link stuff: not obvious how to get contig coords in new API

	#get the information for the show_alignment link
	push(@{$id{$f->id()}->{features}}, $f );
	my ($raw_ctg,$seq_start)=$self->{'container'}
	->raw_contig_position($f->start);
	if($raw_ctg->can('embl_offset')) { #sometimes this $raw_ctg is a 'gapcontig' (?whatever that is) and can't embl_offset
	    my $current_offset=$raw_ctg->embl_offset;
#	    print STDERR "Checking:",$raw_ctg->embl_offset,"\n";
		my (undef,$seq_end)=$VirtualContig->raw_contig_position($f->end);
		#convert the contig coordinates into clone coordinates and the contig id into clone_id
		$id{$f->id()}->{seq_start}=$seq_start+$current_offset if $seq_start+$current_offset< $id{$f->id()}->{seq_start} || not exists( $id{$f->id()}->{seq_start});
		$id{$f->id()}->{seq_end}=$seq_end+$current_offset if $seq_end+$current_offset> $id{$f->id()}->{seq_end} || not exists( $id{$f->id()}->{seq_end});
	    }
	    next unless $raw_ctg->can('id');	#i.e. if 'gapcontig'
	    unless (exists  $id{$f->id()}->{contig}){
	      $id{$f->id()}->{contig}=$raw_ctg->id;
	      $id{$f->id()}->{contig}=~s/\.\d+$//;
	      $id{$f->id()}->{blastdb}=$blastdb if $blastdb;
	    }



	  $current_linkout{alignment}="http://www.gramene.org/perl/blast/show_alignment.pl?clone_id=".$id{$i}->{contig}."&clone_start=".$id{$i}->{seq_start}.
	    "&clone_end=".$id{$i}->{seq_end}."&feature_type=$feature&feature_id=$i&blastdb=$blastdb" if $blastdb;
