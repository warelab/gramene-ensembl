###############################################################################
#
#   Name:           ExportView::GFFExporter;
#
#   Description:    Object for dumping data in GFF
#
#   History:        2001-08-07  jws - original version
#
###############################################################################
package ExportView::GFFExporter;

use strict;

use POSIX qw(strftime);

use ExportView::DBAdaptor;
use ExportView::Helper;
use ExportView::Out;

sub new {
    my ($class,$cgi,$databases,$seqname) = @_;
    my $self = {};
    bless $self, $class;
    $self->cgi($cgi);
    $self->{databases} = $databases;
    $self->seqname($seqname) if $seqname;
    return $self;
}

#######################################################
# export_data - main subroutine called from ExportView
#######################################################
sub export_data {
    my $self = shift;
    my $seqout;
    my $cgi = $self->cgi;
#{
#my @p=$cgi->param;
#for my $p(@p) {
#    my @v=$cgi->param($p);
#    print STDERR "$p=",join(",",@v),"\n";
#}
#}
    my $dbadaptor = ExportView::DBAdaptor->new($cgi,$self->{databases});
    my $helper = ExportView::Helper->new;

    my $output_format = $cgi->param('out') || 'text';
    my $delim = $cgi->param('gff_format') || 'gff';

    unless ($cgi->param('gff_similarity')   || 
            $cgi->param('gff_repeat')       ||
            $cgi->param('gff_genscan')   ||
            $cgi->param('gff_gene')         ||
            $cgi->param('gff_external')     ||
            $cgi->param('gff_marker')       ||
            $cgi->param('gff_variation') 
            ){
      die "Please select at least one type of feature for export in the 'Select export options' section.";
    }

    my %delim_lookup = (
                        gff =>  "\t",
                        tab =>  "\t",
                        csv =>  ",",
                    );
    
    my $glue = $delim_lookup{$delim} || "\t";

    my $vc;
    eval { $vc = $dbadaptor->fetch_vc; };
    $@ && die $@;
    die $helper->missing_field_error if !$vc;

    ##################################################
    # Export requested features in appropriate format
    ##################################################
    my @out;
    tie @out, 'ExportView::Out' if $output_format eq 'direct';

    push @out, 
             "##gff-version 2"
            ,"##Type DNA"
            ,"##date ".strftime("%Y-%m-%d",localtime),
            ,join(" ",'##sequence-region',$vc->seq_region_name,1,$vc->length)  
     ;

    
    my @common_fields = qw( seqname
                            source
                            feature
                            start
                            end
                            score
                            strand
                            frame
                         );
                         
    # build list of additional fields
    my @other_fields;
    if ($cgi->param('gff_similarity') || $cgi->param('gff_repeat') ){
        push @other_fields, qw(hid hstart hend);
    }
    if ($cgi->param('gff_marker')){
        push @other_fields, qw(marker);
    }
    if ($cgi->param('gff_genscan')){
        push @other_fields, qw(genscan);
    }
    if ($cgi->param('gff_gene')){
        push @other_fields, qw(exon_id transcript_id gene_id);
    }

    # accumulate data
    if ($cgi->param('gff_similarity')){
        my %oksource= map {($_,1)} $cgi->param('gff_similarity_source');
#       print STDERR map { "$_ ok=$oksource{$_}\n" } keys %oksource;
        foreach my $feature (@{$vc->get_all_SimilarityFeatures()}){
            if(%oksource) {
                next unless $oksource{ $feature->analysis->gff_source};
            } else {
                next if $feature->analysis->gff_source eq 'Proprietary';
            }
            my $ref = { 'hid'   => $feature->hseqname,
                        'hstart'=> $feature->hstart,
                        'hend'  => $feature->hend,
                        };
            my $fields = join($glue, @{$self->common_fields($feature)});
            $fields .= $glue . $self->other_fields(\@other_fields, $delim, $ref);
            push @out, $fields;
        }    
    }
    if ($cgi->param('gff_marker')){
        my $first=1;
        foreach my $feature (@{$vc->get_all_MarkerFeatures()}){
            my $syn = $feature->marker->display_MarkerSynonym ;
            my ($name,$source);
            if($syn) {
                $name=$syn->name;
                $source=$syn->source;
            } else {
                print STDERR "marker feature id=",$feature->dbID,", marker id=",$feature->marker->dbID," no synonym\n";
                $name="[".$feature->marker->dbID."]";
                undef $source;
            }
            my $ref = { 'marker' => $name };
            my $fields = join($glue, @{$self->common_fields($feature,$source,'locus')});
            $fields .= $glue . $self->other_fields(\@other_fields, $delim, $ref);
            push @out, $fields;
        }    
    }
    if ($cgi->param('gff_repeat')){
        foreach my $feature (@{$vc->get_all_RepeatFeatures()}){
            my $ref = { 'hid'   => $feature->repeat_consensus()->name(),
                        'hstart'=> $feature->hstart,
                        'hend'  => $feature->hend,
                        };
            my $fields = join($glue, @{$self->common_fields($feature,undef,'repeat_region')});
            $fields .= $glue . $self->other_fields(\@other_fields, $delim, $ref);
            push @out, $fields;
        }    
    }
    if ($cgi->param('gff_external')){
        foreach my $feature (@{$vc->get_all_ExternalFeatures()}){
            my $fields = join($glue, @{$self->common_fields($feature)});
            $fields .= $glue . $self->other_fields(\@other_fields, $delim);
            push @out, $fields;
        }    
    }
    if ($cgi->param('gff_variation')){
        foreach my $feature (@{$vc->get_all_SNPs()}){
            my $fields = join($glue, @{$self->common_fields($feature)});
            $fields .= $glue . $self->other_fields(\@other_fields, $delim);
            push @out, $fields;
        }    
    }
    if ($cgi->param('gff_genscan')){
        foreach my $feature (@{$vc->get_all_PredictionTranscripts()}){
            my $feature_name = $feature->stable_id 
                             || 'prediction_id_'.$feature->dbID;
            foreach my $exon (@{$feature->get_all_Exons()}){
                my $fields = join($glue, @{$self->common_fields($exon,undef,'exon')});
                my $ref = {'genscan' => $feature_name };
                $fields .= $glue . $self->other_fields(\@other_fields, $delim, $ref);
                push @out, $fields;
            }
        }    
    }
    if ($cgi->param('gff_gene')){
        foreach my $feature (@{$vc->get_all_Genes()}){
            my $gene_name = $feature->stable_id 
                             || 'gene_id_'.$feature->dbID;
            foreach my $transcript(@{$feature->get_all_Transcripts()}){
                my $transcript_name = $transcript->stable_id 
                                 || 'transcript_id_'.$transcript->dbID;
                foreach my $exon(@{$transcript->get_all_Exons()}){
                    my $exon_name = $exon->stable_id 
                                     || 'exon_id_'.$exon->dbID;
                    my $ref = { 'exon_id' => $exon_name,
                                'transcript_id' => $transcript_name,
                                'gene_id' => $gene_name,
                                'gene_type' => $self->gene_type($feature),
                                };
                    my $fields = join($glue, @{$self->common_fields($exon,undef,'exon')});
                    $fields .= $glue . $self->other_fields(\@other_fields, $delim, $ref);
                    push @out, $fields;
                }
            }
        }    
    }

    unless ($output_format eq 'direct' or scalar (@out) > 0){
        die ("This query produces no results.  Try exporting a different region.");
    }
    
    # Add the field names to the output, if required
    unless ($delim eq 'gff'){
        unshift @out, join($glue, @common_fields, @other_fields);
    }

    if ($output_format eq 'zip'){
        my ($zip_fh,$zip_url) =  $helper->get_new_zip_filename;
        open (ZIPOUT, ">$zip_fh") or die ("Cannot create new zip file.  Please try again later");
    
        foreach my $line(@out){
            print ZIPOUT "$line\n";
        }
        close ZIPOUT;
        my $status = system(" gzip $zip_fh");
        die("Cannot create zipped export.  Please try again later.") unless $status == 0;
        $helper->send_redirect("$zip_url.gz");
    }
    else { 
        $helper->send_header($output_format);
        
        if ($output_format eq 'html'){
            print "<h3>".$vc->id."</h3>";
            print "<pre>";
        }
        foreach my $line(@out){
            print $line."\n";
        }
        
        print "</pre>" if $output_format eq 'html';
    
        $helper->send_footer($output_format);
    }
    untie @out if $output_format eq 'direct';
}


sub cgi{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'cgi'} = $value;
    }
    return $self->{'cgi'};
}


sub seqname{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'seqname'} = $value;
    }
    return $self->{'seqname'};
}


sub common_fields{
    my ($self, $feature, $default_source_tag, $default_primary_tag) = @_;
    my ($source, $tag, $start, $end, $str,$score,$frame,$name,$strand);

    if ($feature->can('score')) {
        $score = $feature->score();
    }
    $score = '.' unless defined $score;

    if ($feature->can('frame')) {
        $frame = $feature->frame();
    }
    $frame = '.' unless defined $frame;
    
    if ($feature->can('strand')) {
        $strand = $feature->strand();
    }
    if (! $strand) {
        $strand = ".";
    } elsif( $strand == 1 ) {
        $strand = '+';
    } elsif ( $feature->strand == -1 ) {
        $strand = '-';
    }
    
    if (defined $self->seqname) {
        $name=$self->seqname;
    } else {
        if($feature->can('entire_seq') && $feature->entire_seq()) {
          $name = $feature->entire_seq()->name();
        }
        if(!$name && $feature->can('seqname')) {
          $name = $feature->seqname();
          $name =~ s/\s/_/g;
        }

        $name ||= 'SEQ';
    }

    if ($feature->can('source_tag')){
        $source = $feature->source_tag;
        $source =~ s/\s/_/g;
    }

    if(!defined($source)) {
        if ($feature->can('source_tag')){
            $source = $feature->source_tag;
            #warn "$source from source_tag\n";
        }
        if( ! $source && $feature->can('analysis')
                                && $feature->analysis){
            $source = $feature->analysis->gff_source;
            #warn "$source from analysis gff_source\n";
        }
        #warn "source already = ensembl\n" if $source eq 'ensembl';
        #warn "no source\n" unless $source;
        $source ||= $default_source_tag || 'ensembl';
    }
    $source =~ s/\s/_/g;
    $source =~ s/:/_/g;

    if ($feature->can('primary_tag')){
        $tag = $feature->primary_tag;
        $tag =~ s/\s/_/g;
    }
    $tag ||= $default_primary_tag || '.';

    if ($feature->can('start')){
        $start = $feature->start;
    }
    
    if ($feature->can('end')){
        $end = $feature->end;
    }

    my @results = ( $name,
                    $source,
                    $tag,
                    $start,
                    $end,
                    $score,
                    $strand,
                    $frame,
                  );
    return \@results;

}


sub other_fields{
    my ($self, $fieldlist, $delim, $dataref) = @_;
    my %fields;
    @fields{@$fieldlist} = ();    # build hash of fields from fieldlist
     
    if ($dataref){  # if passed some field data, load it into the hash
        %fields = (%fields, %{$dataref});
    }
    
    # Munge the fields together appropriately for the delimiter type
    my $munged;
    if ($delim eq 'tab'){
        $munged = join ("\t", @fields{@$fieldlist});
    }
    elsif ($delim eq 'csv'){
        $munged = join (",", @fields{@$fieldlist});
    }
    elsif ($delim eq 'gff'){
        foreach my $field(@$fieldlist){
            if (defined $fields{$field}){
                $munged .= "$field ".$fields{$field}."; ";
            }
        }
    }
    return $munged;
}

=head2 gene_type

  Arg[1]      : Bio::EnsEMBL::Gene $gene - a gene object
  Example     : my $gene_type = $self->gene_type($gene);
  Description : This method returns the type of a gene; directly if gene->type
                returns something meaningful (like in Vega) or by checking
                known-ness otherwise (for Ensembl)
  Return type : String - the gene type
  Exceptions  : none
  Caller      : internal

=cut

sub gene_type {
    my ($self, $gene) = @_;

    # sanity checks
    return unless $gene->isa('Bio::EnsEMBL::Gene');

    my $type;
    if ($gene->can('biotype')) {
        $type = $gene->biotype;
    } elsif ( $gene->can('type') ) {
        $type = $gene->type;
    } else {
        return;
    }

    # find out known-ness of gene for ensembl genes (since gene->type isn't
    # really useful)
    if ($type eq 'ensembl') {
        $gene->is_known ? (return 'Known') : (return 'Novel');
    } else {
        return $type;
    }
}

1;
