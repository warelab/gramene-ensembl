#!/lab/bin/perl 

=head1 NAME

genes_from_tigr.pl - make ensembl genes from TIGR genes

=cut



##
# There maybe some problem with the phase setting
# but it doesn't seem to matter


#for opening database by species 


use lib map { "/usr/local/ensembl-live/$_" } 
        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
             ensembl-draw/modules ensembl-compara/modules );

use strict;
use warnings;
use Data::Dumper qw(Dumper); # For debug

use Getopt::Long;
use Pod::Usage;
use Date::Calc;

#you may not want all of these:
#use Gramene::Config;
#use Bio::EnsEMBL::DBLoader;
#use EnsWeb;
#use EnsEMBL::DB::Core;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::GoXref;

use DBI;


use lib ("/usr/local/gramene/scripts/ensembl/MSU-fetch/TIGRXMLLIB/");
use TIGR_XML_parser;
use Gene_obj;



=head1 SYNOPSIS

genes_from_tigr.pl  [options] 
 
 Options:
    --help              help message
    --man               full documentation
    --registry_file
    --species           species in EnsEMBL Web conf to use for db
    --v                 add to verbosity
    --analysis          analysis logic_name


=head1 OPTIONS

=over 4

=item B<--registry_file> 

the EnsEMBL db config file

=item B<--species> 

A species in the EnsEMBL web config

=item B<--analysis>

logic name of the analysis to be assigned to genes added. Def tigr_gene

=item B<--v>

verbose

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

list of chromosomes xml files in TIGR data to process


=head1 DESCRIPTION

Loads genes from the TIGR XML into the Gramene-Ensembl
database.  The gene xrefs (DBEntries) and display xref are also set.

If there is a problem during the load, the following SQL may help
clear up the database (use with care);

Or use the script to remove all the TIGR genes and their derivatives 
~/scripts/gramene_ensembl/delete_genes_by_source.sql

Note: The following stmt only deletes those fully stored genes
If only part of a gene is stored, it may not get deleted
Also the xref and object_xref records are not deleted

delete ge, tr, ex, ex_tr, pr, ge_sid, tr_sid, ex_sid, pr_sid
#Select ge_sid.stable_id,tr_sid.stable_id,ex_sid.stable_id,pr_sid.stable_id
from   gene ge, 
       transcript tr, 
       exon ex, 
       exon_transcript ex_tr, 
       translation pr,
       gene_stable_id ge_sid, 
       transcript_stable_id tr_sid,
       exon_stable_id ex_sid,  
       translation_stable_id pr_sid
where  ge.gene_id        = tr.gene_id
and    tr.transcript_id  = ex_tr.transcript_id
and    ex.exon_id        = ex_tr.exon_id
and    tr.transcript_id  = pr.transcript_id
and    ge.gene_id        = ge_sid.gene_id
and    tr.transcript_id  = tr_sid.transcript_id
and    ex.exon_id        = ex_sid.exon_id
and    pr.translation_id = pr_sid.translation_id
and    ge.analysis_id    = ?;

Brute force - removes all genes and products.
delete from gene;        
delete from gene_stable_id;
delete from gene_description;
delete from transcript;  
delete from transcript_stable_id;
delete from translation; 
delete from translation_stable_id;
delete from protein_feature;
delete from exon;
delete from exon_stable_id;
delete from exon_transcript;
delete from xref;
delete from object_xref;
delete from identity_xref;
delete from go_xref;


=cut

my($ensembl_species,$verbose,$analysis,$registry_file);

my %gene_stable_id_list; # Screen for duplicates. E.g. in ITGR3 XML DB
                         # global across chromosomes

{  #Argument Processing
  my $help=0;
  my $man=0;
  
  Getopt::Long::Configure("no_ignore_case");
  GetOptions( "help|?"    => \$help,
	      "man"       => \$man,
	      "species=s" => \$ensembl_species,
	      "v+"        => \$verbose,
	      "analysis=s"=>\ $analysis,
	      "registry_file=s" => \$registry_file,
	      
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  $ensembl_species || ( warn( "Need a --species\n" ) &&
			pod2usage(1) );
  $registry_file || ( warn( "Need a --registry_file\n" ) &&
			pod2usage(1) );
  $analysis ||= 'msu_gene';
  
}

#Connect to core Ensembl db for species:


my $reg = "Bio::EnsEMBL::Registry";	#Use this to get adaptors 
$reg->load_all( $registry_file ) or die "load_all failed";
     
my $gene_adaptor  = $reg->get_adaptor($ensembl_species, 'core', 'Gene')
  or die "can't get Gene adaptor for $ensembl_species";

my $slice_adaptor = $reg->get_adaptor($ensembl_species,'core','Slice') 
  or die "can't get Slice adaptor for $ensembl_species";

my $dbe_adaptor   = $reg->get_adaptor($ensembl_species,'core','DBEntry')
  or die "can't get DBEntry adaptor for $ensembl_species";

my $ana_adaptor   = $reg->get_adaptor($ensembl_species,'core','Analysis')
  or die "can't get Analysis adaptor for $ensembl_species";

     
my $dba=$slice_adaptor->db; #DBAdaptor
my $dbc=$dba->dbc;	#DBConnection
warn "user ".$dbc->username  .", db ".$dbc->dbname."\n";


my $analysis_obj = $ana_adaptor->fetch_by_logic_name($analysis) || 
    ( warn( "No analysis corresponding to $analysis\n" ) && 
      pod2usage(1) );

#my $date = sprintf('%4.4i%2.2i%2.2i',Date::Calc::Today);
my $date = time(); 
#perl buildin func returning the number of non-leap seconds since whatever time the system considers to be the epoch, ensembl code would convert it to real date
my $pwd  = `pwd`;
chomp $pwd;

my $chr;

# Loop through each pseudochromosome

for my $file(@ARGV) {

  unless ($file =~ /\//) {
    #full path not specified.
    $file = "$pwd/$file";
  }

  if( #$file =~ m= [^/0-9]* 0* (\d+) \. [^/]*\z =xms || 
      #$file =~ m= [^/]* (mitochondrion | chloroplast) \. [^/]*\z =xmsi ||
      $file =~ m= (mitochondrion | chloroplast) / chr([CM]) \. [^/]*\z =xmsi){
    # tigr xml file name chr01.xml
    $chr = $2;
    $verbose && warn( "Process $file for chromosome $chr\n");
  }else{
    die "Cannot parse $file for chromosome";
  }
  
  
  # Create TIGR XML Parser, we don't want to add genes
  # from different chromosome to one list
  # so create a new parser for each chr
  
  my $TIGRparser = new TIGR_XML_parser();
  $TIGRparser->capture_genes_from_assembly_xml("$file");
  my @genes = $TIGRparser->get_genes();

  $verbose && warn("get ", scalar @genes , " genes\n");

  # Create slice corresponding to pseudochromosome
  my $slice=$slice_adaptor->fetch_by_region('chromosome',$chr)
      || ( warn( "cannot make slice for chromosome $chr") && next );

  # Loop through each gene locus on pseudochromosome
  for my $gene (@genes) {
    
    $gene->refine_gene_object();
      # TODO: process pseudogenes; Ensembl can handle these fine!
    #next unless $gene->{is_pseudogene};

    my $gene_stable_id = $gene->{pub_locus} ||$gene->{TU_feat_name};
    $verbose && warn("gene stable_id is " . $gene_stable_id . "\n");
    
    unless( $gene_stable_id){
      print STDERR "ERROR! -- No gene stable id can be found\n" .
	$gene->toString() . "\n";
     next;      
    }

      # Fetch and store the gene, now the gene becomes the ensembl gene 
      #instead of the tigr gene object
    my $ensembl_gene = &make_gene( $slice, $gene );
    unless( $ensembl_gene){
      print STDERR "ERROR! -- No ensEMBL gene object was build for gene $gene_stable_id, skip\n";
      next;
    }
    
    #For testing#################
    $gene_adaptor->store($ensembl_gene);
    
    $verbose && warn( "Gene ", $ensembl_gene->stable_id, " stored\n\n" );
    
  }


}


exit;
    
# end of main program

########################## subroutines ######################################

#sub max_id_num {
#my($what,$prefix)=@_;
#my $max=$dbh->selectrow_array("select max(stable_id) from ${what}_stable_id where stable_id like '$prefix%'");
#    print "max_id_num($what,$prefix)-->$max\n" if $verbose;
#    return ( $max=~ /(\d+)/) ? $1 :0;
#}


#----------------------------------------------------------------------
# Creates and returns the Bio::EnsEMBL::Gene object based on the TIGR::TU.
# Calls the &_create_<object> routines for Transcript and DBEntry objs
#


sub make_gene{

  my $slice = shift;   # Bio::EnsEMBL::Slice
  my $tigr_tu = shift; # TIGR::TU element reprsent a gene region
                       # including splicing isoforms

  my $gene_stable_id = $tigr_tu->{pub_locus}   # pub_locus is the
    || $tigr_tu->{TU_feat_name};               # preferred
                                               # stable public name
                                               # for this gene locus
                                               # but sometimes there is no 
  #value for it, in that case, use TIGR internal feat name as substitute

  #We have to have gene_stable_id be proceed
  unless( $gene_stable_id){
    print STDERR "ERROR! -- No gene stable id can be found\n" .
      $tigr_tu->toString() . "\n";
    return;
    
  }
    

  if( $gene_stable_id_list{$gene_stable_id} ++ ){ # Screen for dups

    # The LOCUS is not uniq to the gene but FEAT_NAME is 
    # But LOCUS takes precedence over FEAT_NAME

    if( !$gene_stable_id_list{$tigr_tu->{TU_feat_name}} ++){
      $gene_stable_id = $tigr_tu->{TU_feat_name};
    }else{
      print STDERR "ERROR! -- Duplicate for gene $gene_stable_id - skipping\n" ;
      return;
    }
  }
  
  warn("gene stable id = $gene_stable_id\n");

  my $description = $tigr_tu->get_product_names(); #or com_name;
  my $gene_type   = $tigr_tu->{is_pseudogene} ? 
                    'pseudogene': 
		      format_gene_type( $tigr_tu->{gene_type} );

  # Create Ensembl gene object
  my $gene=Bio::EnsEMBL::Gene->new
      (
       -stable_id   => $gene_stable_id,
       -description => $description,
       -analysis    => $analysis_obj,
       -version     => 1,
       -biotype     => $gene_type,
       -source      => 'tigr',
       -created_date  => $date,
       -modified_date => $date,
       );
  
  # Add XREFs to gene
  # The xref primary id cannot be null
  # make gene TU_feat_name as primary id and pub_locus as display_label
  my $gene_xref = [{
		    db          => 'TIGR_LOCUS',
		    primary_acc          => $tigr_tu->{TU_feat_name} || $gene_stable_id,
		    display_label       => $tigr_tu->{pub_locus} || $gene_stable_id,
		    description => $description    || '',
		   },
		   #{
		   # db          => 'TIGR_FN',
		   # id          => $tigr_tu->{TU_feat_name} || $gene_stable_id,
		   # description => $description    || '',
		   #}
                  ];
  
  foreach my $dbentry( make_dbentries($gene_xref) ){
    $gene->add_DBEntry( $dbentry );
    if( $dbentry->dbname eq 'TIGR_LOCUS' ){ # Use as display_xref
      $gene->display_xref( $dbentry );
    }
  }
  
  # Loop through each transcript.
  my @models = ($tigr_tu);
  push @models, $tigr_tu->get_additional_isoforms();
  $verbose && warn( "Gene ", $gene_stable_id ,
                    " has ", scalar(@models), " Transcripts\n" );

  for my $model (@models) {

    my $transcript = make_transcript($slice, $model, $tigr_tu->{is_pseudogene} );
    next unless $transcript;

    # Attach GO xrefs from TIGR::TU to the EnsEMBL::Translation
    my $translation = $transcript->translation;
    foreach my $goxref( make_goxrefs($tigr_tu) ){
      $translation->add_DBEntry( $goxref );
    }

    # Attach the EnsEMBL::Transcript to the EnsEMBL::Gene 
    $gene->add_Transcript($transcript);

  }

  return $gene;
}


#----------------------------------------------------------------------
# Creates and returns the Bio::EnsEMBL::Transcript object based on 
# the TIGR::TU.
# Calls the &_create_<object> routines for DBEntry objs
#
sub make_transcript {
  
  my $slice = shift; # Bio::EnsEMBL::Slice
  my $model = shift; # TIGRXMLLIB::Gene_obj;
  my $is_pseudogene = shift;

  my $transcript_stable_id = $model->{model_pub_locus} #LOC_Os10g40780.1 LOC_Os10g40780.2
    ||         $model->{Model_feat_name}               # 12010.m21966, 12010.m06845
    ||         $model->{pub_locus}                     # LOC_Os10g40780
    ||         $model->{TU_feat_name};                 # 12010.t03317
  
  
  unless( $transcript_stable_id ){
    print STDERR "ERROR! -- No transcript stable id can be found\n" .
      $model->toString() . "\n";
    return;
  }
    
  my $description = $model->get_product_names();
  #print "transcipt is " . $model->toString;
  my @exons = $model->get_exons() ; #get an ordered list of
                                       # mRNA_exon_obj; the first exon
                                       #of the list corresponds to the
                                       #first exon of the spliced gene
  $verbose && warn( "  Transcript $transcript_stable_id has ",
                    scalar( @exons ) , " Exons\n");

  my @ensembl_exons;
  my ($coding_start_exon,
      $coding_end_exon,
      $coding_start,
      $coding_end,
      $coding_start_exon_index, #will determine phase from this one
      $start_phase);
  
    #Don't set phase until have coding start, which may be defaulted 

    #my $strand = $exons[0]->start <= $exons[0]->end ? 1 : -1;

    #Note: where e_start=e_end from TIGR xml, don't know if + or - strand
    #set these flags when see something unambiguous:
    #and correct 1-basepair exons later if necessary

  my ($some_plus_strand,$some_minus_strand); 
  
  for my $exon_index (0..$#exons) { #use index since need it later
    
    my $e = $exons[$exon_index];
    
    my ($exon5, $exon3) = $e->get_mRNA_exon_end5_end3();
    #my ($exon5, $exon3) = ($e->{end5}, $e->{end3});

$verbose && warn( "EXON5-3: $exon5, $exon3\n");
    my ($exonstart, $exonend) = ($exon5, $exon3);
    my $strand                =    1;
    
    if ( $exonstart > $exonend) {
      ($exonstart,$exonend,$strand)=($exonend,$exonstart,-1) ;
      $some_minus_strand=1;
    } else {
      $some_plus_strand ||= $exonstart < $exonend;
    }
    
    my $exon_stable_id = $e->{feat_name} || 
      $transcript_stable_id.".exon".($exon_index+1);
    
    my $exon = Bio::EnsEMBL::Exon->new
      (
       -start     => $exonstart,
       -end       => $exonend,
       -strand    => $strand,
       -slice     => $slice,
       -stable_id => $exon_stable_id,
       -version   => 1,
       -created_date  => $date,
       -modified_date => $date,
      );
    
    $verbose && warn( "    Exon $exon_stable_id at ",
		      $exon->feature_Slice->name, "\n" );
    
        #warn Dumper( $e );
        #my $phase = ( $e->coding_start - $e->e_start ) * $strand;
        #my $end_phase = ( $phase + $exon->length ) %3;
        #warn( "==> $phase, $end_phase" );

    my($cd5, $cd3) = $e->get_CDS_end5_end3();
    if(defined $cd5 ) {
      if(! defined $coding_start_exon) {
	$coding_start_exon = $exon;
	$coding_start= ( ( $cd5 - $exon5 )* $strand )+ 1;
	$coding_start_exon_index = $exon_index;
      }
      $coding_end_exon = $exon;
      $coding_end = ( ( $cd3 - $exon5 )* $strand )+ 1;
    }
    push @ensembl_exons,$exon;
  }
  
    #All 1bp exons defaulted to plus strand - fix them if they
    # should be minus strand
  if($some_minus_strand && ! $some_plus_strand) { 
    for my $exon (@ensembl_exons) {
      $exon->strand(-1);
    }
  }
  
    #construct transcript (do it now for warning message)
  my $gene_type   = $is_pseudogene ? 
                            'pseudogene' : 
			      format_gene_type($model->{gene_type});
  my $transcript=Bio::EnsEMBL::Transcript->new
    (
     -exons     => \@ensembl_exons,
     -stable_id => $transcript_stable_id,
     -version   => 1,
     -slice     => $slice,
     -biotype   => $gene_type,
     -analysis  => $analysis_obj,
     -created_date  => $date,
     -modified_date => $date,
    );
  
  my $transcript_xref = [{
			  db          => 'TIGR_LOCUS_MODEL',
			  primary_acc          => $model->{Model_feat_name} || $transcript_stable_id,
			  display_label        => $transcript_stable_id || '',
			  description => $description    || '',
			 },
			 #{
			  #db          => 'TIGR_FN',
			  #id          => $model->{Model_feat_name} || $transcript_stable_id,
			  #description => $description    || '',
			# }
                        ];
    # Transcript xrefs...
    foreach my $dbentry( make_dbentries($transcript_xref) ){
      $transcript->add_DBEntry( $dbentry );
      if( $dbentry->dbname eq 'TIGR_LOCUS_MODEL' ){ # Use as display_xref
        $transcript->display_xref( $dbentry );
      }
    }
  
  
  
  if(!defined $coding_start_exon) {
        # Qiaoping Yuan <qyuan@tigr.org>, Nov 25, 2003:
        #   "some linkage problems (missed the link between exon and cds) 
        #     ... use the gene model coordinates as the CDS"
      
    $verbose && warn( "  No CDS for Transcript", $transcript_stable_id, 
                        "\n" );
    my( $model5, $model3 ) = $model->get_model_span(); #(end5, end3)

    if( $model5 && $model3 ){
      $verbose && warn ("model5-3: ( $model5, $model3 )\n");
      for my $exon_i (0..$#ensembl_exons) {
	
	my $this_exon=$ensembl_exons[$exon_i];
	
	if( $model5 >= $this_exon->start && $model5 <= $this_exon->end ) {
	  
	  $coding_start_exon = $this_exon;
	  
	  $coding_start      = $this_exon->strand > 0 ? 
	    $model5 - $this_exon->start + 1 :
	      $this_exon->end - $model5 + 1;
	  
	  $coding_start_exon_index = $exon_i;
	  $verbose && warn ("    start $model5 -> $coding_start\n");
	}
	
	if( $model3 >= $this_exon->start && $model3 <= $this_exon->end ) {
	  
	  $coding_end_exon = $this_exon;
	  
	  $coding_end      = $this_exon->strand>0 ? $model3-$this_exon->start+1
	  : $this_exon->end-$model3+1;
	  $verbose && warn ("    end $model3 -> $coding_end\n");
	}
      }
    }
  }
  
  unless(defined $coding_start_exon) { #default to beginning of first exon 
    $coding_start_exon=$ensembl_exons[0];
    $coding_start=1;
    $coding_start_exon_index=0;
  }
  unless(defined $coding_end_exon) {
    $coding_end_exon=$ensembl_exons[$#ensembl_exons];
    $coding_end=$coding_end_exon->end-$coding_end_exon->start+1;
  }

#    #to set phases even if no coding region.  
#Redundant at the moment:
#    defined($coding_start_exon_index) or $coding_start=1, $coding_start_exon_index=0;

    #Set phases:
    #phase is 0 when coding start is 1
  $start_phase= ( 4 - $coding_start %3 ) %3;
    #--No, for some reason just set start phase of 1st exon to 0
    # 1st exon or 1st coding exon?? -try the latter
  $start_phase=0;
    #Go backwards to set phases
  my $end_phase=$start_phase;
  for(my $j=$coding_start_exon_index-1;$j>=0;--$j) {
    $ensembl_exons[$j]->end_phase($end_phase);
        #now subtract length of that exon mod 3
    $end_phase = (3+$end_phase-( $ensembl_exons[$j]->end-$ensembl_exons[$j]->start+1)%3)%3;
    $ensembl_exons[$j]->phase($end_phase);
  }
    #Now go forwards to set phases
  for(my $k=$coding_start_exon_index;$k<=$#ensembl_exons;$k++) {
    $ensembl_exons[$k]->phase($start_phase);
    $start_phase=($start_phase+$ensembl_exons[$k]->end-$ensembl_exons[$k]->start+1)%3;
    $ensembl_exons[$k]->end_phase($start_phase);
  }

    # Create the translation, and add to transcript
  if(defined($coding_start_exon) && defined($coding_end_exon) ) {
      #yes, this is redundant at present
      
    my $translation=Bio::EnsEMBL::Translation->new
      (
       -start_exon => $coding_start_exon,
       -end_exon   => $coding_end_exon,
       -seq_start  => $coding_start,
       -seq_end    => $coding_end,
       -stable_id  => $transcript_stable_id, # Reuse ID
       -version    => 1,
       -created_date  => $date,
       -modified_date => $date,       
      );
    $transcript->translation($translation);
    
  }

    #NB: transcript xrefs, not translation
    #for my $xref (@{$modelxrefs}) {
    #  my $entry=Bio::EnsEMBL::DBEntry->new(
#            -primary_id=>$xref->{id}
#           ,-dbname=>$xref->{db}
#           ,-release=>1        #I guess
#           ,-display_id=>$xref->{id}
#           ,-description=>$xref->{description}
#       );
#       $transcript->add_DBLink($entry);
#       $transcript->display_xref or 
#           $transcript->display_xref($entry);
#           # better than nothing (it'd better be) 
#       print join ( " ", $transcript->stable_id,$entry->dbname
#                    ,$entry->primary_id,$entry->display_id),"\n";
#    }

    return $transcript;

}


#----------------------------------------------------------------------
# Ensembl DBEntry objects are used for gene/transcript etc xrefs 
#
# There are only cDNA evidence xref in the tigrv4 xml
#
# create entry in xref table
# Before making it, need to check whether the same entry already exist, 
# if does, use the existing one instead of making new one

sub make_dbentries{

  my $xref  = shift;
  my @data  = @{ $xref || [] };
  my @xrefs = ();

  foreach my $entry( @data ){
    push @xrefs, Bio::EnsEMBL::DBEntry->new
        (
         -dbname      => $entry->{db} || '',
         -primary_id  => $entry->{primary_acc} || 'noid',
         -display_id  => $entry->{display_label} || '', # Filed required; Ens v28
         -description => $entry->{description} || '',
         -version     => 1,
         -release     => 1,
         );
  }
  return @xrefs;

  
}

#----------------------------------------------------------------------
# Creates EnsEMBL::GoXref objects. Similar to make_dbentries
#
# create entry in xref/go_xref tables
# Before making it, need to check whether the same entry already exist, 
# if does, use the existing one instead of making new one

sub make_goxrefs{
  my $gene_obj = shift;

  my @xrefs = ();

  my @gene_ontology_objs = $gene_obj->get_gene_ontology_objs();
  
  foreach my $entry( @gene_ontology_objs ){
    my $xref =  Bio::EnsEMBL::GoXref->new
      (
       -dbname      => 'GO',
       -primary_id  => $entry->{go_id} || 'noid', #store gene fails w/o primary id
       -display_id  => $entry->{go_id} || '', # Filed required; Ens v28
       -description => $entry->{go_descript} || '',
       -version     => 1,
       -release     => 1,
      );
    
    my @evidences = $entry->get_evidence();
    my $evidence  = $evidences[0] if (scalar @evidences >= 1);
    
    $xref->add_linkage_type( $evidence->ev_code() || 'IEA' );
    push @xrefs, $xref;
  }
  return @xrefs;
}

#----------------------------------------------------------------------
# DEPRECATED: Now handled in MAIN
#
sub store_gene {
  warn( "DEPRECATED store_gene" );
  return;
    my ($gene,$tuxrefs)=@_;
    return unless @{$gene->get_all_Transcripts};
    #$gene_counter++;
    print "store_gene \n" if $verbose;
    $gene->analysis($analysis_obj);
    for my $xref (@{$tuxrefs}) {
      last; # TODO: Reinstate XREFS! 
        my $entry=Bio::EnsEMBL::DBEntry->new(
             -primary_id=>$xref->{id}
            ,-dbname=>$xref->{db}
            ,-release=>1        #I guess
            ,-display_id=>$xref->{id}
            #,-description=>$xref->{description}
        );

      $gene->add_DBLink($entry);
      $gene->display_xref or 
          $gene->display_xref($entry);
            # better than nothing (it'd better be) 
            # if add it after gene:
        #$dbe_adaptor->store($entry,$gene->dbID,"Gene")
        #         or die ("$_");
        print join ( " ", 'gene',$entry->dbname
                     ,$entry->primary_id,$entry->display_id),"\n";
    }

    #$geneadaptor->store($gene);
    warn join( ", ", @{$gene->get_all_Transcripts} );
    print $gene->stable_id," stored\n" if $verbose;


    #print "ready to transform " if $verbose;
    #$gene->transform;
    #print "-- ok\n" if $verbose;
    #$description ||= 'tigr annotaton';
    # TODO: Handle gene description correctly
    #$dbh->db->do("insert into gene_description (gene_id,description) values(?,?)",{},$gene->dbID,$description);
}

sub format_gene_type{

  my $gene_type = shift;
  
  $gene_type =~ m/ protein [^A-Za-z]* coding /xmsi && return "protein_coding";
  $gene_type =~ m/ (\w+) RNA /xmsi && return lc($1) . 'RNA';
  return $gene_type;
}

__END__

=head1 OUTPUT

=item B<Standard Output>

=head1 NOTES
    
Anything flagged as a pseudogene will be skipped

Non-coding genes are not handled.

A Coding gene without CDS will have its coding region determined by its model's coordinates.
Or defaulted to the whole transcript.

=head1 AUTHOR

   Gramene Project (www.gramene.org)
   Stein Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

