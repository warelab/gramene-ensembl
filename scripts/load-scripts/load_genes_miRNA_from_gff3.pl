#!/usr/local/bin/perl -w 

=pod

=head1 NAME

load_genes_withMakerAED_gff3.pl - parse genes out from JGI gff3 file and load them into EnsemblDB
                              use Poplar v2.0 gff3 file as example

=head1 SYNOPSIS

  load_genes_withMakerAED_gff3.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -l|--logic_name       the logic name for the gene track
  --no_insert        Do not make changes to the database. For debug.
  --noncoding	nonCoding genes tRNA, miRNA ...
  -b|--biotype		nonCoding biotyoe such as miRNA, tRNA, lincRNA

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s|--species>
  Use this species entry from the registry file [REQUIRED].

B<-l|--logic_name>
  the logic name for these set of genes, the track name shown on the contigview, default = ensembl


B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

  Parses the miRNA gff3 file and uses the data therein to populate 
  the gene, exon, transcript and transcript_attrib tables 
  in the ensembl DB indicated by the --species and --ensembl_registry.
  It will also load the mirBase IDs into xref, object_to_xref tables

  example of the gff file:
1        maizev4 miRNA_primary_transcript        6412284 6412406 .       +       .       ID=MI0013239;Alias=MI0013239;Name=zma-MIR528a
1        maizev4 miRNA   6412285 6412305 .       +       .       ID=MIMAT0014028;Alias=MIMAT0014028;Name=zma-miR528a-5p;Derives_from=MI0013239
1        maizev4 miRNA   6412386 6412406 .       +       .       ID=MIMAT0015362;Alias=MIMAT0015362;Name=zma-miR528a-3p;Derives_from=MI0013239

http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0013239
external db_name
miRBase
miRBase_trans_name
it OK to omit the Structure attribute and simply load them as single exon genes ( chr1:6412284-6412406), with miRNA transcript_attrib (attrib_type_id=15) values 1-21, 102-122 and external link them to mirBase?

B<The Ensembl Registry>

  The database connection details for both Ensembl and interproscan
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'maize', 
    '-group'   => 'core', 
    '-dbname'  => 'zea_mays_core_30_bac20', );
  ---

B<Restoring the database>

  If the script bombs out half-way through, your database will be in a
  partial loaded state (i.e. in a bit of a mess). Here is some SQL to
  help put it right;

  DELETE FROM gene, transcript, transcript_attrib, exon_transcript, exon,
               analysis, object_xref ox, xref x
  WHERE       gene.gene_id = transcript.gene_id
  AND         transcript.transcript_id=transcript_attrib.transcript_id
  AND         transcript.transcript_id=exon_transcript.transcript_id
  AND         exon_transcript.exon_id=exon.exon_id
  AND         gene.analysis_id=analysis.analysis_id
  AND	      join on ox and x ...
  AND         analysis.logic_name="fgenesh_gene"

  The above could be improved with some left joins so that the
  constraints do not fail for partly-loaded gene models.


Maintained by Sharon Wei <weix@cshl.edu>

=cut

#delete from gene; delete from transcript; delete from translation; delete from exon; delete from exon_transcript; delete from gene_stable_id; delete from transcript_stable_id; delete from translation_stable_id; delete from exon_stable_id;

use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Readonly;
use List::Util qw( first );
use List::MoreUtils;

#Readonly my @NAME_FILEDS => qw(NAME ALIAS ID);
Readonly my @NAME_FILEDS => qw(NAME ID ALIAS);

Readonly my $UTR_REGEX   => qr{UTR}xms;
Readonly my $MIBASE => 'miRBase_gene_name';
Readonly my $MIBASE_TRPT => 'miRBase_trans_name';
Readonly my $noncoding_exon_phase =>  -1;

use vars qw( $BASEDIR );

BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;
use Date::Calc;


use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS $NONCODING $BIOTYPE );

my $date = time(); 

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $logic_name, $no_insert, $non_coding, $biotype );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "logic_name=s"       => \$logic_name,
        "no_insert"          => \$no_insert,
	"noncoding"          => \$non_coding,
	"biotype=s"	     => \$biotype,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args

  my $gff_file = $ARGV[0] || pod2usage("\nNeed the path to a gff file\n");
  $species    || pod2usage("\nNeed a --species\n");
  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  $logic_name ||= 'ensembl';
  
  map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
  } $reg, $gff_file;



  # Put stuff in the database?
  $I= $no_insert ? 0 : 1; 
  $NONCODING = $non_coding ? 1 : 0;
 
  $BIOTYPE = $biotype ? $biotype : undef;

warn ("noncoding flag is $NONCODING\n"); #exit;

  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) &&
                pod2usage() );

  # Create the analysis
  $ANALYSIS = Bio::EnsEMBL::Analysis->new
      ( -logic_name => $logic_name ); #'GeneModel_RiceIndica_BGI' 

  # Create a GFF stream
  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


  print( "Loading genes for $species\n" );
  my $pre_text = "  Target DB: ";
  foreach my $dba( $ENS_DBA ){
    print( $pre_text.$dba->dbc->dbname ."@". $dba->dbc->host .
          ":". $dba->dbc->port ."\n");
  }
}


my $GENES = {};
my %TRPT2GENE = ();
my @priTrpt_lines;
my @miRNA_lines;

while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /\#/ );
  next if ( $line =~ /^\s+/ );
  chomp $line;
  if ( $line =~ /\tmiRNA_primary_transcript\t/i){
      push @priTrpt_lines, $line ;
  }else{
      push @miRNA_lines, $line ;
  }
}

#first pass for priTrpt 
#1        maizev4 miRNA_primary_transcript        6412284 6412406 .       +       .       ID=MI0013239;Alias=MI0013239;Name=zma-MIR528a
for my $data_lines(  \@priTrpt_lines, \@miRNA_lines){
    for my $feature_line (  @{$data_lines} ){
	print "$feature_line\n";  #####
  
	my $line_hash = parse_line($feature_line);
	my $feature   = $line_hash->{FEATURE};
	my $attribs   = $line_hash->{ATTRIBS};
	
	my $gene_name_pre =  get_first_func($attribs, @NAME_FILEDS);
	my $gene_name = $gene_name_pre;
	if( $feature eq 'miRNA'){
	    $gene_name = lc $gene_name_pre;
	    $gene_name =~ s/\-[35]p$//;

	    my $mirna_coord = join '|', map {$line_hash->{$_} } qw(SEQ_NAME START END STRAND);
	    die "failed to calculate relative position on the primary transcript $gene_name for miRNA $gene_name_pre" 
		unless ( calculate_mirna_trpt_pos( $GENES->{$gene_name}, $mirna_coord) );
	    	#ATTRIBS needs deep copy
	    for my $k(keys %{$attribs}){
		warn("attrib copy for $k\n");
		$GENES->{$gene_name}->{ATTRIBS}->{$k} = $attribs->{$k} unless defined $GENES->{$gene_name}->{ATTRIBS}->{$k}; 
	    }
	    
	    next;   
	}
	
	my $gene_id   = $attribs->{ID};
	my $mirbase_id =  $attribs->{ALIAS}; 
    
	$attribs->{DESCRIPTION} ||= "mirBase gene $mirbase_id";
	my $gene_description = $attribs->{DESCRIPTION} || $gene_name;

	warn "Gene $gene_id $gene_name $gene_description\n";#####

	my $gene_hash_key = lc $gene_name;
	unless( $GENES->{$gene_hash_key} ){
	    
	    for my $k ( keys %$line_hash){
		$GENES->{$gene_hash_key}->{$k} = $line_hash->{$k};
	    }
	    $GENES->{$gene_hash_key}->{GENE_NAME} = $gene_name;
	    $GENES->{$gene_hash_key}->{GENE_ID} = $gene_id;
	    $GENES->{$gene_hash_key}->{DESCRIPTION} = $gene_description;
	}

	
    }
}

#2nd pass for priTrpt 
  #1        maizev4 miRNA   6412285 6412305 .       +       .       ID=MIMAT0014028;Alias=MIMAT0014028;Name=zma-miR528a-5p;Derives_from=MI0013239
  #1        maizev4 miRNA   6412386 6412406 .       +       .       ID=MIMAT0015362;Alias=MIMAT0015362;Name=zma-miR528a-3p;Derives_from=MI0013239


# compute more features and store in the hash
#
#foreach my $g( keys %$GENES ){
#  my $genedata = $GENES->{$g};
#}
 
warn ("Done parsing\n");

# Load the data
my $sa = $ENS_DBA->get_adaptor('Slice');
my $n = 1;
my $number = scalar( keys %$GENES );

foreach my $gene_name_key( keys %$GENES ){
  my $agenedata = $GENES->{$gene_name_key};
  my $eGene = Bio::EnsEMBL::Gene->new();
  $eGene->analysis($ANALYSIS);
  $eGene->stable_id( $agenedata->{GENE_NAME} );
  $eGene->description( $agenedata->{DESCRIPTION} );
  $eGene->biotype( $BIOTYPE );
  $eGene->version(1);  
  $eGene->created_date( $date );
  $eGene->modified_date( $date );

  my $seq_region_name = $agenedata->{SEQ_NAME};
  my $slice = $sa->fetch_by_region( undef, $seq_region_name );
  
  # Add XREFs to gene
  #print Dumper ($agenedata->{ATTRIBS});
  add_xrefs( $eGene, $agenedata->{ATTRIBS}, 'GENE');  
  
  my $eTranscript = Bio::EnsEMBL::Transcript->new();
  my $trpt_name = $gene_name_key;
  $eTranscript->stable_id( $trpt_name );
  $eTranscript->analysis( $ANALYSIS );
  $eTranscript->version(1);
  $eTranscript->description( $agenedata->{DESCRIPTION} );
  $eTranscript->biotype( $BIOTYPE );
  $eTranscript->created_date( $date );
  $eTranscript->modified_date( $date );

  for my $relative_pos (@{$agenedata->{ATTRIBS}->{'relative_pos'}} ){
		
	#warn ("relative_pos=$relative_pos\n");
	my $mirna_attrib = Bio::EnsEMBL::Attribute->new
       		(-CODE => 'miRNA',
        	-VALUE => "$relative_pos");
	$eTranscript->add_Attributes($mirna_attrib);

   }

  my $transcript_xref;

    # Transcript xrefs...
    #print "For transcript\n";
    #print Dumper ($atrptdata->{ATTRIBS});
  add_xrefs( $eTranscript, $agenedata->{ATTRIBS}, 'TRANSCRIPT');
      
  my $eExon = new Bio::EnsEMBL::Exon;
  $eExon->start($agenedata->{START});
  $eExon->end($agenedata->{END});
  $eExon->strand($agenedata->{STRAND});	   
  $eExon->slice($slice);
  $eExon->stable_id( $trpt_name.'.exon' );
  $eExon->version( 1 );
      
  $eExon->phase( $noncoding_exon_phase );
  $eExon->end_phase( $noncoding_exon_phase );
  $eTranscript->add_Exon($eExon);


    ##################################
    # EnsEMBL add transcript to gene #
    ##################################
    
  $eGene->add_Transcript($eTranscript);
 
 
  if( $I ){
    my $dbid = $ENS_DBA->get_GeneAdaptor->store($eGene);
    print "// dbID = $dbid for gene = $gene_name_key ($n of $number)\n";
  }
  $n++;
}

#warn Dumper( $GENES );


exit();

sub parse_line{
 # Split gff line,\
  # start always <= end even for - strand genes
  #1        maizev4 miRNA_primary_transcript        6412284 6412406 .       +       .       ID=MI0013239;Alias=MI0013239;Name=zma-MIR528a
  #1        maizev4 miRNA   6412285 6412305 .       +       .       ID=MIMAT0014028;Alias=MIMAT0014028;Name=zma-miR528a-5p;Derives_from=MI0013239
  #1        maizev4 miRNA   6412386 6412406 .       +       .       ID=MIMAT0015362;Alias=MIMAT0015362;Name=zma-miR528a-3p;Derives_from=MI0013239

    my $line = shift;

    my( $seqname, $source, $feature, 
	$start, $end, $score, $strand, $frame, 
	$attribute ) = split( /\s+/, $line, 9 );
    #$feature = uc($feature);
    $seqname =~ s/chr0*//i;
    $strand = $strand eq '-' ? -1 : 1;
    
    my %attribs;
    foreach my $id( split( /;/, $attribute ) ){
	
	if( $id  =~ m/([^=]*)=(.*)/ ){
	    #ID=MI0013239;Alias=MI0013239;Name=zma-MIR528a
	    #ID=MIMAT0014028;Alias=MIMAT0014028;Name=zma-miR528a-5p;Derives_from=MI0013239
	    
	    my $key   =  uc $1;
	    my $value = $2;
	    $value =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
	    $value =~ s/^"//;
	    $value =~ s/"$//;
	    $attribs{$key} = $value;
	    warn ("parse line attrib $key => $value\n");
	    
	}
    }
    
    return { SEQ_NAME => $seqname,
	     SOURCE => $source,
	     FEATURE => $feature,
	     START => $start,
	     END => $end,
	     SCORE => $score,
	     STRAND => $strand,
	     FRAME => $frame,
	     ATTRIBS => \%attribs,
    };
}


sub calculate_mirna_trpt_pos{ 

  my ($priTrpt_info, $miRNA_coord) = @_;
  
  my $pri_seqname = $priTrpt_info->{SEQ_NAME};
  my $pri_start = $priTrpt_info->{START};
  my $pri_end = $priTrpt_info->{END};
  my $pri_strand = $priTrpt_info->{STRAND};

  my ($miRNA_seqname, $miRNA_start, $miRNA_end, $miRNA_strand) = split /\|/, $miRNA_coord;

  #warn ("pristrand is $pri_strand, miRNA string is $miRNA_coord, strand is $miRNA_strand\n");
  warn ("mature and pri miRNA not on same seqname, $miRNA_seqname v.s $pri_seqname") and return 0 
	unless($miRNA_seqname eq $pri_seqname);
  warn ("mature and pri miRNA not on same strand, $miRNA_strand v.s $pri_strand") and return 0 
        unless($miRNA_strand == $pri_strand);

	
  my ($rel_start, $rel_end) = (0,0);
  if( $miRNA_strand == 1){

	$rel_start = $miRNA_start-$pri_start+1;
	$rel_end = $miRNA_end-$pri_start+1;
	
  }elsif( $miRNA_strand == -1){
     
	$rel_start = $pri_end-$miRNA_end+1;
	$rel_end = $pri_end-$miRNA_start+1;

  }else{
	warn ("Unrecognized strand $miRNA_strand");
	return 0;
  }
  #warn("$miRNA_seqname, $miRNA_start, $miRNA_end, $miRNA_strand == $pri_seqname, $pri_start, $pri_end, $pri_strand,$rel_start-$rel_end \n" );
  push @{$priTrpt_info->{ATTRIBS}->{'relative_pos'}}, "$rel_start-$rel_end";
  return 1;
  
}

sub add_xrefs{ # 1.ensembl_object 
               # 2.attribs
               # 3.type (gene or transcript)
  # Add XREFs to gene/transcript
  # The xref primary id cannot be null
  # for Rice make gene TU_feat_name as primary id and pub_locus as display_label

  my $eObj    = shift;
  my $attribs = shift;
  my $type    = uc(shift);
  my $xref;
  
  for my $k (keys %{$attribs}){
      warn ("debug attrib in add_xrefs $k => $attribs->{$k}\n");
  }

  my $db_name = $type eq 'GENE' ? $MIBASE : $MIBASE_TRPT; 
  my $xref_label_key = $type eq 'GENE' ? 'ALIAS' : 'DERIVES_FROM';
  warn ("$type, $db_name xref key is $xref_label_key, accession is ", $attribs->{$xref_label_key} ); 
  $xref =        {
                       db                   => $db_name,
                       primary_acc          => $attribs->{$xref_label_key},
                       display_label        => $attribs->{$xref_label_key},
                       description          => $attribs->{DESCRIPTION}   || '',
                   };

  my $flag = 0;
  foreach my $dbentry( make_dbentries([$xref]) ){
    $eObj->add_DBEntry( $dbentry );
    #print "dbentry->dbname=", $dbentry->dbname, "\n";
    $eObj->display_xref( $dbentry ) unless $flag++;
  }
  

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

  my $axref  = shift;
  my @data  = @{ $axref || [] };
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

sub get_first_func  { 
  
  my $attribs_local = shift;

  first { defined } map { $attribs_local->{$_}  } @_;

}

1;
