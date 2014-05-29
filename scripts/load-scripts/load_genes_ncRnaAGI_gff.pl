##!/usr/local/bin/perl -w 
 
=pod

=head1 NAME

load_genes_ncRNA_gff.pl - load non_coding RNA genes into EnsemblDB

=head1 SYNOPSIS

  load_genes_ncRNA_gff.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -l|--logic_name       analysis.logic_name to use for the genes.
  -n|--no_insert        Do not make changes to the database. For debug.
  -v|--verbose          Print verbose output
  -source|--source	annotation source such as JGI, MIPS
  -correction           self correct on start > stop cases

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
  analysis.logic_name to use for the genes. Default is ensembl.

B<--source>
  annotation source, for xref

B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

B<-v|--verbose>
  Print verbose output

=head1 DESCRIPTION

  Parses ncRNA annotations from AGI gff file , and uses the data
  therein to populate the gene, exon, transcript 
  tables (inc stable_id tables) in the ensembl DB indicated by the
  --species and --ensembl_registry.

    the data format looks like:
Obart_chr6      Infernal        ncRNA   14823665        14823778        61.3    -       .       ID=ncRNA_1;Name=5S_rRNA;family=RF00001
Obart_chr3      Infernal        ncRNA   7404078 7404194 58.2    +       .       ID=ncRNA_2;Name=5S_rRNA;family=RF00001
Obart_chr7      Infernal        ncRNA   869871  869987  54.8    +       .       ID=ncRNA_3;Name=5S_rRNA;family=RF00001
Obart_chr1      Infernal        ncRNA   7274467 7274580 48.5    +       .       ID=ncRNA_4;Name=5S_rRNA;family=RF00001
Obart_chr4      Infernal        ncRNA   22629876        22630036        142.2   -       .       ID=ncRNA_5;Name=U1;family=RF00003

Oglab01_0032    MIPS_Oglab_v2.0 tRNA    75794   75865   61.03   -       .       ID=trna_165;Name=tRNA-Gln
Oglab01_0048    MIPS_Oglab_v2.0 rRNA_fragment   187155  187398  398.00  +       .       ID=rRNA_2;Name=5S_rRNA
Oglab01_0053    MIPS_Oglab_v2.0 tRNA    16782   16853   71.08   +       .       ID=trna_166;Name=tRNA-Gln
...
Oglab01_0259    MIPS_Oglab_v2.0 tRNA    226597  226681  40.47   +       .       ID=trna_190;Name=tRNA-Met
Oglab01_0259    MIPS_Oglab_v2.0 tRNA_intron     226635  226645  .       +       .       ID=trna_190_intron;Parent=trna_190

B<The Ensembl Registry>

  The database connection details for both Ensembl and interproscan
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'Populus_trichocarpa', 
    '-group'   => 'core', 
    '-dbname'  => 'populus_trichocarpa_core_42_jgi2004', );
  ---

B<Restoring the database>

  If the script bombs out half-way through, your database will be in a
  partial loaded state (i.e. in a bit of a mess). Here is some SQL to
  help put it right;

  DELETE FROM gene, transcript, translation, exon_transcript, exon,
              gene_stable_id, transcript_stable_id,
              translation_stable_id, exon_stable_id, analysis
  WHERE       gene.gene_id = transcript.gene_id
  AND         transcript.transcript_id=translation.transcript_id
  AND         transcript.transcript_id=exon_transcript.transcript_id
  AND         exon_transcript.exon_id=exon.exon_id
  AND         gene.gene_id = gene_stable_id.gene_id
  AND         transcript.transcript_id=transcript_stable_id.transcript_id
  AND         translation.translation_id=translation_stable_id.translation_id
  AND         exon.exon_id=exon_stable_id.exon_id
  AND         gene.analysis_id=analysis.analysis_id
  AND         analysis.logic_name="fgenesh_gene"

  The above could be improved with some left joins so that the
  constraints do not fail for partly-loaded gene models.

  Here is the brute-force approach that removes all gene models
  regardless of analysis;

  delete from gene; 
  delete from gene_stable_id; 
  delete from transcript; 
  delete from transcript_stable_id; 
  delete from translation; 
  delete from translation_stable_id;
  delete from exon_transcript; 
  delete from exon; 
  delete from exon_stable_id;


Maintained by Will Spooner <whs@ebi.ac.uk>

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

use vars qw( $BASEDIR );
BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
  unshift @INC, '/usr/local/ensembl-live/ensembl/modules';

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
use Readonly;


Readonly my $RNA_REGEX => qr/^(RRNA_FRAGMENT|NCRNA|TRNA)$/i;
Readonly my $TRNA_INTRON_REGEX => qr/^(TRNA_INTRON)$/i;
Readonly my $biotype => 'non_coding';

use vars qw( $ENS_DBA $I $INSERT $GFF_HANDLE $ANALYSIS $V 
             $GENES  %PARENT $ANNOT_SOURCE $correction);


my $date = time(); 

#print "date is $date\n";

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $logic_name, $no_insert );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        "logic_name=s"       => \$logic_name,
        "no_insert"          => \$no_insert,
        "verbose"            => \$V,
	"source=s"	     => \$ANNOT_SOURCE,
	"correction"         => \$correction,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args
  
  $ANNOT_SOURCE = 'AGI' unless $ANNOT_SOURCE;

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

  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $reg\n" ) &&
                pod2usage() );

  # Create the analysis
  print( "[INFO] Genes will use analysis.logic_name of $logic_name\n" );
  $ANALYSIS = Bio::EnsEMBL::Analysis->new
      ( -logic_name => $logic_name );

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

#----------------------------------------------------------------------
# Main

while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /\#/ ); 
  next if ( $line =~ /^\s+/ );
  chomp $line;
  $line =~ s/\r//g;  #get rid of weird ^M

  print "\n$line\n" if $V;

  # Split gff line
  my( $seqname, $source, $feature, 
      $start, $end, $score, $strand, $frame, 
      $attribute ) = split( /\s+/, $line, 9 );
  $feature = uc($feature);
  $seqname =~ s/.*chr0*(\d+)/$1/i;

 if ($correction && $start > $end ){
	my $temp;
	$temp = $start;
	$start = $end;
	$end = $temp;
 }

  my %attribs;
  foreach my $id( split( /;\s*/, $attribute ) ){
    #print "split attribute $id\n";
    if( $id  =~ m/(\S+)[\s\=]+\"*([^\"]+)/ ){
      $attribs{uc ($1)} = uc($2);
      #print "attrib mapping $1 => $2\n";
    }
  }

  my $id = $attribs{ID} ;
  my $parent = $attribs{PARENT};
  $PARENT{ $id } = $parent if $parent;
  my $gene_id = $parent ? &get_gene_id( $parent ) : $id;
  
  my $name = $attribs{NAME} ? $attribs{NAME} : '';
  print "id=$id, gene_id=$gene_id,,,\n" if $V;

  my $family = $attribs{FAMILY} ? $attribs{FAMILY} : '';
  print "id=$id, gene_id=$gene_id,,,\n" if $V;

#tRNA, rRNA_fragment, tRNA_intron

  if( $feature =~ $RNA_REGEX ){
      #my ($biotype) = split /_/, $feature; 
      print "[INFO] Processing $biotype gene $gene_id, name $name\n" if $V;
      $GENES->{$gene_id}->{GENE_NAME}  = $id;
      $GENES->{$gene_id}->{GENE_NAME_ALIAS}  = $name;
      $GENES->{$gene_id}->{ATTRIB}  = $attribute;
      $GENES->{$gene_id}->{FAMILY}  = $family;
      $GENES->{$gene_id}->{BIOTYPE}  = $biotype;
      $GENES->{$gene_id}->{SEQ_NAME}   = $seqname;
      $GENES->{$gene_id}->{TRPT_COUNT} = 0;
      $GENES->{$gene_id}->{STRAND}     = $strand eq '+' ? 1 : -1;
      $GENES->{$gene_id}->{STRART}     = $start;
      $GENES->{$gene_id}->{END}     = $end;
      
      $GENES->{$gene_id}->{TRPT}->{$id}->{TRPT_NAME}  = $id; 
      $GENES->{$gene_id}->{TRPT}->{$id}->{TRPT_NAME_ALIAS} = $name;
      $GENES->{$gene_id}->{TRPT}->{$id}->{EXON_COUNT} = 0;
      $GENES->{$gene_id}->{TRPT}->{$id}->{CDS_COUNT}  = 0;
      $GENES->{$gene_id}->{TRPT}->{$id}->{START}  = $start;
      $GENES->{$gene_id}->{TRPT}->{$id}->{END}  = $end;
      
  }elsif( $feature =~ $TRNA_INTRON_REGEX ){

      print "\n[INFO] Processing intron for gene $gene_id, parent=$parent\n" if $V;
      $GENES->{$gene_id}->{TRPT}->{$parent}->{INTRON}{$start} = $end;      
  }
  
}

warn ( "[INFO] found " . scalar( keys %$GENES) . " genes in gff file\n" );

foreach my $geneid ( sort keys %$GENES ){

    #print "geneid=$geneid\n\n";
  my $genedata = $GENES->{$geneid};
  $genedata->{TRANSCRIPT_COUNT} = keys %{$genedata->{TRPT}};

  warn ( "[INFO] found $genedata->{TRANSCRIPT_COUNT} transcripts in gene $geneid\n" );

  
  foreach my $trptid ( sort keys %{$genedata->{TRPT}} ){
      

      my $trptdata = $genedata->{TRPT}->{$trptid};
      my $trpt_start =  $trptdata->{START};
      my $trpt_end =  $trptdata->{END};
	        
      if( $trptdata->{INTRON} ){
	  #get exons by parsing out intron boundries
	  my $intron_starts = [sort {$a<=>$b} keys %{$trptdata->{INTRON}}];
	  my $intron_ends   = [sort {$a<=>$b} values %{$trptdata->{INTRON}}];
	  
	  # Get exons starts array and ends array
	  print "Parse out exons for $trptid, $trpt_start, $trpt_end\n";
	  my ($exon_starts, $exon_ends) = &get_exon_starts_ends
	      ($trpt_start, $trpt_end, $intron_starts, $intron_ends) ;
	  
	  $trptdata->{EXON_START} = $exon_starts;
	  $trptdata->{EXON_END} = $exon_ends;
	  $trptdata->{EXON_COUNT} = scalar @{$exon_ends};
	  
	  print "$trptid: exon_starts=", join ",", @$exon_starts, "\n" if $V;
	  print "$trptid: exon_starts=", join ",", @$exon_ends, "\n" if $V;
      }else{

	  $trptdata->{EXON_START} = [$trpt_start];
	  $trptdata->{EXON_END} = [$trpt_end];
	  $trptdata->{EXON_COUNT} = 1;
		  
      }
  }
}


#exit;

# Load the data
my $sa = $ENS_DBA->get_adaptor('Slice');
my $n = 1;
my $number = scalar( keys %$GENES );
warn( "[INFO] Storing $number genes to DB\n" );
foreach my $geneid( sort keys %$GENES ){

  my $genedata = $GENES->{$geneid};
  my $family = $genedata->{FAMILY};

  my $slice = $sa->fetch_by_region( undef, $genedata->{SEQ_NAME} );

  unless( $slice ){
      print STDERR "ERROR: No seq_region found for ", $genedata->{SEQ_NAME},"\n";
      next;
  }

  print "seq name is",$genedata->{SEQ_NAME} , "\n";
  print "slice name is", $slice->name, "\n";
  my $eGene = Bio::EnsEMBL::Gene->new();
  $eGene->analysis($ANALYSIS);
  $eGene->stable_id( $genedata->{GENE_NAME} );
  $eGene->description( $genedata->{ATTRIB}  || $genedata->{GENE_DESCRIPTION} || $genedata->{GENE_NAME_ALIAS} || '');
  $eGene->version(1);  
  $eGene->biotype($genedata->{BIOTYPE});
  $eGene->created_date( $date );
  $eGene->modified_date( $date );
  $eGene->slice($slice);

  # Add XREFs to gene
  # The xref primary id cannot be null
  # make gene TU_feat_name as primary id and pub_locus as display_label

  my $gene_src = "${ANNOT_SOURCE}_GENE";

  my $gene_xref = [{
		    db                   => $gene_src,
		    primary_acc          => $genedata->{GENE_NAME},
		    display_label        => $genedata->{GENE_NAME_ALIAS} || $genedata->{GENE_NAME},
		    description          => $genedata->{ATTRIB} || $genedata->{GENE_DESCRIPTION}   || $genedata->{GENE_NAME_ALIAS} || '',
		   },
                  ];
  
  foreach my $dbentry( make_dbentries($gene_xref) ){
    $eGene->add_DBEntry( $dbentry );
    if( $dbentry->dbname eq $gene_src ){ # Use as display_xref
      $eGene->display_xref( $dbentry );
    }
  }
  
  

  foreach my $trptid( sort keys %{$genedata->{TRPT}} ){
    
    my $trptdata = $genedata->{TRPT}->{$trptid};
    my $trpt_name =  $trptdata->{TRPT_NAME};

    my $eTranscript = Bio::EnsEMBL::Transcript->new();
    $eTranscript->analysis( $ANALYSIS );
    $eTranscript->stable_id( $trpt_name );
    $eTranscript->description( $genedata->{ATTRIB} || $genedata->{GENE_DESCRIPTION} || $trptdata->{TRPT_NAME_ALIAS} || '');
    $eTranscript->version( 1 );    
    $eTranscript->biotype( $trptdata->{BIOTYPE} || $biotype );
    $eTranscript->created_date( $date );
    $eTranscript->modified_date( $date );
    $eTranscript->slice( $slice );

    my $trpt_src = "${ANNOT_SOURCE}_TRANSCRIPT"; 
    my $transcript_xref = [{
			    db            => $trpt_src,
			    primary_acc   => $trpt_name,
			    display_label => $trptdata->{TRPT_NAME_ALIAS} || $trpt_name || '',
			    description   => $genedata->{ATTRIB} || $genedata->{GENE_DESCRIPTION} || $trptdata->{TRPT_NAME_ALIAS} || '',
			   },
			  ];

    # Transcript xrefs...
    foreach my $dbentry( make_dbentries($transcript_xref) ){
      $eTranscript->add_DBEntry( $dbentry );
      if( $dbentry->dbname eq $trpt_src ){ # Use as display_xref
        $eTranscript->display_xref( $dbentry );
      }
    }

    #my $slice = $sa->fetch_by_region( undef, $genedata->{SEQ_NAME} );
    
    my $start_exon;
    my $end_exon;
    my $trans_start_exon;
    my $trans_end_exon;
    my $translation_offset_left;
    my $translation_offset_right;
    my $exon_start_phase = 0;
    my $exon_end_phase   = 0;
    my $phase_diff       = 0;
    my $last_exon_phase  = 0;
    
    
    #########
    # EXONS #
    #########
    
    for ( my $exon = 0; $exon < $trptdata->{EXON_COUNT}; $exon++ ){
      
      my $exon_start  = $trptdata->{EXON_START}[$exon];
      my $exon_end    = $trptdata->{EXON_END}[$exon];
      
      my $eExon = new Bio::EnsEMBL::Exon;
      $eExon->start($exon_start);
      $eExon->end($exon_end);
      $eExon->strand($genedata->{STRAND});	   
      $eExon->slice($slice);
      $eExon->stable_id( $trptid.'.exon'.($exon+1) );
      $eExon->version( 1 );
      
      $eTranscript->add_Exon($eExon);
 
      #phase data is required to create transcirpt
      $phase_diff = 0;
      $exon_start_phase = 0;
      $exon_end_phase   = 0;
	
	
      # set exon phase
      $eExon->phase($exon_start_phase);
      $eExon->end_phase($exon_end_phase);
	
      #start exon
      unless ( $start_exon ) {
	  if ( $exon == 0 ) {
	      $start_exon = $eExon; 
	  }
      }
      
      # Final exon
      unless ( $end_exon ) {
	  if ( $exon == $trptdata->{EXON_COUNT}-1 ) {
	      $end_exon = $eExon;
	  }
      }
      
  } #_ END of exon loop _#
      

   
    $eTranscript->start_Exon($start_exon);
    $eTranscript->end_Exon($end_exon);
    
    #my $eTranslation = new  Bio::EnsEMBL::Translation
     # (
     #  -START_EXON  => $trans_start_exon,
     #  -END_EXON    => $trans_end_exon,
     #  -SEQ_START   => $translation_offset_left,
     #  -SEQ_END     => $translation_offset_right,
     #  -STABLE_ID   => $trpt_name,
     #  -VERSION     => '1',
     # );
    
    #########################################
    # EnsEMBL add translation to transcript #
    #########################################
    #print "before translation\n";
    #$eTranscript->translation($eTranslation);
    
    ##################################
    # EnsEMBL add transcript to gene #
    ##################################
    
    $eGene->add_Transcript($eTranscript);
    
  }
  
  my $dbid = '???';
  if( $I ){
      eval{  $dbid = $ENS_DBA->get_GeneAdaptor->store($eGene); };
      if ($@){
         print STDERR "ERROR adding gene $genedata->{GENE_NAME}, $@";   
	 next;
      }else{
         print "// dbID = $dbid for gene = $genedata->{GENE_NAME} ($n of $number)\n";
      }
  }

  $n++;

}

warn( "[INFO] Stored $number genes to DB\n" );
warn( "[INFO] DONE!\n" );

sub get_gene_id{

  my $id = shift;

  #Find gene source
  
  my $gene = $id;
  while( defined $PARENT{ $id } and $gene = $PARENT{ $id }){
    $id = $gene;
  }    

  return $gene;
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



sub get_exon_starts_ends{

    my ($trpt_start, $trpt_end, $intron_starts, $intron_ends) = @_;

    my @intron_starts_srt = sort {$a<=>$b} @$intron_starts;
    my @intron_ends_srt = sort {$a<=>$b} @$intron_ends;
    
    #check the data integrity, intron_ends need to >= intron_starts
    # trpt_start smaller than intron_starts
    # trpt_end greater than intron_ends

    my $intron_start_cnt = scalar @intron_starts_srt;
    my $intron_end_cnt   = scalar @intron_ends_srt;

    if( $intron_start_cnt != $intron_end_cnt ){
	warn("ERROR: intron starts count ($intron_start_cnt) not equal to intron ends count ($intron_end_cnt)");
	return undef;
    }

    if( $trpt_start >= $intron_starts_srt[0] ){
	warn("ERROR: trpt_start ($trpt_start) not smaller than smallest intron start $intron_starts_srt[0]");
	return undef;
    }

    if( $trpt_end <= $intron_ends_srt[-1] ){
	warn("ERROR: trpt_end ($trpt_end) not greater than biggest intron end $intron_ends_srt[-1]");
	return undef;
    }
    
    my @exon_starts = ($trpt_start);
    my @exon_ends ; #= ($trpt_end);
    for (my $i=0; $i<$intron_end_cnt; $i++){
	if($intron_ends_srt[$i] < $intron_starts_srt[$i]){
	    warn("ERROR: for the index $i intron, end ($intron_ends_srt[$i]) is smaller than start ( $intron_starts_srt[$i])\n");
	    return undef;
	}
    
	push @exon_starts, $intron_ends_srt[$i]+1;
	push @exon_ends, $intron_starts_srt[$i]-1;
	
    }
    
    push @exon_ends, $trpt_end;
    
    
    return (\@exon_starts, \@exon_ends);
}


1;
