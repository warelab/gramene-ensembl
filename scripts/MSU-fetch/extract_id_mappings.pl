#!/usr/local/bin/perl -w 

=pod

=head1 NAME

extract_id_mappings.pl - extract the mapping between MSU internal id to public locua ID

=head1 SYNOPSIS

  load_genes_from_msu_gff3.pl [options] gff_file

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
 
=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

=head1 DESCRIPTION


Maintained by Sharon Wei <weix@cshl.edu>

=cut

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
  #unshift @INC, $BASEDIR.'/bioperl-live';
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

# make it a function for reuse,
# could parse in an anonymous sub from command line 
# "sub {}", then $func->() to invoke
# for grape gff
  #Gene GSVIVG00000001001
  #Trpt GSVIVT00000001001
  #prtn GSVIVP00000001001

sub parse_gene_name {

  my $prefix = "GSVIV";
  if(/${prefix}[A-Z](\d+)/i){
    return "${prefix}G$1";
  }else{
    return undef;
  }
}


use vars qw( $GFF_HANDLE );

my $date = time(); 

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Process args

  my $gff_file = $ARGV[0] || pod2usage("\nNeed the path to a gff file\n");

  map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
  } $gff_file;


  # Create a GFF stream
  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );


}


my $GENES = {};
my %TRPT2GENE = ();


while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /\#/ ); 
  next if ( $line =~ /^\s+/ );
  chomp $line;

#  print "$line\n";  #####
  #warn(Dumper($TRPT2GENE));

  # Split gff line,\
  # start always <= end even for - strand genes
  my( $seqname, $source, $feature, 
      $start, $end, $score, $strand, $frame, 
      $attribute ) = split( /\s+/, $line, 9 );
  $feature = uc($feature);
  $seqname =~ s/chr//i;

  my %attribs;
  foreach my $id( split( /;/, $attribute ) ){

    if( $id  =~ m/([^=]*)=(.*)/ ){
      #ID=13101.t00001;Name=TBC%20domain%20containing%20protein%2C%20expressed;Alias=LOC_Os01g01010
      #ID=11562.t00100;Name="ORF249";Note="pub_locus=LOC_Osp1g01070"

      my $key   = uc $1;
      my $value = $2;
      $value =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
      $value =~ s/^"//;
      $value =~ s/"$//;
      $attribs{$key} = $value;
      
    }
  }


  
  #get gene name for this feature
  my ($gene_id, $gene_name, $gene_description); 
  my ($transcript_id, $transcript_name, $parent_id);
  my $organell = 0; #a flag for mitochondrion/chloroplast gff3,
                    #they were converted from xml using tigr tools

  if( $feature eq 'GENE' ){

    $gene_id   = $attribs{'ID'};
    $gene_name = $attribs{'ALIAS'} || $gene_id;

    if( defined $attribs{'NOTE'} &&
      $attribs{'NOTE'} =~ /pub_locus=([\w.]+)/i ){
      $gene_name = $1;
      $organell = 1;
    }
    
    $gene_description = $attribs{'NAME'};
    
 #   print "Gene $gene_id $gene_name $gene_description\n";#####
    unless( $GENES->{$gene_id} ){
      $GENES->{$gene_id}->{GENE_NAME} = $gene_name;
      $GENES->{$gene_id}->{GENE_DESCRIPTION} = $gene_description;
      $GENES->{$gene_id}->{SEQ_NAME}  = $seqname;
      $GENES->{$gene_id}->{START}     = $start; #always smaller than end
      $GENES->{$gene_id}->{END}       = $end;
      $GENES->{$gene_id}->{STRAND}    = $strand eq '+' ? 1 : -1;
      $GENES->{$gene_id}->{ORGANELL}  = 1 if $organell;
    }
  }
  elsif( $feature eq 'MRNA' ){
   
    $transcript_id   = $attribs{ID};
    $transcript_name = $attribs{ALIAS} || $transcript_id;
    $gene_id         = $attribs{PARENT};

    next if ( $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id} );

    if ( $GENES->{$gene_id}->{ORGANELL} ){
      my $organell_trpt_cnt = ++$GENES->{$gene_id}->{ORGANELL_CNT};
      $transcript_name = $GENES->{$gene_id}->{GENE_NAME} . '.' . $organell_trpt_cnt; 
    };
 #   print "TRANSCRIPT_ID = $transcript_id, transcript_name = $transcript_name\n";
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{TRPT_NAME} 
      = $transcript_name;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{SEQ_NAME} 
      = $seqname;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{START} 
      = $start; 
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{END}  
      = $end;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{EXON_COUNT} 
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_COUNT} 
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_COUNT} 
      = 0;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{STRAND} 
      =  $strand eq '+' ? 1 : -1;

    $TRPT2GENE{$transcript_id} = $gene_id;
    
  }
  elsif( $feature eq 'CDS' ){ # use CDS and UTR as a proxy for exon

    $transcript_id = $attribs{'PARENT'} ;
    $gene_id = $TRPT2GENE{$transcript_id};
 #   print"CDS for $transcript_id, $gene_id\n";#####
    unless ( $gene_id ){
      die("cannot find gene for transcript $transcript_id at " . 
	  $GFF_HANDLE->input_line_number );
      
    }
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_COUNT} ++;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_START} ||= [];
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_END}   ||= [];
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_START}}, 
      $start;
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{CDS_EXON_END}}, 
      $end;

  }

  elsif( $feature eq 'FIVE_PRIME_UTR' || $feature eq 'THREE_PRIME_UTR'  ){
    $transcript_id = $attribs{'PARENT'};
    $gene_id = $TRPT2GENE{$transcript_id};
 #   print"UTR $transcript_id, $gene_id\n";#####
    unless ( $gene_id ){
      die("cannot find gene for transcript $transcript_id at ".
	 $GFF_HANDLE->input_line_number );
      
    }
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_COUNT} ++;
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_START} ||= [];
    $GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_END}   ||= [];
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_START}}, 
      $start;
    push @{$GENES->{$gene_id}->{TRANSCRIPTS}->{$transcript_id}->{UTR_EXON_END}}, 
      $end;
    
  }
  else{
    die( "Unrecognized feature $feature" );
  }

}


# compute more features and store in the hash
#
foreach my $g( keys %$GENES ){
  my $genedata = $GENES->{$g};

  print join "\t", ('gene', $g, "$genedata->{GENE_NAME}\n" );  
  
  foreach my $t( keys %{$genedata->{TRANSCRIPTS}} ){

    my $trptdata = $genedata->{TRANSCRIPTS}->{$t};

    print join "\t", ('transcript', $t, "$trptdata->{TRPT_NAME}\n" );
    
  }
}

