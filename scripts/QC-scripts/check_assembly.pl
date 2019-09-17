#!/lab/bin/perl 

=head1 NAME

check_assembly.pl - Compare the fasta sequence from data providers ftp site with that dumped from the ensembl database, and report the ones that don't match 

=head1 SYNOPSIS

check_assembly.pl  [options] fasta_file1 fasta_file2 ...
 
 Options:
    --help		help message
    --man		full documentation
    --registry_file     Default is $GrameneEnsemblDir/conf/ensembl.registry
    --species           species represent one of the genome assembly in ensembl db
    --acs               the assembly coordinate system of the sequences to be compared
    --ccs               the component coordinate system of the sequences to be compared
    --chrwrong		output fasta file of problematic chromosome pieces 
    --verbose

=head1 OPTIONS

=over 4

    
=item B<--chrwrong> 
 
name for fasta file  to receive asc (for example:chromosome) pieces that don't match the
ccs(for example:clone) part they're supposed to

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

fasta files for the assembly coordinate system


=head1 AUTHOR

   Sharon Wei (weix@cshl.edu)
   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut




use lib map { "/usr/local/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules conf  ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);

#use lib '.';

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;
use File::Glob;

use Bio::SeqIO;
use Bio::PrimarySeq;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;




my $verbose;
my %clone2acc;
my $ensembl_species;
my $registry_file;
my ($outchrwrong, $acs, $ccs);
    
{  #Argument Processing
  my $help=0;
  my $man=0;
  my ($filechrwrong);
  GetOptions( "help|?"=>\$help,
	      "man"=>\$man,
	      "species=s" => \$ensembl_species,
	      "registry_file=s" => \$registry_file,
	      "v|verbose" => \$verbose,
	      "chrwrong=s"=>\$filechrwrong,
	      "acs=s"=>\$acs,
	      "ccs=s"=>\$ccs,
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  if($filechrwrong) {
    $outchrwrong= Bio::SeqIO->new('-format' => 'fasta'
				  , '-file' => ">$filechrwrong")
      or die "can't open $filechrwrong for output:$!";
  }
}

###
### Get Ensembl DB adaptor
###

$ENV{'ENSEMBL_SPECIES'}    = $ensembl_species;

my $reg="Bio::EnsEMBL::Registry";	#Use this to get adaptors 
$registry_file    ||= $ENV{GrameneEnsemblDir}.'/conf/ensembl.registry';
Bio::EnsEMBL::Registry->load_all( $registry_file );
#  or die "load_all($registry_file) failed";

my $slice_adaptor =$reg->get_adaptor($ensembl_species,'core','Slice') 
  or die "can't get Slice adaptor for $ensembl_species";


my $dba=$slice_adaptor->db; #DBAdaptor
my $dbc=$dba->dbc;	#DBConnection
warn "user ".$dbc->username
  .", db ".$dbc->dbname."\n";


my $ens_dbh=$dbc->db_handle; #you can use this as a DBI database handle


for my $tigr_fasta (@ARGV) {
  
  my $inp_fa=Bio::SeqIO->new( -file => $tigr_fasta, '-format'=> 'fasta')
    or die "opening $tigr_fasta:$!" ;
  
  while( my $nextseq= $inp_fa->next_seq  ) {
    
    #create Bio::Seq for chr
    print "next_seq_display_id=", $nextseq->display_id, "\n";
    
    my $src_seq = uc($nextseq->seq);
    my $tigrseq= Bio::Seq->new (-display_id => $nextseq->display_id,
				-seq => $src_seq,
			       );
	  
    #parse out the database chr name
    my $chr;
    if( $acs =~ /chr/i){ 

	if(
	 $tigrseq->display_id =~ m/ chr0*(\d+)$ /ixms || 
	 $tigrseq->display_id =~ m/ ( chloroplast | mitochondrion ) /ixms 
	 ){
		$chr=$1;
	}elsif( $tigrseq->display_id =~ m/ (chr)? \s*(\S+) /ixms) {
      
    		 $chr=$2; 
	}else{ $chr= $tigrseq->display_id;}

	$chr    =~ s/chr(omosome)?_?0*//i;
	print "chr = $chr\n";
    } else {
      $chr=$tigrseq->display_id;

      warn "$acs seq: $chr\n";

    }

	if ($chr =~ /\w+3s$/i){
                $chr = '3s'; #chr3s 
	}

	warn ("ensembl db chr is $acs:$chr\n");

    #create chr slice 
    my $chr_slice=$slice_adaptor->fetch_by_region($acs,$chr);
    $chr_slice or warn "$acs $chr in $tigr_fasta not in ensembl db" 
      and next;

    #compare chr sequence between src and that from ensembl API
    print "$acs:$chr compara the sequences between source and that retrieved from db by ensembl API\n";
    my $ensembl_seq = uc( $chr_slice->seq );
    if( $src_seq ne $ensembl_seq ){
	print "$acs:$chr ERROR!!! sequences mismatch between source and ensembl\n";
	
	$outchrwrong->write_seq($chr_slice) if $outchrwrong;
	#next;
    }else{
	print "$acs:$chr MATCHED!!! sequences between source and ensembl\n";
    }

    #compara the sequences at the assembly path for component coord system $ccs
    print "$acs:$chr | inspect projection from $acs to $ccs\n";
    for my $segment ( @{$chr_slice->project($ccs)} ) {
      #print("In projection for: ", $segment->from_start, ", ",$segment->from_end, "\n");

      my $clone_slice=$segment->to_Slice();
      my $tigr_segseq=$tigrseq->subseq($segment->from_start,$segment->from_end);
      
      #somehow ensembl convert non ATGCN into N
      #$tigr_segseq =~ tr/ATGCN/N/c; #this is what clean_dna.pl does
      
      if( $tigr_segseq ne $clone_slice->seq ) {
	my $tigr_seg_name = join '', ("$chr:", 
				      $segment->from_start, 
				      "-", 
				      $segment->from_end);
	my $clone_name = join '', ($clone_slice->seq_region_name, ":",
				  $clone_slice->start,
				  '-',
				  $clone_slice->end);


	print "mismatch $chr:",$segment->from_start," to "
	  ,$segment->from_end," vs "
	    ,$clone_slice->seq_region_name,":"
	      ,$clone_slice->start," to "
		,$clone_slice->end, " ", $clone_slice->strand
		  ,"\n>$tigr_seg_name\n$tigr_segseq\n>$clone_name\n", $clone_slice->seq, "\n";
	
	if($outchrwrong) {
	  my $seq=$clone_slice->project($acs)->[0]->to_Slice;
	  $outchrwrong->write_seq($seq);
	}

      } else {
	print "$chr:",$segment->from_start," to ",$segment->from_end," == "
	  ,$clone_slice->seq_region_name,":"
	    ,$clone_slice->start," to "
	      ,$clone_slice->end
		,"\n" if $verbose;
      }

      
    } #end of each projection
    #last;
  } #end of while each asm seq

} #end of for each file


    
# end of main program


