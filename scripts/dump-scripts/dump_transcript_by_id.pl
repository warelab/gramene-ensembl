#!/usr/local/bin/perl 

=head1 NAME

dump_transcripts.pl - Make fasta files of transcripts, and 3' and 5'
    regions, for all Genes or for Genes given as arguments.
	

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-50/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
#use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } 
#        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
#             ensembl-draw/modules ensembl-compara/modules );
use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;


#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Gene;
use Getopt::Long;
use Pod::Usage;


=head1 SYNOPSIS

dump_transcripts.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --id_file           the file containing transcript ids

=head1 OPTIONS

=over 4


=item B<--exclude>

    Genes to ignore.  


=item B<--exclude-clone>
    
    Clones to ignore 
    This may be a comma-separated list.
    
=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped


=item B<--id_file>

the file containing transcript ids

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back


=cut

my ($species, $registry, $idfile);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"id_file"=>\$idfile
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  pod2usage("Needs ID file") unless $idfile;
}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $trpt_adaptor=$ENS_DBA->get_TranscriptAdaptor;
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 

open my $fh, $idfile or die "cannot open file $idfile";
my @trptids = grep { /\S+/} map { chomp; } <$fh>;

my %count;


$seqio = new Bio::SeqIO(-format => 'fasta',
			-file => ">${species}.fasta",
			#-fh     => \*STDOUT
		       );


foreach my $id (@trptids) {
  #print "! $geneid\n";
  
  $count{trpt}++;
  
  my($trans, $header) = split /\s+/, $id;
  
  my $cdna_seq=$trans->spliced_seq;
  $seq_obj = Bio::Seq->new(
			   -display_id => $comp_id,
			   -seq => $cdna_seq,
			  );
  
  $seqio->write_seq($seq_obj);
  $count{qualified_transcripts}++;
  
  
  #last;
}

for my $k (sort keys %count){
  print "$k = $count{$k}\n";
}



  
########################## subroutines ######################################

__END__


=head1 OUTPUT


=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

