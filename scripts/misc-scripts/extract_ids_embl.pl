#!/usr/local/bin/perl

=pod


=head1 NAME

extract_id_embl.pl : get fasta sequences from a fasta file by their sequences names

=head1 SYNOPSIS

extract_id_embl.pl

Options:

    --help how brief help and exit
    --man    Show full documentation
    --embl the file with sequences
    --outfile the output file name

=head1 DESCRIPTION

Given a list sequence names and a fasta file with sequences, retrieve the fasta sequences whose sequence names are on the list and print out in the standard output


=cut

use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use Pod::Usage;

my ( $seqfile, $help, $outfile, $listfile);
my $debug;

GetOptions(
	   'help' => \$help,
	   'embl=s' => \$emblfile,
	   'outfile=s' => \$outfile,
	   'debug'    => \$debug,
	   ) or die "failed at GetOptions";

if( $help || !$emblfile ){
    pod2usage(-verbose => 2);
}

my @keys = qw(protein_id note db_xref locus_tag translation);
my $key_reg = join '|', @keys; 
my %gene_hash;
my @primary_tag = qw(gene mRNA CDS);

print "DEBUG: parse EMBL file $emblfile\n" if $debug;
my $seqio  = Bio::SeqIO->new( '-format' => 'EMBL' , -file => $emblfile);


while (my $seqobj = $seqio->next_seq() ) {
   
	my $gene_stable_id = '';
	my $trpt_stable_id = '';

	for my $feat_object ($seqobj->get_SeqFeatures) {
		
    		my $ptag = $feat_object->primary_tag;
		
 	   for my $tag ($feat_object->get_all_tags) {
		next unless $tag =~ /^$key_reg$/i;
		print "DEBUG matched key_reg $key_reg for tag $tag, gene_stable_id=$gene_stable_id\n" if $debug;
		my $value = join ',', ($feat_object->get_tag_values($tag)); 
		print "$ptag :: $tag = $value\n" if $debug;

		if( ! exists $gene_hash{$value} and $ptag eq 'gene' and $tag eq 'locus_tag'  ){	
			$gene_stable_id = $value;
			next;
		}

		if ($ptag eq 'mRNA' and $tag eq 'locus_tag'){
                        $gene_stable_id = $value;
                        next;
                }

		if ($gene_stable_id and $ptag eq 'mRNA' and $tag eq 'note'){
			$value =~ s/transcript *//gi;
			$trpt_stable_id = $value;
			print "DEBUG got transcript ID $trpt_stable_id\n" if $debug;
			next;
		}
		if ($gene_stable_id and $trpt_stable_id and $ptag eq 'CDS' and $tag eq 'note'){
			$value =~ s/.*transcript *//gi;
			$trpt_stable_id = $value;
                        next;
                }
		
		if ($gene_stable_id and $trpt_stable_id and $ptag eq 'CDS' and $tag eq 'protein_id'){
                        $gene_hash{$gene_stable_id}{mRNAs}{$trpt_stable_id}{protein_id} = $value;
                	next;
		}

		if ($gene_stable_id and $trpt_stable_id and $ptag eq 'CDS' and $tag eq 'translation'){
                        $gene_hash{$gene_stable_id}{mRNAs}{$trpt_stable_id}{translation} = $value;
                	next;
		}

		if ($gene_stable_id and $trpt_stable_id and $ptag eq 'CDS' and $tag eq 'db_xref'){
                        $gene_hash{$gene_stable_id}{mRNAs}{$trpt_stable_id}{db_xref} = $value;
                	next;
		}

	    }
	}

}

print "Print mapping to $outfile \n" if $debug;

open (FH, '>', $outfile) or die "Cannot open $outfile to write";

for my $gid( sort keys %gene_hash) {

	print "gene is $gid\n" if $debug;
	for my $tid ( sort keys %{$gene_hash{$gid}{mRNAs}} ){
		my $pid = $gene_hash{$gid}{mRNAs}{$tid}{protein_id};
		my $pseq = $gene_hash{$gid}{mRNAs}{$tid}{translation};
		my $xrefs = $gene_hash{$gid}{mRNAs}{$tid}{db_xref};
		
		print "tid is $tid\n pid=$pid\n xrefs=$xrefs\n" if $debug;
		
		print  FH join "\t", (
				$gid,
				$tid,
				$pid,
				$xrefs,
				"$pseq\n"
		);
	}

#	last;
}


__END__

#       foreach $key ( $ac->get_all_annotation_keys() ) {
#           @values = $ac->get_Annotations($key);
#           foreach $value ( @values ) {
#              # value is an Bio::AnnotationI, and defines a "as_text" method
#              print "Annotation ",$key," stringified value ",$value->as_text,"\n";
#              # also defined hash_tree method, which allows data orientated
#              # access into this object
#              $hash = $value->hash_tree();
#           }
#       }


