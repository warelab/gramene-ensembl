use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

my $name=shift;
my $file=shift;


my $in  = Bio::SeqIO->new(   
			     -file   => $file,
			     -format => 'Fasta',
			     );    


my $out= Bio::SeqIO->new(
			 -format => 'Fasta',
			 -file   => ">$file.$name",
			 );


while (my $seqobj = $in->next_seq()) {
    
    my $id=$seqobj->display_id;
    print "id=$id, name=$name\n";
    if( $id =~ /$name/i){
	print "found";
	$out->write_seq($seqobj);
    }
}
