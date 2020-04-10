#!/lab/bin/perl -w
  

=head1 Name

    split_fasta.pl - split a fasta file into smaller chuncks 

=head1 SYNOPSIS

    split_fasta.pl OPTIONS

=head2 OPTIONS


    -i            input_file
    -d            directory for ouput files
    -clean        remake the directory, clean everything in it 
    -min_lines    minimum lines each splitted file has
    -fn           number of splitted files
    -index_start  the starting index of the splitted files, has to be no smaller than 1, default is 1

    
=head1 Author

Sharon Wei (weix@cshl.edu)


=cut




use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path;

use Bio::SeqIO;
use Bio::Seq;

my($help, $input_file, $out_dir, $min_lines, $fn, $index_start, $clean, $interproscan);

GetOptions(
	   "help"          => \$help,
	   "i=s"           => \$input_file,
	   "d=s"           => \$out_dir,
	   "clean"         => \$clean,
	   "min_lines=i"   => \$min_lines,
	   "fn=i"          => \$fn,
	   "index_start=i" => \$index_start,
	   "interproscan"  => \$interproscan,
	   ) || pod2usage();


if($help || !$input_file || !$out_dir || !($min_lines || $fn) ){
    pod2usage();
}


$index_start ||= 1;


if($index_start < 1){
    
    print "index_start has be greater and equal to 1\n";
    pod2usage();
}


# prepare out_dir
#

if($clean){
    rmtree($out_dir, 0, 1) if -d $out_dir;
    eval { mkpath( $out_dir, 0, 0755 ) };  #rwxr-xr-x
    die "Fail to make dir: '$out_dir': $@\n" if $@; 
}else{
    die "$out_dir not exist" unless -d $out_dir;
}

split_file(
	   out_dir => $out_dir,
	   ori_file => $input_file,
	   min_lines => $min_lines,
	   split_cnt => $fn,
	   index_start => $index_start,
	   );



sub split_file{

    my (%args) = @_;
    
    my $out_dir   = $args{out_dir};
    my $ori_file  = $args{ori_file};
    my $min_lines = $args{min_lines};
    my $split_cnt = $args{split_cnt};
    my $index_start = $args{index_start};
    
    
    warn "Start splitting $ori_file\n";
    
    my ($ori_file_basename, $path, $type) = fileparse( $ori_file, qr{\.[^.]*} );
    
    #Usually, the query dataset split by the number of files its split into
    
    if(  $split_cnt  ){
    
	my $total_sequences = qx[grep '>' $ori_file | wc -l];
	my $seq_per_file = int($total_sequences / $split_cnt);
	my $remainder    = $total_sequences % $split_cnt;
	my $in  = Bio::SeqIO->new(   
				     -file   => $ori_file,
				     -format => 'Fasta',
				     );    
	my $file_counter = $index_start;
	my $outfile = $out_dir . '/' . $ori_file_basename . '.' . $file_counter . ".fa";
	
	my $out= Bio::SeqIO->new(
				 -format => 'Fasta',
				 -file   => ">$outfile",
				 );
	my $current_seq_count = 0;
	my $total_seq_count = 0;
      SQ:
	while (my $seqobj = $in->next_seq()) {
	    $current_seq_count++;
	    $total_seq_count++;

	    if($interproscan){ #interproscan cannot take illegal aa such as . *
	    	my $new_seq = $seqobj->seq;
	       	$new_seq =~ s/[.x*]+$//i;
	
		unless ($new_seq =~ /[.*]/){
		
			$seqobj = new Bio::Seq(
				-display_id => $seqobj->display_id(),
				-seq        => $new_seq,
				);
			$out->write_seq($seqobj);
		}
		
	     }else{
	    	$out->write_seq($seqobj);
	     }	

	    if ($total_seq_count == $total_sequences) { last SQ; }
	    if ( $remainder && $current_seq_count == $seq_per_file + 1 || 
		 !$remainder && $current_seq_count == $seq_per_file ) {
		# reset the seq count and change the out file
		$current_seq_count = 0;
		$file_counter++;
		--$remainder if $remainder > 0;

		$outfile = $out_dir . '/' . 
		    $ori_file_basename . '.' . $file_counter . ".fa";
				
		$out= Bio::SeqIO->new(
				      -format => 'Fasta',
				      -file   => ">$outfile",
				      );
	    }
	    
	}
	
	
	
	
	# the target genome file split by minimum lines
	
    }elsif($min_lines){
	
	open OF, $ori_file or die("Can't open the file: $ori_file. $!\n");
	
	my $line_counter;
	my $file_counter    = $index_start;
    
	my $split_file_name = "$out_dir/${ori_file_basename}.${file_counter}.fa";

	open (OUT_FILE, ">$split_file_name") or
	    die("Can't open the outfile $split_file_name for writing: $!\n");
	
	$file_counter++;

	while(<OF>){
	    
	    $line_counter++;
      
	    if( m/^>/ && $line_counter > $min_lines ){

		$split_file_name = "$out_dir/${ori_file_basename}.${file_counter}.fa";
	    
		open (OUT_FILE, ">$split_file_name")  or
		    die("Can't open the outfile for writing: $!\n");
		
		$file_counter++;
		$line_counter = 0;
	    }
	    
	    print OUT_FILE $_ unless (/^$/s);
	}

	close OUT_FILE or die("Can't close out_file $split_file_name: $!\n");
	close OF or die("Can't close $ori_file: $!\n");
	
    }else{
	die("No splitting parameters defined (split count or min lines)\n");
    }
    
    

    return;
    
    
    
}







__END__

	#use the card dealing algorithm
	my @fh;
	my $line   = 0;
	my $file ;

	while(<OF>) {
	    
	    next if /^\s*$/s;

	    if( m/^>/ ){
		$line = $line%$split_cnt;
		
		unless($fh[$line]){
		    $file = ">$out_dir/${ori_file_basename}." . $line+1 . ".fa";
		    open my  $fh, $file or die( "Can't open $file: $!");
		    $fh[$line] = $fh;
		}
	    }
	    print { $fh[$line] } $_;
	    
	}
