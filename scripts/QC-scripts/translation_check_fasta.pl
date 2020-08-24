#!/lab/bin/perl 

=head1 NAME

translation_check_fasta.pl - check translations given in a fasta file
 	against translations of corresponding genes in the database
	ok=at least one transcript of the gene has a matching translation

=cut


use lib '/usr/local/ensembl-live/ensembl/modules';
use lib '/usr/local/ensembl-live/ensembl-compara/modules';

use strict; 
use warnings;


use Carp qw(cluck);

use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Registry;


=head1 SYNOPSIS

translation_check_fasta.pl  [options] fasta-file [ ... fasta-file]
 
 Options:
    --id		how to make translation, transcript, or gene id from info in FASTA file
    --help		help message
    --man		full documentation
    --species         species in EnsEMBL registry to use for db [required]
    --registry_file   Default is $GrameneEnsemblDir/conf/ensembl.registry


=head1 OPTIONS

=over 4

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

=item B<--id>

3 components separated by /
what:
    gene (stable id) - require that at least one translation matches
      --this is the default for historical reasons
    transcript (stable id)
    translation (stable id)
    xref:I<database> (I<database> is name for xref)
	get translation if exists, else transcript, else gene

whence:
    id(default)
    description

parsing:
    reg expression which should set $1 to the id
    Else it's the whole description

E.g.
TIGR chrNN.pep

--id='xref:TIGR_FN//([^|]+)'

    
    
=back

=head1 ARGUMENTS

 Fasta files to process

 Yes, it should process standard input if no files are
 given, but that's not implemented.
 

=cut

my $rank1=1;
my %grape_id_prefix = (GENE       => 'GSVIVG',
		       TRANSCRIPT => 'GSVIVT',
		       TRANSLATION => 'GSVIVP',
		      );
my $ensembl_species=$ENV{ENSEMBL_SPECIES};
my $registry_file;
my ($idwhat,$idwhence,$idregex,$idxref,$verbose, $idfunc);
    {
	my $id='';
	my $help=0;
	my $man=0;
	GetOptions( "id=s"=>\$id,
		    #"idfunc=s"=>\$idfunc,
		    "species=s" => \$ensembl_species,
		    "registry_file=s" => \$registry_file,
		    "v+"=>\$verbose,
		    "help|?"=>\$help,
		    "man"=>\$man,
		    "rank1=i"=>\$rank1,
		  )
	  or pod2usage(2);
	pod2usage(-verbose => 2) if $man;
        pod2usage(1) if $help;
	($idwhat,$idwhence,$idregex)=split /\//,$id;
	$idwhat ||='gene';
	$idwhence ||='id';
	$idxref=$1 if $idwhat =~ /xref:(.+)/;
	pod2usage(2) unless
	    ($idxref || $idwhat eq 'gene' || $idwhat eq 'transcript' 
		     || $idwhat eq 'translation' )
	    and ( $idwhence eq 'id' || $idwhence eq 'description')
	    and  ( ! $idregex or $idregex =~ /\(.+\)/ )  #need to capture the id with parentheses
	    ;
	$idregex=qr($idregex) if $idregex;
	
	
    }
# a nice idea, never used
#    my $code;
#    if ($geneidcode) {
#        $code='sub geneid {  local($_)=@_;'.$code.'}';
#	eval $code;
#	print STDERR "eval\n $code\n===========\n$@\n" and pod2usage(1) if $@;
#    }

#my $get_id_func = eval "sub { my $aaid=shift; $aaid=$idfunc; return $1; }" if $idfunc;
#print "created function from $idfunc\n" if ( defined $get_id_func);

$verbose||=0;
my $verbose_count= $verbose>1?1e6
                  :$verbose?20
		  :3; # it's decremented


###
### Get Ensembl DB adaptor
###

    $ENV{'ENSEMBL_SPECIES'}    = $ensembl_species;

    my $reg="Bio::EnsEMBL::Registry";	#Use this to get adaptors 
    $registry_file    ||= $ENV{GrameneEnsemblDir}.'/conf/ensembl.registry';
Bio::EnsEMBL::Registry->load_all( $registry_file );
          #or die "load_all($registry_file) failed";
     
  my $slice_adaptor =$reg->get_adaptor($ensembl_species,'core','Slice') 
     or die "can't get Slice adaptor for $ensembl_species";

     
     my $dba=$slice_adaptor->db; #DBAdaptor
     my $dbc=$dba->dbc;	#DBConnection
     warn "user ".$dbc->username
           .", db ".$dbc->dbname."\n";


my $ens_dbh=$dbc->db_handle; #you can use this as a DBI database handle


my $gene_adaptor=$reg->get_adaptor($ensembl_species,'core','Gene');
my $transcript_adaptor=$reg->get_adaptor($ensembl_species,'core','Transcript');
$transcript_adaptor || die "cannot get $transcript_adaptor";
my $translation_adaptor=$reg->get_adaptor($ensembl_species,'core','Translation');
my $dbe_adaptor=$reg->get_adaptor($ensembl_species,'core','DBEntry');
my $coord_system_adaptor=$reg->get_adaptor($ensembl_species,'core','CoordSystem');

print "rank1=$rank1\n";
my $rank1_coord_system=$coord_system_adaptor->fetch_by_rank($rank1)->name;
my $seq_coord_system=$coord_system_adaptor->fetch_by_name('seqlevel')->name;
#print "rank coord_system_ $rank1_coord_system\n $seq_coord_system\n";
my $stdCodonTable=Bio::Tools::CodonTable->new();

my %count;

my $fixes='';

print STDERR "id what=$idwhat, whence=$idwhence\n";

while(my $infile=shift) {
    my $seqin=Bio::SeqIO->new( '-format' => 'Fasta', -file => $infile)
			  or print STDERR "can't open $infile:$!\n"
			     and next;
    FASTA_SEQ:
    while(my $fasta_translation=$seqin->next_seq()) {
#	if($geneidcode) {
#	    $geneid=geneid($fasta_translation);
#	}

	my $id=$fasta_translation->$idwhence;
	#$id =~ s/P(\d+)/T$1/;
	#$id = $get_id_func->($id);


	if($idregex) {
	    if($id=~$idregex) { 
	        $id=$1;
	    } else {
	        print STDERR "No id in $id\n";
		++$count{'not found'};
		next;
	    }
	}elsif( $id =~ /GSVIV[A-Z](\d+)/ ){ #this is the grape ids
	  $id = $grape_id_prefix{uc($idwhat)}.$1;
	}
	 
	  
	
	print "id = $id\n";

	my @GeneScriptLation;   #each [$geneobject, $transcriptobject,
			#    $transcript->translate, translated sequence];
			#last two filled in later because ->translate can fail
			# (if screwed up db loading enough)
	#if find via gene, this will include all transcripts of genes found

	if ($idxref) {
	    my @gene_ids=$dbe_adaptor->list_gene_ids_by_extids($id); 

	    if(@gene_ids) {
		for my $gii (@gene_ids) {
		   my $gene= $gene_adaptor->fetch_by_dbID($gii);
	           push @GeneScriptLation, map { [ $gene,$_ ] }
		      @{$gene->get_all_Transcripts};
	        }
	    } else {
	        print STDERR "xref $id not found\n";
		++$count{'not found'};
		next;
	    }
	} elsif ($idwhat eq 'gene') {
	    my $gene;
	    eval {  $gene=$gene_adaptor->fetch_by_stable_id($id); };
	    print STDERR "$@\n" and ++$count{'not found'} and next if $@;
	    print STDERR "no gene fetched by $id\n" and next unless $gene;
	    @GeneScriptLation= map { [$gene,$_] }
		      @{$gene->get_all_Transcripts};
	} elsif ($idwhat eq 'transcript') {
	    my $transcript; #print "id = $id\n";
	    $id =~ s/\_P/_T/; print "nid = $id\n";#protein id is 'Zm00001e019034_P004', transcript id is 'Zm00001e019034_T004'
	    eval {  $transcript=$transcript_adaptor->fetch_by_stable_id($id); };
	    print STDERR "$@\n" and ++$count{'not found'} and next if $@;
	    print STDERR "no transcript fetched by $id\n" and ++$count{'not found'} and next unless $transcript;
	    @GeneScriptLation= ( [ $gene_adaptor->fetch_by_transcript_id( $transcript->dbID)
				,$transcript ]);
	} else { #must be translation
	    my $transcript;
	    eval {  $transcript=$transcript_adaptor->fetch_by_translation_stable_id($id); };
	    print STDERR "$@\n" and ++$count{'not found'} and next if $@;
	    @GeneScriptLation= ([ $gene_adaptor->fetch_by_transcript_id( $transcript->dbID)
				,$transcript ]);
	}
		
	my (%seq_regions,$strand,%ens_seq);
	
	my $fa_seq=lc($fasta_translation->seq);
	#$fa_seq=~ s/\*$//; #TIGR puts * at end for stop codon always,
			   #ensembl doesn't
	$fa_seq =~ s/[.*x]+$//; #JGI sometimes put X when the aa cannot be determined for
			  #example: a hanning 'A', or . or * for stop codon
	$fa_seq =~ s/\./*/g; #internal stop codon, we want to unify them to be represented by *;
	#print "fa_seq=$fa_seq\n";

	foreach my $gtt (@GeneScriptLation) {
	    my ($gene,$trans) =@$gtt;
	    my $pep;
	    eval { $pep=$trans->translate; };
	    $@ and ++$count{'translate failed'} 
	       and print $trans->stable_id,"=$id translate failed: $@\n" 
	       and die; #and next FASTA_SEQ;   #weix
	    my $trpt_id=$trans->stable_id;
#warn("DEBUG: trpt_id=$trpt_id");
	    $gtt->[2]=$pep;
	    $gtt->[3]=lc($pep->seq());
	    #print "pep_seq\n".$gtt->[3]."\n";
	    ++$count{'ok'.scalar(@GeneScriptLation)} and next FASTA_SEQ
	     if( $fa_seq  eq $gtt->[3] );
	    print "$fa_seq\n  eq\n $gtt->[3]\n";
	   

	    print "$trpt_id=> $gtt->[3] (ensembl_translation)\n";
            print "$trpt_id=> $fa_seq (fasta_seq)\n";
 
	    my ($uniform_fa_seq, $uniform_ensembl) = ($fa_seq, $gtt->[3]);
	    
	    $uniform_ensembl =~ s=[*]+.*==;
	    $uniform_fa_seq =~ s=[*.]+.*==; #some file has * or . to indicate stop codon

	    ++$count{'ok'.scalar(@GeneScriptLation)} and next FASTA_SEQ
	      if( $uniform_fa_seq  eq $uniform_ensembl );

	    print ">$trpt_id (uniform_ensembl)\n$uniform_ensembl\n>$trpt_id (uniform_fa_seq)\n$uniform_fa_seq\n";
	    my $ori_len_ens = length($gtt->[3]) ;
            my $trunc_len_ens = length($uniform_ensembl);
            my $ori_len_fa = length($fa_seq) ;
            my $trunc_len_fa = length($uniform_fa_seq);
            print "Truncate at 1st stop codon to get uniformed sequences, internal stop codon for ensembl: $ori_len_ens -> $trunc_len_ens\n" if ( $ori_len_ens != $trunc_len_ens);
            print "Truncate at 1st stop codon to get uniformed sequences, internal stop codon for fasta: $ori_len_fa -> $trunc_len_fa\n" if ( $ori_len_fa != $trunc_len_fa );
	   ++$count{'mismatch-internal-stop'} if ($ori_len_ens != $trunc_len_ens or $ori_len_fa != $trunc_len_fa) ;
        }

	#no transcript of this gene matches -need to report 
	#error

	my @allstrands=();
	my @fa_seq=unpack('C*',$fa_seq);
	my %fa_mismatch=();

	foreach my $gtt (@GeneScriptLation) {
	    my ($gene,$trans,$pep,$ens_seq)=@$gtt;
	    print "gene ",$gene->stable_id,", transcript ",$trans->stable_id,", translation ",$pep->id,"\n";
	    ## Get a list of exons for this transcript
	    my @exons=@{$trans->get_all_Exons};
	    my $lationstart=$trans->translation->start;

	    my ($forward,$reverse);
	    foreach my $exon (@exons){
		{  #make sure to report top & bottom slices this is on
		 my $slice=$slice_adaptor->fetch_by_Feature($exon);
		     #$exon->slice might be whole chromosome. this is just the exon's footprint
		 $seq_regions{$slice->seq_region_name}=1;
		 for my $cs ($rank1_coord_system, $seq_coord_system) {
		     if($slice->coord_system_name ne $cs) {
		         for my $seg ( @{$slice->project($cs)} ) {
			     $seq_regions{$seg->to_Slice->seq_region_name}=1;
			 }
		     }
		 }
		}

		$forward=1 if $exon->strand==1;
		$reverse=1 if $exon->strand==-1;
	    }
	    $strand = $forward ? $reverse ? 'forward & reverse strands'
					  : 'forward strand'
			       : 'reverse strand';
	    push @allstrands,$strand;
	    my @ens_seq=unpack('C*',$ens_seq);

	    my ($mismatchposn,$mismatchcount,$lastmismatch)=find_mismatch(\@ens_seq,\@fa_seq);
	    $fa_mismatch{$mismatchposn} =$fa_mismatch{$lastmismatch} =1;

	    my $remaining=$mismatchposn*3+$lationstart-1;
	    my $mm_exon='?';
	    my $structure=join (''
			      , map { $_->stable_id." (".$_->length.") phase "
				     .$_->phase." to ".$_->end_phase."\n" 
				     .join("\n",$_->seq->seq =~ /(.{1,60})/g) ."\n"
				    } @exons );
	    foreach my $exon (@exons) {
		$remaining -= $exon->length;
		if($remaining<0) { $mm_exon=$exon->stable_id; $remaining+=$exon->length; last; }
	    }
	    print "transcript ",$trans->stable_id,", on: ",join("+",sort keys %seq_regions),", $strand  $mismatchcount mismatches, the first at $mismatchposn ($mm_exon:$remaining) :\nEnsEMBL (tx start $lationstart)\n".show_mismatch($ens_seq,$mismatchposn,$lastmismatch)."\n$structure\n";
	    if($mismatchposn<=$#ens_seq) {
		print "first mismatch codon ",substr($trans->seq->seq,3*$mismatchposn,3),"\n"; 
		#do_3tx($trans->seq->seq,$lationstart,$fa_seq,$mismatchcount);
	        if($mismatchposn==length($fa_seq)) {
		    #Simple: Ensembl seq is too long
		    #want to output --adjust id=-n for  adjust_translation_end.pl
		    #n should be transcript length-translation start- 3*desired peptide length 
		    $fixes.='--adjust '.$trans->stable_id.'='
			  .(length($fa_seq)*3+$lationstart-$trans->seq->length)."\n";
		}
	    }
	} # end for each transcript
	print "Fasta ",$fasta_translation->id,"\n".show_mismatch($fa_seq,keys %fa_mismatch)."\n\n" ;
	$count{"mismatch ".join("/",sort @allstrands)}++;
    }
}

print "\n",map { "$_: $count{$_}\n" } sort keys %count;

print "\n$fixes";


# --------- Subroutines ------------


sub do_3tx {
    my ($dna_seq,$previous_start,$desired_peptide,$previous_mismatch)=@_;
#    print join(",",length($dna_seq),length($desired_peptide),$previous_mismatch),"\n" if $verbose_count;
    my $pad=' ' x (60-length($dna_seq)%60);
    my @seqs=($dna_seq.$pad);
    my $besttx=undef;
#    my $bestmiscnt=$previous_mismatch;
    my $bestmatches=length($desired_peptide)-$previous_mismatch
                    +abs(length($dna_seq)/3-length($desired_peptide));
    print "bestmatches $bestmatches\n";
    my ($best,$bestslimis);
    for my $i (0..($previous_start+1)) {	# so 0..2 if start is 1
	my $tx=lc($stdCodonTable->translate(substr($dna_seq,$i)));
	#compare $tx and $desired_peptide
	my ($slimis)=@{sliding_mismatch($tx,$desired_peptide)};
	   #--only using one of the best matches below, so why bother 
	my $matches=$slimis->[1];
#	print " $i,$firstmis,$miscnt,",substr($tx,0,50),"\n" if $verbose_count;
	if($matches>$bestmatches) {
	    $best=$i+1; #0-based -> 1-based
	    $bestmatches=$matches;
	    $besttx=$tx;
	    $bestslimis=$slimis;
	}
	push @seqs,(' ' x $i).join('  ',split('', $tx)).$pad if($i<=2); #'if' since only need to show 3 frames
    }
    print "bestmatches $bestmatches\n";
    my $bestmiscnt=length($desired_peptide)-$bestmatches
		+abs(length($dna_seq)/3-length($desired_peptide));
    print "bestmiscnt $bestmiscnt\n";
    if($besttx) {	#this is set only if found something better
	if($best == $previous_start) {
	  print "Translation of transcript from $best to end yields with offset ",$bestslimis->[0],", $bestmiscnt mismatches between ",$bestslimis->[2]," and ",$bestslimis->[3],":\n"
	    ."Adjust translation end!\n";
	} else {
	  print "Translation start at $best instead of $previous_start offset "
	      ,$bestslimis->[0],": $bestmiscnt mismatches, from "
	      ,$bestslimis->[2],"to",$bestslimis->[3],":\n"
	    .show_mismatch($besttx,$bestslimis->[2],$bestslimis->[3])."\n";
	  print "wanted:\n"
	    .show_mismatch($desired_peptide,$bestslimis->[2]+$bestslimis->[0]
				    ,$bestslimis->[3]+$bestslimis->[0])."\n"
	    if $bestmiscnt;
	}
    }
    if($verbose_count && $bestmiscnt>1) { #Don't bother to be verbose if we solved the problem
	print length($dna_seq),"\n" ;
	for my $i (0.. int( (length($dna_seq)-1)/60 )) {
	    print map { substr($_,$i*60,60)."\n" } @seqs;
	}
	print "\n";
	$verbose_count--;
    }
}

sub find_mismatch {
    my ($tx,$peptide)=@_; #really no difference. could be called seq1 and seq2
    ref $tx or $tx=[ unpack('C*',$tx) ];
    ref $peptide or $peptide=[ unpack('C*',$peptide) ];
    my ($mismatchposn,$mismatchcount,$lastmismatch);
    my $end= $#$peptide<$#$tx? $#$peptide:$#$tx;
    for($mismatchposn=0;$mismatchposn<=$end;$mismatchposn++) { 
	last if $tx->[$mismatchposn] ne $peptide->[$mismatchposn];
    }
    for(my $i=$mismatchposn,$mismatchcount=0;$i<=$#$peptide && $i<=$#$tx;$i++) { 
	++$mismatchcount and $lastmismatch=$i if $tx->[$i] ne $peptide->[$i];
    }
    $lastmismatch=$mismatchposn unless defined $lastmismatch;
    $mismatchcount += abs(scalar(@$peptide)-scalar(@$tx));
    return ($mismatchposn,$mismatchcount,$lastmismatch);
}

sub show_mismatch {
    my ($seq,@mm)=@_;
    $seq=join("\n",$seq =~ /(.{1,60})/g)."//";;
    cluck('empty mm') unless @mm;
    cluck('undef in mm') if grep { !defined $_ } @mm;
    
    my $order = 0;
    foreach my $mm (sort { $b<=>$a } @mm) {
	$mm+=int($mm/60); #puts mark at the beginning of a line instead
			  #  of the end of the previous line, ok?
	$mm+=2*$order; #should we count the <> we inserted for prevous mm
	$mm=length($seq) if $mm>length($seq);
	$seq=substr($seq,0,$mm).'<>'.substr($seq,$mm);
	++$order;
    }
    $seq;
}



sub sliding_mismatch { #returns array of arrayref. ref is to offset, number of matches
    my ($seqa,$seqb)=@_;
    ref $seqa or $seqa=[ unpack('C*',$seqa) ];
    ref $seqb or $seqb=[ unpack('C*',$seqb) ];
    #offset=offsets for b relative to a
    my %offsets= ( 0 => 1, 1=> 1, -1 => 1, 2=>1, -2=>1, 3=>1, -3=>1);
    my $offset;
    for $offset (keys %offsets) { $offsets{$offset+$#$seqb-$#$seqa}=1 };
    my $mostmatches=-1;
    my @results;	# list of [ offset,number of matches,first mismatch,last mismatch]
    for $offset(keys %offsets) {
        my($firstmis,$lastmis,$matches)=(-1,-1,0);
	my $start=$offset<0?-$offset:0;
	my $end=$#$seqa<$#$seqb-$offset?$#$seqa:$#$seqb-$offset;
	for (my $i=$start; $i<=$end;$i++) {
	    if($seqa->[$i] eq $seqb->[$i+$offset]) {
	        $matches++;
	    } else {
	        $firstmis ||= $i;
		$lastmis=$i;
	    }
	}
	if($matches>$mostmatches) {
	    @results=([$offset,$matches,$firstmis,$lastmis]);
	    $mostmatches=$matches;
	} elsif($matches==$mostmatches) {
	    push @results,[$offset,$matches,$firstmis,$lastmis];
	}
    }
    print "sliding_mismatch(",scalar(@$seqa),",",scalar(@$seqb),")\n"
	  ,(map { "    ".join(",",@$_)."\n" } @results);
#	  ,(map { "    ".join(",",ref($_)?("[".join(",",@$_)."]"):$_)."\n" } @results);
    return \@results;
}


__END__

=head1 OUTPUT

Info on mismatches


=item B<Standard Output>

Mainly for debugging

=head1 NOTES
    


=head1 AUTHOR

   Steven Schmidt
   Gramene Project (www.gramene.org)
   Stein Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut


