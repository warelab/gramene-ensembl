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

cdna_check_fasta.pl  [options] fasta-file [ ... fasta-file]
 
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
my $dbe_adaptor=$reg->get_adaptor($ensembl_species,'core','DBEntry');
my $coord_system_adaptor=$reg->get_adaptor($ensembl_species,'core','CoordSystem');

print "rank1=$rank1\n";
my $rank1_coord_system=$coord_system_adaptor->fetch_by_rank($rank1)->name;
my $seq_coord_system=$coord_system_adaptor->fetch_by_name('seqlevel')->name;
#print "rank coord_system_ $rank1_coord_system\n $seq_coord_system\n";

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
	    eval {  $transcript=$transcript_adaptor->fetch_by_stable_id($id); };
	    print STDERR "$@\n" and ++$count{'not found'} and next if $@;
	    print STDERR "no transcript fetched by $id\n" and next unless $transcript;
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
	$fa_seq=~ s/\*$//; #TIGR puts * at end for stop codon always,
			   #ensembl doesn't
	$fa_seq=~ s/x$//; #JGI sometimes put X when the aa cannot be determined for
			  #example: a hanning 'A'
	#print "fa_seq=$fa_seq\n";

	foreach my $gtt (@GeneScriptLation) {
	    my ($gene,$trans) =@$gtt;
	    my $cdna;
	    eval { $cdna=$trans->seq; };
	    $@ and ++$count{'cdna failed'} 
	       and print $trans->stable_id,"=$id cdna failed: $@\n" 
	       and next FASTA_SEQ;
	    my $trpt_id=$trans->stable_id;

	    $gtt->[2]=$cdna;
	    $gtt->[3]=lc($cdna->seq);
	    ++$count{'ok'.scalar(@GeneScriptLation)} and next FASTA_SEQ
	     if( $fa_seq  eq $gtt->[3] );
	    #print "$fa_seq\n  eq\n $gtt->[3]\n";
	    
        }

	#no transcript of this gene matches -need to report 
	#error

	foreach my $gtt (@GeneScriptLation) {
	    my ($gene,$trans,$cdna,$cdna_lc)=@$gtt;
	    print ">gene ",$gene->stable_id,", transcript ",$trans->stable_id."\n$cdna_lc\n";
	    print ">Fasta ",$fasta_translation->id,"\n$fa_seq\n";
	    
       }
    }
}

print "\n",map { "$_: $count{$_}\n" } sort keys %count;

print "\n$fixes";


# --------- Subroutines ------------


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


