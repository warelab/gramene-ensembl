#!/usr/local/bin/perl


use lib map {  "/data/ware/weix/gramene_ensembl/$_" }
            qw ( ensembl-live/ensembl/modules
		 bioperl-live );
        

use strict;
use Pod::Usage;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

use Bio::SeqIO;
use Bio::Seq;
use Date::Calc;

=head1 SYNOPSIS

ensembl_template.pl  [options] 
 
 Options:

    --help              help message
    --man               full documentation
    --species           ensembl species, default is oryza_sativa
    --registry_file     Default is $GrameneEnsemblDir/conf/ensembl.registry
    --coord_system      the coord_system of the dumped sequences
    --repeat_mask       soft mask the sequences
    --list_file         the file containing list of seq_region_names
    --list_cs           the coord_system the listed names are on

=head1 OPTIONS

=over 4

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

=item B<--species>

ensembl species, default is oryza_sativa

=item B<--registry_file>

Set the registry file for the ensembl database, default is $GrameneEnsemblDir/conf/ensembl.registry 

=item B<--coord_system>

The coord system the sequences to be dumped are on, for example chromosome, clone, contig 

=item B<--repeat_mask>

Soft Repeat mask the sequnences


=item B<--list_file>

the file contain a list of seq_region_names

=item B<--list_cs>

the coord system the listed names are on such as clone

=back

=cut


my ($registry_file, $help, $man, $coord_system, $species, $repeat_mask);
my ($list_file, $list_cs);

GetOptions( "help|?"          => \$help,
            "man"             => \$man,
            "species=s"       => \$species,
            "registry_file=s" => \$registry_file,
            "coord_system=s"  => \$coord_system,
	    "repeat_mask"     => \$repeat_mask,
	    "list_file=s"     => \$list_file,
	    "list_cs=s"       => \$list_cs,
          )
  or pod2usage(2);

pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;
pod2usage("Need both list_file and list_cs") 
    if ($list_file && !$list_cs ||
	$list_cs && !$list_file) ;

$registry_file ||= "$ENV{GrameneEnsemblDir}/conf/ensembl.registry";
$species ||= "Oryza_sativa";

###
# Load database adaptors

Bio::EnsEMBL::Registry->load_all( $registry_file );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( warn( "No core DB for $species set in $registry_file\n" ) &&
              pod2usage(1) );

###
# 
my $rundate = sprintf( '%4.4i%2.2i%2.2i', Date::Calc::Today );
my $file_basename = "$species-$coord_system-$rundate";

my $seq_io = Bio::SeqIO->new
  ( -format=>'fasta',
    -file=> $repeat_mask ? ">${file_basename}RM.fa" : ">$file_basename.fa",
  );


###
# Get object adaptors
# For more methods, see online API doc at
# http://www.ensembl.org/info/software/Pdoc/ensembl/index.html



my $sla = $ENS_DBA->get_adaptor('slice');

if( $list_file ){

    open my $fh, $list_file or die "Cannot open $list_file";
    my @seq_region_names = <$fh>;
    chomp @seq_region_names;

    for my $name(@seq_region_names){
	
	my $seq_slice=$sla->fetch_by_region($list_cs, $name);
	my $projected_segs = $seq_slice->project( $coord_system);

	for my $seg( @{$projected_segs} ){
	    my $slice = $seg->to_Slice();

	    my $seqstr;

	    if($repeat_mask){ 
		$seqstr = $slice->get_repeatmasked_seq(['RepeatMask'], 1)->seq;
		$seqstr ||= $slice->seq;
	    }else{
		$seqstr = $slice->seq;
	    }

	    my $seq = Bio::Seq->new(
				    -display_id => $slice->seq_region_name,
				    -seq        => $seqstr,
				    );
	    $seq_io->write_seq( $seq );
	}
    }
    
}else{

    foreach my $slice( @{$sla->fetch_all($coord_system)} ){

	my $seqstr;
	if($repeat_mask){ 
	    $seqstr = $slice->get_repeatmasked_seq(['RepeatMask'], 1)->seq;
	    $seqstr ||= $slice->seq;
	}else{
	    $seqstr = $slice->seq;
	}
	
	my $seq = Bio::Seq->new(
				-display_id => $slice->seq_region_name,
				-seq        => $seqstr,
				);
	$seq_io->write_seq( $seq );
	
    }
}


exit; 
