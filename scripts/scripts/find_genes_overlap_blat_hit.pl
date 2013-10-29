#!/usr/local/bin/perl


use lib map { $ENV{HOME} . "/gramene_ensembl/$_" }
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

#For maize, we want the working-set genes
use constant ATTRIB_CODE => 'working-set';  


=head1 NAME
    
find_genes_overlap_blat_hit.pl 
    
=head1 SYNOPSIS

find_genes_overlap_blat_hit.pl  [options] blat_mapping_file1 blat_mapping_file2 ...
 
 Options:

    --help              help message
    --man               full documentation
    --species           ensembl species
    --registry_file     the file store connections info to ensembl databases
    --cs                the coordinate system the sequnece name is on 
    --index             the colum number of the seq name
    --start_index       the colum number of the start coordinate

=head1 OPTIONS

=over 4

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

=item B<--species>

ensembl species in the registry file

=item B<--registry_file>

Set the registry file for the ensembl database

=item B<--cs>

The coord system the mappings or blat hits are on

=item B<--index>

The colum index of the target seq region name, for example, 0,1,2,...

=item B<--start_index>

The colum index of the hit start coordinate on the target seq region, 
for example, 0,1,2,...

=back

=head1 DESCRIPTION

find genes overlaping with the hits of query sequences mapped on some genome coordinate system

=cut



my ($registry_file, $help, $man, $coord_system, 
    $species, $index,  $start_index);

GetOptions( "help|?"          => \$help,
            "man"             => \$man,
            "species=s"       => \$species,
            "registry_file=s" => \$registry_file,
            "cs=s"            => \$coord_system,
	    "index=i"         => \$index,
	    "start_index=i"   => \$start_index,
          )
  or pod2usage(2);

pod2usage(-verbose => 2) if ($man || !$index || !$start_index 
			     || !$species || !$registry_file || !$coord_system);
pod2usage(1) if $help;

###
# Load database adaptors, get the db adaptor for the species of study
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( warn( "No core DB for $species set in $registry_file\n" ) &&
              pod2usage(1) );

###
# Get slice adaptor
my $sla = $ENS_DBA->get_adaptor('slice');


###
#
for my $f (@ARGV){

    
    open my $fh, $f or die "Cannot open $f to read";

    while (my $mapping = <$fh>){
	
	chomp( $mapping );
	next if $mapping =~ /HIT_ID/i;

	#print "$mapping\n";
	my @data = split /\t/, $mapping;
	my $seq_region_name = $data[$index];
	$seq_region_name =~ s/Chr\.?\s*//i;
	my $strand = $data[$index+1];
	my $contig_start = $data[$start_index];
	next unless ($contig_start =~ /^\d+$/ && $contig_start != 0);
	my $contig_end   = $data[$start_index+1];

	#print "$coord_system, $seq_region_name\n";
	
	my $slice = $sla->fetch_by_region($coord_system, $seq_region_name, 
					  $contig_start, $contig_end, $strand);
	
	

	my $gene_listref = $slice->get_all_Genes;

	if( @{$gene_listref} > 0){
	    for my $gene( @{$gene_listref}){
		my $gene_stable_id = $gene->stable_id;
		my $attrib_lstref = $gene->get_all_Attributes( ATTRIB_CODE );
		if (scalar @{$attrib_lstref} > 0){
		    
		    print "$mapping\t$gene_stable_id\tworking-set\n";
		    
		}
		#else{
		#    print "$mapping\t$gene_stable_id\n";
		#}
	    }
	}else{
	
	    print "$mapping\n";
	}
	
    }
  }  

exit;

__END__

