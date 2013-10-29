#!/lab/bin/perl 

=head1 NAME

dump_go_protein.pl 

=cut


BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'EnsemblDir'} ||= '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'EnsemblDir'}."/ensembl-live/$_" } 
        qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-variation/modules ensembl-draw/modules ensembl-compara/modules );

use strict;
use warnings;
use Pod::Usage;
use POSIX qw(strftime);
use Bio::EnsEMBL::Registry;
use Getopt::Long;



=head1 SYNOPSIS

dump_go_protein.pl  [options] 
 
 Options:
    --species		species in EnsEMBL Web conf to use for db [required]
    --registry_file	ensembl registry file,
    --help		help message
    --man		full documentation

=head1 OPTIONS

=over 4

=item B<--species> 

A species in the EnsEMBL web config in $GrameneEnsemblDir/conf


=item B<--registry_file>

registry file

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

NONE

=cut






my ($ensembl_species, $registry_file,$DUMPDIR);
    {  #Argument Processing
	my $help=0;
	my $man=0;
	GetOptions( "help|?"=>\$help,
		    "man"=>\$man,
		    "species=s" => \$ensembl_species,
		    "registry_file=s" => \$registry_file,
		  )
	  or pod2usage(2);
	pod2usage(-verbose => 2) if $man;
        pod2usage(1) if $help;
        $ensembl_species or pod2usage(2);
    }

# Load database adaptors
Bio::EnsEMBL::Registry->load_all( $registry_file );
my $db=Bio::EnsEMBL::Registry->get_DBAdaptor($ensembl_species,'core','slice') or die "can't get  adaptor for $ensembl_species";

#my $db=$sa->db;	#Bio::Ensembl::DBAdaptor

my $out_fh;
my $out_file = $ensembl_species.'_go_protein.txt';
open($out_fh,">$out_file") or die "$!";




my $dbc = $db->dbc;

my $sql = 
    q[
	select trs.stable_id,x.description,x.dbprimary_acc,gx.linkage_type
	from   xref x,
	       object_xref ox,
	       go_xref gx,
	       translation_stable_id trs,
	       translation tr,
	       external_db db
	where  
	       trs.translation_id = tr.translation_id
	and    ox.xref_id = x.xref_id
	and    ox.object_xref_id = gx.object_xref_id
	and    x.external_db_id = db.external_db_id
	and    db.db_name        = 'GO'
    ];

my $sql_translation = q[ 
	and    ox.ensembl_id     = tr.translation_id
	and    ox.ensembl_object_type = 'Translation'
    ];



print $out_fh join("\t",qw(!db db_obj_id db_obj_symbol qualifier term_acc dbxref evn_code with aspect obj_name obj_synonym obj_type species date assign)),"\n";



my $sth = $dbc->prepare($sql.$sql_translation);
$sth->execute() or die $dbc->errstr;

my %seen_assocs;
while(my ($gene,$go_name,$go_acc,$evn) = $sth->fetchrow_array()){
    my $assoc = join("\t",($gene,$go_name,$go_acc,$evn));
    unless($seen_assocs{$assoc}){
	&print_GOC_file($out_fh,$gene,$go_name,$go_acc,$evn);
	$seen_assocs{$assoc}=1;
    }
}



close($out_fh);
$sth->finish;

sub print_GOC_file{

    my ($out_fh,$gene,$go_name,$term_acc,$evn_code) = @_;

    ## fields for GOC format requirements
    my $db = 'GR_Ensembl';
    my $with = '';
    my $qualifier = '';
    my $date = strftime "%Y%m%d", localtime;
    my $assign='GR';
    my $obj_type='Ensembl gene';
    my $dbxref = 'GR_REF:8396';
    my $species = 'taxon:xxxx';  # dummy param.
    #my %term_type_to_aspect = ('2'=>'P','3'=>'C','4'=>'F');
    my $aspect='NA'; # dummy param.
    ###


    my ($db_obj_id,$db_obj_symbol,$obj_name) = $gene;
    my $obj_synonym = '';

    print $out_fh join("\t",(
			    $db,
			    $db_obj_id,
			    $db_obj_symbol,
			    $qualifier,
			    $term_acc,
			    $dbxref,
			    $evn_code,
			    $with,
			    $aspect,
			    $obj_name,
			    $obj_synonym,
			    $obj_type,
			    $species,
			    $date,
			    $assign
			)),"\n";


}

exit;
    
# end of main program

########################## subroutines ######################################

__END__

=head1 OUTPUT

NONE

=head1 NOTES
    


=head1 AUTHOR

   (@cshl.edu)
   Gramene Project (www.gramene.org)
   Stein Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

