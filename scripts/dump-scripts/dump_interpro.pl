#!/usr/local/bin/perl 


=head1 NAME

dump_interpro.pl 

=cut


BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'EnsemblDir'} ||= '/usr/local/ensembl-live/'; 
}



use lib map { $ENV{'EnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-variation/modules ensembl-draw/modules ensembl-compara/modules );



use strict;
use warnings;
use Pod::Usage;
use POSIX qw(strftime);
use Bio::EnsEMBL::Registry;
use Getopt::Long;



=head1 SYNOPSIS

dump_interpro.pl  [options] 
 
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






my ($ensembl_species, $registry_file,$interpro_file,$DUMPDIR);
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
my $out_file = $ensembl_species.'_interpro.txt';
open($out_fh,">$out_file") or die "$!";

my $dbc = $db->dbc;


my %interpros;
&get_interpros($dbc,\%interpros);

my $sql = 
    q[
	select    ts.stable_id,
              i.interpro_ac,
              pf.hit_name,
              pf.evalue,
              pf.perc_ident,
              pf.seq_start, 
              pf.seq_end
	from      protein_feature pf,
              translation_stable_id ts,
              interpro i
        where  pf.translation_id = ts.translation_id
         and  pf.hit_name = i.id
    ];

print $out_fh join("\t",qw(!protein_id Interpro_id Interpro_name domain evalue
perc_ident seq_start seq_end)),"\n";



my $sth = $dbc->prepare($sql);
$sth->execute() or die $dbc->errstr;

my %seen_assocs;
while(my @data = $sth->fetchrow_array()){
    my $protein = shift @data;
    my $ip_acc  = shift @data;
    my $assoc = join("\t",($protein,$ip_acc));
    unless($seen_assocs{$assoc}){
        my $ip_name = $interpros{$ip_acc};
        if($ip_name){
           print $out_fh join("\t",($protein,$ip_acc,$ip_name,@data)),"\n";
           $seen_assocs{$assoc}=1;
        }else{
           print "$ip_acc no Name\n";
       }
    }
}



close($out_fh);
$sth->finish;



sub get_interpros{
  my ($dbc,$interpros)=@_;
  my $sth_ip = $dbc->prepare(q[
    SELECT dbprimary_acc, x.description 
    FROM xref x, external_db e
    where x.external_db_id = e.external_db_id
    and e.db_name = 'Interpro';
  ]);
   
  $sth_ip->execute();
  while(my ($ipr,$name) = $sth_ip->fetchrow_array()){
       $interpros->{$ipr}=$name;
  }
  $sth_ip->finish;
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

