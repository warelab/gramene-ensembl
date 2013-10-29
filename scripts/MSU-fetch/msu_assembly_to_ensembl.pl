#!/usr/local/bin/perl


=head1 NAME

tigr_assembly_to_ensembl.pl - 

=cut


BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'EnsemblDir'} ||= '/usr/local/ensembl-live/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/gramene/scripts/ensembl/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
use lib map { $ENV{'EnsemblDir'}."/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);


use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;

#you may not want all of these:
use Gramene::Config;
use DBI;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Registry;
use Text::RecordParser::Tab;

=head1 SYNOPSIS

tigr_assembly_to_ensembl.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry_file
    --species		species in EnsEMBL Web conf to use for db
    
=head1 OPTIONS

=over 4

=item B<--registry_file>

The registry file contain DB connection information for species core databases

=item B<--species>

The species name in the conf file, species.ini

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

list of MSU agp files
    
Chr   Chr_start      Chr_end  clone BAC_start BAC_end  ori

Chr1    1       31687   OSJNOa264G09    1       31687   +
Chr1    31688   175927  P0672D08        29490   173729  +
Chr1    175928  302974  P0436E04        17598   144644  +
Chr1    302975  438996  P0005A05        54813   190834  +
Chr1    438997  542578  P0482C06        79198   182779  +
Chr1    542579  661386  P0439B06        17184   135991  +


=cut

my %clone2acc;
my $ensembl_species=$ENV{ENSEMBL_SPECIES};
my $registry_file;
my $clone2acc_file;
    
{  #Argument Processing
  
  my $help=0;
  my $man=0;
  GetOptions( "help|?"          =>\$help,
	      "man"             =>\$man,
	      "species=s"       => \$ensembl_species,
	      "registry_file=s" => \$registry_file,
	      "id_file=s"       => \$clone2acc_file,
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  pod2usage(-verbose => 2) if !$registry_file;
  
}

pod2usage(-verbose => 2) unless @ARGV ;
print "@ARGV\n";
#exit;

###
### Get Ensembl DB adaptor
###

$ENV{'ENSEMBL_SPECIES'} = $ensembl_species;

my $reg = "Bio::EnsEMBL::Registry";	#Use this to get adaptors 
$reg->load_all( $registry_file ) or die "load_all failed";
     
my $slice_adaptor =$reg->get_adaptor($ensembl_species,'core','Slice') 
  or die "can't get Slice adaptor for $ensembl_species";

my $attribute_adaptor =$reg->get_adaptor($ensembl_species,'core','Attribute')
     or die "can't get Attribute adaptor for $ensembl_species";

my $coord_system_adaptor = $reg->get_adaptor($ensembl_species,'core'
                                     ,'CoordSystem');

     
my $dba=$slice_adaptor->db; #DBAdaptor
my $dbc=$dba->dbc;	#DBConnection
warn "user ".$dbc->username
  .", db ".$dbc->dbname."\n";


#my $ens_dbh=$dbc->db_handle; #you can use this as a DBI database handle

my $assembly_sth = $dbc->prepare("insert into assembly 
      (asm_seq_region_id ,cmp_seq_region_id ,asm_start ,asm_end ,cmp_start ,cmp_end ,ori)
      values  (?, ?,?,?,?,?,?)");

#Get clone - accession mapping
#
my %clonename2acession;
if ( -s $clone2acc_file > 0 ){

  print "found $clone2acc_file\n";
  open my $fh, $clone2acc_file or die "cannot open $clone2acc_file";
  
 
  while( my $idmap = <$fh> ){

    chomp$idmap;
    my @ids = split /\t/, $idmap;
    
    $clonename2acession{ uc($ids[0]) } = $ids[1];
    
    
  }
}

###
### Get TIGR ASSEMBLY data
###
#

my %chrs;
my $p = new Text::RecordParser::Tab;


for my $file (@ARGV){
  
 # if( $file =~ m= [^/0-9]* 0* (\d+) \. [^/]*\z =xms ){
 #   $chr = $1;
 #   print STDERR "Process $file for chromosome $chr\n";
 # }else{
 #   die "Cannot parse $file for chromosome";
 # }
  
  print "file=$file\n";
  
  $p->filename($file);
  $p->bind_fields( qw[chr chr_start chr_end clone_name bac_start bac_end ori] );

  my $records = $p->fetchall_arrayref( { Columns => {} } );
  
 
  for my $r( @{$records} ){
    
    my $clone_name = uc ($r->{'clone_name'});
    my $accession;

    print "$clone_name => $clonename2acession{ $clone_name }\n";

    if( $clone_name =~ /^Syng/i){
      #these are Syngenta contigs
      $accession = $clone_name;
    }else{
      $accession = $clonename2acession{ $clone_name } ?
                 $clonename2acession{ $clone_name }: $clone_name;
    }

    #print "$clone_name , $accession\n";

    my $chr = $r->{'chr'};
    $chr =~ s/chr(omosome)?0*//i;

    $chrs{$chr} = $slice_adaptor->fetch_by_region('chromosome',$chr) unless $chrs{$chr};

    my $chr_slice = $chrs{$chr} ;
    $chr_slice or die "Cannot get slice for chr $chr";

    
    if( $r->{'chr_start'} < 1 || $r->{'chr_start'} > $r->{'chr_end'}
	|| $r->{'chr_end'} > $chr_slice->end) {
      
      print STDERR "bad coordinates ",$r->{'chr_start'}," to ",$r->{'chr_end'}," accession=$accession'}", " clone=",$r->{'clone_name'}," on chrom=$chr\n";
      next;
    }
    

    my $clone_slice = $slice_adaptor->fetch_by_region('clone', $accession);
    
    die "Can't found the clone using $accession" unless $clone_slice;
    
    
    my ( $cmp_start, $cmp_end )=( $r->{'bac_start'}, $r->{'bac_end'} );
    die "Component coordinates cannot be null $chr, $accession\n" unless ($cmp_start && $cmp_end);

    my $ori = $r->{'ori'} eq '+' ? 1:-1;

    if( $cmp_start > $cmp_end ) {
      die "Wrong assumptions about orientations";

    }
    
    print "chr_seq_region_id = " . $chr_slice->get_seq_region_id . "\n";
    $assembly_sth->execute($chr_slice->get_seq_region_id, $clone_slice->get_seq_region_id
			   ,$r->{'chr_start'}, $r->{'chr_end'}
			   ,$cmp_start ,$cmp_end ,$ori);

    
  }

}

exit;
    
# end of main program


=head1 OUTPUT

=item B<Standard Output>

=head1 NOTES
    
Assumes TIGR xml files have been loaded in same server as
core ensembl db for species
and are accessible by the same user.

=head1 AUTHOR


   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

