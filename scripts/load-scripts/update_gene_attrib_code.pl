#!/usr/local/bin/perl -w


BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/';
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/';
}


use lib map { $ENV{GrameneEnsemblDir} . "/$_" }
            qw ( ensembl/modules
		 ensembl-compara/modules
		 bioperl-live );

use strict;
use Bio::EnsEMBL::Registry;
use Bio::AlignIO;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Data::Dumper qw(Dumper);
use Pod::Usage;

=head1 NAME

    update_gene_attrib_code.pl set gene attrib code for a list of genes, this originated from the need to identify only the high confident genes for sorghum db provided we have a list of high confident gene names, you could also apply it to transcript by specifying -type transcript

=head SYNOPSIS

    update_gene_attrib_code.pl [OPTIONS] files

=head1 OPTIONS


=over 4

=item B<--reg_conf>

 the registry file for connecting to ensembl db 

=item B<--species>
  
 the species name in the registry file

=item B<--attrib_code>

  the gene attrib code you give this set of genes

=item B<--type>

  gene or transcipt, default is gene

=back


=head1 DESCRIPTION


  perl  /usr/local/gramene/scripts/ensembl/load_scripts/update_gene_attrib_code.pl -reg_conf ~/registries/bhsqldw.registry -species sorghum_bicolor ~/data/compara_stats/sorghum_high_confidence.gid 


=head1 AUTHORS
    Sharon Wei

=head1 COPYRIGHT





=cut


my $help;
my $reg_conf;

my $species;
my ($attrib_code, $attrib_value);
my $type='gene';

GetOptions(
	   "help"                 => \$help,
	   "reg_conf=s"           => \$reg_conf,
	   "species=s"            => \$species,
	   "attrib_code=s"        => \$attrib_code,
	   "attrib_value=s"       => \$attrib_value,
	   "type=s"               => \$type,
	   );

$type = lc $type;
pod2usage( verbose => 2) if $help;
pod2usage( verbose => 2) if ($type ne 'gene' and $type ne 'transcript' );

my $reg = "Bio::EnsEMBL::Registry";

$reg->no_version_check(1);
$reg->load_all($reg_conf);

my $DBAdaptor = $reg->get_DBAdaptor($species, 'core');
throw("Registry configuration file has no data for connecting to core db <$species>")
    if ( !$DBAdaptor );

my $dbh = $DBAdaptor->dbc->db_handle;

my %sql;


$sql{get_attrib_type} = qq[
select attrib_type_id from attrib_type where code = ?
];

$sql{insert_attrib_type} = qq[
insert into attrib_type (code, name) values(?, ?)
];


$sql{get_gene_id} = qq[
SELECT   ${type}_id
From     ${type}_stable_id
WHERE    stable_id = ?
];

$sql{insert_gene_attrib} = qq[
insert into ${type}_attrib (${type}_id, attrib_type_id, value) values(?, ?, ?)
];



my @gene_ids;
for my $file (@ARGV){
  print "process file $file\n";
  open my $fh, $file or die "can not open $file to read";
    
  my $gsth = $dbh->prepare($sql{get_gene_id});

  my ($line, $gsid, $tag_value);
  while ( $line = <$fh> ){
    chomp $line;
    ($gsid, $tag_value) = split ' ', $line;
    $tag_value ||= $attrib_value;
    $gsth->execute( $gsid );

    my ($gid) = $gsth->fetchrow_array;
    push @gene_ids, [$gid, $tag_value];

    #print "get $gid for $gsid\n";
    
  }
}


my $attrib_type_id;
my $attrib_type_sth = $dbh->prepare( $sql{get_attrib_type} ) or 
  die "cannot create sth for  $sql{get_attrib_type}";
$attrib_type_sth->execute( $attrib_code );
($attrib_type_id) = $attrib_type_sth->fetchrow_array ;

unless( $attrib_type_id ){
  my $attrib_code_insert_sth = $dbh->prepare($sql{insert_attrib_type});
  my $rv=$attrib_code_insert_sth->execute( $attrib_code,  $attrib_code);
  die "failed insertion into attrib_type" unless $rv;
}
$attrib_type_sth->execute( $attrib_code );
($attrib_type_id) = $attrib_type_sth->fetchrow_array ;


my $insert_gene_attrib_sth = $dbh->prepare($sql{insert_gene_attrib}) or 
die "cannot prepare $sql{insert_gene_attrib}";

for my $gidref( @gene_ids){
  my ($gid, $tag_value) = @$gidref;  
  print "insert $gid, $attrib_type_id, $tag_value\n";
  my $rv = $insert_gene_attrib_sth->execute( $gid, $attrib_type_id,  $tag_value);
  die "failed insertion into ${type}_attrib" unless $rv;
  
}

__END__
