#!/usr/bin/perl -w

#transfer/dump genes. The source genes/transcripts can be all 
#transcripts in query db or use IDs from an input file

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';
use vars qw[ $VERSION ];
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaDnaAlignFeature;



$VERSION = sprintf "%d.%02d", q$Revision: 1.10 $ =~ /(\d+)\.(\d+)/;
pod2usage(0) unless $ARGV[0];

my ($_help, $_version, $sdbhost, $sdbport, $sdbuser, $sdbpass, $sdbname, 
    $sdnadb, $tdbhost, $tdbport, $tdbuser, $tdbpass, $tdbname, $idfile, $type, 
    $istrans, $coord_system, $addid, $idbase, $keepid, $addscore, $newtype,
    $glen, $plen, $write, $usechr, $baselen, $logicname);
GetOptions(
	   'h|help'        => \$_help,
	   'v|version'     => \$_version,
	   'dbhost=s'      => \$sdbhost,
	   'dbport=s'      => \$sdbport,
	   'dbuser=s'      => \$sdbuser,
	   'dbpass=s'      => \$sdbpass,
	   'dbname=s'     => \$sdbname,
	   'dnadb=s'       => \$sdnadb,
	   'tdbhost=s'      => \$tdbhost,
	   'tdbport=s'      => \$tdbport,
	   'tdbuser=s'      => \$tdbuser,
	   'tdbpass=s'      => \$tdbpass,
	   'tdbname=s'     => \$tdbname,
	   'logic_name=s'   => \$logicname,

	   ) or die;
pod2usage(2) if $_help;
if ( $_version ) {
    my $prog = basename( $0 );
    print "$prog v$VERSION\n";
    exit(0);
}

$sdbport ||= '3306';
$tdbhost ||= $sdbhost;
$tdbuser ||= $sdbuser;
$tdbport ||= $sdbport;
$tdbpass ||= $sdbpass;
$tdbpass ||= $sdbname;

print STDERR "Target: $tdbhost $tdbport $tdbname\n";

my $sdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					       -user   => $sdbuser,
					       -pass   => $sdbpass,
					       -dbname => $sdbname,
					       -host   => $sdbhost,
					       -port => $sdbport,
					       -driver => 'mysql',
					       );
my $tdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					       -user   => $tdbuser,
					       -pass   => $tdbpass,
					       -dbname => $tdbname,
					       -host   => $tdbhost,
					       -port => $tdbport,
					       -driver => 'mysql',
					       );

$coord_system = 'toplevel';
print "coord_system is $coord_system, logicname is $logicname\n";

my $slice_adaptor = $tdba->get_SliceAdaptor;

my $transcripts;

my $repeat_features = get_all_repeats($sdba,  $coord_system, $logicname);

print "repeat_features to copy: ", scalar @$repeat_features, "\n";


my $rcdba = $tdba->get_RepeatConsensusAdaptor();
my $rfdba = $tdba->get_RepeatFeatureAdaptor();
my $analdba = $tdba->get_AnalysisAdaptor();

my %repeat_adaptor_hash = (repeat_consensus => $rcdba,
			   repeat_features => $rfdba);
my %repeat_consensus_hash;
my %analysis_hash;

my %repeat_load_hash;
my $num=0;

foreach my $repeat_f ( @$repeat_features){
    
    my $repeat_c = $repeat_f->repeat_consensus;
    
    my $rcid = $repeat_c->dbID;
    my $new_repeat_c = $repeat_consensus_hash{$rcid};
    unless($new_repeat_c){
	$rcdba->store($repeat_c);
	$repeat_consensus_hash{$rcid} = $repeat_c;	 
    }
    
    $repeat_f->repeat_consensus($repeat_c);
    
    

    my $anal = $repeat_f->analysis;
    my $analid = $anal->dbID;
    my $new_anal = $analysis_hash{$analid};
    unless( $new_anal ){
	$analdba->store($anal);
	$analysis_hash{$analid} = $anal;
    }
    
    $repeat_f->analysis( $anal );

    $rfdba->store( $repeat_f);
    $num++;

    
    print "$num repeat_features copied\n";
    #last;
    
}


exit;



###########################################################


# retrieve all transcripts on rice all chromosomes
sub get_all_repeats {
    my ($dba, $coord_system, $logicname) = @_;
    
    my $slice_adaptor = $dba->get_SliceAdaptor;
    my @repeat_features;
    my @slices = @{$slice_adaptor->fetch_all($coord_system)};
    print STDERR scalar @slices, " $coord_system(s)\n";
    
    for my $slice (@slices) {
	my $repeats = $slice->get_all_RepeatFeatures($logicname);
#	print scalar @$repeats, " repeat regions\n";

	push @repeat_features,  @$repeats;
    }
    
       
    return \@repeat_features;
}


# ----------------------------------------------------

=head1 NAME

repeat-copy.pl

=head1 SYNOPSIS

  repeat-copy.pl -sdbhost host -sdbport dbport -sdbuser user -sdbname name 
    -sdnadb dnadb -tdbhost host -tdbport port -tdbuser user -tdbname name 
    -logicname logic_name -coord_system coord

Options:

  -h|--help        Show brief help and exit
  -v|--version     Show version and exit
  -dbhost=s         source db host
  -dbport=s        source db port
  -dbuser=s         db user 
  -dbname=s        db for query genes
  -dnadb=s          db for dna seq 
  -tdbhost=s        target  db host
  -sdbport=s        target db port
  -tdbuser=s         db user 
  -tdbname=s        db for target genes
  -logic_name        logic name for the repeat analysis
  
=head1 DESCRIPTION

To copy genes from one db to another 

=head1 SEE ALSO

perl.

=head1 AUTHOR

Sharon Wei E<lt>weix@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2005 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

