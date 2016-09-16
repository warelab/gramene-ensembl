#!/lab/bin/perl -w

BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/'; 
}

# The first shall be last...
#use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

#use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
#        qw ( bioperl-live modules ensembl/modules conf
#	     ensembl-external/modules ensembl-draw/modules
#	     ensembl-compara/modules );

use lib map { "$ENV{'GrameneEnsemblDir'}/ensembl-live/$_" } 
        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
             ensembl-draw/modules ensembl-compara/modules );

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;

#use Gramene::Config;
use DBI;
use DBI qw(:sql_types);

use Data::Dumper qw(Dumper); # For debug

use FindBin qw( $Bin );
use File::Basename qw( dirname );


# Import EnsEMBL modules
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

=head1 SYNOPSIS

load_dna_align_feature.pl  [options] 
 
 Options:
    --species		species in EnsEMBL registry to use for db [required]
    --v 		makes more verbose
    --help		help message
    --man		full documentation
    --coord_system	coord system to use
    --registry_file	Default is $GrameneEnsemblDir/conf/ensembl.registry
    --logic_name	Allowed logic name
    --replace		Remove all hits for each track being loaded


=head1 OPTIONS

=over 4

=item B<--species> 

A species in the EnsEMBL registry

=item B<--coord_system>

Coordinate system name (e.g. 'clone', 'chromosome' )

Defaults to the sequence level coordinate system 

Can be over-ridden by data file: if region name is chr_NUMBER
perhaps followed by _NUMBER then chromosome coordinates are used
and the 1st NUMBER is the chromosome name.

=item B<--logic_name>

Allowed logic name (can be repeated).

If none given then all are allowed.

=item B<--replace>

Before loading a track, remove all existing hits

If you manage to load several files into the same track,
the remove happens only once.


=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

files to process (produced by Gramene Pipeline)

filename is prefix+logic_name+'.whatever'

=cut

use vars qw($ENS_DBA $species $logic_name  $coord_system
		$replace $verbose);
$verbose=0;

{  #Argument Processing
  my $help=0;
  my $man=0;
  my $registry_file;
  GetOptions
      ( 
        "help|?"          => \$help,
        "man"             => \$man,
        "species=s"       => \$species,
	"coord_system=s"  => \$coord_system,
        "registry_file=s" => \$registry_file,
        "logic_name=s"    => \$logic_name,
        "replace"         => \$replace,
        "v+"		  => \$verbose,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Validate the input files
  $registry_file    ||= $ENV{GrameneEnsemblDir}.'/conf/ensembl.registry';
  @ARGV || ( warn( "Need some input files to process\n" ) && pod2usage(1) );
  map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
  } $registry_file, @ARGV;
  warn( "Found ",scalar @ARGV," files to process...\n" );

  # Load the ensembl file
  $species || ( warn( "Need a --species\n" ) && pod2usage(1) );
  Bio::EnsEMBL::Registry->load_all( $registry_file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( warn( "No core DB for $species set in $registry_file\n" ) &&
                pod2usage(1) );

}


my $slice_adaptor    = $ENS_DBA->get_adaptor('Slice');
my $analysis_adaptor = $ENS_DBA->get_adaptor('Analysis');
my $feature_adaptor  = $ENS_DBA->get_adaptor('DnaAlignFeature');

unless ($coord_system) {
    my $coord_system_adaptor=  $ENS_DBA->get_adaptor('CoordSystem');
    $coord_system=$coord_system_adaptor->fetch_sequence_level->name;
    warn("$coord_system coordinates\n");
}


my %replaced;

foreach my $file( @ARGV ){
  warn( "Processing $file\n" );
  warn ($logic_name) if $verbose;

  my $analysis = &fetch_analysis($logic_name); 
  if($replace && !$replaced{$logic_name}++ && $analysis->dbID) {
      print "replace is $replace, analysisID is ", $analysis->dbID if $verbose;
      $feature_adaptor->db->dbc->do(
       "delete from dna_align_feature where analysis_id=".$analysis->dbID)
       or die("delete $logic_name = ".$analysis->dbID.": $DBI::errstr");
  } 
  
  open( FEATURES, $file ) || die ("Cannot open $file: $!");
  while( my $line = <FEATURES> ){
    next if ($line =~ /^\#/ || $line =~ /^\s*$/);
    chomp $line;
    
    warn ($line) if $verbose;

    #1       genomic_organellar_insert_annotation    nuclear_mitochondrial   3820191 3820421 100     +       .       ID=32_cluster_24;Target=mitochondria genome 352597 352827;evalue=2E-31;length_orgDNA=85
    my( $chr, 
        $src, 
        $type, 
        $chr_start, 
        $chr_end, 
	$score, 
	$chr_strand,,
        $phase, 
        $attribs ,
	 ) = split /\t/, $line;
    
    my %attributes = map{$_=~ s/^\s+//; $_=~s/\s$//; $_ =~ s/\s+/_/;  split /\=/; }(split /;/, $attribs);
#print "score is $score\n";
#next;
    for (keys %attributes){
	
	print "$_ => $attributes{$_}\n";
    }

    $chr =~ s/^chr(omosome)?_?//;
    warn "$coord_system:$chr\n" if $verbose;

    my $slice = $slice_adaptor->fetch_by_region($coord_system,$chr) ||
    ( warn( "  Cannot fetch $coord_system $chr - skipping\n" ) && next );

    #my $cigar = ( $seq_region_end-$seq_region_start+1 ) . "M";

    
    my ($hit_name, $hit_start, $hit_end);

    if($attributes{Target} =~ /(\w+)\s(\d+)\s(\d+)/){
	($hit_name, $hit_start, $hit_end) = ($1, $2, $3);
	delete $attributes{Target};
    }elsif( defined $attributes{ID} ){
	$hit_name =  $attributes{ID};
	$hit_start = 1;
	$hit_end = $chr_end-$chr_start+1;
	delete $attributes{ID};
    }else{
	warn("ERROR: No hit name can be determined from attribute 'Target' or 'ID' \n");
	next;
    }

    my $evalue = undef;
    if( defined $attributes{evalue} ){  
	my $e = $attributes{evalue};
	$evalue = $e =~ /E/i ? sprintf("%.10g", $e) : sprintf ("%.10f", $e);
    }
    delete $attributes{evalue};

    my $cigar_string = $chr_end-$chr_start+1;

    my $external_data = join ';', map{ "$_=$attributes{$_}" } (keys %attributes);
    print "exteral_data=$external_data\n" if $verbose;

    my $strand = $chr_strand =~ /\+/ ? 1 : -1;
    my $feature = Bio::EnsEMBL::DnaDnaAlignFeature->new
        ( 
          -start        => $chr_start,
          -end          => $chr_end,
          -strand       => $strand,
          -hstart       => $hit_start,
          -hend         => $hit_end,
          -hstrand      => 1,
          -hseqname     => $hit_name,
          -analysis     => $analysis,
          -slice        => $slice,
          -cigar_string => $cigar_string,
          -p_value 	=> $evalue,        #key cannot be evalue
          -score        => $score,
	  -extra_data => $external_data,   #key cannot be external_data
          -percent_id   => $score, 
          );

    $feature_adaptor->store($feature);
    warn( "    Created feature $hit_name at ",
          $feature->feature_Slice->name , "\n") if $verbose;

  }


}

exit;

    
# end of main program

########################## subroutines ######################################

#======================================================================
# Returns an Analysis object; either fetched from the DB if one exists,
# or created fresh, in which case it is stored to the DB.
sub fetch_analysis{
  my $logic_name = shift || die("Need a logic_name" );
  my $db_file    = shift || '';
  my $adaptor = $ENS_DBA->get_adaptor('Analysis'); 

  my %args = ( -logic_name=>$logic_name, 
               $db_file ? (-db_file=>$db_file) : () );

  # Hard-coded nastyness to make analyses correspond to Ensembl
  if( $logic_name eq 'MyAnalysis' ){
    $args{-logic_name} = 'MyLogic';
    $args{-db}         = 'MyDB';
    $args{-program}    = 'MyProgram';
  }

  my $analysis;
  if( $analysis = $adaptor->fetch_by_logic_name($args{-logic_name}) ){
    # Analysis found in database already; use this.
    return $analysis;
  }

  # No analysis - create one from scratch
  $analysis = Bio::EnsEMBL::Analysis->new(%args);
  warn("created analysis for $logic_name.\n");
  return $analysis;
}

#======================================================================
1;

__END__

=head1 OUTPUT

Explanatory only

=head1 NOTES
    
Creates entry in analysis table if doesn't exist


=head1 AUTHOR

   Will Spooner, Kiran Kumar, Steven Schmidt
   Gramene Project (www.gramene.org)
   Stein Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

