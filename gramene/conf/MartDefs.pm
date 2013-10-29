
#!/usr/local/bin/perl -w
#======================================================================
#   
#   Name:        MartDefs.pm
#   
#   Description: module to create/store/retrieve martview config data 
#                using the EnsemblMart::MartInfo module
#
#======================================================================

package MartDefs;
use strict;
use Storable qw( freeze thaw );

use Data::Dumper;

use EnsemblMart::MartInfo;
use EnsemblMart::LocationAdaptor;
use SiteDefs;
use SpeciesDefs;

use vars qw( @ISA @EXPORT );
use Exporter;
@ISA     = qw(Exporter);
@EXPORT  = qw( is_filter_available 
	       all_species 
	       all_foci
	       species_by_focus
	       foci_by_species 
	       assembly_info
	       strains_available
	       traits_available
	       encode_available);

#----------------------------------------------------------------------
# Some class variables
my $CONFFILE = "${SiteDefs::ENSEMBL_SERVERROOT}/conf/martconf.packed";
my $CONF = undef();

#----------------------------------------------------------------------
# Prepare the config

=head2 new

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub new{

  # If config is stored, just deserialise 
  if( -e $CONFFILE ){
    open FH, $CONFFILE || die ("Can't open $CONFFILE: $!" );
    local $/ = undef;
    $CONF = thaw(<FH>);
    
    warn( ( '-' x 78 ) ."\n" );
    warn("[MARTCONF][INFO] Retrieving conf from $CONFFILE\n" );
    warn( ( '-' x 78 ) ."\n" );
    #warn Dumper( $CONF );
    return 1;
  }

  # No config found; create new
  
  # Get database handle
  my $dbh = SpeciesDefs::db_connect_multi_species( 'ENSEMBL_MART' ) || 
    warn( "[MARTCONF][WARN] Could not connect to mart database \n" ) && return 0;

  # Get MartInfo instance, and use to get all species+foci
  my $mart_info     = EnsemblMart::MartInfo->new( $dbh );
  my $mart_location = EnsemblMart::LocationAdaptor->new( $dbh );
  my @species   = @{$mart_info->species_available()};
  my @foci      = @{$mart_info->foci_available() };

  $CONF->{all_species} = \@species;
  $CONF->{all_foci}    = \@foci;
  $CONF->{all_filters}    = [];
  $CONF->{all_attributes} = [];
  $CONF->{foci_by_species} = {};
  $CONF->{species_by_focus} = {};
  $CONF->{chromosomes_by_species} = {};
  $CONF->{bands_by_species} = {};
  $CONF->{filters} = {};
  $CONF->{attributes} = {};
  $CONF->{strains} = {};
  $CONF->{traits} = {};
  $CONF->{trait_regions_by_species} = {};
  $CONF->{encode} = {};
  $CONF->{encode_regions_by_species} = {};

  # Loop for each species
  foreach my $sp( @species ){
    
    # Update mart_info with the species
    $mart_info->species($sp);
    $mart_location->species($sp);
    
    my $assembly = $mart_info->retrieve_assembly($sp);
    $CONF->{assembly} ->{$sp} = $assembly;
    
    #my @strains = $mart_info->strains_available($sp);
    $CONF->{strains} ->{$sp} = $mart_info->strains_available($sp);
    
    $CONF->{traits} ->{$sp} = $mart_info->traits_available($sp);
    $CONF->{encode} ->{$sp} = $mart_info->encode_available($sp);

    # Loop for each focus
    foreach my $fo( @{$mart_info->foci_available} ){
      
      # Get a list of all available filters for the species+focus
      my $filters = $mart_info->filters_available($fo);
      my $attribs = $mart_info->attributes_available($fo);
      
      # Add each filter to the master config
      if( @$filters or @$attribs ){
        $CONF->{foci_by_species} ->{$sp} ||= {}; 
        $CONF->{species_by_focus}->{$fo} ||= {};
        $CONF->{foci_by_species} ->{$sp}->{$fo} = 1;
        $CONF->{species_by_focus}->{$fo}->{$sp} = 1;
      }

      # Save typing!
      my $fconf = $CONF->{filters};
      foreach( @$filters ){
	$fconf->{$_} ||= {};
	$fconf->{$_}->{$sp} ||= {};
	$fconf->{$_}->{$sp}->{$fo} = 1;
      }
      my $aconf = $CONF->{attributes};
      foreach( @$attribs ){
	$aconf->{$_} ||= {};
	$aconf->{$_}->{$sp} ||= {};
	$aconf->{$_}->{$sp}->{$fo} = 1;
      }
    }

    # Update chromosomes by species
    $CONF->{chromosomes_by_species}->{$sp} = $mart_info->chrs_available;
    $CONF->{bands_by_species}->{$sp} = {$mart_location->gen_chr_band_hash};
    $CONF->{trait_regions_by_species}->{$sp} = {$mart_location->gen_trait_regions_hash};
    $CONF->{encode_regions_by_species}->{$sp} = {$mart_location->gen_encode_regions_hash};
  }
  # Get a sorted list of all filters
  $CONF->{all_filters}    = [ sort keys %{$CONF->{filters}} ];
  $CONF->{all_attributes} = [ sort keys %{$CONF->{attributes}} ];

  warn("[MARTCONF][INFO] Writing Configuration to $CONFFILE");
  warn( ( '-' x 78 ) ."\n" );

  # Save away the conf
  open( FH, "> $CONFFILE" ) || die ("Can't open $CONFFILE for writing: $!" );
  print FH ( freeze($CONF) );
}

#----------------------------------------------------------------------

=head2 all_species

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub all_species {
  if( ! $CONF ){ new() };
  return @{$CONF->{all_species}};
}
#----------------------------------------------------------------------

=head2 assembly_info

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub assembly_info {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{assembly}->{$species} || return ();
  return $CONF->{assembly}->{$species};
}
#----------------------------------------------------------------------

=head2 strains_available

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub strains_available {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{strains}->{$species} || return [];
  return @{$CONF->{strains}->{$species}};
}
#----------------------------------------------------------------------

=head2 traits_available

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub traits_available {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{traits}->{$species} || return [];
  return @{$CONF->{traits}->{$species}};
}
#----------------------------------------------------------------------

=head2 encode_available

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub encode_available {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{encode}->{$species} || return [];
  return @{$CONF->{encode}->{$species}};
}
#---------------------------------------------------------------------
=head2 all_foci

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub all_foci {
  if( ! $CONF ){ new() };
  return @{$CONF->{all_foci}};
}

#----------------------------------------------------------------------

=head2 all_filters

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub all_filters {
  if( ! $CONF ){ new() };
  return @{$CONF->{all_filters}};
}

#----------------------------------------------------------------------

=head2 all_attributes

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub all_attributes {
  if( ! $CONF ){ new() };
  return @{$CONF->{all_attributes}};
}

#----------------------------------------------------------------------

=head2 all_names

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub all_names {
  if( ! $CONF ){ new() };
  return( @{$CONF->{all_filters}}, @{$CONF->{all_attributes}} );
}




#----------------------------------------------------------------------

=head2 species_by_focus

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub species_by_focus {
  my $focus = shift;
  if( ! $CONF ){ new() };
  $CONF->{species_by_focus}->{$focus} || return ();
  return sort keys %{$CONF->{species_by_focus}->{$focus}};
}

#----------------------------------------------------------------------

=head2 foci_by_species

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub foci_by_species {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{foci_by_species}->{$species} || return ();
  return sort keys %{$CONF->{foci_by_species}->{$species}};
}

#----------------------------------------------------------------------

=head2 chromosomes_by_species

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub chromosomes_by_species {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{chromosomes_by_species}->{$species} || return ();
  return @{$CONF->{chromosomes_by_species}->{$species}};
}

#----------------------------------------------------------------------


=head2 bands_by_species

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub bands_by_species {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{bands_by_species}->{$species} || return ();
  return %{$CONF->{bands_by_species}->{$species}};
}

#----------------------------------------------------------------------


=head2 trait_regions_by_species

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub trait_regions_by_species {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{trait_regions_by_species}->{$species} || return ();
  return %{$CONF->{trait_regions_by_species}->{$species}};
}

#----------------------------------------------------------------------


=head2 encode_regions_by_species

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub encode_regions_by_species {
  my $species = shift;
  if( ! $CONF ){ new() };
  $CONF->{encode_regions_by_species}->{$species} || return ();
  return %{$CONF->{encode_regions_by_species}->{$species}};
}
#----------------------------------------------------------------------

=head2 is_focus_available

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub is_focus_available {
  if( ! $CONF ){ new };
  my $species = shift;
  my $focus   = shift;
  my $species_conf = $CONF->{foci_by_species}->{$species} || return 0;
  return $species_conf->{$focus} ? 1 : 0;
}

#----------------------------------------------------------------------

=head2 is_filter_available

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub is_filter_available {
  if( ! $CONF ){ new }
  my $name    = shift;
  my $species = shift;
  my $focus   = shift;
  my $name_conf = $CONF->{filters}->{$name} || return 0;
  my $species_conf = $name_conf->{$species} || return 0;
  return $species_conf->{$focus} ? 1 : 0;
}

#----------------------------------------------------------------------

=head2 is_attribute_available

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub is_attribute_available {
  if( ! $CONF ){ new }
  my $name    = shift;
  my $species = shift;
  my $focus   = shift;
  my $name_conf = $CONF->{attributes}->{$name} || return 0;
  my $species_conf = $name_conf->{$species} || return 0;
  return $species_conf->{$focus} ? 1 : 0;
}





#----------------------------------------------------------------------
1;
