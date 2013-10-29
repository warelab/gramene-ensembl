package EnsEMBL::Web::UserConfig::Vmapview;
use strict;
use EnsEMBL::Web::UserConfig;
use vars qw(@ISA);
@ISA = qw(EnsEMBL::Web::UserConfig);

sub init {
  my ($self) = @_;
  $self->{'_label'}           = 'above',
  $self->{'_band_labels'}     = 'on',
  $self->{'_image_height'}    = 450,
  $self->{'_top_margin'}      = 40,
  $self->{'_band_links'}      = 'yes',
  $self->{'_userdatatype_ID'} = 109;
  
  $self->{'general'}->{'Vmapview'} = {
      '_artefacts'   => [qw(Vmarkers Videogram Vaccbacs Vssrs Vclones Vsupercontigs)],# Vsupercontigs)],
#      '_artefacts'   => [qw(Vtigr Vgenes Vtos17 Vmarkers Videogram Vfpcmap)], # gramene artefacts
      '_options'   => [],
      
      '_settings' => {
	  'width'       => 500, # really height <g>
	  'bgcolor'     => 'background1',
	  'bgcolour1'   => 'background1',
	  'bgcolour2'   => 'background1',
      },
      
      # Maize specific mapview density configs

      'Vmarkers' => { # Maize-specific                                             
	  'on'          => 'on',
	  'pos'         => '95',
	  'width'       => 60,
	  'logicname'   => 'Density_cornsensus',
	  'track_labels'=> {'Density_cornsensus'         => 'Maize Overgos'},
	  'colours'     => {'Density_cornsensus'         => 'red'},
      },
      
      'Vaccbacs' => {
	  'on'           => 'on',
	  'pos'          => '97',
	  'width'        => 60,
	  'logicname'    => 'Density_acc_bac_map',
	  'track_labels' => {'Density_acc_bac_map' => 'Accessioned BACs'},
	  'colours'      => {'Density_acc_bac_map' => 'blue'},
      },
      
      'Vclones' => {
	  'on'           => 'on',
	  'pos'          => '96',
	  'width'        => 60,
	  'logicname'    => 'Density_bac_map',
	  'track_labels' => {'Density_bac_map' => 'FPC Clones' },
	  'colours'      => {'Density_bac_map' => 'black' },
      },
      
      
      
      
      'Vssrs' => {
	  'on'           => 'on',
	  'pos'          => '94',
	  'width'        => 60,
	  'logicname'    => 'Density_ssr_marker',
	  'track_labels' => {'Density_ssr_marker' => 'Electronic SSRs'},
	  'colours'      => {'Density_ssr_marker' => 'purple'},
      },
      
      
      
      'Videogram' => {
	  'on'          => "on",
	  'pos'         => '1000',
	  'width'       => 24,
	  'bandlabels'  => 'on',
	  'totalwidth'  => 100,
	  'col'         => 'g',
	  'padding'     => 6,
      },
      
      'Vsupercontigs' => {
	  'on'          => 'on',
	  'pos'         => '400',
	  'width'       =>  20,
	  'totalwidth'  =>  100,
	  'padding'     =>  6,
	  'col'         => 'blue',
	  'col_ctgs1'    => 'blue',
	  'col_ctgs2'    => 'darkgreen',
	  'lab'         => 'black',
	  'include_labelling' => 1,
	  'available'   => 'features mapset_core_bins',
      }
  }
}

### The following configurations are all gramene specific.
### Kept for reference.

#    'Vtigr' => { # Maize-specific                                                               
#      'on'          => 'off',                                                                   
#      'pos'         => '90 ',                                                                   
#      'width'       => 60,                                                                      
#      'col_genes'   => 'black',                                                                 
#      'col_gc'      => 'red',                                                                   
#      'logicname'   => 'PercentGC Density_tigr_gene',                                           
#      'track_labels'=> {'PercentGC'         => '% CG',                                          
#                        'Density_tigr_gene' => 'TIGR Genes' },                                  
#      'colours'     => {'PercentGC'         => 'red',                                           
#                        'Density_tigr_gene' => 'black' },                                       
#    },                                                                                          
#    'Vmarkers' => { # Maize-specific                                                            
#      'on'          => 'off',                                                                   
#      'pos'         => '95',                                                                    
#      'width'       => 60,                                                                      
#      'logicname'   => 'Density_Rice_SSR Density_Rice_rflp_marker',                             
#      'track_labels'=> {'Density_Rice_SSR'         => 'SSR Markers',                            
#                        'Density_Rice_rflp_marker' => 'RFLP Markers' },                         
#      'colours'     => {'Density_Rice_SSR'         => 'red',                                    
#                        'Density_Rice_rflp_marker' => 'green' },                                
#    },                                                                                 

#    'Vtos17' => { # Maize-specific
#      'on'          => 'off',
#      'pos'         => '97',
#      'width'       => 60,
#      'logicname'   => 'Density_Rice_tos17_insert Density_Rice_T_DNA_Insert',
#      'track_labels'=> {'Density_Rice_tos17_insert' => 'TOS17 Inserts',
#                        'Density_Rice_T_DNA_Insert' => 'T-DNA Inserts'},
#      'colours'     => {'Density_Rice_tos17_insert' => 'purple',
#                        'Density_Rice_T_DNA_Insert' => 'magenta'},
#    },
#    'Vfpcmap' => { # Maize-specific 
#      'on'          => 'off',
#      'pos'         => '98',
#      'width'       => 60,
#      'logicname'   => 'Density_bac_map Density_acc_bac_map',
#      'track_labels'=> {'Density_bac_map'     => 'BACs',
#                        'Density_acc_bac_map' => 'Acc BACs'},
#      'colours'     => {'Density_bac_map'     => 'black',
#                        'Density_acc_bac_map' => 'red'},
#    },

#    'Vgenes' => {
#      'on'          => 'off',
#      'pos'         => '100',
#      'width'       => 60,
#      'col_genes'   => 'black',
#      'col_xref'    => 'rust',
#      'col_pred'    => 'black',
#      'col_known'   => 'rust',
#      'logicname' => 'knownGeneDensity geneDensity'
#    },
#    'Vrefseqs' => {
#      'on'          => 'off',
#      'pos'         => '110',
#      'width'       => 60,
#      'col'         => 'blue',
#      'logicname' => 'refseqs'

#    },        
#    'Vpercents' => {
#      'on'          => 'off',
#      'pos'         => '200',
#      'width'       => 60,
#      'col_gc'      => 'red',
#      'col_repeat'  => 'black',
#      'logicname' => 'PercentageRepeat PercentGC'
#    },    
#    'Vsnps' => {
#      'on'          => 'off',
#      'pos'         => '300',
#      'width'       => 60,
#      'col'         => 'blue',
#      'logicname' => 'snpDensity'
#    },        

1;
