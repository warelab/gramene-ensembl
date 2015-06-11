package BlastIndex;
use base qw/Bio::EnsEMBL::Production::Pipeline::PipeConfig::FASTA_conf/;
sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },
    #Override of options

        'registry' => '/usr/local/ensembl-79/gramene-live/grmblast/conf/Reg.pm', # default option to refer to Reg.pm, should be full path
        'base_path' => '/usr/local/ensembl-79/blastdb/', #where do you want your files
        

        # Species to run; use this to restrict to subset of species but be aware you
        # are still open to the reuse checks. If you do not want this then use
        # force_species
        #species => ['Arabidopsis_thaliana', 'Oryza_sativa'],
        
        dump_types => ['dna', 'cdna'],
 
        ### Defaults

        pipeline_name => 'fasta_dump_'.$self->o('release'),

        wublast_exe => 'xdformat',
        ncbiblast_exe => 'makeblastdb',
        blat_exe => 'faToTwoBit',
        port_offset => 30000,

        email => $self->o('ENV', 'USER').'@cshl.edu',

 
  };
}
1;
