package Bio::EnsEMBL::GlyphSet::tigr_transcript;
use strict;
use vars qw(@ISA);
use Carp qw(cluck);
use EnsWeb;
use Bio::EnsEMBL::GlyphSet_transcript;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end eprof_dump);

@ISA = qw(Bio::EnsEMBL::GlyphSet_transcript);

sub my_label {
    my $self = shift;
warn("Returning tigr_transcript label");
    return "TIGR Genes";
    #return $self->{'config'}->{'_draw_single_Transcript'} || ( EnsWeb::species_defs->AUTHORITY.' trans.');

}

sub colours {
    my $self = shift;
    my $Config = $self->{'config'};
    return $Config->get('transcript_lite','colours');
}

sub features {
  my ($self) = @_;
  
  return $self->{'container'}->get_all_Genes('tigr_gene');	#lc(EnsWeb::species_defs->AUTHORITY));
}


sub colour {
    my ($self, $gene, $transcript, $colours, %highlights) = @_;

    my $genecol = $colours->{ "_".$transcript->external_status };

    if(exists $highlights{$transcript->stable_id()}) {
        print STDERR "\n\n\a================\nhas transcript stable id returning $colours->{'superhi'}\n";
      return ($genecol, $colours->{'superhi'});
    } elsif(exists $highlights{$transcript->external_name()}) {
        print STDERR "\n\n\a================\nhas transcript external name returning $colours->{'superhi'}\n";
      return ($genecol, $colours->{'superhi'});
    } elsif(exists $highlights{$gene->stable_id()}) {
        print STDERR "\n\n\a================\nhas gene stable id returning $colours->{'hi'}\n";
      return ($genecol, $colours->{'hi'});
    }
      
    return ($genecol, undef);
}

sub href {
    my ($self, $gene, $transcript, %highlights ) = @_;

    my $gid = $gene->stable_id();
    my $tid = $transcript->stable_id();
    
    return ( $self->{'config'}->get('transcript_lite','_href_only') eq '#tid' && exists $highlights{$gene->stable_id()} ) ?
        "#$tid" : 
        qq(/$ENV{'ENSEMBL_SPECIES'}/geneview?gene=$gid);

}

sub zmenu {
    my ($self, $gene, $transcript) = @_;
    my $tid = $transcript->stable_id();
    my $pid = $transcript->translation->stable_id(),
    my $gid = $gene->stable_id();
    my $id   = $transcript->external_name() eq '' ? $tid : ( $transcript->external_db.": ".$transcript->external_name() );
    my $zmenu = {
        'caption'                       => EnsWeb::species_defs->AUTHORITY." Gene",
        "00:$id"			=> "",
	"01:Gene:$gid"                  => "/$ENV{'ENSEMBL_SPECIES'}/geneview?gene=$gid&db=core",
        "02:Transcript:$tid"    	        => "/$ENV{'ENSEMBL_SPECIES'}/transview?transcript=$tid&db=core",                	
        '04:Export cDNA'                => "/$ENV{'ENSEMBL_SPECIES'}/exportview?tab=fasta&type=feature&ftype=cdna&id=$tid",
        
    };
    
    if($pid) {
    $zmenu->{"03:Peptide:$pid"}=
    	qq(/$ENV{'ENSEMBL_SPECIES'}/protview?peptide=$pid&db=core);
    $zmenu->{'05:Export Peptide'}=
    	qq(/$ENV{'ENSEMBL_SPECIES'}/exportview?tab=fasta&type=feature&ftype=peptide&id=$pid);	
    }
    
    my $DB = EnsWeb::species_defs->databases;

    if($DB->{'ENSEMBL_EXPRESSION'}) {
      $zmenu->{'06:Expression information'} = 
	"/$ENV{'ENSEMBL_SPECIES'}/sageview?alias=$gid";
    }

    return $zmenu;
}

sub text_label {
    my ($self, $gene, $transcript) = @_;
    my $tid = $transcript->stable_id();
    my $id  = ($transcript->external_name() eq '') ? 
      $tid : $transcript->external_name();

    if( $self->{'config'}->{'_both_names_'} eq 'yes') {
        return $tid.(($transcript->external_name() eq '') ? '' : " ($id)" );
    }

    return $self->{'config'}->{'_transcript_names_'} eq 'yes' ? ($transcript->external_name || 'NOVEL') : $tid;    
  }

sub legend {
    my ($self, $colours) = @_;
    return ('genes', 900, 
        [
            EnsWeb::species_defs->AUTHORITY.' predicted genes (known)' => $colours->{'_KNOWN'},
            EnsWeb::species_defs->AUTHORITY.' predicted genes (novel)' => $colours->{'_'}
        ]
    );
}

sub error_track_name { return EnsWeb::species_defs->AUTHORITY.' transcripts'; }

1;
