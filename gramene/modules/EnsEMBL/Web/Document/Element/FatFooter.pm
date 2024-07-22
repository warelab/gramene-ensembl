=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package EnsEMBL::Web::Document::Element::FatFooter;

### Optional fat footer - site-specific, so see plugins 

use strict;

use base qw(EnsEMBL::Web::Document::Element);

sub content {
  my $species_defs = shift->species_defs;
  #my $sister_sites = qq(<p><a href="http://oge.gramene.org">OGE Browser</a></p>
  #			<p><a href="http://maizev4.gramene.org">Maize B73_RefGen_v4 PreRelease Browser</a></p>);
  
  my $sister_sites = qq(<p><a href="https://vitis.gramene.org">Grapevine Gramene Browser</a></p>
                        <p><a href="https://maize-pangenome.gramene.org/">Maize Gramene Browser</a></p>
			<p><a href="https://www.sorghumbase.org/">Sorghumbase Browser</a></p>);
  my $html = '<hr /><div id="fat-footer">';

  $html .= qq(
              <div class="column-three left">
                <h3>About Us</h3>
                <p><a href="http://www.gramene.org/about-gramene">About Gramene</a></p>
	        <p><a href="http://oryza.gramene.org/feedback">Contact us</a></p>
                <p><a href="http://www.gramene.org/cite">Citing Gramene</a></p>
                <!--<p><a href="http://www.ebi.ac.uk/about/privacy">Privacy policy</a></p>-->
                <!--<p><a href="http://www.ensemblgenomes.org/info/about/cookies">Cookies</a></p>-->
                <!--<p><a href="http://www.ebi.ac.uk/Information/termsofuse.html">EMBL-EBI Terms of use</a></p>-->
                <!--<p><a href="http://ensemblgenomes.org/info/about/legal">Disclaimer</a></p>-->
              </div>
  );


 #$html .= qq(
 #             <div class="column-four left">
 #               <h3>Get help</h3>
 #               <p><a href="/info/website/">Using this website</a></p>
 #               <p><a href="http://ensemblgenomes.org/info">Documentation</a></p>
 #               <p><a href="/info/website/upload">Adding custom tracks</a></p>
 #               <p><a href="/info/website/ftp/index.html">Downloading data</a></p>
 #             </div>
 # );

  #foreach("bacteria","fungi","plants","protists","metazoa"){
    #$sister_sites .= qq(<p><a href="http://$_.ensembl.org">Ensembl ${\ucfirst($_)}</a></p>) if $species_defs->EG_DIVISION ne $_;
  #}

  $html .= qq(
              <div class="column-three left">
                <h3>Our sister sites</h3>
                $sister_sites
              </div>
  );


  $html .= qq(
              <div class="column-three left">
                <h3>Follow us</h3>
                <p><a class="media-icon" href="http://www.gramene.org/blog">
                  <img alt="[RSS logo]" title="Gramene blog" src="/i/rss_icon_16.png"></a>
                  <a href="http://www.gramene.org/blog">Blog</a></p>
                <p><a class="media-icon" href="https://twitter.com/ensemblgenomes">
                  <img alt="[twitter logo]" title="Follow us on Twitter!" src="/i/twitter.png"></a>
                    <a href="https://twitter.com/intent/follow?original_referer=http%3A%2F%2Fwww.gramene.org%2F&ref_src=twsrc%5Etfw&screen_name=GrameneDatabase&tw_p=followbutton">Twitter</a></p>
              </div>
  );

  $html .= '</div>';

  return $html;
}

1;
