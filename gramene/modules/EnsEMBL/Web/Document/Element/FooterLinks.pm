package EnsEMBL::Web::Document::Element::FooterLinks;

### Replacement footer links for www.ensembl.org

use strict;

sub content {

  return qq(
    <script type="text/javascript">

    var _gaq = _gaq || [];
    _gaq.push(['_setAccount', 'UA-1624628-5']);
    _gaq.push(['_trackPageview']);

    (function() {
        var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
        ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
        var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
    })();

    </script>
    <div class="twocol-right right unpadded">
      <a href="http://www.ensemblgenomes.org">About&nbsp;EnsemblGenomes</a> | 
      <a href="/info/about/contact/index.html">Contact&nbsp;Us</a> | 
      <a href="/info/website/help/index.html">Help</a>
    </div>) 
  ;
}

1;

