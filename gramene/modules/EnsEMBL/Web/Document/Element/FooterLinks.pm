package EnsEMBL::Web::Document::Element::FooterLinks;

### Replacement footer links for www.ensembl.org

use strict;

sub content {

  return qq(
    <script>
      (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

      ga('create', 'UA-1624628-14', 'sorghumbase.org');
      ga('send', 'pageview');

    </script>
    <div class="twocol-right right unpadded">
      <a href="https://www.sorghumbase.org/">About&nbsp;Sorghumbase</a> | 
      <a href="https://www.sorghumbase.org/contact">Contact&nbsp;Us</a> | 
      <a href="/info/website/index.html">Help</a>
    </div>) 
  ;
}

1;

