// Gramene: the cookie banner links a Privacy Policy and Terms of Use (GDPR_POLICY_URL, GDPR_TERMS_URL).
// Gramene only has a privacy policy (CSHL's), so with no terms URL show just the privacy link.
Ensembl.LayoutManager = Ensembl.LayoutManager.extend({

  showGDPRCookieBanner: function () {
    var shown = this.base();

    if (shown && !Ensembl.gdpr_terms_url) {
      $('.cookie-message .msg').html(
        "This website requires cookies, and the limited processing of your personal data in order to function. By using the site you are agreeing to this as outlined in our " +
        "<a target='_blank' href='" + Ensembl.gdpr_policy_url + "'>Privacy Policy</a>."
      );
    }

    return shown;
  }
});
