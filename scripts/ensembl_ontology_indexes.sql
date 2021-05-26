# Andrew added two index child_parent_idx, is_root_idx to speed up the query for ensembl browser gene page
# He added them on cabot:ensembl_ontology_88

CREATE TABLE `closure` (
  `closure_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `child_term_id` int(10) unsigned NOT NULL,
  `parent_term_id` int(10) unsigned NOT NULL,
  `subparent_term_id` int(10) unsigned DEFAULT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  `ontology_id` int(10) unsigned NOT NULL,
  `confident_relationship` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`closure_id`),
  KEY `parent_subparent_idx` (`parent_term_id`,`subparent_term_id`),
  KEY `child_parent_idx` (`child_term_id`,`parent_term_id`,`ontology_id`)
) ENGINE=MyISAM AUTO_INCREMENT=8432534 DEFAULT CHARSET=latin1;

CREATE TABLE `term` (
  `term_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ontology_id` int(10) unsigned NOT NULL,
  `subsets` text,
  `accession` varchar(64) NOT NULL,
  `name` varchar(255) NOT NULL,
  `definition` text,
  `is_root` int(11) NOT NULL DEFAULT '0',
  `is_obsolete` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`term_id`),
  UNIQUE KEY `accession_idx` (`accession`),
  UNIQUE KEY `ontology_acc_idx` (`ontology_id`,`accession`),
  KEY `name_idx` (`name`),
  KEY `is_root_idx` (`is_root`)
) ENGINE=MyISAM AUTO_INCREMENT=188801 DEFAULT CHARSET=latin1;
