-- MySQL dump 10.13  Distrib 5.1.61, for redhat-linux-gnu (x86_64)
-- 
-- Host: mysql-eg-live-web    Database: ensembl_blast
-- ------------------------------------------------------
-- Server version    5.1.49-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

-- 
-- Table structure for table `alignment`
-- 

DROP TABLE IF EXISTS `alignment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `alignment` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `job_id` varchar(48) NOT NULL,
  `source` varchar(16) DEFAULT NULL,
  `species` varchar(48) DEFAULT NULL,
  `qset` varchar(48) DEFAULT NULL,
  `location` varchar(255) DEFAULT NULL,
  `qstart` bigint(11) DEFAULT NULL,
  `qend` bigint(11) DEFAULT NULL,
  `identity` smallint(6) DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  `result` longtext,
  `region` varchar(32) DEFAULT NULL,
  `tstart` bigint(20) DEFAULT NULL,
  `tend` bigint(20) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=MyISAM AUTO_INCREMENT=91485 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `alignment_copy`
-- 

DROP TABLE IF EXISTS `alignment_copy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `alignment_copy` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `job_id` varchar(48) NOT NULL,
  `source` varchar(16) DEFAULT NULL,
  `species` varchar(48) DEFAULT NULL,
  `qset` varchar(48) DEFAULT NULL,
  `location` varchar(255) DEFAULT NULL,
  `qstart` bigint(11) DEFAULT NULL,
  `qend` bigint(11) DEFAULT NULL,
  `identity` smallint(6) DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  `result` longtext,
  `region` varchar(32) DEFAULT NULL,
  `tstart` bigint(20) DEFAULT NULL,
  `tend` bigint(20) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=MyISAM AUTO_INCREMENT=28240 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_hit`
-- 

DROP TABLE IF EXISTS `blast_hit`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_hit` (
  `hit_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  PRIMARY KEY (`hit_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=604160 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_hit20111215`
-- 

DROP TABLE IF EXISTS `blast_hit20111215`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_hit20111215` (
  `hit_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  PRIMARY KEY (`hit_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=15 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_hit20111216`
-- 

DROP TABLE IF EXISTS `blast_hit20111216`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_hit20111216` (
  `hit_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  PRIMARY KEY (`hit_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=64637 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_hsp`
-- 

DROP TABLE IF EXISTS `blast_hsp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_hsp` (
  `hsp_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  `chr_name` varchar(32) DEFAULT NULL,
  `chr_start` int(10) unsigned DEFAULT NULL,
  `chr_end` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`hsp_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=7957405 DEFAULT CHARSET=latin1 MAX_ROWS=705032704 AVG_ROW_LENGTH=4000;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_hsp20111215`
-- 

DROP TABLE IF EXISTS `blast_hsp20111215`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_hsp20111215` (
  `hsp_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  `chr_name` varchar(32) DEFAULT NULL,
  `chr_start` int(10) unsigned DEFAULT NULL,
  `chr_end` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`hsp_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=1610 DEFAULT CHARSET=latin1 MAX_ROWS=705032704 AVG_ROW_LENGTH=4000;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_hsp20111216`
-- 

DROP TABLE IF EXISTS `blast_hsp20111216`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_hsp20111216` (
  `hsp_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  `chr_name` varchar(32) DEFAULT NULL,
  `chr_start` int(10) unsigned DEFAULT NULL,
  `chr_end` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`hsp_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=521051 DEFAULT CHARSET=latin1 MAX_ROWS=705032704 AVG_ROW_LENGTH=4000;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_result`
-- 

DROP TABLE IF EXISTS `blast_result`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_result` (
  `result_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  PRIMARY KEY (`result_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=86937 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_result20111215`
-- 

DROP TABLE IF EXISTS `blast_result20111215`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_result20111215` (
  `result_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  PRIMARY KEY (`result_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=13 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_result20111216`
-- 

DROP TABLE IF EXISTS `blast_result20111216`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_result20111216` (
  `result_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ticket` varchar(32) DEFAULT NULL,
  `object` longblob,
  PRIMARY KEY (`result_id`),
  KEY `ticket` (`ticket`)
) ENGINE=MyISAM AUTO_INCREMENT=13467 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_table_log`
-- 

DROP TABLE IF EXISTS `blast_table_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_table_log` (
  `table_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `table_name` varchar(32) DEFAULT NULL,
  `table_type` enum('TICKET','RESULT','HIT','HSP') DEFAULT NULL,
  `table_status` enum('CURRENT','FILLED','DELETED') DEFAULT NULL,
  `use_date` date DEFAULT NULL,
  `create_time` datetime DEFAULT NULL,
  `delete_time` datetime DEFAULT NULL,
  `num_objects` int(10) DEFAULT NULL,
  PRIMARY KEY (`table_id`),
  KEY `table_name` (`table_name`),
  KEY `table_type` (`table_type`),
  KEY `use_date` (`use_date`),
  KEY `table_status` (`table_status`)
) ENGINE=MyISAM AUTO_INCREMENT=7 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `blast_ticket`
-- 

DROP TABLE IF EXISTS `blast_ticket`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_ticket` (
  `ticket_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `create_time` datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  `update_time` datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  `ticket` varchar(32) NOT NULL DEFAULT '',
  `object` longblob,
  PRIMARY KEY (`ticket_id`),
  UNIQUE KEY `ticket` (`ticket`),
  KEY `create_time` (`create_time`),
  KEY `update_time` (`update_time`)
) ENGINE=MyISAM AUTO_INCREMENT=38437 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `ena_alignment`
-- 

DROP TABLE IF EXISTS `ena_alignment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ena_alignment` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `job_id` varchar(48) NOT NULL,
  `source` varchar(16) DEFAULT NULL,
  `species` varchar(48) DEFAULT NULL,
  `qset` varchar(48) DEFAULT NULL,
  `location` varchar(255) DEFAULT NULL,
  `qstart` bigint(11) DEFAULT NULL,
  `qend` bigint(11) DEFAULT NULL,
  `identity` smallint(6) DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  `result` longtext,
  `region` varchar(32) DEFAULT NULL,
  `tstart` bigint(20) DEFAULT NULL,
  `tend` bigint(20) DEFAULT NULL,
  `fcount` int(11) DEFAULT NULL,
  `ftext` text,
  PRIMARY KEY (`id`)
) ENGINE=MyISAM AUTO_INCREMENT=28295 DEFAULT CHARSET=latin1 ROW_FORMAT=DYNAMIC;
/*!40101 SET character_set_client = @saved_cs_client */;

-- 
-- Table structure for table `journal`
-- 

DROP TABLE IF EXISTS `journal`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `journal` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `create_time` datetime DEFAULT NULL,
  `job_id` varchar(48) NOT NULL DEFAULT '',
  `progress` smallint(6) DEFAULT NULL,
  `status` varchar(16) DEFAULT NULL,
  `counter` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`,`job_id`)
) ENGINE=MyISAM AUTO_INCREMENT=617 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-06-11 16:03:04 
