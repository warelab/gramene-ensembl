/* 2006-09-08
 * 
 * This file is intended as a temporary patch to restore the
 * performance of the Ensembl database. It does so by undoing
 * the changes of ensembl/sql/patch_37_38.sql that modified the
 * type of the "value" column from "varchar()" to "text".
 */
ALTER TABLE misc_attrib CHANGE value value VARCHAR(255) NOT NULL DEFAULT '';
ALTER TABLE seq_region_attrib  value value VARCHAR(255) NOT NULL DEFAULT '';
ALTER TABLE gene_attrib        value value VARCHAR(255) NOT NULL DEFAULT '';
ALTER TABLE transcript_attrib  value value VARCHAR(255) NOT NULL DEFAULT '';
ALTER TABLE translation_attrib value value VARCHAR(255) NOT NULL DEFAULT '';