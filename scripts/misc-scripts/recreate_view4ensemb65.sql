DROP VIEW if EXISTS exon_stable_id;

DROP VIEW if EXISTS gene_stable_id;

DROP VIEW if EXISTS operon_stable_id;

DROP VIEW if EXISTS operon_transcript_stable_id;

DROP VIEW if EXISTS translation_stable_id;

DROP VIEW if EXISTS transcript_stable_id;

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW exon_stable_id (exon_id, stable_id, version, created_date, modified_date) AS (SELECT exon_id, stable_id, version, created_date, modified_date FROM exon);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW gene_stable_id (gene_id, stable_id, version, created_date, modified_date) AS (SELECT gene_id, stable_id, version, created_date, modified_date FROM gene);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW operon_stable_id (operon_id, stable_id, version, created_date, modified_date) AS (SELECT operon_id, stable_id, version, created_date, modified_date FROM operon);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW operon_transcript_stable_id (operon_transcript_id, stable_id, version, created_date, modified_date) AS (SELECT operon_transcript_id, stable_id, version, created_date, modified_date FROM operon_transcript);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW translation_stable_id (translation_id, stable_id, version, created_date, modified_date) AS (SELECT translation_id, stable_id, version, created_date, modified_date FROM translation);

CREATE DEFINER = CURRENT_USER SQL SECURITY INVOKER VIEW transcript_stable_id (transcript_id, stable_id, version, created_date, modified_date) AS (SELECT transcript_id, stable_id, version, created_date, modified_date FROM transcript);


