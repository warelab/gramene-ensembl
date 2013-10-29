/*
This SQL assigns the repeat_consensus.repeat_type fields in the
Ensembl core database based on the MIPS repeat_consensus.repeat_class
values. This is possible by the way that the repeat_class values
represent the reCat ontology

See http://www.warelab.org/bugs/view.php?id=1974 for more details.
*/

UPDATE repeat_consensus
SET    repeat_type = 'Transposons/Unclassified'
WHERE  repeat_class LIKE '02%';

UPDATE repeat_consensus
SET    repeat_type = 'Type I Transposons/LTR' 
WHERE  repeat_class LIKE '02.01.01%';

UPDATE repeat_consensus
SET    repeat_type = 'Type I Transposons/LINE' 
WHERE  repeat_class LIKE '02.01.02.05%';

UPDATE repeat_consensus
SET    repeat_type = 'Type I Transposons/SINE' 
WHERE  repeat_class LIKE '02.01.02.10%';

UPDATE repeat_consensus
SET    repeat_type = 'Type II Transposons' 
WHERE  repeat_class LIKE '02.05%';

UPDATE repeat_consensus
SET    repeat_type = 'Type III Transposons' 
WHERE  repeat_class LIKE '02.10%';

UPDATE repeat_consensus
SET    repeat_type = 'RNA repeats' 
WHERE  repeat_class LIKE '10.01%';

UPDATE repeat_consensus
SET    repeat_type = 'Simple Sequence Repeat' 
WHERE  repeat_class LIKE '01%';

UPDATE repeat_consensus
SET    repeat_type = 'Unknown' 
WHERE  repeat_class LIKE '90%' 
    OR repeat_class LIKE '98%'
    OR repeat_class LIKE '99%';

