#
# Table structure for table 'quality_score'
#

CREATE TABLE quality_score (
  seq_region_id          int(10) unsigned not null,
  window_size            smallint unsigned not null,
  position	         int unsigned not null,
  score                  blob,

  PRIMARY KEY (seq_region_id, window_size, position)
) MAX_ROWS = 15000000 AVG_ROW_LENGTH = 841 COLLATE=latin1_swedish_ci;
