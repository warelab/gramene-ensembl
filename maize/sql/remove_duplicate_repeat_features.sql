create temporary table extra_repeat_features as
select max(repeat_feature_id) as repeat_feature_id
  from repeat_feature
 group by seq_region_id, seq_region_start, seq_region_end, seq_region_strand, repeat_consensus_id
 having count(*) >= 2;

delete
  from repeat_feature
 using repeat_feature, extra_repeat_features
 where repeat_feature.repeat_feature_id = extra_repeat_features.repeat_feature_id;
