create temporary table good_repeat_features as
select min(repeat_feature_id) as repeat_feature_id
  from repeat_feature
 group by seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id;
 
create index good_repeat_features_index
       using btree on good_repeat_features (repeat_feature_id);

delete repeat_feature
  from repeat_feature
       left join good_repeat_features using (repeat_feature_id)
 where good_repeat_features.repeat_feature_id is null;
