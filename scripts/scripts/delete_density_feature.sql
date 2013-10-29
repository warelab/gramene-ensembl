delete a.*, dt.*, df.* from analysis a join density_type dt using (analysis_id) join density_feature df using (density_type_id);
