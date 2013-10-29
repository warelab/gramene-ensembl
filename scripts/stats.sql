#get repeat feature counts

select rc.repeat_type, count(*) from repeat_consensus rc join repeat_feature rf using (repeat_consensus_id) group by rc.repeat_type;

