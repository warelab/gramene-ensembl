#!/bin/sh

FROM=$1
TO=$2

echo "copy analysis description from $FROM to $TO"

mysql -u weix -pwarelab -h cabot -e \
"insert into ${TO}.analysis_description select b3a.analysis_id, b1ad.description, b1ad.display_label, b1ad.displayable, b1ad.web_data  from ${TO}.analysis b3a left join ${TO}.analysis_description b3ad on (b3a.analysis_id=b3ad.analysis_id)      join ${FROM}.analysis b1a on (b3a.logic_name=b1a.logic_name) join ${FROM}.analysis_description b1ad on (b1a.analysis_id=b1ad.analysis_id) where b3ad.analysis_id is null"



