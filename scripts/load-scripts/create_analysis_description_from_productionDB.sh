#!/bin/sh


LG_FILE=$1
echo "logic_name fileis $LG_FILE"

LGS=`sed -n "s/.*/'&',/p" $LG_FILE | tr '\n' ' '`

echo "logic_names are $LGS"

echo "select 1, ad.logic_name, ad.description, ad.display_label, 1, wd.data from analysis_description ad join analysis_web_data awd using (analysis_description_id) join web_data wd using (web_data_id) where logic_name in (${LGS}'')";

mysql -Nq -u weix -pwarelab -h cabot ensembl_production -e "select 1, ad.logic_name, ad.description, ad.display_label, 1, wd.data from analysis_description ad join analysis_web_data awd using (analysis_description_id) join web_data wd using (web_data_id) where logic_name in (${LGS}'')";

echo "if nothing returned from above query, try the following"


echo "select 1, ad.logic_name, ad.description, ad.display_label, 1, null from analysis_description ad join analysis_web_data awd using (analysis_description_id) where logic_name in (${LGS}'')"
mysql -Nq -u weix -pwarelab -h cabot ensembl_production -e "select 1, ad.logic_name, ad.description, ad.display_label, 1, null from analysis_description ad join analysis_web_data awd using (analysis_description_id) where logic_name in (${LGS}'')";
