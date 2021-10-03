#!/bin/bash

just_name=$1
tabl_out="TABLE_Isolate_Typing_results.txt"
bin_out="BIN_Isolate_Typing_results.txt"
#contamination_level=10
printf "$just_name\t" >> "$tabl_out"
printf "$just_name," >> "$bin_out"
###EMM TYPE OUTPUT###
emm_out="NF"
while read -r line
do
    if [[ -n "$line" ]]
    then
        justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
        if [[ "$emm_out" == "NF" ]]
        then
            emm_out="$justTarget"
        else
            emm_out="$emm_out;$justTarget"
        fi
    fi
done <<< "$(sed 1d *__emm-Type__Results.txt)"
printf "$emm_out\t" >> "$tabl_out"
printf "$emm_out," >> "$bin_out"
###MLST OUTPUT###
sed 1d "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt | while read -r line
do
    MLST_tabl=$(echo "$line" | cut -f2-9)
    echo "MLST line: $MLST_tabl\n";
    printf "$MLST_tabl\t" >> "$tabl_out"
    MLST_val=$(echo "$line" | awk -F" " '{print $2}')
    printf "$MLST_val," >> "$bin_out"
done #< "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt
#tail -n+2 "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt | cut -f2-9 >> "$tabl_out"

###Features Targets###
while read -r line
do
    FEAT_targ=$(echo "$line" | cut -f2)
    printf "$FEAT_targ\t" >> "$tabl_out"
done < TEMP_protein_Results.txt

###PBP_ID Output###
justPBPs="NF"
sed 1d TEMP_pbpID_Results.txt | while read -r line
do
    if [[ -n "$line" ]]
    then
	justPBPs=$(echo "$line" | awk -F"\t" '{print $2}')
    fi
    printf "$justPBPs\t" >> "$tabl_out"
done

###Resistance Targets###
while read -r line
do
    #RES_targ=$(echo "$line" | cut -f2)
    #printf "$RES_targ\t" >> "$tabl_out"
    printf "$line\t" | tr ',' '\t' >> "$tabl_out"
done < RES-MIC_"$just_name"
printf "\n" >> "$tabl_out"

cat BIN_Features_Results.txt | sed 's/$/,/g' >> "$bin_out"
cat BIN_Res_Results.txt >> "$bin_out"
printf "\n" >> "$bin_out"
