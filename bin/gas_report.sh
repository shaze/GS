#!/bin/bash

base=$1
tabl_out=$2
bin_out=$3
emm=$4


printf "$base\t" >> "$tabl_out"
printf "$base," >> "$bin_out"
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
done <<< "$(sed 1d $emm)"
printf "$emm_out\t" >> "$tabl_out"
printf "$emm_out," >> "$bin_out"

nlen=`wc -l MLST_${base}__mlst__Streptococcus_pyogenes__results.txt | awk '{print $1}'`
if [[ "$nlen" -eq "1" ]]; then
    MLST_tabl=".       .      .      .      .      .      .      ."
    printf "$MLST_tabl\t" >>  "$tabl_out"
    printf ".,">> $bin_out
else
  sed 1d "MLST_${base}__mlst__Streptococcus_pyogenes__results.txt" | while read -r line
  do
    MLST_tabl=$(echo "$line" | cut -f2-9)
    echo "MLST line: $MLST_tabl\n";
    printf "$MLST_tabl\t" >> "$tabl_out"
    MLST_val=$(echo "$line" | awk -F" " '{print $2}')
    printf "$MLST_val," >> "$bin_out"
  done 
fi


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
done < "RES-MIC_$base"

if [[ -e $(echo ./velvet_output/*_Logfile.txt) ]]
then
    vel_metrics=$(echo ./velvet_output/*_Logfile.txt)
    print "velvet metrics file: $vel_metrics\n";
    velvetMetrics.pl -i "$vel_metrics";
    line=$(cat velvet_qual_metrics.txt | tr ',' '\t')
    printf "$line\t" >> "$tabl_out"
    #pwd | xargs -I{} echo {}"/velvet_output/contigs.fa" >> "$tabl_out"
else
    printf "NA\tNA\tNA\tNA\t$readPair_1\tNA\n" >> "$tabl_out"
fi

#printf "*\t*\t*\t*\t*\t*" >> $tabl_out
printf "\n" >> "$tabl_out"

cat BIN_Features_Results.txt | sed 's/$/,/g' >> "$bin_out"
cat BIN_Res_Results.txt >> "$bin_out"
printf "\n" >> "$bin_out"
