#!/bin/bash


###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
base=$1
tabl_out="TABLE_Isolate_Typing_results.txt"
bin_out="BIN_Isolate_Typing_results.txt"
printf "$base\t" >> "$tabl_out"
printf "$base," >> "$bin_out"


###Serotype Output###
sero_out="NF"
while read -r line
do
    if [[ -n "$line" ]]
    then
        justTarget=$(echo "$line" | awk -F"\t" '{print $3}')
        if [[ "$sero_out" == "NF" ]]
        then
            sero_out="$justTarget"
        else
            sero_out="$sero_out;$justTarget"
        fi
    fi
done <<< "$(sed 1d TEMP_SeroType_Results.txt)"
printf "$sero_out\t" >> "$tabl_out"
printf "$sero_out," >> "$bin_out"



###MLST OUTPUT###
sed 1d "MLST_${base}__mlst__Streptococcus_agalactiae__results.txt" | while read -r line
do
    MLST_tabl=$(echo "$line" | cut -f2-9)
    printf "$MLST_tabl\t" >> "$tabl_out"
    MLST_val=$(echo "$line" | awk -F" " '{print $2}')
    printf "$MLST_val," >> "$bin_out"
done

###PBP_ID Output###
justPBPs="NF"
sed 1d TEMP_pbpID_Results.txt | while read -r line
do
    if [[ -n "$line" ]]
    then
        justPBPs=$(echo "$line" | awk -F"\t" '{print $2}' | tr ':' '\t')
	justPBP_BIN=$(echo "$line" | awk -F"\t" '{print $2}' | tr ':' ',')
    fi
    printf "$justPBPs\t" >> "$tabl_out"
    printf "$justPBP_BIN," >> "$bin_out"
done

###Resistance Targets###
#while read -r line
#do
#    RES_targ=$(echo "$line" | cut -f2)
#    printf "$RES_targ\t" >> "$tabl_out"
#done < TEMP_Res_Results.txt
#cat BIN_Res_Results.txt | sed 's/$/,/g' >> "$bin_out"

###Resistance Targets###
while read -r line
do
    #RES_targ=$(echo "$line" | cut -f2)
    #printf "$RES_targ\t" >> "$tabl_out"
    printf "$line\t" | tr ',' '\t' >> "$tabl_out"
done < "RES-MIC_$base"

###Surface Targets###
while read -r line
do
    SURF_targ=$(echo "$line" | cut -f2)
    printf "$SURF_targ\t" >> "$tabl_out"
done < TEMP_Surface_Results.txt
cat BIN_Surface_Results.txt >> "$bin_out"

if [[ -e $(echo ./velvet_output/*_Logfile.txt) ]]
then
    vel_metrics=$(echo ./velvet_output/*_Logfile.txt)
    print "velvet metrics file: $vel_metrics\n";
    velvetMetrics.pl -i "$vel_metrics";
    line=$(cat velvet_qual_metrics.txt | tr ',' '\t')
    printf "$line\t" >> "$tabl_out"

    printf "$readPair_1\t" >> "$tabl_out";
    pwd | xargs -I{} echo {}"/velvet_output/contigs.fa" >> "$tabl_out"
else
    printf "NA\tNA\tNA\tNA\t$readPair_1\tNA\n" >> "$tabl_out"
fi
#printf "\n" >> "$tabl_out"

#printf "\n" >> "$tabl_out"
#printf "\n" >> "$bin_out"
