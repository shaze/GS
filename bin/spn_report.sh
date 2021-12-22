#!/bin/bash

base=$1
sampl_out=$2
tabl_out=$3

bin_out="BIN_Isolate_Typing_results.txt"

bin_out="BIN_Isolate_Typing_results.txt"
printf "$base\t" >> "$tabl_out"
printf "$base," >> "$bin_out"

###Serotype Output###
sero_out="NF"
pili_out="neg"
while read -r line
do
    if [[ -n "$line" ]]
    then
        justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
	if [[ "$justTarget" == "PI-1" ]]
	then
            if [[ "$pili_out" == "neg" ]]
            then
		pili_out="1"
            elif [[ "$pili_out" == "2" ]]
	    then
		pili_out="1:2"
            fi
	elif [[ "$justTarget" == "PI-2" ]]
	then
            if [[ "$pili_out" == "neg" ]]
            then
                pili_out="2"
            elif [[ "$pili_out" == "1" ]]
	    then
                pili_out="1:2"
            fi
        else
            if [[ "$sero_out" == "NF" ]]
            then
		sero_out="$justTarget"
            else
		sero_out="$sero_out;$justTarget"
            fi
	fi
    fi
done <<< "$(sed 1d OUT_SeroType_Results.txt)"
printf "$sero_out\t$pili_out\t" >> "$tabl_out"
printf "$sero_out,$pili_out\t" >> "$bin_out"

###MLST OUTPUT###
sed 1d MLST_${base}__mlst__Streptococcus_pneumoniae__results.txt | while read -r line
do
    MLST_tabl=$(echo "$line" | cut -f2-9)
    echo "MLST line: $MLST_tabl\n";
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

###bLactam Predictions###
#sed 1d "BLACTAM_MIC_RF.txt" | while read -r line
#sed 1d "BLACTAM_MIC_RF_with_SIR.txt" | while read -r line
#do
#    pbpID=$(tail -n1 "TEMP_pbpID_Results.txt" | awk -F"\t" '{print $2}')
#    if [[ ! "$pbpID" =~ .*NF.* ]] && [[ ! "$pbpID" =~ .*NEW.* ]]
#    then
#	echo "No NF or NEW outputs for PBP Type"
#	bLacTab=$(echo "$line" | tr ' ' '\t')
#	printf "$bLacTab\t" >> "$tabl_out"
#	bLacCom=$(echo "$line" | tr ' ' ',')
#	printf "$bLacCom," >> "$bin_out"
#    else
#	echo "One of the PBP types has an NF or NEW"
#	printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\t" >> "$tabl_out"
#	printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF," >> "$bin_out"
#    fi
#done

pbpID=$(tail -n1 "TEMP_pbpID_Results.txt" | awk -F"\t" '{print $2}')
if [[ ! "$pbpID" =~ .*NF.* ]] #&& [[ ! "$pbpID" =~ .*NEW.* ]]
then
    echo "No NF outputs for PBP Type"
    bLacTab=$(tail -n1 "BLACTAM_MIC_RF_with_SIR.txt" | tr ' ' '\t')
    printf "$bLacTab\t" >> "$tabl_out"
    #bLacCom=$(echo "$line" | tr ' ' ',')
    #printf "$bLacCom," >> "$bin_out"
else
    echo "One of the PBP types has an NF"
    printf "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> "$tabl_out"
    #printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF," >> "$bin_out"
fi


###Resistance Targets###
while read -r line
do
    #RES_targ=$(echo "$line" | cut -f2)
    #printf "$RES_targ\t" >> "$tabl_out"
    printf "$line\t" | tr ',' '\t' >> "$tabl_out"
done < RES-MIC_"$base"

if [[ -e $(echo ./velvet_output/*_Logfile.txt) ]]
then
    vel_metrics=$(echo ./velvet_output/*_Logfile.txt)
    echo "velvet metrics file: $vel_metrics\n";
    velvetMetrics.pl -i "$vel_metrics";
    line=$(cat velvet_qual_metrics.txt | tr ',' '\t')
    printf "$line\t" >> "$tabl_out"

    printf "$readPair_1\t" >> "$tabl_out";
    pwd | xargs -I{} echo {}"/velvet_output/contigs.fa" >> "$tabl_out"
else
    printf "NA\tNA\tNA\tNA\t$readPair_1\tNA\n" >> "$tabl_out"
fi

