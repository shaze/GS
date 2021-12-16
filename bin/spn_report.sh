#!/bin/bash

base=$1
tabl_out=$2
sampl_out=$3




printf "$base\t" >> "$tabl_out"
printf "\tSerotype:\n" >> "$sampl_out"
lineNum=$(cat TEMP_SeroType_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no serotype was found
    firstLine=$(head -n1 TEMP_SeroType_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Serotype\n" >> "$sampl_out"
    printf "No_Serotype" >> "$tabl_out"
else
    count=0
    misc_target=()
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -eq 1 ]]
        then
	    #print the header line to the 'SAMPL' output file
            printf "\t\t$line\n" >> "$sampl_out"
        else
            #If misc. resistance target is greater than the contamination threshold then add that
            #misc. resistance target to the output array 'misc_target'
            printf "\t\t$line\n" >> "$sampl_out"
            justTarget=$(echo "$line" | awk -F"\t" '{print $3}')
	    Depth=$(echo "$line" | awk -F"\t" '{print $4}')
	    if [[ $Depth = *"/"* ]]
	    then
		depth1=$(echo $Depth | cut -f1 -d/)
		depth2=$(echo $Depth | cut -f2 -d/)
		justDepth=$([ $depth1 '<' $depth2 ] && echo "$depth1" || echo "$depth2")
	        #echo "Sero duel target: $depth1 | $depth2 | $justDepth" > TEST_sero_table_out.txt
		#min=$([ $var1 '<' $var2 ] && echo "$var1" || echo "$var2")
	    else
		justDepth=$Depth
	    fi
	    justMatchType=$(echo "$line" | awk -F"\t" '{print $2}')
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
		echo "Target $justTarget is a match"
		#misc_target+=("$justTarget($justDepth|$justMatchType)")
		#printf "$justTarget;" >> TEMP_table_results.txt
		misc_target+=("$justTarget")
	    fi
        fi
    done < TEMP_SeroType_Results.txt
    #if the output array 'misc_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if [ ${#misc_target[@]} -eq 0 ];
    then
        printf "No_Serotype" >> "$tabl_out"
    else
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g'
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g' >> "$tabl_out"
    fi
fi
printf "\t" >> "$tabl_out"
EOF

: <<EOF
###MLST OUTPUT###
printf "\tMLST:\n" >> "$sampl_out"
count=0
while read -r line
do
    count=$(( $count + 1 ))
    if [[ "$count" -eq 1 ]]
    then
        printf "\t\t$line\n" >> "$sampl_out"
    else
        printf "\t\t$line\n" >> "$sampl_out"
        MLST_tabl=$(echo "$line" | cut -f2-9)
        printf "$MLST_tabl\t" >> "$tabl_out"
    fi
done < "$out_nameMLST"__mlst__Streptococcus_pneumoniae__results.txt
EOF

###PBP_ID OUTPUT###
printf "\tPBP_ID Code:\n" >> "$sampl_out"
lineNum=$(cat TEMP_pbpID_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no PBP results were found
    firstLine=$(head -n1 TEMP_pbpID_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_PBP_Type\n" >> "$sampl_out"
    #printf "$firstLine\t" >> TEMP_table_title.txt
    printf "No_PBP_Type\t" >> "$tabl_out"
else
    count=0
    while read -r line
    do
	count=$(( $count + 1 ))
        #justPBPs=$(echo "$line" | cut -f2-4)
	justPBPs=$(echo "$line" | awk -F"\t" '{print $2}')
	if [[ "$count" -eq 1 ]]
	then
            printf "\t\t$justPBPs\n" >> "$sampl_out"
            #printf "$justPBPs\t" >> TEMP_table_title.txt
	else
            printf "\t\t$justPBPs\n" >> "$sampl_out"
            printf "$justPBPs\t" >> "$tabl_out"
	fi
    done < TEMP_pbpID_Results.txt
fi
#EOF

###MISC. RESISTANCE###
printf "\tMisc. GBS Resistance:\n" >> "$sampl_out"
lineNum=$(cat TEMP_miscR_Results.txt | wc -l)
misc_contamination_level=7
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no misc. resistance were found
    firstLine=$(head -n1 TEMP_miscR_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Resistance\n" >> "$sampl_out"
    #printf "$firstLine\t" >> TEMP_table_title.txt
    if grep -q failed MISC_*.log;
    then
        printf "**Failed Read Mapping - Results not accurate**" >> "$tabl_out"
    else
	printf "No_Resistance" >> "$tabl_out"
    fi
else
    count=0
    misc_target=()
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -eq 1 ]]
        then
	    #print the header line to the 'SAMPL' output file
            printf "\t\t$line\n" >> "$sampl_out"
            #printf "$line\t" >> TEMP_table_title.txt
        else
            #If misc. resistance target is greater than the contamination threshold then add that
            #misc. resistance target to the output array 'misc_target'
            printf "\t\t$line\n" >> "$sampl_out"
            justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
            justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
            justMatchType=$(echo "$line" | awk -F"\t" '{print $2}')
            if [[ $(echo "$justDepth > $misc_contamination_level" | bc) -eq 1 ]] || [[ "$justDepth" == "NA" ]]
            then
		#this condition checks if the target is one of the extraction alleles. If so, it will append (extract) to output
		#if [[ "$justMatchType" == "imperfect" && " ${extract_arr[@]} " =~ " ${justTarget} " ]]
		if [[ "$line" =~ "Extract" ]]
		then
                    echo "Target $justTarget is a match but needs extraction"
		    #misc_target+=("$justTarget($justDepth|$justMatchType)")
		    #printf "$justTarget;" >> TEMP_table_results.txt
		    misc_target+=("$justTarget(extract)")
		else
		    echo "Target $justTarget is a match"
		    misc_target+=("$justTarget($justMatchType)")
		    #misc_target+=("$justTarget")
		fi
	    fi
        fi
    done < TEMP_miscR_Results.txt
    #if the output array 'misc_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if grep -q failed MISC_*.log;
    then
        printf "**Failed Read Mapping - Results not accurate**" >> "$tabl_out"
    else
	if [ ${#misc_target[@]} -eq 0 ];
	then
            printf "No_Resistance" >> "$tabl_out"
	else
            printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g'
            printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g' >> "$tabl_out"
	fi
    fi
fi
printf "\t" >> "$tabl_out"

###Surface / Secretory Protein Output (Not including T-Antigens)###
#printf "\tProtein Targets:\n" >> "$sampl_out"
#lineNum=$(cat TEMP_protein_Results.txt | wc -l)
#if [[ "$lineNum" -eq 1 ]]
#then
#    #if the file only contains the header line then no surface/secretory targets were typed
#    firstLine=$(head -n1 TEMP_protein_Results.txt)
#    printf "\t\t$firstLine\n\t\tNo_Protein_Targets\n" >> "$sampl_out"
#    printf "No_Protein_Targets" >> "$tabl_out"
#else
#    count=0
#    prot_target=()
#    while read -r line
#    do
#        count=$(( $count + 1 ))
#        if [[ "$count" -eq 1 ]]
#        then
#            #print the header line to the 'SAMPL' output file
#            printf "\t\t$line\n" >> "$sampl_out"
#        else
#            #If surface/secretory target is greater than the contamination threshold then add that
#            #target to the output array 'prot_target'
#            printf "\t\t$line\n" >> "$sampl_out"
#            justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
#            justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
#            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
#            then
#                echo "Target $justTarget is a match"
#                #printf "$justTarget-($justDepth);" >> TEMP_table_results.txt
#                #prot_target+=("$justTarget($justDepth)")
#		prot_target+=("$justTarget")
#            fi
#        fi
#    done < TEMP_protein_Results.txt
#    #done < TEMP_no-Tantigen_Results.txt
#    #if the output array 'prot_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
#    if [ ${#prot_target[@]} -eq 0 ];
#    then
#        printf "No_Protein_Targets" >> "$tabl_out"
#    else
#        printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g'
#        printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g' >> "$tabl_out"
#    fi
#fi
#printf "\t" >> "$tabl_out"

: <<EOF
###ARG-ANNOT and ResFinder Resistance Gene Typing Output###
printf "\tGeneral Resistance Targets:\n\t\tDB_Target\tMatch_Type\tDepth\n" >> "$sampl_out"
genRes_target=()
if [[ -s "$out_nameARG"__fullgenes__ARGannot_r1__results.txt ]]
then
    count=0
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -ne 1 ]]
        then
	    isIdentical="identical"
	    justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
	    if [[ -n "$justDiffs" ]]
	    then
		isIdentical="imperfect"
	    fi
	    justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
	    justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
	    #printf "\t\tARGannot_r1_$justTarget\t$isIdentical\t$justDepth\n" >> "$sampl_out"
	    printf "\t\tARG_$justTarget\t$isIdentical\t$justDepth\n" >> "$sampl_out"
	    if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
	    then
		echo "Target $justTarget is a match"
		genRes_target+=("ARG_$justTarget")
	    fi
	fi
    done < "$out_nameARG"__fullgenes__ARGannot_r1__results.txt
fi
if [[ -s "$out_nameRES"__fullgenes__ResFinder__results.txt ]]
then
    count=0
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -ne 1 ]]
        then
	    isIdentical="identical"
	    justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
	    if [[ -n "$justDiffs" ]]
	    then
		isIdentical="imperfect"
	    fi
	    justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
	    justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
	    #printf "\t\tResFinder_$justTarget\t$isIdentical\t$justDepth\n" >> "$sampl_out"
	    printf "\t\tRF_$justTarget\t$isIdentical\t$justDepth\n" >> "$sampl_out"
	    if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
	    then
		echo "Target $justTarget is a match"
		genRes_target+=("RF_$justTarget")
	    fi
	fi
    done < "$out_nameRES"__fullgenes__ResFinder__results.txt
fi
#if the output array 'genRes_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
if [ ${#genRes_target[@]} -eq 0 ];
then
    printf "No_Gen_Resistance_Targets" >> "$tabl_out"
    printf "\t\tNo_Gen_Resistance_Targets\n" >> "$sampl_out"
else
    printf '%s\n' "${genRes_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g'
    printf '%s\n' "${genRes_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g' >> "$tabl_out"
fi
printf "\t" >> "$tabl_out"

###PlasmidFinder Plasmid Typing Output###
printf "\tPlasmid Prediction Targets:\n\t\tTarget\tMatch_Type\tDepth\n" >> "$sampl_out"
if [[ -s "$out_namePLAS"__fullgenes__PlasmidFinder__results.txt ]]
then 
    count=0
    while read -r line
    do
	count=$(( $count + 1 ))
	if [[ "$count" -ne 1 ]]
	then
	    isIdentical="identical"
	    justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
	    if [[ -n "$justDiffs" ]]
	    then
		isIdentical="imperfect"
	    fi
	    justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
	    justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
	    printf "\t\t$justTarget\t$isIdentical\t$justDepth\n" >> "$sampl_out"
	    if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
	    then
		echo "Target $justTarget is a match"
		plas_target+=("$justTarget")
	    fi
	fi
    done < "$out_namePLAS"__fullgenes__PlasmidFinder__results.txt
fi
#if the output array 'genRes_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
if [ ${#plas_target[@]} -eq 0 ];
then
    printf "No_Plasmid_Targets" >> "$tabl_out"
    printf "\t\tNo_Plasmid_Targets" >> "$sampl_out"
else
    printf '%s\n' "${plas_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g'
    printf '%s\n' "${plas_target[@]}" | sort | tr '\n' ';' | sed 's/;$//g' >> "$tabl_out"
fi
EOF
printf "\n\n" >> "$sampl_out"
printf "\n" >> "$tabl_out"

rm *.fastq.gz
rm cutadapt_*.fastq
rm *.sam
rm *.bam
rm *.pileup


