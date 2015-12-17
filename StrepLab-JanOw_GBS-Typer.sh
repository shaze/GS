#!/bin/bash -l

#export PATH=/scicomp/home/ycm6/TEMP_GBS-Typing:$PATH

## -- begin embedded SGE options --
read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1/job-control.txt)
## -- end embedded SGE options --

###Load Modules###
#. /usr/share/Modules/init/bash
module load perl/5.12.3
module load ncbi-blast+/2.2.29
module load BEDTools/2.17.0
module load Python/2.7
module load samtools/0.1.18
module load bowtie2/2.1.0
module load freebayes/0.9.21
module load prodigal/2.60

###This script is called for each job in the qsub array. The purpose of this code is to read in and parse a line of the job-control.txt file
###created by 'StrepLab-JanOw_GBS-wrapr.sh' and pass that information, as arguments, to other programs responsible for various parts of strain
###characterization (MLST, serotype and antibiotic drug resistance prediction).

readPair_1=${PARAM[0]}
readPair_2=${PARAM[1]}
allDB_dir=${PARAM[2]}
batch_out=${PARAM[3]}
sampl_out=${PARAM[4]}
#mlst_ref=${PARAM[5]}
#mlst_def=${PARAM[6]}





###Start Doing Stuff###
cd "$sampl_out"
batch_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)}')
out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')  ###Use This For Batches off the MiSeq###
#out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-1)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')   ###Otherwise Use This###
just_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')
out_nameMLST=MLST_"$just_name"
out_nameSERO=SERO_"$just_name"
out_nameMISC=MISC_"$just_name"
out_namePBP=PBP_"$just_name"
out_namePROT=PROT_"$just_name"

###Call MLST###
mod-srst2.py --mlst_delimiter '_' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMLST" --save_scores --mlst_db "$allDB_dir/Streptococcus_agalactiae.fasta" --mlst_definitions "$allDB_dir/sagalactiae.txt"
###Check and extract new MLST alleles###
/scicomp/home/ycm6/TEMP_GBS-Typing/MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_agalactiae__results.txt "$out_nameMLST"__*.Streptococcus_agalactiae.sorted.bam "$allDB_dir/Streptococcus_agalactiae.fasta"

###Call GBS Serotype###
/scicomp/home/ycm6/TEMP_GBS-Typing/GBS_serotyper.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/GBS_seroT_Gene-DB_Final.fasta" -n "$out_nameSERO"

###Call GBS bLactam Resistances###
/scicomp/home/ycm6/TEMP_GBS-Typing/GBS_PBP-Gene_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/GBS_bLactam_Ref.fasta" -n "$out_namePBP"

###Call GBS Misc. Resistances###
/scicomp/home/ycm6/TEMP_GBS-Typing/GBS_miscRes_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -m GBS_miscR_Gene-DB_Final.fasta -v GBS_vancR_Gene-DB_Final.fasta -n "$out_nameMISC"

###Type Surface adn Secretory Proteins###
/scicomp/home/ycm6/TEMP_GBS-Typing/GBS_surface-secretory_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -p GBS_protein_Gene-DB_Final.fasta -n "$out_namePROT"

#Add in perl script to find contamination threshold here
contamination_level=10




###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
sampl_out="TEMP_GBS_Typing_Results.txt"

printf "$just_name\n" >> "$sampl_out"
printf "$just_name\t" >> TEMP_table_results.txt
###SEROTYPE OUTPUT###
printf "\tSerotype:\n" >> "$sampl_out"
lineNum=$(cat TEMP_SeroType_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no serotype was found
    firstLine=$(head -n1 TEMP_SeroType_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Serotype\n" >> "$sampl_out"
    printf "No_Serotype" >> TEMP_table_results.txt
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
	    justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
	    justMatchType=$(echo "$line" | awk -F"\t" '{print $2}')
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
		echo "Target $justTarget is a match"
		misc_target+=("$justTarget($justDepth|$justMatchType)")
		#printf "$justTarget;" >> TEMP_table_results.txt
	    fi
        fi
    done < TEMP_SeroType_Results.txt
    #if the output array 'misc_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if [ ${#misc_target[@]} -eq 0 ];
    then
        printf "No_Serotype" >> TEMP_table_results.txt
    else
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';'
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' >> TEMP_table_results.txt
    fi
fi
printf "\t" >> TEMP_table_results.txt

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
        printf "$MLST_tabl\t" >> TEMP_table_results.txt
    fi
done < "$out_nameMLST"__mlst__Streptococcus_agalactiae__results.txt

###PBP_ID OUTPUT###
printf "\tPBP_ID Code:\n" >> "$sampl_out"
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
        printf "$justPBPs\t" >> TEMP_table_results.txt
    fi
done < TEMP_pbpID_Results.txt

###MISC. RESISTANCE###
printf "\tMisc. Resistance:\n" >> "$sampl_out"
lineNum=$(cat TEMP_miscR_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no misc. resistance were found
    firstLine=$(head -n1 TEMP_miscR_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Resistance\n" >> "$sampl_out"
    #printf "$firstLine\t" >> TEMP_table_title.txt
    printf "No_Resistance" >> TEMP_table_results.txt
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
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
                echo "Target $justTarget is a match"
		misc_target+=("$justTarget($justDepth|$justMatchType)")
		#printf "$justTarget;" >> TEMP_table_results.txt
	    fi
        fi
    done < TEMP_miscR_Results.txt
    #if the output array 'misc_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if [ ${#misc_target[@]} -eq 0 ];
    then
        printf "No_Resistance" >> TEMP_table_results.txt
    else
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';'
        printf '%s\n' "${misc_target[@]}" | sort | tr '\n' ';' >> TEMP_table_results.txt
    fi
fi
printf "\t" >> TEMP_table_results.txt

###Surface / Secretory Protein Output (Not including T-Antigens)###
printf "\tProtein Targets:\n" >> "$sampl_out"
lineNum=$(cat TEMP_protein_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no surface/secretory targets were typed
    firstLine=$(head -n1 TEMP_protein_Results.txt)
    printf "\t\t$firstLine\n\t\tNo_Protein_Targets\n" >> "$sampl_out"
    printf "No_Protein_Targets" >> TEMP_table_results.txt
else
    count=0
    prot_target=()
    while read -r line
    do
        count=$(( $count + 1 ))
        if [[ "$count" -eq 1 ]]
        then
            #print the header line to the 'SAMPL' output file
            printf "\t\t$line\n" >> "$sampl_out"
        else
            #If surface/secretory target is greater than the contamination threshold then add that
            #target to the output array 'prot_target'
            printf "\t\t$line\n" >> "$sampl_out"
            justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
            justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
            if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
            then
                echo "Target $justTarget is a match"
                #printf "$justTarget-($justDepth);" >> TEMP_table_results.txt
                prot_target+=("$justTarget($justDepth)")
            fi
        fi
    done < TEMP_protein_Results.txt
    #done < TEMP_no-Tantigen_Results.txt
    #if the output array 'prot_target' is not empty, print out the sorted types to the 'TEMP_table_results.txt' file
    if [ ${#prot_target[@]} -eq 0 ];
    then
        printf "No_Protein_Targets" >> TEMP_table_results.txt
    else
        printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';'
        printf '%s\n' "${prot_target[@]}" | sort | tr '\n' ';' >> TEMP_table_results.txt
    fi
fi
printf "\n" >> "$sampl_out"
printf "\n" >> TEMP_table_results.txt

#cat TEMP_table_results.txt
#cat "$sampl_out" >> "$batch_out"/SAMPL_GBS_"$batch_name"_Typing_Results.txt
#cat TEMP_table_results.txt >> "$batch_out"/TABLE_GBS_"$batch_name"_Typing_Results.txt
#if [[ -e TEMP_newPBP_allele_info.txt ]]
#then
#    cat TEMP_newPBP_allele_info.txt >> "$batch_out"/UPDATR_GBS_"$batch_name"_PBP-DB_updater.txt
#fi


###Unload Modules###
module unload perl/5.12.3
module unload ncbi-blast+/2.2.29
module unload BEDTools/2.17.0
module unload Python/2.7
module unload samtools/0.1.18
module unload bowtie2/2.1.0
module unload freebayes/0.9.21
module unload prodigal/2.60
