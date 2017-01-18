#!/bin/bash -l

#export PATH=/scicomp/home/ycm6/TEMP_GBS-Typing:$PATH
temp_path=$(pwd)
export PATH=$PATH:$temp_path:/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/JanOw_Dependencies

## -- begin embedded SGE options --
read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1/job-control.txt)
## -- end embedded SGE options --

###Load Modules###
#. /usr/share/Modules/init/bash
module load perl/5.22.1
module load ncbi-blast+/2.2.29
module load BEDTools/2.17.0
module load freebayes/0.9.21
module load prodigal/2.60
module load cutadapt/1.8
module load srst2/0.1.7

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
out_nameARG=ARG_"$just_name"
out_nameRES=RES_"$just_name"
out_namePLAS=PLAS_"$just_name"

###Call MLST###
srst2 --samtools_args "\-A" --mlst_delimiter '_' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMLST" --save_scores --mlst_db "$allDB_dir/Streptococcus_agalactiae.fasta" --mlst_definitions "$allDB_dir/sagalactiae.txt" --min_coverage 99.999
###Check and extract new MLST alleles###
MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_agalactiae__results.txt "$out_nameMLST"__*.Streptococcus_agalactiae.sorted.bam "$allDB_dir/Streptococcus_agalactiae.fasta"

###Call GBS Serotype###
GBS_serotyper.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/GBS_seroT_Gene-DB_Final.fasta" -n "$out_nameSERO"

###Call GBS bLactam Resistances###
module unload perl/5.22.1
module load perl/5.16.1-MT
PBP-Gene_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/GBS_bLactam_Ref.fasta" -n "$out_namePBP" -s GBS -p 1A,2X
module unload perl/5.16.1-MT
module load perl/5.22.1

###Call GBS Misc. Resistances###
GBS_miscRes_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -m GBS_miscR_Gene-DB_Final.fasta -v GBS_vancR_Gene-DB_Final.fasta -n "$out_nameMISC"

###Type Surface adn Secretory Proteins###
perl /scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/GBS_Scripts_Reference/GBS_surface-secretory_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir" -p GBS_protein_Gene-DB_Final.fasta -n "$out_namePROT"

###Type ARG-ANNOT Resistance Genes###
srst2 --samtools_args '\\-A' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameARG" --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db "$allDB_dir/ARGannot_r1.fasta"

###Type ResFinder Resistance Gene###
srst2 --samtools_args '\\-A' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameRES" --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db "$allDB_dir/ResFinder.fasta"

###Type PlasmidFinder Resistance Gene###
srst2 --samtools_args '\\-A' --input_pe "$readPair_1" "$readPair_2" --output "$out_namePLAS" --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db "$allDB_dir/PlasmidFinder.fasta"

#Add in perl script to find contamination threshold here
contamination_level=10





###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
tabl_out="TABLE_Isolate_Typing_results.txt"
sampl_out="SAMPLE_Isolate__Typing_Results.txt"
res_targ=()
surf_targ=()

printf "$just_name\n" >> "$sampl_out"
printf "$just_name\t" >> "$tabl_out"
###SEROTYPE OUTPUT###
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
	    justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
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
done < "$out_nameMLST"__mlst__Streptococcus_agalactiae__results.txt

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

###MISC. RESISTANCE###
lineNum=$(cat TEMP_miscR_Results.txt | wc -l)
misc_contamination_level=7
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then no misc. resistance were found
    #or there was a error in the read mapping
    firstLine=$(head -n1 TEMP_miscR_Results.txt)
    printf "\t\t$firstLine\n" >> "$sampl_out"
    if grep -q failed MISC_*.log;
    then
        printf "**Failed Read Mapping - Results not accurate**" >> "$tabl_out"
        printf "\t\t**Failed Read Mapping - Results not accurate**\n" >> "$sampl_out"
    else
        printf "No_Resistance" >> "$tabl_out"
        printf "\t\tNo_Resistance\n" >> "$sampl_out"
    fi
else
    extract_arr=(PARCGBS-1 GYRAGBS-1 23SWT-1 23SWT-3 RPLDGBS-1 RPLDGBS-2 RPLVGBS-1 RPLVGBS-2 RPOBgbs-1 RPOBgbs-2 RPOBgbs-3 RPOBgbs-4)
    while read -r line
    do
        #If misc. resistance target is greater than the contamination threshold then add that
        #misc. resistance target to the output array 'misc_target'
	printf "\t\t$line\n" >> "$sampl_out"
        justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
        justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
        justMatchType=$(echo "$line" | awk -F"\t" '{print $2}')
        if [[ $(echo "$justDepth > $misc_contamination_level" | bc) -eq 1 ]]
        then
            #this condition checks if the target is one of the extraction alleles. If so, it will append (extract) to output
            if [[ "$justMatchType" == "imperfect" && " ${extract_arr[@]} " =~ " ${justTarget} " ]]
            then
                #echo "Target $justTarget is a match but needs extraction"
                res_targ+=("$justTarget(extract)")
            elif [[ "$justMatchType" == "identical" ]]
            then
                #echo "Target $justTarget is an identical match"
                res_targ+=("$justTarget")
            else
                #echo "Target $justTarget is a match"
                res_targ+=("$justTarget($justMatchType)")
            fi
        fi
    done < TEMP_miscR_Results.txt
fi

###ARG-ANNOT and ResFinder Resistance Gene Typing Output###
if [[ -s "$out_nameARG"__fullgenes__ARGannot_r1__results.txt ]]
then
    firstLine=$(head -n1 "$out_nameARG"__fullgenes__ARGannot_r1__results.txt)
    printf "\t\t$firstLine\n" >> "$sampl_out"
    i=1
    while read -r line
    do
        test $i -eq 1 && ((i=i+1)) && continue
	printf "\t\t$line\n" >> "$sampl_out"
        #isIdentical="identical"
        #justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
        #if [[ -n "$justDiffs" ]]
        #then
        #    isIdentical="imperfect"
        #fi
        justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
        justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
        if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
        then
            #echo "Target $justTarget is a match"
            res_targ+=("ARG_$justTarget")
        fi
    done < "$out_nameARG"__fullgenes__ARGannot_r1__results.txt
fi
if [[ -s "$out_nameRES"__fullgenes__ResFinder__results.txt ]]
then
    firstLine=$(head -n1 "$out_nameRES"__fullgenes__ResFinder__results.txt)
    printf "\t\t$firstLine\n" >> "$sampl_out"
    i=1
    while read -r line
    do
        test $i -eq 1 && ((i=i+1)) && continue
        printf "\t\t$line\n" >> "$sampl_out"
        #isIdentical="identical"
        #justDiffs=$(echo "$line" | awk -F"\t" '{print $7}')
        #if [[ -n "$justDiffs" ]]
        #then
        #    isIdentical="imperfect"
        #fi
        justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
        justDepth=$(echo "$line" | awk -F"\t" '{print $6}')
        if [[ $(echo "$justDepth > $contamination_level" | bc) -eq 1 ]]
        then
            #echo "Target $justTarget is a match"
            res_targ+=("RF_$justTarget")
        fi
    done < "$out_nameRES"__fullgenes__ResFinder__results.txt
fi

###Surface / Secretory Protein Output (Not including T-Antigens)###
lineNum=$(cat TEMP_protein_Results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    #if the file only contains the header line then surface proteins were found
    #or there was a error in the read mapping
    firstLine=$(head -n1 TEMP_protein_Results.txt)
    printf "\t\t$firstLine\n" >> "$sampl_out"
    if grep -q failed PROT_*.log;
    then
        printf "**Failed Read Mapping - Results not accurate**" >> "$tabl_out"
        printf "\t\t**Failed Read Mapping - Results not accurate**\n" >> "$sampl_out"
    else
        printf "No_Surface_Targets" >> "$tabl_out"
        printf "\t\tNo_Resistance\n" >> "$sampl_out"
    fi
else
    while read -r line
    do
        #If misc. resistance target is greater than the contamination threshold then add that
        #misc. resistance target to the output array 'misc_target'
        printf "\t\t$line\n" >> "$sampl_out"
        justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
        justDepth=$(echo "$line" | awk -F"\t" '{print $4}')
        justMatchType=$(echo "$line" | awk -F"\t" '{print $2}')
        if [[ $(echo "$justDepth > $misc_contamination_level" | bc) -eq 1 ]]
        then
            if [[ "$justMatchType" == "identical" || "$justMatchType" == "imperfect" ]]
            then
                #echo "Target $justTarget is an identical match"
                surf_targ+=("$justTarget")
            else
                #echo "Target $justTarget is a match"
                surf_targ+=("$justTarget($justMatchType)")
            fi
        fi
    done < TEMP_protein_Results.txt
fi




###Print non-bLactam resistance and surface protein matches###
shopt -s nocasematch
###Print non-bLactam resistance matches###
tet_targ="neg"
EC_targ="neg"
FQ_targ="neg"
misc_othr="neg"
for i in "${res_targ[@]}"
do
    echo "input resistance target: $i"
    if [[ "$i" =~ tet ]] ##Extract tet targets##
    then
        if [[ "$tet_targ" == "neg" ]]
        then
            #echo "first tet match"
            tet_targ="$i"
        else
            #echo "new tet match"
            tet_targ+=":$i"
        fi
    elif [[ "$i" =~ erm|mef|lsa|lnu ]] ##Extract EC targets##
    then
        if [[ "$EC_targ" == "neg" ]]
        then
            #echo "first EC match"
            EC_targ="$i"
        else
            #echo "new EC match"
            EC_targ+=":$i"
        fi
    elif [[ "$i" =~ parc|gyra ]] ##Extract FQ target##
    then
        if [[ "$FQ_targ" == "neg" ]]
        then
            #echo "first FQ match"
            FQ_targ="$i"
        else
            #echo "new FQ match"
            FQ_targ+=":$i"
        fi
    else
        if [[ "$misc_othr" == "neg" ]]
        then
            #echo "first other match"
            misc_othr="$i"
        else
            #echo "new other match"
            misc_othr+=":$i"
        fi
    fi
done

###Print surface target matches###
alp_targ="neg"
srr_targ="neg"
pili_targ="neg"
hvga_targ="neg"
for i in "${surf_targ[@]}"
do
    echo "input surface target: $i"
    if [[ "$i" =~ alp|rib ]] ##Extract alp targets##
    then
        if [[ "$alp_targ" == "neg" ]]
        then
            #echo "first alp match"
            alp_targ="$i"
        else
            #echo "new alp match"
            alp_targ+=":$i"
        fi
    elif [[ "$i" =~ srr ]] ##Extract srr target##
    then
        if [[ "$srr_targ" == "neg" ]]
        then
            #echo "first srr match"
            srr_targ="$i"
        else
            #echo "new srr match"
            srr_targ+=":$i"
        fi
    elif [[ "$i" =~ PI ]] ##Extract PI target##
    then
        if [[ "$pili_targ" == "neg" ]]
        then
            #echo "first pili match"
            pili_targ="$i"
        else
            #echo "new pili match"
            pili_targ+=":$i"
        fi
    elif [[ "$i" =~ hvga ]] ##Extract hvga target##
    then
        if [[ "$hvga_targ" == "neg" ]]
        then
            #echo "first hvga match"
            hvga_targ="$i"
        else
            #echo "new hvga match"
            hvga_targ+=":$i"
        fi
    fi
done
shopt -u nocasematch

printf "TET\tEC\tFQ\tOTHER\tALP\tSRR\tPILI\tHVGA\n"
printf "$tet_targ\t$EC_targ\t$FQ_targ\t$misc_othr\t$alp_targ\t$srr_targ\t$pili_targ\t$hvga_targ\n"
printf "$tet_targ\t$EC_targ\t$FQ_targ\t$misc_othr\t$alp_targ\t$srr_targ\t$pili_targ\t$hvga_targ\n" >> "$tabl_out"


###Unload Modules###
module unload perl/5.22.1
module unload ncbi-blast+/2.2.29
module unload BEDTools/2.17.0
module unload freebayes/0.9.21
module unload prodigal/2.60
module unload cutadapt/1.8
module unload srst2/0.1.7
