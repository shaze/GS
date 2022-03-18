#!/bin/bash -l


//  (C) University of the Witwatersrand, Johannesburg 2021
//  Scott Hazelhurst
//  MIT License as specified in https://github.com/shaze/GS/blob/master/LICENSE.md

params.max_forks=10
max_forks = params.max_forks

if (params.batch_dir == "0") {
   println "No input batch directory was given"
   System.exit(12);
}

params.out_dir = "spn_output"


params.cutadapt_cores=4

println workflow.profile

if (workflow.profile.contains("legacy")) {
  cores  = " "
} else {
  cores = " --cores ${params.cutadapt_cores}"
}


params.spn_DB="/dataC/CRDM/SPN_Reference_DB"
db_dir = params.spn_DB
db = file(db_dir)
suffix=params.suffix

params.contamination_level=10
contamination_level=params.contamination_level
mod_blactam="MOD_bLactam_resistance.fasta"

strep_pneum=file("${db_dir}/Streptococcus_pneumoniae.fasta")
strep_pneum_txt=file("${db_dir}/spneumoniae.txt")

sero_gene_db = file("${db_dir}/SPN_Sero_Gene-DB_Final.fasta")

adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

newDB = file("${db_dir}/newDB")

ref3=file("${db_dir}/newDB/Ref_PBP_3.faa")

Channel.fromFilePairs("${params.batch_dir}/*_R{1,2}_001.fastq.gz").
      ifEmpty { println "No files  match pattern: ${params.batch_dir}/*_R{1,2}_001.fastq.gz";
                System.exit(10) }.map { it -> [it[0], it[1][0], it[1][1] }  
		    into { fqPairs1; fqPairs2; fqPairs3; fqPairs4; fqPairs5; fqPairs6 }



//  Used for the older workflow
//process sampleSubSeq {
//   input:
//     tuple val(base), file(pair) from fqPairs1
//   output:
//     tuple val(base), file("${out}_1.fq.gz"), file("${out}_2.fq.gz")\
//            into sample_pairs_1_ch, sample_pairs_2_ch, sample_pairs_3_ch
//   script:
//     out="${base}_s"
//   """
//      seqtk sample ${pair[0]} 600000 | gzip > ${out}_1.fq.gz
//      seqtk sample ${pair[1]} 600000 | gzip > ${out}_2.fq.gz
//   """
// }



pbp_types = ["1A","2B","2X"]

process makeBlastBlactam {
    input:
    file(db)
    each pbp_type from pbp_types
    output:
       file("${bl_out}*") into blast_db_ch
    storeDir db_dir
    script:
       bl_out="Blast_bLactam_${pbp_type}_prot_DB"
     """
     makeblastdb -in $db/$mod_blactam -dbtype prot -out $bl_out
     """
}



process cutAdapt1 {
  maxForks max_forks
  cpus params.cutadapt_cores
  errorStrategy 'finish'
  input:
   set val(base), file(r1), file(r2) from fqPairs1
  output:
   set val(base), file("temp1.fastq"), file("temp2.fastq") into trim1_ch
  script:
   """
     cutadapt $cores -b $adapter1 -q 20 --minimum-length 50 \
         --paired-output temp2.fastq -o temp1.fastq $r2 $r1
   """
}

process cutAdapt2 {
  errorStrategy 'finish'
  cpus params.cutadapt_cores    
  input:
   set val(base), file(r1), file(r2) from trim1_ch
  output:
  set val(base), file(trim1), file(trim2) into \
      trimmed_ch1,trimmed_ch2,trimmed_ch3,trimmed_ch4,trimmed_ch5
  script:
   trim1="cutadapt_${base}_R1_001.fastq"
   trim2="cutadapt_${base}_R2_001.fastq"
   """
   cutadapt $cores -b $adapter2 -q 20 --minimum-length 50 \
             --paired-output $trim1 -o $trim2  $r2 $r1
   """
}




// No point in splitting since this is cheap compared to the rest
process fastQC {
    errorStrategy 'finish'
    input:
      set val(base), file(f1), file(f2) from trimmed_ch1
     output:
      set val(base), file("*/*html") into qc_ch
      file ("*/*{zip,html}") into fastqc_results_ch
    publishDir "${params.out_dir}/${base}/qc_reports", mode: 'copy', overwrite: true    
    script:
    """
     mkdir ./${base}_R1_qc ./${base}_R2_qc
     fastqc $f1 --outdir=./${base}_R1_qc
     fastqc $f2 --outdir=./${base}_R2_qc
    """
}



process srstP {
    maxForks 2*max_forks+2 
    cpus 4
    input:
      set val(base), file(pair) from fqPairs2
      file(strep_pneum)
      file(strep_pneum_txt)
    output:
      set val(base), file(results), file("${name}*bam") into bam_src_ch
      set val(base),  file(results) into mlst_report_ch
    script:
      name = "MLST_$base" 
      results = "${name}__mlst__Streptococcus_pneumoniae__results.txt"
      """
        srst2 --threads 4  --samtools_args '\\-A' --mlst_delimiter '_' --input_pe $pair --output $name --save_scores \
              --mlst_db $strep_pneum --mlst_definitions $strep_pneum_txt --min_coverage 99.999
      """
}


process MLSTalleleChecker {
   input:
      set val(base), file(results), file(bam) from bam_src_ch
      file(strep_pneum)
   output:
      file ("${base}_new_mlst.txt") optional true into new_mlst_ch
      publishDir "${params.out_dir}/new_mlst/", mode:params.publish, overwrite:true
   script:
   """
    echo $bam
    MLST_allele_checkr.pl $results $bam $strep_pneum
    if [ -e Check_Target_Sequence.txt ]; then
       cp Check_Target_Sequence.txt  ${base}_new_mlst.txt
    fi
   """
}



process seroTyper {
   maxForks max_forks
   input:
     set val(base),  file(pair) from fqPairs3
     file(sero_gene_db)
   output:
     set val(base), file(result) into sero_res_ch
   script:
   result="OUT_SeroType_Results.txt"
   """
   SPN_Serotyper.pl -1 ${pair[0]} -2 ${pair[1]} -r $sero_gene_db -n $base
   """
}




process getVelvetK {
  input:
    set val(base), file(f1), file(f2) from trimmed_ch2  
  output:
    tuple(base), stdout into velvet_k_ch
  """
     velvetk.pl --best --size 2.2M  $f1 $f2
  """
}


process velvet {
  cpus params.max_velvet_cpus
  input:
    tuple val(base), file(f1), file(f2), val(vk) from trimmed_ch3.join(velvet_k_ch)
  output:
    tuple val(base), file(velvet_output) optional true into velvet_ch, velvet_report_ch, velvet_gs_res
  tuple val(base), val("LOW_COVERAGE"), file("LOW_COVERAGE") optional true into low_coverage_ch
  publishDir "${params.out_dir}/${base}/", mode: params.publish, overwrite: true, pattern: "velvet_output"
  script:
     k = vk.trim()
  """
   export OMP_NUM_THREADS=${params.max_velvet_cpus}
   VelvetOptimiser_strepLabWits.pl -s $k -e $k -o "-scaffolding no"	\
                    -f "-shortPaired -separate -fastq $f1 $f2" -d velvet_output;
  """
}
	 


fqPairs4.join(trimmed_ch4)
         .map { base, pair, trim1, trim2 -> [base, pair[0], pair[1], trim1, trim2] }
	 .set { pbp_inputs_ch }




process LoTrac {
   input:
      set val(base),  file(raw1), file(raw2), file(trim1), file(trim2),\
          file(velvet_output)  from pbp_inputs_ch.join(velvet_ch)
      file(db)
    output:
      tuple val(base), file("*fasta") optional true into lotrac_output_ch
      tuple val(base), val("query_seq_bad"), file("BAD_SEQ") optional true into query_length_prob_ch    
	  tuple val(base), file("EXTRACT*") into mic_process_prep
   script:
   """     
   LoTrac_target.pl -1 $raw1 -2 $raw2 -q $db/$mod_blactam \
                     -S 2.2M -f -n $base -o ./
   """
}


process pbpGeneTyper {
   input:
     tuple  val(base),  file(fastas) from lotrac_output_ch
     file(db)
   output:
      set val(base), file("TEMP_pbpID_Results.txt") into pbp_res_ch
      set val(base), file("TEMP_pbpID_Results.txt") into pbp_res1_ch
      file("${base}_newPBP_allele_info.txt")  optional true into newpbp_ch
   script:
   """     
   PBP-Gene_TyperWits.pl -1 XX -2 YY -r $db/$mod_blactam -n $base  -s SPN -p 1A,2B,2X
   if [ -e newPBP_allele_info.txt ]; then
      cp newPBP_allele_info.txt ${base}_newPBP_allele_info.txt
   fi
   """
}




g_res_typer="SPN_Res_Typer.pl"
g_res_gene_db="SPN_Res_Gene-DB_Final.fasta"





process prepMICAA { 
  input:
      tuple val(base), file(extracts) from mic_process_prep
   output:
      tuple val(base), file("Sample*") into mic_aa
   script:
      """
      PBP_AA_sampledir_to_MIC.sh
      """
}

process buildPBPAAtable { 
  input:
      tuple val(base), file(samples) from mic_aa
      file newDB
   output:
      tuple val(base), file("*PBP_AA_table.csv") into aa_tables 
   script:
      """
         Build_PBP_AA_table.R
      """
}


process getMIC_MMRF { 
  input:
      tuple val(base), file(samples) from aa_tables
      file newDB
  output:
      tuple val(base), file("Sample_PBPtype_MIC2_Prediction.csv") into mic2_prediction
  script:
    """      
    AAtable_to_MIC_MM_RF_EN.R
    """
}

process micFormat {
   input:
     tuple val(base), file(prediction) from mic2_prediction
   output:
      tuple val(base), file(rf_file) into mic2_formatted_ch
   script:
      fout="MIC_formatted_with_SIR.csv"
      rf_file="BLACTAM_MIC_RF_with_SIR.txt"
      """
      MIC_format_with_SIR.R Sample_PBPtype_MIC2_Prediction.csv
      spn_mic_format.sh $fout $rf_file
      """
}


process gsResTyper {
   input:
   set val(base),  file(pair), file(cut1), file(cut2), file(velvet_output) \
        from fqPairs6.join(trimmed_ch5).join(velvet_gs_res)
     file(db)
   output:
     tuple val(base), file("OUT_Res_Results.txt"), file("BIN_Res_Results.txt") into gs_res_ch
     tuple val(base), file("BIN_Res_Results.txt") into gs_res1_ch
   script:
   """
     $g_res_typer -1 ${pair[0]} -2 ${pair[1]} -d $db  -r $g_res_gene_db -n $base
   """
}


process target2MIC {
    input:
      tuple val(base), file(out_res), file(bin_res) from gs_res_ch
    output:
      tuple val(base), file("RES-MIC*")  into mic_res_ch
    script:
    """
      SPN_Target2MIC.pl $out_res $base
    """
}

    
// process miscResistance {

//   input:
//      tuple val(base), file(f1), file(f2) from sample_pairs_2_ch
//      file (db)
//   output:
//      tuple val(base), file(out) into misc_res_ch
//   script:
//     out="misc_$base"
//   """
//     SPN_miscRes_Typer.pl -1 $f1 -2 $f2 -r $db -m miscDrug_Gene-DB_Final.fasta -v vanDrug_Gene-DB_Final.fasta -n $outC
//   """
// }






process multiqc {
    input:
    file ('fastqc/*') from fastqc_results_ch.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"
    publishDir "${params.out_dir}", mode: params.publish, overwrite: true
    script:
    """
    multiqc .
    """

}


reports_ch=
  sero_res_ch.
   join(mlst_report_ch).
join(pbp_res_ch).
join(mic2_formatted_ch).
join(mic_res_ch).
join(velvet_report_ch)

process reportSample {

    input:
       set val(base), \
           file(out_sero), \
  	   file(srst_res),\
           file(pbp_res),\
	   file(mic2_fmt),\
	   file(mic_res),\
	   file(velvet_output) \
    	      from reports_ch
      errorStrategy 'finish'
      output:
        file(tabl_out) into reports mode 'flatten'
     publishDir "${params.out_dir}/", mode: 'copy', overwrite: true    
     script:
      tabl_out="${base}_TABLE_Isolate_Typing_results.txt"
      """
      spn_report.sh  $base $base $tabl_out 
      """   
 }



// low_coverage_ch = Channel.empty() // we may get from future impls of velvet

process reportGlobal {
   input: 
      file(reps) from reports.flatten().toList()
      file(newpbps) from newpbp_ch.toList()
      file(db)
       val low_cov from low_coverage_ch.mix (query_length_prob_ch).map { it[0]+"-"+it[1] }.ifEmpty("None").toList()
    output:
      file(rep_name)
   publishDir "${params.out_dir}", mode: 'copy', overwrite: true
   script:
       lc=(low_cov).join(",")
       rep_name="spn-"+new java.text.SimpleDateFormat("yyyy-MM-dd-HHmmss").format(new Date())+".xlsx"
       """
       combined_report_generic.py SPN $suffix ${params.batch_dir} ${params.out_dir}  \
                       $rep_name $lc
       echo $low_cov > see
       """
 }

