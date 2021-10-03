#!/bin/bash -l



max_forks = params.max_forks
db_dir = params.allDB_dir

db=file(db_dir)

staged_batch_dir = file (params.batch_dir)



strep_agal=file("$db_dir/Streptococcus_agalactiae.fasta")
strep_agal_txt=file("$db_dir/sagalactiae.txt")
sero_gene_db=file("$db_dir/GBS_seroT_Gene-DB_Final.fasta")
lactam_db=file("$db_dir/GBS_bLactam_Ref.fasta")
output_dir="gbs_typing"
suffix=params.suffix


Channel.fromFilePairs("${params.batch_dir}/*_R{1,2}_001*"+suffix)
       .into { fqPairs1; fqPairs2; fqPairs3; fqPairs4;  fqPairs5 ; fqPairs6}


process cutAdapt1 {
  maxForks max_forks
  cpus 4
  errorStrategy 'finish'
  input:
   set val(base), file(pair) from fqPairs1
  output:
   set val(base), file("temp1.fastq"), file("temp2.fastq") into trim1_ch
  script:
   """
     cutadapt --cores 4 -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output temp2.fastq -o temp1.fastq $pair
   """
}

process cutAdapt2 {
   errorStrategy 'finish'
  input:
   set val(base), file(r1), file(r2) from trim1_ch
  output:
  set val(base), file(trim1), file(trim2) into \
         trimmed_ch1,trimmed_ch2,trimmed_ch3,trimmed_ch4
  script:
   trim1="cutadapt_${base}_S1_L001_R1_001.fastq"
   trim2="cutadapt_${base}_S1_L001_R2_001.fastq"
   """
   cutadapt --cores 4 -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output $trim1 -o $trim2  $r1 $r2
   """
}


process fastQC {
    errorStrategy 'finish'
    input:
      set val(base), file(f1), file(f2) from trimmed_ch1
    output:
      set val(base), file("*/*html") into qc_ch
      file ("*/*{zip,html}") into fastqc_results_ch
    publishDir "${params.out_dir}/${base}/qc_reports", mode: params.publish, overwrite: true
    script:
    """
     mkdir ./${base}_R1_qc ./${base}_R2_qc
     fastqc $f1 --outdir=./${base}_R1_qc
     fastqc $f2 --outdir=./${base}_R2_qc
    """
}


process srstP {
    maxForks max_forks
    cpus 4
    input:
      set val(base), file(pair) from fqPairs2
      file(strep_agal)
      file(strep_agal_txt)
    output:
      set val(base), file(results), file("${name}*bam") into bam_src_ch
      set val(base),  file(results) into bam_src1_ch
    script:
      name = "MLST_$base"
      results = "${name}__mlst__Streptococcus_agalactiae__results.txt"
      """
        srst2 --threads 4 --samtools_args '\\-A'  --mlst_delimiter '_' --input_pe $pair --output $name --save_scores --mlst_db  $strep_agal --mlst_definitions $strep_agal_txt --min_coverage 99.999 
      """
}

process MLSTalleleChecker {
   input:
      set val(base), file(results), file(bam) from bam_src_ch
      file(strep_agal)
   output:
      file ("${base}_new_mlst.txt") optional true into new_mlst_ch
      publishDir "${params.out_dir}/new_mlst/", mode:params.publish, overwrite:true
   script:
   """
    echo $bam
    MLST_allele_checkr.pl $results $bam $strep_agal
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
   result="TEMP_SeroType_Results.txt"
   """
   GBS_Serotyper.pl -1 ${pair[0]} -2 ${pair[1]} -r $sero_gene_db -n $base
   """
}


fqPairs4.join(trimmed_ch4)
         .map { base, pair, trim1, trim2 -> [base, pair[0], pair[1], trim1, trim2] }
	 .set { pbp_inputs_ch }




process getVelvetK {
  input:
    set val(base), file(f1), file(f2) from trimmed_ch2  
  output:
    tuple(base), stdout into velvet_k_ch
  """
     velvetk.pl --best --size 1.8M  $f1 $f2
  """
}


process velvet {
  cpus params.max_velvet_cpus
  input:
    tuple val(base), file(f1), file(f2), val(vk) from trimmed_ch3.join(velvet_k_ch)
  output:
  tuple val(base), file(velvet_output) optional true into velvet_ch, velvet_report_ch
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
	 


	     
	 

process LoTrac {
   input:
      set val(base),  file(raw1), file(raw2), file(trim1), file(trim2),\
          file(velvet_output)  from pbp_inputs_ch.join(velvet_ch)
      file(db)
    output:
      tuple val(base), file("*fasta") optional true into lotrac_output_ch
      tuple val(base), val("query_seq_bad"), file("BAD_SEQ") optional true into query_length_prob_ch    
   script:
   """     
   LoTrac_target.pl -1 $raw1 -2 $raw2 -q $db/GBS_bLactam_Ref.fasta \
                     -S 2.2M -L 0.95 -f -n $base -o ./
   """
}

process pbpGeneTyper {
   input:
     tuple val(base),  file(fastas) from lotrac_output_ch
     file(db)
   output:
      set val(base), file("TEMP_pbpID_Results.txt") into pbp_res_ch
      set val(base), file("TEMP_pbpID_Results.txt") into pbp_res1_ch
      file("${base}_newPBP_allele_info.txt") optional true into newpbp_ch
  errorStrategy "finish"
  script:
   """     
   PBP-Gene_TyperWits.pl -1 XX -2 YY -r $db/GBS_bLactam_Ref.fasta \
                         -n $base  -s GBS -p 1A,2B,2X
   if [ -e newPBP_allele_info.txt ]; then
      cp newPBP_allele_info.txt ${base}_newPBP_allele_info.txt
   fi
   """
}

process gbsResTyper {
   input:
     set val(base),  file(pair) from fqPairs5
     file(db)
   output:
     set val(base), file("TEMP_Res_Results.txt") into gbs_res_ch
   script:
   """
   GBS_Res_Typer.pl -1 ${pair[0]} -2 ${pair[1]} -d $db \
                    -r GBS_Res_Gene-DB_Final.fasta -n $base
   """
}




process gbsTarget2Mic {
  input:
     set val(base), file(gbs), file(pbp) from gbs_res_ch.join(pbp_res_ch)
  output:
     set val(base), file("RES-MIC*") into gbs_target_ch
  script:
  """
    GBS_Target2MIC.pl $gbs $base $pbp
  """
}



process gbsSurfaceTyper {
   maxForks max_forks
   input:
     set val(base),  file(pair) from fqPairs6
     file(db)
   output:
     set val(base), file("TEMP_Surface_Results.txt"), file("BIN_Surface_Results.txt") into surface_res_ch
  script:
  """
   GBS_Surface_Typer.pl -1 ${pair[0]} -2 ${pair[1]} \
               -r $db -p GBS_Surface_Gene-DB_Final.fasta -n $base
  """
}



reports_ch = sero_res_ch.
                 join(bam_src1_ch).
		 join(velvet_report_ch).
		 join(pbp_res1_ch).
		 join(gbs_target_ch).
		 join(surface_res_ch)




process reportSample {

   input:
      set val(base), \
	  file(sero_result), \
          file(srst_res), \
	  file(pbp_res), file(vo),\
	  file(gbs_target), file(temp_surface), file(bin_surface)  \
	  from reports_ch
   output:
      file("${base}_*") into reports mode 'flatten'
   script:
     tabl_out="TABLE_Isolate_Typing_results.txt"
     bin_out="BIN_Isolate_Typing_results.txt"
     """
     report.sh  $base 
     cp $tabl_out ${base}_${tabl_out}
     cp $bin_out ${base}_${bin_out}
     """   
}


/* These paths are required because they appear in the output report */


process reportGlobal {
  input: 
     file(reps) from reports.toList()
     file(newpbps) from newpbp_ch.toList()
     file(newmlst) from new_mlst_ch.toList()
     file(staged_batch_dir)
     val low_cov from low_coverage_ch.mix (query_length_prob_ch).map { it[0]+"-"+it[1] }.ifEmpty("None").toList()
  output:
     file(rep_name)
  publishDir "${params.out_dir}", mode: 'copy', overwrite: true
  script:
      lc=(low_cov).join(",")
      rep_name="gbs-"+new java.text.SimpleDateFormat("yyyy-MM-dd-HHmmss").format(new Date())+".xlsx"
      """
      combined_report.py $suffix ${params.batch_dir} ${params.out_dir}  \
                      $rep_name $lc
      echo $low_cov > see
      """
}


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

