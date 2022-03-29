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

   
params.out_dir = "gas_output"

staged_batch_dir = file (params.batch_dir)
params.strepA_DB="/dataC/CRDM/GAS_Reference_DB"
db_dir = params.strepA_DB

bl_out="blast_frwd_primr-nucl_DB"
nhr=file("${db_dir}/${bl_out}.nhr")
nin=file("${db_dir}/${bl_out}.nin")
nsq=file("${db_dir}/${bl_out}.nsq")

suffix=params.suffix



params.cutadapt_cores=4

println workflow.profile

if (workflow.profile.contains("legacy")) {
  cores  = " "
} else {
  cores = " --cores ${params.cutadapt_cores}"
}



db=file(db_dir)




strep_pyog=file("$db_dir/Streptococcus_pyogenes.fasta")
strep_pyog_txt=file("$db_dir/spyogenes.txt")

sero_gene_db=file("$db_dir/GAS_seroT_Gene-DB_Final.fasta")
lactam_db=file("$db_dir/GAS_bLactam_Ref.fasta")
output_dir=params.out_dir

adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"


Channel.fromFilePairs("${params.batch_dir}/*_R{1,2}_001.fastq.gz").
      ifEmpty { println "No files  match pattern: ${params.batch_dir}/*_R{1,2}_001.fastq.gz";
                System.exit(10) }.
      into { fqPairs1; fqPairs2; fqPairs3; fqPairs4;  fqPairs5 ; fqPairs6}


primer_db=file("$db_dir/frwd_primr-DB_Final.fasta")
emm_db=file("$db_dir/emm_Gene-DB_Final.fasta")


bowtie_sources = Channel.fromPath(["$db_dir/GAS_features_Gene-DB_Final.fasta","$db_dir/GAS_Res_Gene-DB_Final.fasta"])

process makeBowtieIndices {
    input:
    file f from bowtie_sources
    output:
      file("$f.*") into bowtie_indices
      file("sync") into bowtie_sync
    storeDir db_dir
    script:
    """
    bowtie2-build $f $f
    touch sync
    """
}

process makeBlastDBprimer {
    input:
      file(primer_db)
      file(bowtie_sync) // this is just a hack to have a barrier
    output:
      file("${bl_out}.nhr")
      file("${bl_out}.nin")
      file("${bl_out}.nsq")        
    storeDir db_dir
      
    script:
     """
     makeblastdb -in $primer_db -dbtype nucl -out $bl_out
     """
}



process makeBlastEMMr {
    input:
      file(emm_db)
    output:
     file("emm_db") into blast_emm_db_ch
    script:
       bl_out="blast_emm_Gene-nucl_DB"
     """
     makeblastdb -in $emm_db -dbtype nucl  -out emm_db/$bl_out
     """
}

pbp_types = ["2X"]

process makeBlactamDB {
  input:
    file(db)
  each pbp_type from pbp_types
  output:
    file("$db/${blastDB_name}*p*") into blast_blactam_ch
  storeDir db_dir  
  script:
    blastDB_name = "Blast_bLactam_${pbp_type}_prot_DB"
    blast_seq = "GAS_bLactam_${pbp_type}-DB.faa"
  """
   makeblastdb -in $db/$blast_seq -dbtype prot -out $db/$blastDB_name
  """
}  

process cutAdapt1 {
  maxForks max_forks
  cpus params.cutadapt_cores
  errorStrategy 'finish'
  input:
   set val(base), file(pair) from fqPairs1
  output:
   set val(base), file("temp1.fastq"), file("temp2.fastq") into trim1_ch
  script:
   """
     cutadapt  $cores  -b $adapter1 -q 20 --minimum-length 50 \
         --paired-output temp2.fastq -o temp1.fastq $pair
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
    maxForks max_forks
    input:
      set val(base), file(pair) from fqPairs2
      file(strep_pyog)
      file(strep_pyog_txt)
    output:
      set val(base), file(results), file("${name}*bam") into bam_src_ch
      set val(base),  file(results) into bam_src1_ch
    script:
      name = "MLST_$base"
      results = "${name}__mlst__Streptococcus_pyogenes__results.txt"
      """
        hostname
        srst2 --samtools_args '\\-A'  --mlst_delimiter '_' --input_pe $pair \
	      --output $name --save_scores \
	      --mlst_db  $strep_pyog --mlst_definitions $strep_pyog_txt --min_coverage 99.999 
      """
}

process MLSTalleleChecker {
   input:
      set val(base), file(results), file(bam) from bam_src_ch
      file(strep_pyog)
   output:
      stdout ch
      file ("${base}_new_mlst.txt") optional true into new_mlst_ch
   errorStrategy 'finish'
   script:
   """
    echo $bam
    MLST_allele_checkr.pl $results $bam $strep_pyog
    if [ -e Check_Target_Sequence.txt ]; then
       cp Check_Target_Sequence.txt  ${base}_new_mlst.txt
    fi

   """
}



process getVelvetK {
  input:
    set val(base), file(f1), file(f2) from trimmed_ch2  
  output:
    tuple(base), stdout into velvet_k_ch
  """
     velvetk.pl --best --size 1.8M  $f1 $f2 > Kval
     cat Kval
  """
}


process velvet {
  cpus 2
  input:
    tuple val(base), file(f1), file(f2), val(vk) from trimmed_ch3.join(velvet_k_ch)
  output:
  tuple val(base), file("velvet_output") into \
	velvet_bl_ch, velvet_lo_ch, velvet_emm_ch, velvet_gsres_ch, velvet_report_ch
  script:
     k = vk.trim()
  """
   hostname
   export OMP_NUM_THREADS=${params.max_velvet_cpus}
   echo $k > Kval
   VelvetOptimiser.pl -s $k -e $k -o "-scaffolding no" \
                    -f "-shortPaired -separate -fastq $f1 $f2" -d velvet_output;
  """
}


																																												
process blast {
    memory '4.GB'
    maxForks max_forks
  input:
     tuple val(base), file(velvet_output) from velvet_bl_ch
     file(nhr)
     file(nin)
     file(nsq)
  output:
  tuple val(base), file(res) into blast_ch
  script:
    res="contig-vs-frwd_nucl.txt"
    """
     sleep 1
     ls 
     blastn -db blast_frwd_primr-nucl_DB -query velvet_output/contigs.fa -outfmt 6 -word_size 4 \
          -out $res
     """
}


process emmTyper {
   input:
     tuple val(base),  file(contig_v_frwd), file("velvet_output"),file(emm_db_dir) \
  	   from blast_ch.join(velvet_emm_ch).combine(blast_emm_db_ch)
   output:
     set val(base), file(out) into emm_res_ch
   script:
     out= "${base}_emm_type_results.txt"
   """
   core_emm_typer.pl $emm_db_dir/blast_emm_Gene-nucl_DB  $base 
   """
}


fqPairs4.join(trimmed_ch4)
         .map { base, pair, trim1, trim2 -> [base, pair[0], pair[1], trim1, trim2] }
	 .set { pbp_inputs_ch }


process LoTrac {
   errorStrategy 'finish'
   input:
      set val(base),  file(raw1), file(raw2), file(trim1), file(trim2),\
          file(velvet_output)  from pbp_inputs_ch.join(velvet_lo_ch)
      file(db)
    output:
      tuple val(base), file("*fasta") optional true into lotrac_output_ch
      tuple val(base), val("query_seq_bad"), file("BAD_SEQ") optional true into query_length_prob_ch    
   script:
   """    
   printenv PATH 
   LoTrac_target.pl -1 $raw1 -2 $raw2 -q $db/GAS_bLactam_Ref.fasta \
                     -S 2.2M -L 0.95 -f -n $base -o ./
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
   PBP-Gene_TyperWits.pl -1 XX -2 YY -r $db/GAS_bLactam_Ref.fasta -n $base  -s GAS -p 2X
   if [ -e newPBP_allele_info.txt ]; then
      cp newPBP_allele_info.txt ${base}_newPBP_allele_info.txt
   fi
   """
}


g_res_typer="GAS_Res_Typer.pl"
g_res_gene_db="GAS_Res_Gene-DB_Final.fasta"

process gsResTyper {
   input:
   set val(base),  file(pair), file(cut1), file(cut2), file(velvet_output) \
        from fqPairs5.join(trimmed_ch5).join(velvet_gsres_ch)
     file(db)
   output:
     set val(base), file("TEMP_Res_Results.txt"), file("BIN_Res_Results.txt") into gs_res_ch
     tuple val(base), file("BIN_Res_Results.txt") into gs_res1_ch
   script:
   """
     $g_res_typer -1 ${pair[0]} -2 ${pair[1]} -d $db  -r $g_res_gene_db -n $base
   """
}


g_target2MIC = "GAS_Target2MIC.pl"

process gsTarget2Mic {
  input:
    set val(base), file(gas_temp), file(gbs_bin), file(pbp) from \
        gs_res_ch.join(pbp_res_ch)
  output:
     set val(base), file("RES-MIC*") into gs_target_ch
  script:
  """
     $g_target2MIC $gas_temp $base $pbp
  """
}



g_sf_typer="GAS_Features_Typer.pl"
g_sf_typer_parm="-f GAS_features_Gene-DB_Final.fasta"



process gsFeatureTyper {
   input:
     set val(base),  file(pair) from fqPairs6
     file(db)
   errorStrategy 'finish'
   output:
     set val(base), file("TEMP_protein_Results.txt"), file("BIN_Features_Results.txt") \
       into surface_res_ch
  script:
  """
   $g_sf_typer -1 ${pair[0]} -2 ${pair[1]} -d $db $g_sf_typer_parm -n $base
  """
}


reports_ch = emm_res_ch
               .join(bam_src1_ch)
	       .join(pbp_res1_ch)
	       .join(velvet_report_ch)
	       .join(gs_target_ch)
	       .join(surface_res_ch)
	       .join(gs_res1_ch)

process reportSample {

   input:
      set val(base), \
	  file(emm_result), \
	  file(srst_res),
	  file(pbp_res), file(velvet_output),\
	  file(gs_target), file(temp_surface), file(bin_surface),  \
	  file(bin_res)\
	  from reports_ch
   errorStrategy 'finish'
   output:
      tuple file(tabl_out), file(bin_out) into reports mode 'flatten'
   script:
     tabl_out="${base}_TABLE_Isolate_Typing_results.txt"
     bin_out="${base}_BIN_Isolate_Typing_results.txt"
     """
     gas_report.sh  $base $tabl_out $bin_out $emm_result
     ls
     """   
}



low_coverage_ch = Channel.empty() // we may get from future impls of velvet

process reportGlobal {
  if (System.getenv("SHELL") == "/bin/zsh") {
    beforeScript 'source /usr/share/Modules/init/zsh'
  }
  module 'python/3.9'
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
      rep_name="gas-"+new java.text.SimpleDateFormat("yyyy-MM-dd-HHmmss").format(new Date())+".xlsx"
      """
      python3 --version  
      combined_report_generic.py A $suffix ${params.batch_dir} ${params.out_dir}  \
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
    publishDir "${params.out_dir}", mode: 'copy', overwrite: true
    script:
    """
    export PATH=/opt/exp_soft/python39/bin/:${PATH}
    multiqc .
    """
}

process versions {
    publishDir "${params.out_dir}", mode: 'copy', overwrite: true
    script:
      file("software_version.txt")
      """
      echo cutadapt `cutadapt --version` > software_version.txt
      blastn -version   >> software_version.txt
      echo velvet `velvetg | tail -n+2 | head -n1 `  >> software_version.txt
      srst2 --version   >> software_version.txt
      echo freebayes `freebayes --version`   >> software_version.txt
      prodigal -v   >> software_version.txt
      bowtie2 --version | head -n 1   >> software_version.txt
      samtools |& head -n 3   >> software_version.txt
      """
}
