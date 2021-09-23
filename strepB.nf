#!/bin/bash -l



max_forks = params.max_forks
db_dir = params.allDB_dir

db=file(db_dir)
strep_agal=file("$db_dir/Streptococcus_agalactiae.fasta")
strep_agal_txt=file("$db_dir/sagalactiae.txt")
sero_gene_db=file("$db_dir/GBS_seroT_Gene-DB_Final.fasta")
lactam_db=file("$db_dir/GBS_bLactam_Ref.fasta")
output_dir="output"

Channel.fromFilePairs("${params.batch_dir}/*_R{1,2}_001.fastq.gz").into { fqPairs1; fqPairs2; fqPairs3; fqPairs4;  fqPairs5 ; fqPairs6}


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
   set val(base), file(trim1), file(trim2) into trimmed_ch1,trimmed_ch2
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
    publishDir "${params.out_dir}/${base}/qc_reports", mode: 'copy', overwrite: true, pattern: "*html"    
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
      file(strep_agal)
      file(strep_agal_txt)
    output:
      set val(base), file(results), file("${name}*bam") into bam_src_ch
      set val(base),  file(results) into bam_src1_ch
    script:
      name = "MLST_$base"
      results = "${name}__mlst__Streptococcus_agalactiae__results.txt"
      """
        srst2 --samtools_args '\\-A'  --mlst_delimiter '_' --input_pe $pair --output $name --save_scores --mlst_db  $strep_agal --mlst_definitions $strep_agal_txt --min_coverage 99.999 
      """
}

process MLSTalleleChecker {
   input:
      set val(base), file(results), file(bam) from bam_src_ch
      file(strep_agal)
   output:
      file ("${base}_new_mlst.txt") optional true into new_mlst_ch
      publishDir "${params.out_dir}/new_mlst/", copy:true, overwrite:true
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


fqPairs4.join(trimmed_ch2)
         .map { base, pair, trim1, trim2 -> [base, pair[0], pair[1], trim1, trim2] }
	 .set { pbp_inputs_ch }


process pbpGeneTyper {
   input:
     set val(base),  file(raw1), file(raw2), file(trim1), file(trim2) from pbp_inputs_ch
     file(db)
   output:
      set val(base), file("TEMP_pbpID_Results.txt") into pbp_res_ch
      set val(base), file("TEMP_pbpID_Results.txt"), file(velvet_output) into pbp_res1_ch
      file("${base}_newPBP_allele_info.txt") optional true into newpbp_ch
      file("velvet_output")
   publishDir "${params.out_dir}/${base}/", overwrite: true, pattern: "velvet_output"
   script:
   /* The trim files are actually used but PBP-Gene_Typer constructs the name from the raw file names */
   """     
   which perl
   which makeblastdb
   PBP-Gene_Typer.pl -1 $raw1 -2 $raw2 -r $db/GBS_bLactam_Ref.fasta -n $base  -s GBS -p 1A,2B,2X
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
   GBS_Res_Typer.pl -1 ${pair[0]} -2 ${pair[1]} -d $db -r GBS_Res_Gene-DB_Final.fasta -n $base
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
   input:
     set val(base),  file(pair) from fqPairs6
     file(db)
   output:
     set val(base), file("TEMP_Surface_Results.txt"), file("BIN_Surface_Results.txt") into surface_res_ch
  script:
  """
   GBS_Surface_Typer.pl -1 ${pair[0]} -2 ${pair[1]} -r $db -p GBS_Surface_Gene-DB_Final.fasta -n $base
  """
}


reports_ch = sero_res_ch.join(bam_src1_ch).join(pbp_res1_ch).join(gbs_target_ch).join(surface_res_ch)




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



abs_output = (new File(params.out_dir)).getCanonicalPath()

process reportGlobal {
  input: 
     file(reps) from reports.toList()
     file(newpbps) from newpbp_ch.toList()
     file(newmlst) from new_mlst_ch.toList()
  output:
     file(rep_name)
  publishDir "${params.out_dir}"
  script:
      rep_name=new java.text.SimpleDateFormat("yyyy-MM-dd-HHmmss").format(new Date())+".xlsx"
      """
      combined_report.py ${params.batch_dir}  $abs_output  $rep_name
      """
}


process multiqc {
    input:
    file ('fastqc/*') from fastqc_results_ch.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"
    publishDir "${params.out_dir}"
    script:
    """
    multiqc .
    """
}
