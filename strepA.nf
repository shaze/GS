#!/bin/bash -l

params.max_forks=10
max_forks = params.max_forks

params.out_dir = "gas_output"
params.batch_dir="/dataC/CRDM/gas_test_small/"
params.allDB_dir="/dataC/GAS_Reference_DB"
db_dir = params.allDB_dir
suffix=params.suffix




db=file(db_dir)

strep_pyog=file("$db_dir/Streptococcus_pyogenes.fasta")
strep_pyog_txt=file("$db_dir/spyogenes.txt")

sero_gene_db=file("$db_dir/GAS_seroT_Gene-DB_Final.fasta")
lactam_db=file("$db_dir/GAS_bLactam_Ref.fasta")
output_dir=params.out_dir

adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"


Channel.fromFilePairs("${params.batch_dir}/*_R{1,2}_001.fastq.gz").into { fqPairs1; fqPairs2; fqPairs3; fqPairs4;  fqPairs5 ; fqPairs6}


primer_db=file("$db_dir/frwd_primr-DB_Final.fasta")
emm_db=file("$db_dir/emm_Gene-DB_Final.fasta")


bowtie_sources = Channel.fromPath(["$db_dir/GAS_features_Gene-DB_Final.fasta","$db_dir/GAS_Res_Gene-DB_Final.fasta"])

process makeBowtieIndices {
    input:
    file f from bowtie_sources
    output:
      file("$f.*") into bowtie_indices
    storeDir db_dir
    script:
    """
    bowtie2-build $f $f
    """
}

process makeBlastDBprimer {
    input:
      file(primer_db)
      file(x) from bowtie_indices.toList() // this is just a hack to have a barrier
    output:
     file("${bl_out}*") into blast_db_ch
    //storeDir db_dir
    script:
       bl_out="blast_frwd_primr-nucl_DB"
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
    file("$db/${blastDB_name}*pto") into blast_blactam_ch
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
  cpus 4
  errorStrategy 'finish'
  input:
   set val(base), file(pair) from fqPairs1
  output:
   set val(base), file("temp1.fastq"), file("temp2.fastq") into trim1_ch
  script:
   """
     cutadapt --cores 4 -b $adapter1 -q 20 --minimum-length 50 \
         --paired-output temp2.fastq -o temp1.fastq $pair
   """
}

process cutAdapt2 {
   errorStrategy 'finish'
  input:
   set val(base), file(r1), file(r2) from trim1_ch
  output:
  set val(base), file(trim1), file(trim2) into \
      trimmed_ch1,trimmed_ch2,trimmed_ch3,trimmed_ch4,trimmed_ch5
  script:
   trim1="cutadapt_${base}_R1_001.fastq"
   trim2="cutadapt_${base}_R2_001.fastq"
   """
   cutadapt --cores 4 -b $adapter2 -q 20 --minimum-length 50 \
             --paired-output $trim1 -o $trim2  $r1 $r2
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
   errorStrategy 'finish'
   script:
   """
    echo $bam
    MLST_allele_checkr.pl $results $bam $strep_pyog
   """
}



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
  cpus 2
  input:
    tuple val(base), file(f1), file(f2), val(vk) from trimmed_ch3.join(velvet_k_ch)
  output:
  tuple val(base), file(velvet_output) into \
	velvet_bl_ch, velvet_lo_ch, velvet_emm_ch, velvet_gsres_ch, velvet_report_ch
  script:
     k = vk.trim()
  """
   export OMP_NUM_THREADS=${params.max_velvet_cpus}
   VelvetOptimiser.pl -s $k -e $k -o "-scaffolding no" \
                    -f "-shortPaired -separate -fastq $f1 $f2" -d velvet_output;
  """
}

process blast {
  input:
    tuple val(base), file(contig), \
  	  file(bl1), file(bl2), file(bl3), file(bl4), file(bl5), file(bl6), file(bl7) \
	  from velvet_bl_ch.combine(blast_db_ch)
  output:
  tuple val(base), file(res) into blast_ch
  script:
    res="contig-vs-frwd_nucl.txt"
    """
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
   input:
      set val(base),  file(raw1), file(raw2), file(trim1), file(trim2),\
          file(velvet_output)  from pbp_inputs_ch.join(velvet_lo_ch)
      file(db)
    output:
      tuple val(base), file("*fasta") optional true into lotrac_output_ch
      tuple val(base), val("query_seq_bad"), file("BAD_SEQ") optional true into query_length_prob_ch    
   script:
   """     
   LoTrac_target.pl -1 $raw1 -2 $raw2 -q $db/GAS_bLactam_Ref.fasta \
                     -S 2.2M -L 0.95 -f -n $base -o ./
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

