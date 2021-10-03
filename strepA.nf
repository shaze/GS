#!/bin/bash -l

params.max_forks=10
max_forks = params.max_forks

params.out_dir = "gas_output"
params.batch_dir="/dataC/CRDM/gas_test_small/"
params.allDB_dir="/dataC/GAS_Reference_DB"
db_dir = params.allDB_dir




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
      file(x) from bowtie_indices.toList() // this is just a hack to hack 
    output:
     file("${bl_out}*") into blast_db_ch
    storeDir db_dir
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
     file("${bl_out}*") into blast_emm_db_ch
    storeDir db_dir
    script:
       bl_out="blast_emm_Gene-nucl_DB"
     """
     makeblastdb -in $emm_db -dbtype nucl  -out $bl_out
     """
}

pbp_types = ["2X"]

process makeBlactamDB {
  input:
    file(db)
  each pbp_type from pbp_types
  output:
    file("$db/${blastDB_name}*pto") into blast_blactam_ch
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
      trimmed_ch1,trimmed_ch2,trimmed_ch3,trimmed_ch4
  script:
   trim1="cutadapt_${base}_R1_001.fastq"
   trim2="cutadapt_${base}_R2_001.fastq"
   """
   cutadapt --cores 4 -b $adapter2 -q 20 --minimum-length 50 \
             --paired-output $trim1 -o $trim2  $r1 $r2
   """
}


process fastQC {
    errorStrategy 'finish'
    input:
      set val(base), file(f1), file(f2) from trimmed_ch1
    output:
      set val(base), file("*/*html") into qc_ch
    publishDir "${params.out_dir}/qc_reports", mode: 'copy', overwrite: true    
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
        srst2 --samtools_args '\\-A'  --mlst_delimiter '_' --input_pe $pair --output $name --save_scores --mlst_db  $strep_pyog --mlst_definitions $strep_pyog_txt --min_coverage 99.999 
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
  tuple val(base), file(velvet_output) into velvet_ch, velvet_report_ch
  script:
     k = vk.trim()
  """
    VelvetOptimiser.pl -s $k -e $k -o "-scaffolding no" \
                    -f "-shortPaired -separate -fastq $f1 $f2" -d velvet_output;
  """
}

process blast {
  input:
    tuple val(base), file(contig), file(db) from velvet_ch.combine(blast_db_ch)
  output:
  tuple file(base), file(res) into blast_ch
  script:
    res="contig-vs-frwd_nucl.txt"
    """
     blastn -db $blastdb -query velvet_output/contigs.fa -outfmt 6 -word_size 4 \
          -out $res
     """
}


process emmTyper {
   input:
     tuple val(base),  file(contig_v_frwd) from blast_ch.combine(blast_emm_db_ch)
   output:
     set val(base), file(out) into emm_res_ch
   script:
     out= "${base}__emm-Type__Results.txt"
   """
   core_emm_typer.pl $emm_db_dir  $base 
   """
}


fqPairs4.join(trimmed_ch4)
         .map { base, pair, trim1, trim2 -> [base, pair[0], pair[1], trim1, trim2] }
	 .set { pbp_inputs_ch }


process pbpGeneTyper {
   input:
     tuple  val(base),  file(raw1), file(raw2), file(trim1), file(trim2),	\
            file(blactam)  from pbp_inputs_ch.combine(blast_blactam_ch)
     file(db)
   output:
      set val(base), file("TEMP_pbpID_Results.txt") into pbp_res_ch
      set val(base), file("TEMP_pbpID_Results.txt"),
          file("newPBP_allele_info.txt") into pbp_res1_ch
   script:
   /* The trim files are actually used but PBP-Gene_Typer constructs the name from the raw file names */
   """     
   which perl
   which makeblastdb
   PBP-Gene_Typer.pl -1 $raw1 -2 $raw2 -r $db/GAS_bLactam_Ref.fasta -n $base  -s GAS -p 2X
   touch newPBP_allele_info.txt
   """
}


g_res_typer="GAS_Res_Typer.pl"
g_res_gene_db="GAS_Res_Gene-DB_Final.fasta"

process gsResTyper {
   input:
     set val(base),  file(pair) from fqPairs5
     file(db)
   output:
     set val(base), file("TEMP_Res_Results.txt") into gs_res_ch
   script:
   """
     $g_res_typer -1 ${pair[0]} -2 ${pair[1]} -d $db  -r $g_res_gene_db -n $base
   """
}


g_target2MIC = "GAS_Target2MIC.pl"

process gsTarget2Mic {
  input:
     set val(base), file(gbs), file(pbp) from gs_res_ch.join(pbp_res_ch)
  output:
     set val(base), file("RES-MIC*") into gs_target_ch
  script:
  """
     $_target2MIC $base $pbp
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
     set val(base), file("TEMP_protein_Results.txt"), file("BIN_Features_Results.txt") into surface_res_ch
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

process reportSample {

   input:
      set val(base), \
	  file(emm_result), \
	  file(srst_res),
	  file(pbp_res), file(new_pbp), file(velvet_output),\
	  file(gs_target), file(temp_surface), file(bin_surface)  \
	  from reports_ch
   output:
      set file("${base}_*"), file("${base}_${bin_out}") into reports
   publishDir "${params.out_dir}/", mode: 'copy', overwrite: true    
   script:
     tabl_out="TABLE_Isolate_Typing_results.txt"
     bin_out="BIN_Isolate_Typing_results.txt"
     """
     gas_report.sh  $base 
     cp $tabl_out ${base}_${tabl_out}
     cp $bin_out ${base}_${bin_out}
     if [ ` /usr/bin/stat -c '%b' $new_pbp` -gt 0 ]; then mv newPBP_allele_info.txt ${base}_newPBP_allele_info.txt; fi
     """   
}

