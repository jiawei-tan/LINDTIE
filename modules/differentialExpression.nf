/*
Module      : differentialExpression
Description : Performs differential expression analysis, novel contig detection, and filtering.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

/*
  Process: build_transcript_matrix
    - Gathers all .quant files (case + controls) per case sample, and runs 
      custom Python script to create the transcript count matrix.
    
    Container:
      - conda-forge::pandas=2.3.0
*/
process build_transcript_matrix {

    tag "${sample_id}"
    label 'process_medium'

    container 'oras://community.wave.seqera.io/library/pandas:2.3.0--ca05b8eafb33dbea'

    input:
      tuple val(sample_id), path(case_quant), path(control_quants), val(sample_names)

    output:
      tuple val(sample_id), path("transcript_counts_matrix.tsv"), emit: transcript_matrix

    script:
    def sample_names_str = sample_names.join(',')
    // Handle the case where control_quants might be a single file or a list
    def control_files_str = control_quants instanceof List ? control_quants.join(' ') : control_quants.toString()
    def all_quant_files = "${case_quant} ${control_files_str}"

    """
    python ${params.code_base}/DE/LINDTIE_build_transcript_matrix.py \\
    ${all_quant_files} "${sample_names_str}" transcript_counts_matrix.tsv
    """
}

//------------------------------------------------------------------------------
/*
  Process: compare_transcript_oarfish (differential expression / novel contig detection)
    - If RUN_DE is true (default), runs custom R script to perform DE analysis. 
      Otherwise calls custom Python script to perform novel contig detection.

    Container:
      - conda-forge::r-dplyr=1.1.4
      - conda-forge::r-data.table=1.17.6
      - bioconda::bioconductor-edger=4.4.0
      - conda-forge:r-jsonlite=2.0.0 
      - conda-forge:r-readr=2.1.5 
      - conda-forge::r-arrow=19.0.1
      - conda-forge::pandas=2.3.0 (added for Python script)
*/
process compare_transcript_oarfish {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${sample_id}_output/03-DifferentialExpression", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/bioconductor-edger_r-arrow_r-data.table_r-dplyr_pruned:5cb8cbf534468a05'

    input:
      tuple val(sample_id), path(case_quant), path(control_quants), path(case_meta), path(control_metas), path(case_parquet), path(control_parquets), path(trans_fasta), path(transcript_matrix)
    
    output:
      tuple val(sample_id), path("DE.log"), emit: de_log, optional: true
      tuple val(sample_id), path("DE_transcript_significant.txt"), emit: de_results, optional: true
      tuple val(sample_id), path("DE_transcript_full_results.txt"), emit: de_full_results, optional: true
      tuple val(sample_id), path("DE_MDS_plot.png"), emit: mds_plot, optional: true
      tuple val(sample_id), path("DE_MD_plot.png"), emit: md_plot, optional: true
      tuple val(sample_id), path("DE_QLDisp_plot.png"), emit: ql_disp_plot, optional: true

    script:
    // Use proper boolean check - defaults to true if RUN_DE is not set
    if (params.RUN_DE ?: true) {
      """
      # Create oarfish_output directory structure
      mkdir -p oarfish_output
      
      # Copy case files
      cp ${case_quant} oarfish_output/
      cp ${case_meta} oarfish_output/
      cp ${case_parquet} oarfish_output/
      
      # Copy control files
      for file in ${control_quants}; do
          cp \$file oarfish_output/
      done
      
      for file in ${control_metas}; do
          cp \$file oarfish_output/
      done
      
      for file in ${control_parquets}; do
          cp \$file oarfish_output/
      done

      # List files for debugging
      echo "Files in oarfish_output directory:"
      ls -la oarfish_output/

      cp ${params.code_base}/DE/LINDTIE_de_methods.R ./

      Rscript ${params.code_base}/DE/LINDTIE_compare_transcript.R \\
        ${sample_id} \\
        oarfish_output \\
        ${trans_fasta} \\
        DE_transcript_significant.txt \\
        --FDR=${params.fdr} \\
        --minCPM=${params.min_cpm} \\
        --minLogFC=${params.min_logfc}
      """
    } else {
      // This block runs when params.RUN_DE is false
      """
      echo "Running novel contig detection instead of differential expression analysis"
      
      python ${params.code_base}/DE/LINDTIE_get_novel_contigs.py \\
        ${transcript_matrix} ${trans_fasta} \\
        > DE_transcript_significant.txt
      
      # Create empty files for optional outputs that won't be generated
      touch DE.log
      touch DE_transcript_full_results.txt
      touch DE_MDS_plot.png
      touch DE_MD_plot.png
      touch DE_QLDisp_plot.png
      """
    }
}

//------------------------------------------------------------------------------
/*
  Process: filter_de_contigs
    - Filters significant contigs (DE contigs) from the assembled transcripts using 
      custom Python script.

  Container:
      - conda-forge::pandas=2.3.0
      - bioconda::bio=1.8.0
*/
process filter_de_contigs {

    tag "${sample_id}"
    label 'process_short'

    publishDir "${sample_id}_output/03-DifferentialExpression", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/bio_pandas:7fb76c90a9b5a085'

    input:
      tuple val(sample_id), path(assembled_fa), path(de_results)

    output:
      tuple val(sample_id), path("DE_contigs.fasta"), emit: de_contigs

    script:
    """
    python ${params.code_base}/util/filter_fasta.py \\
      ${assembled_fa} \\
      ${de_results} --col_id transcript_id \\
      > DE_contigs.fasta
    """
}

//------------------------------------------------------------------------------
/*
  Process: align_contigs_to_genome
    - Aligns DE contigs to the hg38 genome, sorts, and indexes the BAM file.

  Container:
  - bioconda::minimap2=2.30
  - bioconda::samtools=1.22
*/
process align_contigs_to_genome {

    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output/03-DifferentialExpression", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/minimap2_samtools:e98addfcfd60e8e7'
		
    input:
      tuple val(sample_id), path(de_contigs), path(hg38_fasta), path(hg38_splice)

    output:
      tuple val(sample_id), path("DE_contigs_mapped_to_hg38.bam"), emit: hg38_with_md_bam
      tuple val(sample_id), path("DE_contigs_mapped_to_hg38.bam.bai"), emit: hg38_with_md_bai

    script:
    """
    set -e

    ## 1) Align with minimap2 (splice)
    minimap2 -t ${task.cpus} -ax splice --secondary=no --junc-bed ${hg38_splice} ${hg38_fasta} ${de_contigs} | \\
      samtools view -@ ${task.cpus} -bS - | \\
      samtools sort -@ ${task.cpus} -o DE_contigs_mapped_to_hg38_sorted.bam

    ## 2) Add MD tags and index
    samtools calmd -b DE_contigs_mapped_to_hg38_sorted.bam ${hg38_fasta} \\
      > DE_contigs_mapped_to_hg38.bam

    samtools index DE_contigs_mapped_to_hg38.bam
    """
}
