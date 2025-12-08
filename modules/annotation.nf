/*
Module      : annotation
Description : Performs contig annotation, refinement, VAF estimation, and post-processing.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

/*
  Process: annotate_contigs
    - Performs contig annotation using custom Python script.

    Container:
      - conda-forge::pandas=2.3.0
      - conda-forge::numpy=2.3.0
      - bioconda::pysam=0.23.3
      - conda-forge::intervaltree=3.1.0
*/
process annotate_contigs {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${sample_id}_output/04-Annotation", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/pysam_intervaltree_numpy_pandas:af06dba6fd8ba217'

    input:
      tuple val(sample_id), path(hg38_with_md_bam), path(hg38_with_md_bai), path(ann_info), path(tx_annotation)

    output:
		  tuple val(sample_id), path("annotated_contigs.bam"), emit: anno_bam
		  tuple val(sample_id), path("annotated_contigs.bam.bai"), emit: anno_bai
		  tuple val(sample_id), path("annotated_contigs_info.tsv"), emit: anno_info
      tuple val(sample_id), path("annotation.log"), emit: anno_log
		  tuple val(sample_id), path("annotated_contigs.vcf"), emit: anno_vcf

    script:
    """
    python ${params.code_base}/annotate/LINDTIE_contigs_annotation.py \\
      ${sample_id} \\
      DE_contigs_mapped_to_hg38.bam \\
      ${ann_info} \\
      ${tx_annotation} \\
      annotated_contigs.bam \\
      annotated_contigs_info.tsv \\
      --minClip=${params.min_clip} \\
      --minGap=${params.min_gap} \\
      --minMatch=${params.min_match} \\
      --log annotation.log \\
      > annotated_contigs.vcf
    """
}

//------------------------------------------------------------------------------
/*
  Process: refine_annotation
    - Refines annotated contigs using custom Python script.
  
  Container:
  - conda-forge::pandas=2.3.0
  - conda-forge::numpy=2.3.0
  - bioconda::pysam=0.23.3
  - conda-forge::intervaltree=3.1.0
  - bioconda::bio=1.8.0
  - bioconda::pybedtools=0.12.0
  - bioconda::samtools=1.22
*/
process refine_annotation {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${sample_id}_output/FinalOutput", mode: 'copy', pattern: '*.{vcf,bam,bai}'
    publishDir "${sample_id}_output/FinalOutput/log", mode: 'copy', pattern: '*.log'

    container 'oras://community.wave.seqera.io/library/bio_pybedtools_pysam_samtools_pruned:704dfffe82e57426'

    input:
      tuple val(sample_id), path(anno_info), path(anno_vcf), path(anno_bam), path(anno_bai), path(tx_annotation), path(hg38_fasta)

    output:
		  tuple val(sample_id), path("refined_annotated_contigs.bam"), emit: refined_annotated_contigs_bam
      tuple val(sample_id), path("refined_annotated_contigs.bam.bai"), emit: refined_annotated_contigs_bai
		  tuple val(sample_id), path("refined_annotated_contigs_info.tsv"), emit: refined_annotated_contigs_info
		  tuple val(sample_id), path("refined_annotated_contigs.vcf"), emit: refined_annotated_contigs_vcf
		  tuple val(sample_id), path("refine_annotation.log"), emit: refine_annotation_log
    script:
    """
    python ${params.code_base}/annotate/LINDTIE_refine_annotation.py \\
      ${anno_info} \\
      ${anno_vcf} \\
      ${anno_bam} \\
      ${tx_annotation} \\
      ${hg38_fasta} \\
      refined_annotated_contigs \\
      --minClip ${params.min_clip} \\
      --minGap ${params.min_gap} \\
      --mismatches ${params.splice_motif_mismatch} \\
      --log refine_annotation.log \\
      > refined_annotated_contigs.vcf

    ## Index the refined_annotated_contigs BAM
    samtools index refined_annotated_contigs.bam
    """
}

//------------------------------------------------------------------------------
/*
  Process: filter_refined_annotated_contigs_fasta
    - Filters the de_contigs FASTA to extract only novel contigs based on annotated_info.tsv using custom Python script.

  Container:
    - conda-forge::pandas=2.3.0
    - bioconda::bio=1.8.0
*/
process filter_refined_annotated_contigs_fasta {

    tag "${sample_id}"
    label 'process_short'

    publishDir "${sample_id}_output/FinalOutput", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/bio_pandas:7fb76c90a9b5a085'

    input:
      tuple val(sample_id), path(de_contigs), path(anno_info)

    output:
		  tuple val(sample_id), path("refined_annotated_contigs.fasta"), emit: refined_annotated_contigs_fasta

    script:
    """
    python ${params.code_base}/util/filter_fasta.py \\
      ${de_contigs} \\
      ${anno_info} --col_id contig_id \\
      > refined_annotated_contigs.fasta
    """
}

//------------------------------------------------------------------------------
/*
  Process: estimate_vaf
    - Estimates VAF using custom R script.

    Container:
    - conda-forge::r-dplyr=1.1.4
    - conda-forge::r-data.table=1.17.6
    - bioconda::bioconductor-tximport=1.34.0
*/
process estimate_vaf {

    tag "${sample_id}"
    label 'process_short'

    container 'oras://community.wave.seqera.io/library/bioconductor-tximport_r-data.table_r-dplyr:e41e296cf4f11db2'

    input:
      tuple val(sample_id), path(transcript_matrix), path(case_quant), path(refined_annotated_contigs_info), path(trans_fasta), path(tx2gene)

    output:
      tuple val(sample_id), path("vaf_estimates.txt"), emit: vaf_estimates

    script:
    """
    Rscript ${params.code_base}/annotate/LINDTIE_estimate_VAF.R \\
      ${transcript_matrix} \\
      ${case_quant} \\
      ${refined_annotated_contigs_info} \\
      ${trans_fasta} \\
      ${tx2gene} \\
      vaf_estimates.txt
    """
}

//------------------------------------------------------------------------------
/*
  Process: post_process
    - Filter and collate refined annotated contigs info to create the final TSV using custom Python script.

    Container:
  - conda-forge::pandas=2.3.0
  - conda-forge::numpy=2.3.0
  - bioconda::pysam=0.23.3
  - bioconda::pybedtools=0.12.0
  - conda-forge::intervaltree=3.1.0
  - bioconda::bio=1.8.0
*/
process post_process {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${sample_id}_output/FinalOutput", mode: 'copy', pattern: '*.tsv'
    publishDir "${sample_id}_output/FinalOutput/log", mode: 'copy', pattern: '*.log'

    container 'oras://community.wave.seqera.io/library/bio_pybedtools_pysam_intervaltree_pruned:b44948a3d69b5ec0'

    input:
      tuple val(sample_id), path(refined_annotated_contigs_info), path(refined_annotated_contigs_fasta), path(transcript_de), path(vaf_estimates)

    output:
		  tuple val(sample_id), path("${sample_id}_results.tsv"), emit: results, optional: true
      tuple val(sample_id), path("${sample_id}_all_variants_ranked_results.tsv"), emit: all_variants_ranked_results, optional: true
      tuple val(sample_id), path("${sample_id}_discarded_results.tsv"), emit: discarded_results, optional: true
      tuple val(sample_id), path("final_post_process.log"), emit: final_post_process_log

    script:
    """
    python ${params.code_base}/annotate/LINDTIE_post_process.py \
      ${sample_id} \
      ${refined_annotated_contigs_info} \
      ${refined_annotated_contigs_fasta} \
      ${transcript_de} \
      ${vaf_estimates} \
      --minClip ${params.min_clip} \
      --minGap ${params.min_gap} \
      ${ params.gene_filter ? "--gene_filter ${params.gene_filter}" : "" } \
      ${ params.var_filter ? "--var_filter ${params.var_filter}" : "" } \
      --log final_post_process.log \
      --all_variants_out ${sample_id}_all_variants_ranked_results.tsv \
      --discard_out ${sample_id}_discarded_results.tsv \
      > ${sample_id}_results.tsv
    """
}
