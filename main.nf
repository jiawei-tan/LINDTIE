#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************
                            LINDTIE
 A Nextflow pipeline to identify aberrant transcripts in cancer using 
 long-read RNA-seq data.

 Authors: Jia Wei Tan
 Contact: tan.j@wehi.edu.au
***************************************************************************/

// Parameter declarations
params.cases_fastq_dir = null
params.controls_fastq_dir = null
params.trans_fasta = null
params.hg38_fasta = null
params.hg38_splice_junctions = null
params.ann_info = null
params.tx_annotation = null
params.tx2gene = null
params.fdr = 0.05
params.min_cpm = 0.1
params.min_logfc = 2
params.min_clip = 20
params.min_gap = 7
params.min_match = "30,0.3"
params.splice_motif_mismatch = 1
params.oarfish_num_bootstraps = 10
params.help = false

if (params.help) {
    log.info """
========================================
            LINDTIE v0.1.0
========================================

General usage:
    nextflow run LINDTIE/main.nf -params-file LINDTIE/params.yaml -profile singularity

Example usage:
    nextflow run LINDTIE/main.nf \\
    -params-file LINDTIE/params.yaml \\
    -profile singularity \\
    --fdr 0.05 \\    
    --min_cpm 0.1 \\
    --min_logfc 2 \\
    --min_clip 20 \\
    --min_gap 7 \\
    --min_match "30,0.3" \\
    --splice_motif_mismatch 1 \\
    --oarfish_num_bootstraps 10 \\
    --gene_filter "TP53,BRCA2" \\
    --var_filter "DEL,INS"

Optional parameters (set via "--params" at runtime or edit the params.yaml file):
--fdr                       : FDR threshold (e.g. 0.05)
--min_cpm                   : Minimum CPM (e.g. 0.1)
--min_logfc                 : Minimum logFC (e.g. 2)
--min_clip                  : Minimum clip length (e.g. 20)
--min_gap                   : Minimum gap (e.g. 7)
--min_match                 : Minimum match (e.g. "30,0.3")
--splice_motif_mismatch     : Splice motif mismatch (e.g. 1)
--oarfish_num_bootstraps    : Number of bootstraps for Oarfish (e.g. 10)
--gene_filter               : List of genes to filter (e.g. "TP53,BRCA2")
--var_filter                : List of variant types to filter (e.g. "DEL,INS")
--help                      : Show this help message
    """
    exit 0
}

/************************** 
* MODULES
**************************/

include { decompress_case_reads } from './modules/decompress'
include { decompress_control_reads } from './modules/decompress'
include { assembly } from './modules/assembly'
include { merge_refTrans_assembly } from './modules/assembly'
include { case_align_quant } from './modules/quantification'
include { control_align_quant } from './modules/quantification'
include { build_transcript_matrix } from './modules/differentialExpression'
include { compare_transcript_oarfish } from './modules/differentialExpression'
include { filter_de_contigs } from './modules/differentialExpression'
include { align_contigs_to_genome } from './modules/differentialExpression'
include { annotate_contigs } from './modules/annotation'
include { refine_annotation } from './modules/annotation'
include { filter_refined_annotated_contigs_fasta } from './modules/annotation'
include { estimate_vaf } from './modules/annotation'
include { post_process } from './modules/annotation'

/************************** 
* WORKFLOW
**************************/
workflow {

    // Define input channels
    ch_case_reads = Channel
        .fromPath("${params.cases_fastq_dir}/*.{fasta,fasta.gz,fa,fa.gz,fastq,fastq.gz,fq,fq.gz}")
        .map { file -> 
            def name      = file.getName()
            def sample_id = name.replaceFirst(/\.(fastq|fasta|fq|fa)(\.gz)?$/, '')
            tuple(sample_id, file)
        }

    ch_control_reads = Channel
        .fromPath("${params.controls_fastq_dir}/*.{fasta,fasta.gz,fa,fa.gz,fastq,fastq.gz,fq,fq.gz}")
        .map { file -> 
            def name       = file.getName()
            def control_id = name.replaceFirst(/\.(fastq|fasta|fq|fa)(\.gz)?$/, '')
            tuple(control_id, file)
        }

    // Process each case sample
    // Decompress reads if needed
    ch_decompressed_case_reads = decompress_case_reads(ch_case_reads)
    ch_decompressed_control_reads = decompress_control_reads(ch_control_reads)
    
    // Assemble
    ch_assembled = assembly(ch_decompressed_case_reads)
    
    // Merge ref
    ch_merged_ref = merge_refTrans_assembly(
        ch_assembled.assembled_fa.combine(Channel.fromPath(params.trans_fasta))
    )

    // Case: Align + Quant
    ch_case_align_quant_result = case_align_quant(
        ch_merged_ref.join(ch_decompressed_case_reads, by: 0)
    )

    // Controls: Align + Quant
    ch_control_align_quant_result = control_align_quant(
        ch_merged_ref.combine(ch_decompressed_control_reads)
    )

    // Group controls by case sample_id
    ch_controls_by_case = ch_control_align_quant_result.control_quant
        .groupTuple(by: 0)  // Group by sample_id (first element)
        .map { sample_id, control_ids, control_quants ->
            tuple(sample_id, control_ids, control_quants)
        }

    // Combine with case data
    ch_matrix_ready = ch_case_align_quant_result.case_quant
        .join(ch_controls_by_case, by: 0)
        .map { sample_id, case_quant, control_ids, control_quants ->
            def all_sample_names = [sample_id] + control_ids
            tuple(sample_id, case_quant, control_quants, all_sample_names)
        }

    // Build transcript matrix
    ch_transcript_matrix_result = build_transcript_matrix(ch_matrix_ready)

    // Group all control outputs by sample_id
    ch_controls_by_case_quant = ch_control_align_quant_result.control_quant
        .groupTuple(by: 0)
        .map { sample_id, control_ids, control_quants ->
            tuple(sample_id, control_quants)
        }

    ch_controls_by_case_meta = ch_control_align_quant_result.control_quant_meta
        .groupTuple(by: 0)
        .map { sample_id, control_ids, control_metas ->
            tuple(sample_id, control_metas)
        }

    ch_controls_by_case_parquet = ch_control_align_quant_result.control_quant_parquet
        .groupTuple(by: 0)
        .map { sample_id, control_ids, control_parquets ->
            tuple(sample_id, control_parquets)
        }

    // Combine all oarfish outputs for DE analysis
    ch_de_input = ch_case_align_quant_result.case_quant
        .join(ch_controls_by_case_quant, by: 0)
        .join(ch_case_align_quant_result.case_quant_meta, by: 0)
        .join(ch_controls_by_case_meta, by: 0)
        .join(ch_case_align_quant_result.case_quant_parquet, by: 0)
        .join(ch_controls_by_case_parquet, by: 0)
        .join(ch_transcript_matrix_result.transcript_matrix, by: 0)
        .combine(Channel.fromPath(params.trans_fasta))
        .map { sample_id, case_quant, control_quants, case_meta, control_metas, case_parquet, control_parquets, transcript_matrix, trans_fasta ->
            tuple(sample_id, case_quant, control_quants, case_meta, control_metas, case_parquet, control_parquets, trans_fasta, transcript_matrix)
        }

    // DE analysis
    ch_de_result = compare_transcript_oarfish(ch_de_input)
    
    // Filter contigs
    ch_filtered = filter_de_contigs(ch_assembled.assembled_fa.join(ch_de_result.de_results, by: 0))

    // Genome alignment
    ch_genome_align = align_contigs_to_genome(
        ch_filtered.de_contigs
            .combine(Channel.fromPath(params.hg38_fasta))
            .combine(Channel.fromPath(params.hg38_splice_junctions))
    )

    // Annotation
    ch_annotated = annotate_contigs(
        ch_genome_align.hg38_with_md_bam.join(ch_genome_align.hg38_with_md_bai, by: 0)
            .combine(Channel.fromPath(params.ann_info))
            .combine(Channel.fromPath(params.tx_annotation))
    )
    
    // Refinement
    ch_refined = refine_annotation(
        ch_annotated.anno_info.join(ch_annotated.anno_vcf, by: 0)
            .join(ch_annotated.anno_bam, by: 0)
            .join(ch_annotated.anno_bai, by: 0)
            .combine(Channel.fromPath(params.tx_annotation))
            .combine(Channel.fromPath(params.hg38_fasta))
    )

    // Filter FASTA
    ch_refined_annotated_contigs_fasta = filter_refined_annotated_contigs_fasta(
        ch_filtered.de_contigs.join(ch_annotated.anno_info, by: 0)
    )

    // VAF
    ch_vaf = estimate_vaf(
        ch_transcript_matrix_result.transcript_matrix
            .join(ch_case_align_quant_result.case_quant, by: 0)
            .join(ch_refined.refined_annotated_contigs_info, by: 0)
            .combine(Channel.fromPath(params.trans_fasta))
            .combine(Channel.fromPath(params.tx2gene))
    )

    // Final
    post_process(
        ch_refined.refined_annotated_contigs_info
            .join(ch_refined_annotated_contigs_fasta.refined_annotated_contigs_fasta, by: 0)
            .join(ch_de_result.de_results, by: 0)
            .join(ch_vaf.vaf_estimates, by: 0)
    )
}