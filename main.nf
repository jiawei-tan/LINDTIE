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

// Defaults
params.assembly_mode = 'hybrid'     // 'hybrid', 'denovo', or 'ref_guided'
params.minimap2_preset = 'map-ont'  // 'map-ont', 'map-pb', 'map-hifi', 'lr:hq'
params.rnabloom2_preset = ''        // '', '-lrpb'
params.RUN_DE = true                // true or false
params.fdr = 0.05
params.min_cpm = 0.5
params.min_logfc = 2
params.min_clip = 20
params.min_gap = 7
params.min_match = "30,0.3"
params.splice_motif_mismatch = 1
params.oarfish_num_bootstraps = 10
params.gene_filter = null
params.var_filter = null
params.help = false
params.version = false

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
    --minimap2_preset "map-pb" \\
    --rnabloom2_preset "lrpb" \\
    --assembly_mode "ref_guided"

Optional parameters:
--assembly_mode             : Strategy: 'hybrid', 'denovo', 'ref_guided' (default: hybrid)
--minimap2_preset           : Minimap2 preset (passed to -ax). Default: 'map-ont'.
                              Options:
                                - 'map-ont' : Oxford Nanopore genomic reads (default)
                                - 'map-pb'  : PacBio CLR genomic reads
                                - 'map-hifi': PacBio HiFi/CCS genomic reads
                                - 'lr:hq'   : Nanopore Q20 genomic reads
--rnabloom2_preset          : RNABloom2 preset. Leave empty for ONT (default). Use '-lrpb' for PacBio.
--RUN_DE                    : Run differential expression analysis. Options: true (default) or false.
--fdr                       : False discovery rate (FDR) threshold (default: 0.05)
--min_cpm                   : Minimum counts per million (CPM) (default: 0.1)
--min_logfc                 : Minimum log fold change (default: 2)
--min_clip                  : Minimum clip length (default: 20)
--min_gap                   : Minimum gap (default: 7)
--min_match                 : Minimum match (default: "30,0.3")
--splice_motif_mismatch     : Splice motif mismatch (default: 1)
--oarfish_num_bootstraps    : Number of bootstraps for Oarfish (default: 10)
--gene_filter               : List of genes to filter (default: NULL)
--var_filter                : List of variant types to filter (default: NULL)
--help                      : Show this help message
--version                   : Show the version of LINDTIE
    """
    exit 0
}

if (params.version) {
    log.info "LINDTIE v0.1.0"
    exit 0
}

// Print parameters to Console for immediate verification
log.info """
================================================================================
                            LINDTIE PARAMETER LOG
================================================================================
"""
params.sort().each { k, v ->
    log.info "\$k".padRight(30) + ": \$v"
}
log.info "================================================================================"


/*************************** MODULES **************************/

include { decompress_case_reads } from './modules/decompress'
include { decompress_control_reads } from './modules/decompress'
include { align_raw_reads_to_hg38 } from './modules/assembly'
include { ref_guided_assembly } from './modules/assembly'
include { de_novo_assembly } from './modules/assembly'
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

/*************************** LOCAL PROCESSES **************************/
/*
  Process: save_params
  Description: Writes the run parameters to a log file.
*/
process save_params {
    tag "${sample_id}"
    label 'process_short'
    
    publishDir "${sample_id}_output", mode: 'copy'

    input:
    val sample_id

    output:
    path "run_parameters.log"

    script:
    """
    cat <<EOF > run_parameters.log
Sample ID             : ${sample_id}

# Parameters for the workflow
assembly_mode         : ${params.assembly_mode}
RUN_DE                : ${params.RUN_DE}

# Tool Presets
minimap2_preset       : ${params.minimap2_preset}
rnabloom2_preset      : ${params.rnabloom2_preset}

fdr                   : ${params.fdr}
min_cpm               : ${params.min_cpm}
min_logfc             : ${params.min_logfc}
min_clip              : ${params.min_clip}
min_gap               : ${params.min_gap}
min_match             : '${params.min_match}'
splice_motif_mismatch : ${params.splice_motif_mismatch}
oarfish_num_bootstraps: ${params.oarfish_num_bootstraps}
gene_filter           : ${params.gene_filter}
var_filter            : ${params.var_filter}
EOF
    """
}

/*************************** WORKFLOW **************************/
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

    // Extract sample IDs and trigger the save_params process
    save_params(ch_case_reads.map { sid, file -> sid })

    // Process each case sample
    // Decompress reads if needed
    ch_decompressed_case_reads = decompress_case_reads(ch_case_reads)
    ch_decompressed_control_reads = decompress_control_reads(ch_control_reads)

    // =========================================================================
    // ASSEMBLY STRATEGY
    // =========================================================================
    
    // Initialize assembly channels as empty
    ch_stringtie_assembly = Channel.empty()
    ch_rnabloom_assembly  = Channel.empty()
    
    // --- Strategy 1: Genome Alignment Required (Hybrid OR Ref-Guided) ---
    if (params.assembly_mode == 'hybrid' || params.assembly_mode == 'ref_guided') {
        
        // 1. Align to Genome
        ch_aligned = align_raw_reads_to_hg38(
            ch_decompressed_case_reads, 
            Channel.fromPath(params.hg38_fasta)
        )
        
        // 2. Select BAM for StringTie2
        // IF Hybrid: use "confident_mapped_bam" (filtered)
        // IF Ref-Guided: use "all_mapped_bam" (standard) OR "confident" depending on preference
        def ch_bam_for_stringtie = (params.assembly_mode == 'ref_guided') 
                                    ? ch_aligned.all_mapped_bam 
                                    : ch_aligned.confident_mapped_bam 

        // 3. Run Ref-Guided Assembly
        ch_ref_guided_assembled = ref_guided_assembly(
            ch_bam_for_stringtie,
            Channel.fromPath(params.gencode_annotation),
            Channel.fromPath(params.hg38_fasta)
        )
        
        ch_stringtie_assembly = ch_ref_guided_assembled.stringtie2_assembled_fa
    }

    // --- Strategy 2: De Novo Assembly (Hybrid OR De Novo) ---
    if (params.assembly_mode == 'hybrid' || params.assembly_mode == 'denovo') {
        
        // Determine input reads
        def ch_reads_for_bloom = (params.assembly_mode == 'hybrid')
                                 ? ch_aligned.rescued_fastq
                                 : ch_decompressed_case_reads

        // Run De Novo Assembly
        ch_de_novo_assembled = de_novo_assembly(ch_reads_for_bloom)
        
        ch_rnabloom_assembly = ch_de_novo_assembled.rnabloom_assembled_fa
    }

    // --- Merge Assemblies ---
    ch_assemblies_to_merge = ch_stringtie_assembly.mix(ch_rnabloom_assembly)
        .groupTuple() 

    ch_merged_ref = merge_refTrans_assembly(
        ch_assemblies_to_merge,
        Channel.fromPath(params.trans_fasta)
    )

    // =========================================================================
    // QUANTIFICATION & DOWNSTREAM ANALYSIS
    // =========================================================================

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
        .groupTuple(by: 0)  
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
        .map { sample_id, case_quant, control_quants, case_meta, control_metas, case_parquet, control_parquets, trans_fasta, transcript_matrix ->
            tuple(sample_id, case_quant, control_quants, case_meta, control_metas, case_parquet, control_parquets, trans_fasta, transcript_matrix)
        }

    // DE analysis
    ch_de_result = compare_transcript_oarfish(ch_de_input)
     
    // Filter contigs
    ch_filtered = filter_de_contigs(ch_merged_ref.merged_ref.join(ch_de_result.de_results, by: 0))

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