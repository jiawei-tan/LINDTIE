#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************
                            LINDTIE
 A Nextflow pipeline to identify aberrant transcripts in cancer using 
 long-read RNA-seq data.
 Authors: Jia Wei Tan
 Contact: tan.j@wehi.edu.au
***************************************************************************/

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
--assembly_mode             : Strategy: 'hybrid', 'denovo', 'ref_guided', 'denovo_subset' (default: hybrid)
--subset_count              : Number of reads to subset to (default: NULL)
--minimap2_preset           : Minimap2 preset (passed to -ax). Default: 'map-ont'.
--rnabloom2_preset          : RNABloom2 preset. Leave empty for ONT (default). Use '-lrpb' for PacBio.
--RUN_DE                    : Run differential expression analysis. Options: true (default) or false.
--detect_viral_integration  : Detect viral integration. Options: true or false (default: false).
--fdr                       : False discovery rate (FDR) threshold (default: 0.05)
--min_cpm                   : Minimum counts per million (CPM) (default: 0.5)
--min_logfc                 : Minimum log fold change (default: 2)
--min_clip                  : Minimum clip length (default: 20)
--min_gap                   : Minimum gap (default: 7)
--min_match                 : Minimum match (default: "30,0.3")
--splice_motif_mismatch     : Splice motif mismatch (default: 1)
--oarfish_num_bootstraps    : Number of bootstraps for Oarfish (default: 10)
--oarfish_growth_rate       : Growth rate for Oarfish (passed to --growth-rate) (default: 0.5)
--gene_filter               : List of genes to filter (default: NULL)
--var_filter                : List of variant types to filter (default: NULL)
--single_sample_min_vaf     : Minimum VAF_sum_WT_TPM to keep a variant when RUN_DE is false (default: 0.1)
--single_sample_cosmic_filter: Require COSMIC tier 1/2 to keep a variant when RUN_DE is false (default: true)
--max_fisher_p_val          : Maximum Fisher exact test p-value to keep; variants with larger p are dropped (default: 0.05)
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
include { controls_exist; effective_run_de; effective_max_fisher_p_val } from './modules/decompress'
include { align_raw_reads_to_hg38 } from './modules/assembly'
include { ref_guided_assembly } from './modules/assembly'
include { subset_reads } from './modules/assembly'
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
include { count_supporting_reads } from './modules/annotation'
include { post_process } from './modules/annotation'

/*************************** LOCAL PROCESSES **************************/
process save_params {
    tag "${sample_id}"
    label 'process_short'
    
    publishDir "${sample_id}_output", mode: 'copy'

    input:
    val sample_id

    output:
    path "run_parameters.log"

    script:
    // Report the values the run actually uses, not the requested ones -- they differ when
    // no controls are present (see modules/decompress.nf).
    def eff_run_de = effective_run_de()
    def eff_max_fisher_p_val = effective_max_fisher_p_val()
    def run_de_note = (eff_run_de.toString() == params.RUN_DE.toString()) ? '' : "  (requested ${params.RUN_DE}; no control reads found)"
    def fisher_note = (eff_max_fisher_p_val == params.max_fisher_p_val) ? '' : "  (requested ${params.max_fisher_p_val}; filter disabled without controls)"
    """
    cat <<EOF > run_parameters.log
Sample ID               : ${sample_id}
Timestamp               : ${new Date().format('yyyy-MM-dd HH:mm:ss')}
assembly_mode           : ${params.assembly_mode}
subset_count            : ${params.subset_count}
RUN_DE                  : ${eff_run_de}${run_de_note}
detect_viral_integration: ${params.detect_viral_integration}
minimap2_preset         : ${params.minimap2_preset}
rnabloom2_preset        : ${params.rnabloom2_preset}
fdr                     : ${params.fdr}
min_cpm                 : ${params.min_cpm}
min_logfc               : ${params.min_logfc}
min_clip                : ${params.min_clip}
min_gap                 : ${params.min_gap}
min_match               : '${params.min_match}'
splice_motif_mismatch   : ${params.splice_motif_mismatch}
oarfish_num_bootstraps  : ${params.oarfish_num_bootstraps}
oarfish_growth_rate     : ${params.oarfish_growth_rate}
gene_filter             : ${params.gene_filter}
var_filter              : ${params.var_filter}
single_sample_min_vaf   : ${params.single_sample_min_vaf}
single_sample_cosmic_filter: ${params.single_sample_cosmic_filter}
max_fisher_p_val        : ${eff_max_fisher_p_val}${fisher_note}
EOF
    """
}

/*************************** WORKFLOW **************************/
workflow {

    // -------------------------------------------------------
    // CONFIGURATION CHECK
    // -------------------------------------------------------
    // Effective run mode, derived from the control reads actually present.
    // NOTE: this is derived rather than assigned back onto `params` -- Nextflow silently
    // ignores reassignment of `params` inside a workflow, which is why the old
    // `params.RUN_DE = false` fallback here never took effect. See modules/decompress.nf.
    def run_de = effective_run_de()

    // If no controls found, fall back to single sample mode and warn the user
    if (!controls_exist() && params.RUN_DE) {
        log.warn "================================================================================"
        log.warn "  WARNING: No control samples found in '${params.controls_fastq_dir}'."
        log.warn "  'RUN_DE' is treated as 'false' for this run."
        log.warn "  Pipeline will run in SINGLE SAMPLE MODE (Novel Contig Detection only)."
        log.warn "================================================================================"
    }

    // Fisher's exact test compares case vs control counts. Control counts are only
    // produced when DE runs, so without them the contingency table is degenerate and
    // every variant gets p = 1 -- any max_fisher_p_val < 1 would discard everything.
    if (!run_de && params.max_fisher_p_val < 1) {
        log.warn "================================================================================"
        log.warn "  WARNING: Running without controls (single sample mode)."
        log.warn "  Fisher's exact test is not meaningful without controls (all p-values = 1)."
        log.warn "  Treating 'max_fisher_p_val' ${params.max_fisher_p_val} as 1 (filter disabled)."
        log.warn "================================================================================"
    }


    // -------------------------------------------------------
    // INPUT CHANNELS
    // -------------------------------------------------------
    // Define input channels
    ch_case_reads = Channel
        .fromPath("${params.cases_fastq_dir}/*.{fasta,fasta.gz,fa,fa.gz,fastq,fastq.gz,fq,fq.gz}")
        .map { file -> 
            def name      = file.getName()
            def sample_id = name.replaceFirst(/\.(fastq|fasta|fq|fa)(\.gz)?$/, '')
            tuple(sample_id, file)
        }

    // Extract sample IDs and trigger the save_params process
    save_params(ch_case_reads.map { sid, file -> sid })

    // Process case samples
    ch_decompressed_case_reads = decompress_case_reads(ch_case_reads)

    // =========================================================================
    // OPTIONAL CONTROL PROCESSING
    // =========================================================================
    
    // Initialize empty channels for controls
    ch_controls_by_case         = Channel.empty()
    ch_controls_by_case_quant   = Channel.empty()
    ch_controls_by_case_meta    = Channel.empty()
    ch_controls_by_case_parquet = Channel.empty()

    if (run_de) {
        ch_control_reads = Channel
            .fromPath("${params.controls_fastq_dir}/*.{fasta,fasta.gz,fa,fa.gz,fastq,fastq.gz,fq,fq.gz}")
            .map { file -> 
                def name       = file.getName()
                def control_id = name.replaceFirst(/\.(fastq|fasta|fq|fa)(\.gz)?$/, '')
                tuple(control_id, file)
            }

        ch_decompressed_control_reads = decompress_control_reads(ch_control_reads)
    }

// =========================================================================
    // ASSEMBLY STRATEGY
    // =========================================================================
    
    // Initialize Channels
    ch_stringtie_assembly = Channel.empty()
    ch_rnabloom_assembly  = Channel.empty()
    
    ch_reads_for_alignment = Channel.empty()
    ch_reads_for_denovo    = Channel.empty()

    // -------------------------------------------------------
    // 1. Route Reads Based on Mode
    // -------------------------------------------------------
    
    if (params.assembly_mode == 'denovo_subset') {
        // --- SPLIT MODE (denovo_subset) ---
        // 1. Split the reads using BBMap
        ch_split_result = subset_reads(ch_decompressed_case_reads, params.subset_count)
        
        // 2. Subset goes DIRECTLY to De Novo
        ch_reads_for_denovo = ch_split_result.subset_fastq
        
        // 3. Remainder goes to Alignment
        ch_reads_for_alignment = ch_split_result.remainder_fastq

    } else if (params.assembly_mode == 'hybrid' || params.assembly_mode == 'ref_guided') {
        // --- STANDARD ALIGNMENT MODES ---
        // Full reads go to alignment
        ch_reads_for_alignment = ch_decompressed_case_reads
    
    } else if (params.assembly_mode == 'denovo') {
        // --- PURE DE NOVO ---
        // Full reads go to de novo
        ch_reads_for_denovo = ch_decompressed_case_reads
    }

    // -------------------------------------------------------
    // 2. Run Alignment (If needed)
    // -------------------------------------------------------
    
    ch_aligned_bam_for_stringtie = Channel.empty()
    
    // Logic: Alignment runs for everything EXCEPT pure denovo
    if (params.assembly_mode != 'denovo') { 
        
        ch_aligned = align_raw_reads_to_hg38(
            ch_reads_for_alignment, 
            Channel.fromPath(params.hg38_fasta),
            Channel.fromPath(params.hg38_splice_junctions)
        )

        // Determine which BAM goes to reference-guided assembly
        if (params.assembly_mode == 'ref_guided' || params.assembly_mode == 'denovo_subset') {
            ch_aligned_bam_for_stringtie = ch_aligned.all_mapped_bam
        } else {
            // hybrid
            ch_aligned_bam_for_stringtie = ch_aligned.confident_mapped_bam
        }

        // IF HYBRID: Rescued reads also go to De Novo
        if (params.assembly_mode == 'hybrid') {
             ch_reads_for_denovo = ch_aligned.rescued_fastq
        }
    }

    // -------------------------------------------------------
    // 3. Run Assemblies
    // -------------------------------------------------------

    // A. Reference Guided Assembly (StringTie2)
    if (params.assembly_mode == 'hybrid' || params.assembly_mode == 'ref_guided' || params.assembly_mode == 'denovo_subset') {
        
        ch_ref_guided_assembled = ref_guided_assembly(
            ch_aligned_bam_for_stringtie,
            Channel.fromPath(params.tx_annotation),
            Channel.fromPath(params.hg38_fasta)
        )
        ch_stringtie_assembly = ch_ref_guided_assembled.stringtie2_assembled_fa
    }

    // B. De Novo Assembly (RNABloom2)
    if (params.assembly_mode == 'hybrid' || params.assembly_mode == 'denovo' || params.assembly_mode == 'denovo_subset') {
        
        ch_de_novo_assembled = de_novo_assembly(ch_reads_for_denovo)
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
        ch_merged_ref.merged_ref.join(ch_decompressed_case_reads, by: 0)
    )

    // Controls: Align + Quant (ONLY IF RUN_DE = TRUE)
    if (run_de) {
        ch_control_align_quant_result = control_align_quant(
            ch_merged_ref.merged_ref.combine(ch_decompressed_control_reads)
        )

        // 1. Prepare controls for Matrix Build (Tuple: sample_id, [ctrl_ids], [ctrl_quants])
        ch_controls_by_case = ch_control_align_quant_result.control_quant
            .groupTuple(by: 0)
            .map { sample_id, control_ids, control_quants ->
                tuple(sample_id, control_ids, control_quants)
            }

        // 2. Prepare controls for DE (Individual channels grouped by sample_id)
        ch_controls_by_case_quant = ch_control_align_quant_result.control_quant
            .groupTuple(by: 0)
            .map { sid, c_ids, c_quants -> tuple(sid, c_quants) }

        ch_controls_by_case_meta = ch_control_align_quant_result.control_quant_meta
            .groupTuple(by: 0)
            .map { sid, c_ids, c_metas -> tuple(sid, c_metas) }

        ch_controls_by_case_parquet = ch_control_align_quant_result.control_quant_parquet
            .groupTuple(by: 0)
            .map { sid, c_ids, c_parquets -> tuple(sid, c_parquets) }
    }

    // =========================================================================
    // CONDITIONAL JOIN LOGIC
    // =========================================================================

    // A. Prepare input for Transcript Matrix
    if (run_de) {
        // Standard join waiting for controls
        ch_matrix_ready = ch_case_align_quant_result.case_quant
            .join(ch_controls_by_case, by: 0)
            .map { sample_id, case_quant, control_ids, control_quants ->
                def all_sample_names = [sample_id] + control_ids
                tuple(sample_id, case_quant, control_quants, all_sample_names)
            }
    } else {
        // Bypass: Just the case, empty lists for controls
        ch_matrix_ready = ch_case_align_quant_result.case_quant
            .map { sample_id, case_quant ->
                // control_quants = [], all_sample_names = [sample_id]
                tuple(sample_id, case_quant, [], [sample_id])
            }
    }

    // Build transcript matrix
    ch_transcript_matrix_result = build_transcript_matrix(ch_matrix_ready)

    // B. Prepare input for DE / Novel Contig Detection
    if (run_de) {
        ch_de_input = ch_case_align_quant_result.case_quant
            .join(ch_controls_by_case_quant, by: 0)
            .join(ch_case_align_quant_result.case_quant_meta, by: 0)
            .join(ch_controls_by_case_meta, by: 0)
            .join(ch_case_align_quant_result.case_quant_parquet, by: 0)
            .join(ch_controls_by_case_parquet, by: 0)
            .join(ch_transcript_matrix_result.transcript_matrix, by: 0)
            // JOIN THE NEW NOVEL_ONLY CHANNEL
            .join(ch_merged_ref.novel_only, by: 0) 
            .combine(Channel.fromPath(params.trans_fasta))
            // Update map to include 'novel_only'
            .map { sample_id, case_quant, c_quants, case_meta, c_metas, case_pq, c_pqs, tx_mat, novel_only, trans_fa ->
                tuple(sample_id, case_quant, c_quants, case_meta, c_metas, case_pq, c_pqs, trans_fa, tx_mat, novel_only)
            }
    } else {
        // Bypass Logic for No Controls
        ch_de_input = ch_case_align_quant_result.case_quant
            .join(ch_case_align_quant_result.case_quant_meta, by: 0)
            .join(ch_case_align_quant_result.case_quant_parquet, by: 0)
            .join(ch_transcript_matrix_result.transcript_matrix, by: 0)
            // JOIN THE NEW NOVEL_ONLY CHANNEL
            .join(ch_merged_ref.novel_only, by: 0) 
            .combine(Channel.fromPath(params.trans_fasta))
            // Update map to include 'novel_only' and pass empty lists for controls
            .map { sample_id, case_quant, case_meta, case_pq, tx_mat, novel_only, trans_fa ->
                tuple(sample_id, case_quant, [], case_meta, [], case_pq, [], trans_fa, tx_mat, novel_only)
            }
    }

    // Run DE (or novel contig detection if !RUN_DE)
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

    // Sequence-verified supporting/spanning read counts per variant. Needs the case
    // genome alignment, so it is skipped in pure denovo mode and an empty table is
    // substituted -- post_process then leaves the columns blank rather than failing.
    if (params.assembly_mode != 'denovo') {
        ch_supporting = count_supporting_reads(
            ch_refined.refined_annotated_contigs_info
                .join(ch_aligned.all_mapped_bam, by: 0)
                .join(ch_aligned.all_mapped_bam_bai, by: 0)
                .combine(Channel.fromPath(params.tx_annotation))
        ).supporting_reads
    } else {
        // Synthesize the header-only table rather than shipping an asset file, so the
        // columns line up with LINDTIE_supporting_reads.py's OUT_COLUMNS.
        def empty_supp = Channel.of(
            'variant_id\tcontig_id\tvariant_type\tsupporting_read_count\tspanning_reads\t' +
            'junction_VAF\tcandidate_depth\tcounted_locus\tsupport_signature\t' +
            'supporting_read_count_reliability'
        ).collectFile(name: 'empty_supporting_reads.tsv', newLine: true)

        ch_supporting = ch_refined.refined_annotated_contigs_info
            .map { sid, info -> sid }
            .combine(empty_supp)
    }

    // Final
    post_process(
        ch_refined.refined_annotated_contigs_info
            .join(ch_refined_annotated_contigs_fasta.refined_annotated_contigs_fasta, by: 0)
            .join(ch_de_result.de_results, by: 0)
            .join(ch_vaf.vaf_estimates, by: 0)
            .join(ch_supporting, by: 0)
            .combine(Channel.fromPath(params.cosmic_tier_data))
    )
}