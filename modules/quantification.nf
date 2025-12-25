/*
Module      : quantification
Description : Aligns reads (case and control) against the merged reference and quantifies transcript expression.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

/*
  Process: case_align_quant
    - Aligns case reads against the merged reference (minimap2; platform-specific),
      then quantifies transcript expression with Oarfish.
      
    Container:
      - bioconda::minimap2=2.30
      - bioconda::samtools=1.22
      - bioconda::qualimap=2.3
      - bioconda::oarfish=0.8.1
*/
process case_align_quant {

    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output/02-Quantification/cases", mode: 'copy', pattern: "${sample_id}.{quant,infreps.pq,meta_info.json}"
    
    container 'oras://community.wave.seqera.io/library/minimap2_oarfish_qualimap_samtools:90c13d0a92fc925a'
    
    input:
      tuple val(sample_id), path(merged_ref), path(reads)

    output:
      tuple val(sample_id), path("minimap_against_merged_refTrans_assembly_unsorted.bam"), emit: case_unsorted_bam
      tuple val(sample_id), path("oarfish_output"), emit: case_oarfish_output
      tuple val(sample_id), path("${sample_id}.quant"), emit: case_quant
      tuple val(sample_id), path("${sample_id}.infreps.pq"), emit: case_quant_parquet
      tuple val(sample_id), path("${sample_id}.meta_info.json"), emit: case_quant_meta
      
    script:
    """
    set -e

    ## Create directories
    mkdir -p oarfish_output

    ## 1) Align with minimap2 using parameterized preset
    minimap2 -t ${task.cpus} -ax ${params.minimap2_preset} --eqx -N 100 \\
        ${merged_ref} ${reads} | \\
        samtools view -@ ${task.cpus} -bS - \\
          > minimap_against_merged_refTrans_assembly_unsorted.bam

    ## 2) Oarfish quant
    oarfish \\
      -j ${task.cpus} \\
      -a minimap_against_merged_refTrans_assembly_unsorted.bam \\
      -o oarfish_output/${sample_id} \\
      --num-bootstraps ${params.oarfish_num_bootstraps} \\
      --filter-group no-filters \\
      --model-coverage

    ## 3) Prepare declared outputs at process root so Nextflow tracks them
    cp -f oarfish_output/${sample_id}.quant ${sample_id}.quant
    cp -f oarfish_output/${sample_id}.infreps.pq ${sample_id}.infreps.pq
    cp -f oarfish_output/${sample_id}.meta_info.json ${sample_id}.meta_info.json
    """
}

//------------------------------------------------------------------------------
/*
  Process: control_align_quant
    - Aligns each control against the same merged reference for a given case (minimap2; platform-specific),
      then quantifies transcript expression with Oarfish.

      Container:
    - bioconda::minimap2=2.30
    - bioconda::samtools=1.22
    - bioconda::oarfish=0.8.1
*/
process control_align_quant {

    tag "${sample_id}/${control_id}"
    label 'process_long'

    publishDir "${sample_id}_output/02-Quantification/controls", mode: 'copy', pattern: "${control_id}.{quant,infreps.pq,meta_info.json}"

    container 'oras://community.wave.seqera.io/library/minimap2_oarfish_samtools:7a3485aa349f407d'

    input:
      tuple val(sample_id), path(merged_ref), val(control_id), path(control_fastq)

    output:
      tuple val(sample_id), val(control_id), path("${control_id}.quant"), emit: control_quant
      tuple val(sample_id), val(control_id), path("${control_id}.infreps.pq"), emit: control_quant_parquet
      tuple val(sample_id), val(control_id), path("${control_id}.meta_info.json"), emit: control_quant_meta
      tuple val(sample_id), val(control_id), path("${control_id}_minimap_unsorted.bam"), emit: control_bam
    script:
    """
    echo "Processing control ${control_id} for sample ${sample_id}"

    ## 1) Align control FASTQ to merged ref
    minimap2 -t ${task.cpus} -ax ${params.minimap2_preset} --eqx -N 100 \\
        ${merged_ref} ${control_fastq} | \\
        samtools view -@ ${task.cpus} -bS - \\
          > ${control_id}_minimap_unsorted.bam

    ## 2) Oarfish quant
    oarfish \\
      -j ${task.cpus} \\
      -a ${control_id}_minimap_unsorted.bam \\
      -o ${control_id} \\
      --num-bootstraps ${params.oarfish_num_bootstraps} \\
      --filter-group no-filters \\
      --model-coverage

    echo "Generated files:"
    ls -la ${control_id}.*
    """
}
