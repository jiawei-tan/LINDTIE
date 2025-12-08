/*
Module      : assembly
Description : Assembles transcripts for each case sample and merges with reference transcriptome.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

/*
  Process: assembly
    - Runs RNA‐Bloom2 to assemble transcripts for each case sample.
    
    Container:
      - bioconda::rnabloom=2.0.1
*/
process assembly {

    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/rnabloom:2.0.1--1a308388e7330445'

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("01-Assembly/rnabloom.transcripts.fa"), emit: assembled_fa

    script:
    """
    ## Create directories
    mkdir -p 01-Assembly
    
    ## Run RNABloom based on platform
    if [ "${params.platform}" == "ONT" ]; then
        rnabloom -long ${reads} -t ${task.cpus} -outdir 01-Assembly
    elif [ "${params.platform}" == "PacBio" ]; then
        rnabloom -lrpb -long ${reads} -t ${task.cpus} -outdir 01-Assembly
    else
        echo "Invalid platform: ${params.platform}"
        exit 1
    fi
    """
}

//------------------------------------------------------------------------------
/*
  Process: merge_refTrans_assembly
    - Concatenates the reference transcriptome with the assembled transcripts.
*/
process merge_refTrans_assembly {
		
    tag "${sample_id}"
    label 'process_short'

    input:
      tuple val(sample_id), path(assembled_fa), path(trans_fasta)

    output:
      tuple val(sample_id), path("merged_refTrans_assembly.fa"), emit: merged_ref

    script:
    """
    cat ${trans_fasta} ${assembled_fa} \\
        > merged_refTrans_assembly.fa
    """
}
