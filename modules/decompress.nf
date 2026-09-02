/*
Module      : decompress
Description : Decompresses .gz files for case and control reads.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

/*
  Shared configuration helpers.

  These derive the pipeline's *effective* run mode from the inputs that are actually
  present. They exist because Nextflow silently ignores reassignment of `params` inside
  a workflow (it warns "defined multiple times"), so the old `params.RUN_DE = false`
  fallback in main.nf never took effect. Deriving the value instead of mutating it works
  both in the workflow body and inside process script blocks.
*/

// The files the control channel would actually read (same glob as main.nf).
def control_read_files() {
    params.controls_fastq_dir
        ? file("${params.controls_fastq_dir}/*.{fasta,fasta.gz,fa,fa.gz,fastq,fastq.gz,fq,fq.gz}")
        : []
}

// A controls directory holding only, say, a stray README does not count as having controls.
def controls_exist() {
    control_read_files().size() > 0
}

// DE requires controls. Coerce explicitly: the string "false" is truthy in Groovy.
def effective_run_de() {
    (params.RUN_DE.toString().toLowerCase() == 'true') && controls_exist()
}

// Control counts are only produced when DE runs. Without them Fisher's exact test is
// degenerate and returns p = 1 for every variant, so any threshold < 1 would discard the
// entire result set -- disable the filter rather than filter everything away.
def effective_max_fisher_p_val() {
    effective_run_de() ? params.max_fisher_p_val : 1
}

/*
  Process: decompress_case_reads
    - Decompresses .gz files for case reads if they are compressed
*/
process decompress_case_reads {

    tag "${sample_id}"
    label 'process_short'

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("decompressed_*"), emit: decompressed_reads

    script:
    def filename = reads.toString()
    def is_gzipped = filename.endsWith('.gz')
    def base_name = filename.replaceAll(/\.gz$/, '')
    def extension = base_name.replaceAll(/.*\./, '')
    def output_name = "decompressed_${sample_id}.${extension}"
    
    if (is_gzipped) {
        """
        echo "Decompressing ${reads} to ${output_name}"
        gunzip -c ${reads} > ${output_name}
        """
    } else {
        """
        echo "File ${reads} is not compressed, copying to ${output_name}"
        cp ${reads} ${output_name}
        """
    }
}

//------------------------------------------------------------------------------
/*
  Process: decompress_control_reads
    - Decompresses .gz files for control reads if they are compressed
*/
process decompress_control_reads {

    tag "${control_id}"
    label 'process_short'

    input:
      tuple val(control_id), path(reads)

    output:
      tuple val(control_id), path("decompressed_*"), emit: decompressed_reads

    script:
    def filename = reads.toString()
    def is_gzipped = filename.endsWith('.gz')
    def base_name = filename.replaceAll(/\.gz$/, '')
    def extension = base_name.replaceAll(/.*\./, '')
    def output_name = "decompressed_${control_id}.${extension}"
    
    if (is_gzipped) {
        """
        echo "Decompressing ${reads} to ${output_name}"
        gunzip -c ${reads} > ${output_name}
        """
    } else {
        """
        echo "File ${reads} is not compressed, copying to ${output_name}"
        cp ${reads} ${output_name}
        """
    }
}