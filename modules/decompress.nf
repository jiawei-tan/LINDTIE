/*
Module      : decompress
Description : Decompresses .gz files for case and control reads.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

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