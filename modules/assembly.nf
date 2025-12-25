/*
Module      : assembly
Description : Assembles transcripts for each case sample and merges with reference transcriptome.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
*/

//------------------------------------------------------------------------------
/*
  Process: align_raw_reads_to_hg38
    - Aligns the raw reads to the reference genome and separates the mapped and unmapped reads.
Container:
      - bioconda::minimap2=2.30
      - bioconda::samtools=1.22
*/
process align_raw_reads_to_hg38 {
    
    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output/01-Assembly/", mode: 'copy', pattern: '*.log'

    container 'oras://community.wave.seqera.io/library/minimap2_samtools:e98addfcfd60e8e7'

    input:
      tuple val(sample_id), path(reads)
      path hg38_fasta

    output:
      tuple val(sample_id), path("confident_mapped.bam"), emit: confident_mapped_bam
      tuple val(sample_id), path("reads_all_sorted.bam"), emit: all_mapped_bam
      tuple val(sample_id), path("rescued.fastq"), emit: rescued_fastq
      tuple val(sample_id), path("read_counts_summary.log"), emit: read_counts_summary, optional: true

    script:
    """
    # 1. Align reads to reference genome
    minimap2 -t ${task.cpus} -ax splice --secondary=no \\
      ${hg38_fasta} ${reads} | \\
      samtools sort -@ ${task.cpus} -o reads_all_sorted.bam

    samtools index reads_all_sorted.bam

    # 2. Single-Pass Routing
    # Pass 'mode' variable to AWK
    samtools view -h reads_all_sorted.bam |
    awk -v mode="${params.assembly_mode}" '
      BEGIN {
        n_total=0; n_confident=0; n_rescued=0;
        r_unmapped=0; r_supp=0; r_low_mapq=0; r_indel=0; r_clip=0; r_gap=0;
      }
      
      /^@/ { 
        print > "confident.sam"
        print > "rescued.sam"
        next 
      }
      
      {
        n_total++
        flag = \$2
        mapq = \$5
        cigar = \$6
        
        is_unmapped = int(flag / 4) % 2
        is_supp = int(flag / 2048) % 2
        is_low_mapq = (mapq < 20)
        
        has_indel = 0
        has_large_gap = 0
        is_clipped = 0
        
        if (!is_unmapped) {
            tmp = cigar
            clipped_bases = 0
            while (match(tmp, /[0-9]+[IDNSM]/)) {
                op = substr(tmp, RLENGTH, 1)
                len = substr(tmp, RSTART, RLENGTH - 1) + 0
                if ((op == "I" || op == "D") && len > 15) has_indel = 1
                if (op == "N" && len > 1000) has_large_gap = 1
                if (op == "S") clipped_bases += len
                tmp = substr(tmp, RSTART + RLENGTH)
            }
            if (clipped_bases > 50) is_clipped = 1
        }

        if (is_unmapped || is_supp || is_low_mapq || has_indel || has_large_gap || is_clipped) {
            print > "rescued.sam"
            n_rescued++
            if (is_unmapped) r_unmapped++
            else if (is_supp) r_supp++
            else if (is_low_mapq) r_low_mapq++
            else if (has_indel) r_indel++
            else if (has_large_gap) r_gap++
            else if (is_clipped) r_clip++
        } else {
            print > "confident.sam"
            n_confident++
        }
      }
      
      END {
        # Only write the log file if mode is hybrid
        if (mode == "hybrid") {
            print "Total Reads Processed: " n_total > "read_counts_summary.log"
            print "--------------------------------" >> "read_counts_summary.log"
            print "Sent to StringTie2 (Confident): " n_confident " (" (n_confident/n_total)*100 "%)" >> "read_counts_summary.log"
            print "Sent to RNA-Bloom2 (Rescued):   " n_rescued " (" (n_rescued/n_total)*100 "%)" >> "read_counts_summary.log"
            print "--------------------------------" >> "read_counts_summary.log"
            print "Rescue Reasons (Hierarchical):" >> "read_counts_summary.log"
            print "  Unmapped:        " r_unmapped >> "read_counts_summary.log"
            print "  Supplementary:   " r_supp >> "read_counts_summary.log"
            print "  Low MAPQ (<20):  " r_low_mapq >> "read_counts_summary.log"
            print "  Large Indel:     " r_indel >> "read_counts_summary.log"
            print "  Large Gap:       " r_gap >> "read_counts_summary.log"
            print "  Soft Clipped:    " r_clip >> "read_counts_summary.log"
        }
      }
    '

    # 3. Convert outputs
    samtools view -b -@ ${task.cpus} confident.sam > confident_mapped.bam
    samtools index confident_mapped.bam
    rm confident.sam

    samtools view -b -@ ${task.cpus} rescued.sam > rescued.bam
    samtools fastq -@ ${task.cpus} -n rescued.bam > rescued.fastq
    rm rescued.sam rescued.bam
    """
}

//------------------------------------------------------------------------------
/*
  Process: ref_guided_assembly
    - Use aligned reads (either all or confident) to perform reference-guided assembly using StringTie2.
Container:
        - bioconda::gffread=0.12.7
        - bioconda::stringtie=2.2.3
*/
process ref_guided_assembly {
    
    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output/01-Assembly/", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/gffread_stringtie:12cc4a646c48604f'

    input:
      tuple val(sample_id), path(mapped_bam)
      path gencode_annotation
      path hg38_fasta

    output:
      tuple val(sample_id), path("stringtie2_novel_transcripts.gtf"), emit: stringtie2_gtf
      tuple val(sample_id), path("stringtie2_assembly.fa"), emit: stringtie2_assembled_fa

    script:
    """
    # Assemble transcripts
    stringtie ${mapped_bam} -G ${gencode_annotation} -p ${task.cpus} -L -o stringtie2.gtf

    # Extract novel transcripts & all exons belonging to those novel transcripts
    awk '\$3=="transcript" && \$0 !~ /reference_id/ {
            print
            tid=""
            for(i=1;i<=NF;i++){
                if(\$i=="transcript_id"){
                    tid=\$(i+1)
                    gsub(/"/,"",tid)
                    gsub(/;/,"",tid)
                    keep[tid]=1
                }
            }
        }
        \$3=="exon" {
            tid=""
            for(i=1;i<=NF;i++){
                if(\$i=="transcript_id"){
                    tid=\$(i+1)
                    gsub(/"/,"",tid)
                    gsub(/;/,"",tid)
                    if(keep[tid]==1) print
                }
            }
        }' stringtie2.gtf > stringtie2_novel_transcripts.gtf

    # Convert GTF to FASTA
    gffread stringtie2_novel_transcripts.gtf -g ${hg38_fasta} -w stringtie2_assembly.fa
    """
}

/*
  Process: de_novo_assembly
    - Use reads (raw or rescued) to perform de novo assembly using RNA-Bloom2.
Container:
      - bioconda::rnabloom=2.0.1
*/
process de_novo_assembly {

    tag "${sample_id}"
    label 'process_long'

    publishDir "${sample_id}_output", mode: 'copy'

    container 'oras://community.wave.seqera.io/library/rnabloom:2.0.1--1a308388e7330445'

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("01-Assembly/rnabloom.transcripts.fa"), emit: rnabloom_assembled_fa

    script:
    // Automatically removes existing dash (if any) and adds a fresh one
    // This makes inputs 'lrpb' and '-lrpb' both valid
    def bloom_mode = "-" + params.rnabloom2_preset.replaceFirst(/^-/, '')
    """
    ## Create directories
    mkdir -p 01-Assembly

    ## Run RNABloom2 - params.rnabloom2_preset is optional only for '-lrpb' for PacBio data.
    rnabloom ${bloom_mode} -long ${reads} -t ${task.cpus} -outdir 01-Assembly
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
      tuple val(sample_id), path(assemblies)
      path trans_fasta

    output:
      tuple val(sample_id), path("merged_refTrans_assembly.fa"), emit: merged_ref

    script:
    """
    cat ${trans_fasta} ${assemblies} > merged_refTrans_assembly.fa
    """
}