'''
Module      : LINDTIE_contigs_annotation
Description : Performs contig annotation.
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
Portability : POSIX
Processes aligned contigs from a BAM file, filters out those matching 
known reference transcript annotations, and annotates the remaining 
contigs for novel structural, splicing, and fusion events.

Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
'''

import pandas as pd
import numpy as np
import pysam
import os
import re
import logging
import sys
import string
import block_helper as bh
import constants
from argparse import ArgumentParser
from intervaltree import IntervalTree
from cv_vcf import CrypticVariant, VCF
from utils import cached, init_logging

PROGRAM_NAME = "LINDTIE_contigs_annotation"
ASCII_LETTERS = np.array(list(string.ascii_letters))

# Global record to track processed query names
record = {}

# Pre-compile regex patterns for GTF parsing
RE_GENE_ID = re.compile(r'gene_id "([\w\-\.\/]+)"')
RE_GENE_NAME = re.compile(r'gene_name "([\w\-\.\/]+)"')

def parse_args():
    description = 'Annotate contigs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log', metavar='LOG_FILE', type=str, help='record program progress in LOG_FILE')
    parser.add_argument(dest='sample', metavar='SAMPLE', type=str, help='Sample name')
    parser.add_argument(dest='bam_file', metavar='BAM_FILE', type=str, help='SAM/BAM file containing contig alignments')
    parser.add_argument(dest='junc_file', metavar='JUNC_FILE', type=str, help='Reference file containing splice junctions')
    parser.add_argument(dest='tx_ref_file', metavar='TX_REF_FILE', type=str, help='Transcriptome GTF reference file')
    parser.add_argument(dest='output_bam', metavar='OUTPUT_BAM', type=str, help='BAM file to write contigs which pass filtering')
    parser.add_argument(dest='contig_info_output', metavar='CONTIG_INFO_OUT', type=str, help='Contig info output file')
    parser.add_argument('--minClip', metavar='MIN_CLIP', type=int, help='Minimum novel block or softclip size')
    parser.add_argument('--minGap', metavar='MIN_GAP', type=int, help='Minimum gap (deletion or insertion) size')
    parser.add_argument('--minMatch', metavar='MIN_MATCH', type=str, help='Comma separated: <min bp>, <min perc>')
    return parser.parse_args()

def set_globals(args):
    global MIN_CLIP, MIN_GAP, MIN_MATCH_BP, MIN_MATCH_PERC

    MIN_CLIP = args.minClip if args.minClip else constants.DEFAULT_MIN_CLIP
    MIN_GAP = args.minGap if args.minGap else constants.DEFAULT_MIN_GAP

    if args.minMatch:
        try:
            min_bp, min_perc = args.minMatch.split(',')
            MIN_MATCH_BP = int(min_bp)
            MIN_MATCH_PERC = float(min_perc)
            if MIN_MATCH_BP < 0:
                raise ValueError("minMatch BP must be >= 0")
            if not (0 <= MIN_MATCH_PERC <= 1):
                raise ValueError("minMatch percentage must be between 0 and 1")
        except ValueError as e:
            logging.error(f"Invalid format for --minMatch: {e}")
            sys.exit(1)
    else:
        MIN_MATCH_BP = constants.DEFAULT_MIN_MATCH_BP
        MIN_MATCH_PERC = constants.DEFAULT_MIN_MATCH_PERC

#=====================================================================================================
# Utility functions
#=====================================================================================================

def get_next_letter(last_letter):
    try:
        next_letter_pos = np.where(ASCII_LETTERS == last_letter[-1])[0][0] + 1
        next_letter = list(string.ascii_letters)[next_letter_pos]
        return last_letter[:-1] + next_letter if len(last_letter) > 1 else next_letter
    except IndexError:
        return last_letter + 'a'

def get_next_id(qname):
    if qname not in record:
        vid = qname + 'a'
        record[qname] = [vid]
        return vid
    else:
        last_letter = re.search('[a-zA-Z]*$', record[qname][-1]).group(0)
        vid = qname + get_next_letter(last_letter)
        record[qname].append(vid)
        return vid

def get_attribute(attributes, pattern):
    re_attr = pattern.search(attributes)
    return re_attr.group(1) if re_attr else ''

def get_gene_lookup(tx_ref_file):
    ref_trees, ex_trees, ex_ref_out = None, None, None
    if tx_ref_file == '':
        return ref_trees, ex_trees, ex_ref_out

    logging.info('Generating lookup for genes...')
    tx_ref = pd.read_csv(tx_ref_file, comment='#', sep='\t', header=None, low_memory=False, usecols=[0, 2, 3, 4, 8])
    
    tx_ref['gene_id'] = tx_ref[8].apply(lambda x: get_attribute(x, RE_GENE_ID))
    tx_ref['gene'] = tx_ref[8].apply(lambda x: get_attribute(x, RE_GENE_NAME))

    gn_ref = tx_ref[[0, 3, 4, 'gene_id', 'gene']].copy()
    gn_ref.columns = ['chrom', 'start', 'end', 'gene_id', 'gene']
    adj_check = (gn_ref.gene_id != gn_ref.gene_id.shift()).cumsum()
    # Using string aggregations to avoid FutureWarning
    gn_ref = gn_ref.groupby(['chrom', 'gene_id', 'gene', adj_check],
                            as_index=False, sort=False).agg({'start': 'min', 'end': 'max'})
    gn_ref = gn_ref.drop_duplicates()

    ref_trees = {}
    chroms = np.unique(gn_ref.chrom.values)
    for chrom in chroms:
        chr_ref = gn_ref[gn_ref.chrom == chrom].drop_duplicates()
        ref_tree = IntervalTree()
        for s, e, g in zip(chr_ref['start'].values, chr_ref['end'].values, chr_ref['gene'].values):
            if g != '':
                ref_tree.addi(s-1, e, g)
        ref_trees[chrom] = ref_tree

    ex_ref = tx_ref[tx_ref[2] == 'exon']
    ex_tree_dfs = []
    ex_trees = {}
    
    for chrom in chroms:
        chr_ref = ex_ref[ex_ref[0] == chrom].drop_duplicates()
        ex_tree = IntervalTree()
        for s, e in zip(chr_ref[3].values, chr_ref[4].values):
            ex_tree.addi(s-1, e)
        ex_tree.merge_overlaps()
        
        tmp = pd.DataFrame([(chrom, tree.begin, tree.end) for tree in ex_tree],
                           columns=['chrom', 'start', 'end'])
        ex_tree_dfs.append(tmp)
        ex_trees[chrom] = ex_tree

    if ex_tree_dfs:
        ex_ref_out = pd.concat(ex_tree_dfs, ignore_index=True)
    else:
        ex_ref_out = pd.DataFrame(columns=['chrom', 'start', 'end'])

    return ref_trees, ex_trees, ex_ref_out

def get_chr_ref(read, ex_ref):
    chr_ref = ex_ref[ex_ref.chrom == read.reference_name]
    if len(chr_ref) == 0:
       logging.debug('WARNING: reference chromosome %s (from read %s) not found in reference', 
                    read.reference_name, read.query_name)
    return chr_ref

def get_juncs(tx):
    starts = tx['exonStarts'].split(',')[1:]
    ends = tx['exonEnds'].split(',')[:-1]
    chroms = [tx['chrom']] * len(starts)
    return list(zip(chroms, ends, starts))

def get_junc_lookup(junc_file):
    logging.info('Generating lookup for known splice junctions...')
    genref = pd.read_csv(junc_file, sep='\t', low_memory=False)
    junc_info = genref.apply(lambda tx: get_juncs(tx), axis=1)

    junc_list = ['%s:%s-%s' % (c, s, e) for jv in junc_info.values for c, s, e in jv]
    juncs = set(junc_list)

    loc_list = [['%s:%s' % (c, s), '%s:%s' % (c, e)] for jv in junc_info.values for c, s, e in jv]
    flat_locs = [l for loc in loc_list for l in loc]
    locs = set(flat_locs)
    
    logging.info('Finished generating splice junction lookup')
    return (juncs, locs)

def get_tx_juncs(read):
    starts, ends = zip(*read.get_blocks())
    blocks = IntervalTree()
    for s, e in zip(starts, ends):
        blocks.addi(s, e)
    blocks.merge_overlaps(strict=False)

    starts = np.sort([block.begin for block in blocks])
    ends = np.sort([block.end for block in blocks])

    chroms = [read.reference_name] * (len(starts)-1)
    tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
    tx_juncs = [junc for junc in tx_juncs if (junc[2] - junc[1]) >= MIN_GAP]
    return tx_juncs

def get_chrom_ref_tree(chrom, ref_trees):
    if chrom in ref_trees: return ref_trees[chrom]
    
    alt_chrom = 'chr%s' % chrom if chrom != 'MT' else 'chrM'
    if alt_chrom in ref_trees: return ref_trees[alt_chrom]

    alt_chrom_2 = chrom.split('chr')[-1]
    alt_chrom_2 = 'MT' if alt_chrom_2 == 'M' else alt_chrom_2
    if alt_chrom_2 in ref_trees: return ref_trees[alt_chrom_2]

    if chrom not in getattr(get_chrom_ref_tree, "logged_missing", set()):
        logging.info('WARNING: reference chromosome %s not found in supplied reference.', chrom)
        if not hasattr(get_chrom_ref_tree, "logged_missing"):
            get_chrom_ref_tree.logged_missing = set()
        get_chrom_ref_tree.logged_missing.add(chrom)
    return None

def do_any_read_blocks_overlap_exons(read, ex_trees, bam_idx):
    for r in bam_idx.find(read.query_name):
        chr_ex = get_chrom_ref_tree(r.reference_name, ex_trees)
        if not chr_ex: continue
        for s, e in r.get_blocks():
            if chr_ex.overlaps(s, e):
                return True
    return False

def get_overlapping_genes(read, ref_trees):
    ref_tree = get_chrom_ref_tree(read.reference_name, ref_trees)
    if not ref_tree:
        return ''
    blocks = read.get_blocks()
    genes = set()
    for block in blocks:
        gns = ref_tree.overlap(block[0], block[1])
        for gn in gns:
            if gn.data != '':
                genes.add(gn.data)
    return '|'.join(sorted(list(genes)))

#=====================================================================================================
# Annotation functions
#=====================================================================================================

def annotate_gaps(cv, read, ci_file):
    gap_idxs = [idx for idx, gap in enumerate(read.cigar) if gap[0] in constants.GAPS and gap[1] >= MIN_GAP]
    for gap_idx in gap_idxs:
        cigar = read.cigar[gap_idx]
        cv.vsize = int(cigar[1])

        block_idx = 0 if gap_idx == 0 else np.max(np.where(np.array([b[0] for b in cv.blocks]) < gap_idx)[0])
        block = cv.blocks[block_idx][1]
        cv.pos = int(block[1])
        cv.cpos = sum([v for c,v in read.cigar[:gap_idx] if c in constants.AFFECT_CONTIG])
        cv.cvsize = cv.vsize if read.cigar[gap_idx][0] == constants.CIGAR['insertion'] else 0

        if cigar[0] == constants.CIGAR['insertion']:
            cv.cvtype = 'INS'
            seq_pos1 = sum([v for c,v in read.cigar[:gap_idx] if c in constants.AFFECT_CONTIG \
                                                                 and c != constants.CIGAR['hard-clip']])
            seq_pos2 = seq_pos1 + cv.cvsize
            seq = read.query_sequence[(seq_pos1-1):seq_pos2]
            cv.ref, cv.alt = seq[:1], seq
        else:
            cv.cvtype = 'DEL'
            seq_pos1 = sum([v for c,v in read.cigar[:gap_idx] if c in constants.AFFECT_REF])
            seq_pos2 = seq_pos1 + cv.vsize
            seq = read.get_reference_sequence()[(seq_pos1-1):seq_pos2]
            cv.ref, cv.alt  = seq, seq[:1]

        cv.vid = get_next_id(read.query_name)
        CrypticVariant.write_contig_info(ci_file, cv)
        print(cv.vcf_output())

def annotate_softclips(cv, read, ci_file):
    sc_idxs = [idx for idx, clip in enumerate(read.cigar) if clip[0] == constants.CIGAR['soft-clip'] and clip[1] >= MIN_CLIP]
    for sc_idx in sc_idxs:
        cigar = read.cigar[sc_idx]
        sc_left = sc_idx == 0
        if 'FUS' not in cv.cvtype and cv.cvtype == '.':
            cv.cvtype = 'UN'
        cv.vsize, cv.cvsize = int(cigar[1]), int(cigar[1])
        
        block_idx = 0 if sc_idx == 0 else np.max(np.where(np.array([b[0] for b in cv.blocks]) < sc_idx)[0])
        block = cv.blocks[block_idx][1]
        cv.pos = int(block[0]) + 1 if sc_left else int(block[1])
        
        rcigar = read.cigar[::-1] if cv.cstrand == '-' else read.cigar
        idx = len(rcigar) - sc_idx - 1 if cv.cstrand == '-' else sc_idx + 1
        cv.cpos = sum([v for c,v in rcigar[:idx] if c in constants.AFFECT_CONTIG])
        
        varseq = read.query_sequence[cv.cpos:(cv.cpos+cv.cvsize)]
        refseq = read.get_reference_sequence()
        cv.ref = refseq[0] if sc_left else refseq[-1]
        cv.alt = '%s]%s:%d]%s' % (varseq, cv.chrom, cv.pos-1, cv.ref) if sc_left else \
                 '%s[%s:%d[%s' % (cv.ref, cv.chrom, cv.pos+1, varseq)
        
        cv.vid = get_next_id(read.query_name)
        CrypticVariant.write_contig_info(ci_file, cv)
        print(cv.vcf_output())

def annotate_block_right(cv, read, cpos, olapping, block, block_idx):
    qseq, rseq = bh.get_block_sequence(read, block_idx)
    seq_left_pos = block[1] - max(olapping.end)
    cv.ref, cv.alt = rseq[(-seq_left_pos):], ']' + qseq[(-seq_left_pos):]
    cv.cpos, cv.pos = cpos, max(olapping.end) + 1
    cv.vsize, cv.cvsize = abs(len(cv.alt)-1 - len(cv.ref)), len(cv.alt)-1
    return cv

def annotate_block_left(cv, read, cpos, olapping, block, block_idx):
    qseq, rseq = bh.get_block_sequence(read, block_idx)
    seq_right_pos = min(olapping.start) - block[0]
    cv.ref, cv.alt = rseq[:seq_right_pos], qseq[:seq_right_pos] + '['
    cv.cpos, cv.pos = cpos, min(olapping.start) - len(cv.ref) + 1
    cv.vsize, cv.cvsize = abs(len(cv.alt)-1 - len(cv.ref)), len(cv.alt)-1
    return cv

def annotate_blocks(cv, read, chr_ref, ci_file):
    cv.parid = '.' 
    novel_blocks = [(idx, block) for idx, block in cv.blocks if bh.is_novel_block(block, chr_ref, MIN_CLIP)]
    rcigar = read.cigar[::-1] if cv.cstrand == '-' else read.cigar
    
    for block_idx, block in novel_blocks:
        idx = len(rcigar) - block_idx - 1 if cv.cstrand == '-' else block_idx
        cpos1 = sum([v for c,v in rcigar[:idx] if c in constants.AFFECT_CONTIG])
        cpos2 = sum([v for c,v in rcigar[:idx+1] if c in constants.AFFECT_CONTIG])

        olapping = chr_ref[np.logical_and(block[0] < chr_ref.start, block[1] > chr_ref.end)]
        left = chr_ref[np.logical_and(block[1] > chr_ref.start, block[1] <= chr_ref.end)]
        right = chr_ref[np.logical_and(block[0] >= chr_ref.start, block[0] < chr_ref.end)]

        if len(left) > 0 and len(right) > 0: # RI
            cv.cvtype = 'RI'
            qseq, rseq = bh.get_block_sequence(read, block_idx)
            seq_right_pos = block[1] - min(left.start)
            seq_left_pos = max(right.end) - block[0]

            cv.pos = block[0] + seq_left_pos + 1
            cv.ref = rseq[seq_left_pos:(-seq_right_pos)]
            cv.alt = ']' + qseq[seq_left_pos:(-seq_right_pos)] + '['
            cv.cpos = cpos1 + seq_left_pos
            cv.vsize, cv.cvsize = abs(len(cv.alt)-2 - len(cv.ref)), len(cv.alt)-2
            cv.vid = get_next_id(read.query_name)
            
        elif len(olapping) > 0: # EE both sides
            cv.cvtype = 'EE'
            cv = annotate_block_left(cv, read, cpos2, olapping, block, block_idx)
            cv.vid = get_next_id(read.query_name)
            print(cv.vcf_output())
            CrypticVariant.write_contig_info(ci_file, cv)

            cv = annotate_block_right(cv, read, cpos1, olapping, block, block_idx)
            cv.vid = get_next_id(read.query_name)
            
        elif len(left) > 0: # EE Left
            cv.cvtype = 'EE'
            cv = annotate_block_left(cv, read, cpos2, left, block, block_idx)
            cv.vid = get_next_id(read.query_name)
            
        elif len(right) > 0: # EE Right
            cv.cvtype = 'EE'
            cv = annotate_block_right(cv, read, cpos1, right, block, block_idx)
            cv.vid = get_next_id(read.query_name)
            
        else: # NE
            qseq, rseq = bh.get_block_sequence(read, block_idx)
            cv.ref, cv.alt = rseq, '[' + qseq + ']'
            cv.pos, cv.cvtype = block[0] + 1, 'NE'
            cv.cpos = cpos1
            cv.vid = get_next_id(read.query_name)
            cv.vsize, cv.cvsize = abs(len(cv.alt)-2 - len(cv.ref)), len(cv.alt)-2

        print(cv.vcf_output())
        CrypticVariant.write_contig_info(ci_file, cv)

def annotate_fusion(args, read, juncs, bam_idx, ex_ref, ref_trees, outbam):
    logging.info("Starting annotate_fusion for read: %s", read.query_name)
    try:
        all_reads = [r for r in bam_idx.find(read.query_name)]
        if len(all_reads) <= 1:
            if len(all_reads) == 1:
                # return outcome of single read annotation
                return annotate_single_read(args, all_reads[0], juncs, ex_ref, ref_trees, outbam)
            return False
 
        def pair_score(r1, r2):
            try:
                if r1.reference_name != r2.reference_name: return (2, 0)
                dist = abs(r1.reference_start - r2.reference_start)
                return (1, dist)
            except Exception:
                return (0, 0)
 
        best_pair = None
        best_key = (-1, -1)
        for i in range(len(all_reads)):
            for j in range(i + 1, len(all_reads)):
                key = pair_score(all_reads[i], all_reads[j])
                if key > best_key:
                    best_key, best_pair = key, (all_reads[i], all_reads[j])
 
        if best_pair is None: return False
        r1, r2 = best_pair
    except Exception as e:
        logging.warning("Exception in annotate_fusion: %s", e)
        return False

    ci_file = args.contig_info_output
    cv1 = CrypticVariant().from_read(r1)
    cv2 = CrypticVariant().from_read(r2)
    cv1.cvtype, cv2.cvtype = 'FUS', 'FUS'
    cv1.genes = get_overlapping_genes(r1, ref_trees)
    cv2.genes = get_overlapping_genes(r2, ref_trees)
    
    if cv1.genes == cv2.genes and cv1.genes != '':
        cv1.cvtype, cv2.cvtype = 'IGR', 'IGR'
        
    if cv1.genes == '' and cv2.genes == '':
        record[read.query_name] = []
        return False

    def process_clip(read, cv):
        candidates = [(idx, clip[0]) for idx, clip in enumerate(read.cigar)
                      if clip[0] in (constants.CIGAR['hard-clip'], constants.CIGAR['soft-clip'])]
        if not candidates: return None, None, None
        
        clip_idx, clip_type = next(((i, t) for i, t in candidates if t == constants.CIGAR['hard-clip']), candidates[0])
        clip_left = (clip_idx == 0)
        
        block_idx = 0 if clip_left else np.max(np.where(np.array([b[0] for b in cv.blocks]) < clip_idx)[0])
        block = cv.blocks[block_idx][1]
        cv.pos = int(block[0]) if clip_left else int(block[1])
        return clip_idx, clip_left, True

    clip_idx1, clip_left1, ok1 = process_clip(r1, cv1)
    clip_idx2, clip_left2, ok2 = process_clip(r2, cv2)
    if not ok1 or not ok2: return False

    try:
        refseq1 = r1.get_reference_sequence()
        refseq2 = r2.get_reference_sequence()
    except Exception:
        return False

    varseq1, varseq2 = '', ''
    cv1.ref = refseq1[0] if clip_left1 else refseq1[-1]
    cv2.ref = refseq2[0] if clip_left2 else refseq2[-1]

    if r1.is_reverse == r2.is_reverse:
        cv1.alt = ('%s]%s:%d]%s' % (varseq1, cv2.chrom, cv2.pos-1, cv1.ref)) if clip_left1 else \
                  ('%s[%s:%d[%s' % (cv1.ref, cv2.chrom, cv2.pos+1, varseq1))
        cv2.alt = ('%s]%s:%d]%s' % (varseq2, cv1.chrom, cv1.pos-1, cv2.ref)) if clip_left2 else \
                  ('%s[%s:%d[%s' % (cv2.ref, cv1.chrom, cv1.pos+1, varseq2))
    else:
        cv1.alt = ('%s[%s:%d[%s' % (varseq1, cv2.chrom, cv2.pos-1, cv1.ref)) if clip_left1 else \
                  ('%s]%s:%d]%s' % (cv1.ref, cv2.chrom, cv2.pos+1, varseq1))
        cv2.alt = ('%s[%s:%d[%s' % (varseq2, cv1.chrom, cv1.pos-1, cv2.ref)) if clip_left2 else \
                  ('%s]%s:%d]%s' % (cv2.ref, cv1.chrom, cv1.pos+1, varseq2))

    cv1.vid = get_next_id(r1.query_name)
    cv2.vid = get_next_id(r2.query_name)
    cv1.parid, cv2.parid = cv2.vid, cv1.vid

    rcigar = r1.cigar[::-1] if cv1.cstrand == '-' else r1.cigar
    idx = len(rcigar) - clip_idx1 - 1 if cv1.cstrand == '-' else clip_idx1 + 1
    cv1.cpos = sum([v for c, v in rcigar[:idx] if c in constants.AFFECT_CONTIG])

    print(cv1.vcf_output())
    print(cv2.vcf_output())
    
    outbam.write(r1)
    outbam.write(r2)
    CrypticVariant.write_contig_info(ci_file, cv1, cv2)

    annotate_single_read(args, r1, juncs, ex_ref, ref_trees, genes=cv1.genes, fusion_processed=True)
    annotate_single_read(args, r2, juncs, ex_ref, ref_trees, genes=cv2.genes, fusion_processed=True)
    return True

def annotate_juncs(cv, read, locs, novel_juncs, ci_file):
    cv.cvsize = 0
    for junc in novel_juncs:
        pos1, pos2 = int(junc[1]), int(junc[2])
        junc_idx = [idx for idx, block in cv.blocks if block[1] == pos1][0]
        junc_type = read.cigar[junc_idx+1][0]
        if junc_type in constants.GAPS or junc_type == constants.CIGAR['soft-clip']:
            continue

        cp = CrypticVariant().from_read(read)
        cp.genes = cv.genes
        varseq, refseq = '', read.get_reference_sequence()

        rcigar = read.cigar[::-1] if cv.cstrand == '-' else read.cigar
        idx = len(rcigar) - junc_idx - 1 if cv.cstrand == '-' else junc_idx + 1
        cpos = sum([v for c,v in rcigar[:idx] if c in constants.AFFECT_CONTIG])
        cv.cpos, cp.cpos = cpos, cpos
        cv.pos, cp.pos = pos1, pos2+1

        rpos = sum([v for c,v in read.cigar[:(junc_idx+1)] if c in constants.AFFECT_REF])
        cv.ref, cp.ref = refseq[rpos-1], refseq[rpos]
        cv.alt = '%s[%s:%d[%s' % (cv.ref, cv.chrom, cv.pos, varseq)
        cp.alt = '%s]%s:%d]%s' % (varseq, cp.chrom, cp.pos, cp.ref)

        cv.vid = get_next_id(read.query_name)
        cp.vid = get_next_id(read.query_name)
        cv.parid, cp.parid = cp.vid, cv.vid

        loc_left = '%s:%d' % (cv.chrom, pos1)
        loc_right = '%s:%d' % (cv.chrom, pos2)
        in_left = loc_left in locs
        in_right = loc_right in locs
        
        if not in_left and not in_right:
            cv.cvtype, cp.cvtype = 'NEJ', 'NEJ'
        elif not (in_left and in_right):
            cv.cvtype, cp.cvtype = 'PNJ', 'PNJ'
        else:
            cv.cvtype, cp.cvtype = 'AS', 'AS'

        print(cv.vcf_output())
        print(cp.vcf_output())
        CrypticVariant.write_contig_info(ci_file, cv, cp)

def annotate_single_read(args, read, juncs, ex_ref, ref_trees, outbam=None, genes='', fusion_processed=False):
    ci_file = args.contig_info_output
    genes = get_overlapping_genes(read, ref_trees) if genes == '' else genes
    fusion = any([op in [constants.CIGAR['hard-clip'], constants.CIGAR['soft-clip']] and val >= MIN_CLIP for op, val in read.cigar])
    if genes == '' and not fusion:
        logging.debug('No gene(s) intersecting read %s; skipping' % read.query_name)
        return False

    has_gaps = any([op in constants.GAPS and val >= MIN_GAP for op, val in read.cigar])
    has_scs = any([op == constants.CIGAR['soft-clip'] and val >= MIN_CLIP for op, val in read.cigar]) and not fusion_processed
    
    tx_juncs = get_tx_juncs(read)
    unknown_juncs = ['%s:%s-%s' % (c, s, e) not in juncs[0] for c, s, e in tx_juncs]
    has_novel_juncs = any(unknown_juncs)

    chr_ref = get_chr_ref(read, ex_ref)
    has_novel_blocks = any([bh.is_novel_block(block, chr_ref, MIN_CLIP) for block in read.get_blocks()])

    if has_gaps or has_scs or has_novel_juncs or has_novel_blocks:
        cv = CrypticVariant().from_read(read)
        cv.genes = genes
        if has_gaps: annotate_gaps(cv, read, ci_file)
        if has_scs: annotate_softclips(cv, read, ci_file)
        if has_novel_juncs:
            novel_juncs_list = [list(x) for x in np.array(tx_juncs)[unknown_juncs]]
            annotate_juncs(cv, read, juncs[1], novel_juncs_list, ci_file)
        if has_novel_blocks: annotate_blocks(cv, read, chr_ref, ci_file)
        if outbam: outbam.write(read)
        return True
    else:
        logging.debug('Nothing to annotate for read %s' % read.query_name)
        return False

def annotate_contigs(args):
    ref_trees, ex_tree, ex_ref = get_gene_lookup(args.tx_ref_file)
    juncs = get_junc_lookup(args.junc_file)

    bam = pysam.AlignmentFile(args.bam_file, 'rc')
    bam_idx = pysam.IndexedReads(bam, multiple_iterators=True)
    bam_idx.build()
    
    outbam_file_unsort = '%s_unsorted.bam' % os.path.splitext(args.output_bam)[0]
    outbam = pysam.AlignmentFile(outbam_file_unsort, 'wb', template=bam)
    ci_file = args.contig_info_output

    logging.info('Checking contigs for non-reference content...')
    
    # Counter for variants found
    found_variants_count = 0
    
    for read in bam.fetch(multiple_iterators=True):
        if read.reference_id < 0: continue

        if not do_any_read_blocks_overlap_exons(read, ex_tree, bam_idx): continue

        if read.query_name in record: continue

        rlen = read.reference_length
        qlen = float(read.query_length)
        if qlen == 0: continue

        if (rlen < MIN_MATCH_BP) or (rlen / qlen) < MIN_MATCH_PERC:
            logging.debug('Skipping contig %s: not enough bases match reference' % read.query_name)
            continue

        allmatch = all([op == constants.CIGAR['match'] for op, val in read.cigar])
        if len(read.get_blocks()) == 1 and allmatch:
            chr_ex = get_chrom_ref_tree(read.reference_name, ex_tree)
            s, e = read.get_blocks()[0]
            if not (chr_ex.overlaps(s, s + 1) and chr_ex.overlaps(e - 1, e)):
                continue

        is_clipped = any([op in [constants.CIGAR['hard-clip'], constants.CIGAR['soft-clip']] and val >= MIN_CLIP for op, val in read.cigar])
        
        annotated = False
        if is_clipped:
            annotated = annotate_fusion(args, read, juncs, bam_idx, ex_ref, ref_trees, outbam)
        else:
            annotated = annotate_single_read(args, read, juncs, ex_ref, ref_trees, outbam)
            
        if annotated:
            found_variants_count += 1

    bam.close()
    outbam.close()

    pysam.sort('-o', args.output_bam, outbam_file_unsort)
    pysam.index(args.output_bam)
    os.remove(outbam_file_unsort)
    
    logging.info(f"Total annotated contigs found: {found_variants_count}")
    if found_variants_count == 0:
        logging.warning("No variants found in any contig. Generating empty output files.")

def main():
    args = parse_args()
    set_globals(args)
    init_logging(args.log)
    print(VCF.get_header(args.sample))
    CrypticVariant.write_contig_header(args.contig_info_output)
    
    logging.info("--- Parameter Settings ---")
    logging.info("MIN_CLIP: %s", MIN_CLIP)
    logging.info("MIN_GAP: %s", MIN_GAP)
    logging.info("MIN_MATCH: %s bp, %s percent", MIN_MATCH_BP, MIN_MATCH_PERC)

    annotate_contigs(args)
    logging.info("Contig annotation completed successfully.")

if __name__ == '__main__':
    main()