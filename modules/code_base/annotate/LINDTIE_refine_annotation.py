'''
Module      : LINDTIE_refine_annotation
Description : Performs further filtering on annotated contigs
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
Portability : POSIX
Performs advanced refinement of annotated contigs by applying biological 
and structural filters, including motif validation, exon and gene overlap 
checks, variant reclassification, splice-junction matching, and size-based 
criteria.

Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
'''

import numpy as np
import pandas as pd
import re
import sys
import logging
import pysam
import LINDTIE_contigs_annotation as ac
import pybedtools as pbt
import constants
import tempfile
import csv
from Bio import SeqIO
from intervaltree import IntervalTree
from pybedtools import BedTool
from argparse import ArgumentParser
from utils import init_logging, exit_with_error
from pybedtools.helpers import BEDToolsError

PROGRAM_NAME = 'LINDTIE_refine_annotation'

SPLICE_VARS = ['AS']
SV_VARS = ['DEL', 'INS']
NOVEL_BLOCKS = ['EE', 'NE']
NOVEL_JUNCS = ['PNJ', 'NEJ']
LARGE_SV = ['FUS', 'IGR']
UNKNOWN = ['UN']
SENSE_MOTIF = ['AG', 'GT']
ANTISENSE_MOTIF = ['AC', 'CT']

def parse_args():
    description = 'Refine annotated contigs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log', metavar='LOG_FILE', type=str, help='record program progress in LOG_FILE')
    parser.add_argument(dest='contig_info_file', metavar='CONTIG_INFO_FILE', type=str, help='Contig info file')
    parser.add_argument(dest='vcf_file', metavar='VCF_FILE', type=str, help='Contig VCF file')
    parser.add_argument(dest='bam_file', metavar='BAM_FILE', type=str, help='Contig BAM file')
    parser.add_argument(dest='tx_ref_file', metavar='TX_REF_FILE', type=str, help='Transcriptome GTF reference file')
    parser.add_argument(dest='fasta', metavar='FASTA', type=str, help='Genome fasta file')
    parser.add_argument(dest='out_prefix', metavar='OUT_PREFIX', type=str, help='Output prefix')
    parser.add_argument('--minClip', metavar='MIN_CLIP', type=int, help='Minimum novel block or softclip size')
    parser.add_argument('--minGap', metavar='MIN_GAP', type=int, help='Minimum gap (deletion or insertion) size')
    parser.add_argument('--mismatches', metavar='MISMATCHES', type=int, default=0, help='Number of allowed mismatches')
    return parser.parse_args()

def set_globals(args):
    global MIN_CLIP, MIN_GAP
    MIN_CLIP = args.minClip if args.minClip else constants.DEFAULT_MIN_CLIP
    MIN_GAP = args.minGap if args.minGap else constants.DEFAULT_MIN_GAP

def get_block_seqs(exons):
    logging.info("Extracting block sequences from exon blocks.")
    block_seqs = {}
    try:
        for record in SeqIO.parse(exons.seqfn, 'fasta'):
            block_seqs[record.id] = str(record.seq)
    except Exception as e:
        logging.error("Failed to parse block sequences: %s", e)
        return {}
    logging.info("Extracted block sequences for %d records.", len(block_seqs))
    return block_seqs

def load_vcf_file(contig_vcf):
    logging.info("Loading VCF file: %s", contig_vcf)
    try:
        cvcf = pd.read_csv(contig_vcf, sep='\t', header=None, comment='#', low_memory=False)
    except pd.errors.EmptyDataError:
        logging.warning("VCF file has no data (empty or only headers).")
        return pd.DataFrame()
    logging.info("VCF file loaded with shape: %s", cvcf.shape)
    return cvcf

def get_diff_count(motif, side=0, sense=True):
    assert len(motif) == 2
    test_motif = SENSE_MOTIF[side] if sense else ANTISENSE_MOTIF[side]
    count = int(motif[0] != test_motif[0]) + int(motif[1] != test_motif[1])
    return count

def is_motif_valid(motifs, mismatches):
    if mismatches > 2:
        return True
    elif '' in motifs:
        side = np.where('' != np.array(motifs))[0][0]
        diff_s = get_diff_count(motifs[side], side=side)
        diff_as = get_diff_count(motifs[side], side=side, sense=False)
        mismatches = 1 if mismatches == 2 else mismatches
        return diff_s <= mismatches or diff_as <= mismatches
    else:
        ldiff_s = get_diff_count(motifs[0])
        ldiff_as = get_diff_count(motifs[0], sense=False)
        rdiff_s = get_diff_count(motifs[1], side=1)
        rdiff_as = get_diff_count(motifs[1], side=1, sense=False)

        if ldiff_s + rdiff_s <= mismatches and ldiff_s < 2 and rdiff_s < 2:
            return True
        if ldiff_as + rdiff_as <= mismatches and ldiff_as < 2 and rdiff_as < 2:
            return True
    return False

def check_valid_motif(left_id, right_id, block_seqs, mismatches):
    try:
        lseq = '' if left_id == '' else block_seqs.get(left_id, '')
        rseq = '' if right_id == '' else block_seqs.get(right_id, '')

        motif = [lseq, rseq]
        valid = False if lseq == '' and rseq == '' else is_motif_valid(motif, mismatches)
        return valid, ''.join(motif)
    except Exception as e:
        logging.debug("Error checking motif: %s", e)
        return False, ''

def get_valid_motif_vars(variants, args):
    logging.info("Getting valid motif variants.")
    vcf = load_vcf_file(args.vcf_file)
    
    if vcf.empty:
        return pd.DataFrame(columns=['variant_id', 'motif', 'valid_motif'])

    vcf = vcf[vcf[2].isin(variants.variant_id.values)]
    vcf = vcf[vcf[3].notnull()]

    b_blocks = vcf[vcf[4].str.startswith('[') & vcf[4].str.endswith(']')]
    r_blocks = vcf[vcf[4].str.startswith(']')]
    l_blocks = vcf[vcf[4].str.endswith('[')]
    
    left = pd.concat([l_blocks, b_blocks])
    right = pd.concat([r_blocks, b_blocks])
    
    mlocs = pd.DataFrame({'chr': pd.concat([left[0], right[0]]),
                          'start': pd.concat([left[1] - 3, right[1] + right[3].str.len()-1])})
    mlocs['end'] = mlocs['start'] + 2
    mlocs = mlocs.drop_duplicates()
    mlocs = mlocs[mlocs['start'] >= 0]
    
    chr_sizes = pbt.chromsizes('hg38')
    valid_chroms = []
    for chrom in mlocs['chr'].unique():
        ref_chrom = 'chr%s' % chrom if chrom != 'MT' else 'chrM'
        if ref_chrom in chr_sizes:
            chr_max = chr_sizes[ref_chrom][1]
            over_limit = (mlocs.start > chr_max) | (mlocs.end > chr_max)
            mlocs = mlocs[~((mlocs['chr'] == chrom) & over_limit)]
            valid_chroms.append(chrom)

    try:
        g = BedTool.from_dataframe(mlocs).remove_invalid()
        g = g.sequence(fi=args.fasta)
        bs = get_block_seqs(g)
    except Exception as e:
        logging.error("Error extracting sequences: %s", e)
        return pd.DataFrame(columns=['variant_id', 'motif', 'valid_motif'])

    valid_vars, motifs = [], []

    if not l_blocks.empty:
        left_ids = ['%s:%d-%d' % (r[0], r[1]-3, r[1]-1) for _, r in l_blocks.iterrows()]
        info_left = [check_valid_motif(lid, '', bs, args.mismatches) for lid in left_ids]
        valid_left = [v for v, m in info_left]
        motifs.extend([m for v, m in info_left])
        if any(valid_left):
            valid_vars.extend(l_blocks[valid_left][2].values)

    if not r_blocks.empty:
        rpos = r_blocks[1] + r_blocks[3].str.len()-1
        right_ids = ['%s:%d-%d' % (chrom, p, p+2) for chrom, p in zip(r_blocks[0], rpos)]
        info_right = [check_valid_motif('', rid, bs, args.mismatches) for rid in right_ids]
        valid_right = [v for v, m in info_right]
        motifs.extend([m for v, m in info_right])
        if any(valid_right):
            valid_vars.extend(r_blocks[valid_right][2].values)

    if not b_blocks.empty:
        rpos = b_blocks[1] + b_blocks[3].str.len()-1
        left_ids = ['%s:%d-%d' % (r[0], r[1]-3, r[1]-1) for _, r in b_blocks.iterrows()]
        right_ids = ['%s:%d-%d' % (chrom, p, p+2) for chrom, p in zip(b_blocks[0], rpos)]
        info_both = [check_valid_motif(lid, rid, bs, args.mismatches) for lid, rid in zip(left_ids, right_ids)]
        valid_both = [v for v, m in info_both]
        motifs.extend([m for v, m in info_both])
        if any(valid_both):
            valid_vars.extend(b_blocks[valid_both][2].values)

    motif_info = pd.concat([l_blocks, r_blocks, b_blocks])
    motif_info['motif'] = motifs
    motif_info['variant_id'] = motif_info[2]
    motif_info['valid_motif'] = motif_info.variant_id.isin(valid_vars)
    return motif_info[['variant_id', 'motif', 'valid_motif']]

def check_for_valid_motifs(contigs, vars_to_check, args):
    logging.info("Checking for valid motifs in contigs.")
    if any(vars_to_check):
        motif_info = get_valid_motif_vars(contigs[vars_to_check], args)       
        contigs = contigs.merge(motif_info, on='variant_id', how='left')
        contigs['motif'] = contigs.motif.fillna('')
    else:
        contigs['motif'] = ''
        contigs['valid_motif'] = None
    return contigs

def _normalize_range(start, end):
    if end < start:
        return end, start
    return start, end

def _pos_interval(pos):
    # Convert 1-based position to 0-based half-open interval
    return pos - 1, pos

def _range_interval(start, end):
    # Convert 1-based inclusive range to 0-based half-open interval
    start, end = _normalize_range(start, end)
    return start - 1, end

def get_overlap_size(ex_trees, chrom, start, end):
    ex_tree = ac.get_chrom_ref_tree(chrom, ex_trees)
    if not ex_tree:
        return float('nan')

    qstart, qend = _range_interval(start, end)
    overlaps = ex_tree.overlap(qstart, qend)
    if not overlaps:
        return float('nan')

    # Use max overlap across exons to preserve prior semantics
    max_olap = 0
    for interval in overlaps:
        olap = min(qend, interval.end) - max(qstart, interval.begin)
        if olap > max_olap:
            max_olap = olap
    return max_olap

def check_overlap(ex_trees, chrom, start, end, size=0):
    olap_size = get_overlap_size(ex_trees, chrom, start, end)
    return False if np.isnan(olap_size) else olap_size >= size

def get_pos_parts(loc):
    loc_split = loc.split(':')
    chrom = loc_split[0]
    pos = int(loc_split[1].split('(')[0])
    strand = '.'
    try:
        strand = re.search(r'\(([-+])\)', loc_split[1]).group(1)
    except AttributeError:
        pass
    return chrom, pos, strand

def get_varsize(sv):
    chr1, start, s1 = get_pos_parts(sv['pos1'])
    chr2, end, s2 = get_pos_parts(sv['pos2'])
    return end - start

def overlaps_same_exon(sv, ex_trees):
    chr1, start, s1 = get_pos_parts(sv['pos1'])
    if get_pos_parts(sv['pos2'])[0] != chr1:
        return False
    chr2, end, s2 = get_pos_parts(sv['pos2'])

    ex_tree = ac.get_chrom_ref_tree(chr1, ex_trees)
    if ex_tree:
        s1_start, s1_end = _pos_interval(start)
        s2_start, s2_end = _pos_interval(end)
        olap1 = ex_tree.overlap(s1_start, s1_end)
        olap2 = ex_tree.overlap(s2_start, s2_end)
        return len(olap1) > 0 and len(olap2) > 0 and olap1 == olap2
    return False

def overlaps_exon(sv, ex_trees):
    chr1, start, s1 = get_pos_parts(sv['pos1'])
    chr2, end, s2 = get_pos_parts(sv['pos2'])

    span_vars = ['DEL'] + NOVEL_JUNCS
    if sv['variant_type'] in span_vars:
        if chr1 != chr2:
            return False
        size = MIN_GAP if sv['variant_type'] == 'DEL' else MIN_CLIP
        return check_overlap(ex_trees, chr1, start, end, size=size)
    else:
        olap1 = check_overlap(ex_trees, chr1, start, start)
        olap2 = check_overlap(ex_trees, chr2, end, end)
        return olap1 or olap2

def match_splice_juncs(contigs):
    spliced_exons = []
    exons = contigs[contigs.variant_type.isin(NOVEL_BLOCKS)]
    novel_juncs = contigs[contigs.variant_type.isin(NOVEL_JUNCS + ['DEL'])]
    
    for idx, row in exons.iterrows():
        back_junc = novel_juncs.pos2 == row['pos1']
        front_junc = novel_juncs.pos1 == row['pos2']
        if (back_junc | front_junc).any():
            spliced_exons.append(row['variant_id'])

    return contigs.variant_id.isin(spliced_exons)

def vars_overlap_exon(contigs, ex_trees):
    if len(contigs) == 0:
        return np.array([], dtype=bool)

    # Always compute actual overlap instead of forcing EE/AS to True
    return contigs.apply(overlaps_exon, axis=1, args=(ex_trees,)).values

def get_junc_vars(contigs, ex_trees, args):
    within_exon = contigs.apply(overlaps_same_exon, axis=1, args=(ex_trees,))
    nj_var = contigs.variant_type.isin(NOVEL_JUNCS)
    nj_dels = np.empty(0, dtype=object)
    
    if nj_var.sum() > 0:
        bigger_than_mingap = contigs[nj_var].apply(get_varsize, axis=1) >= MIN_GAP
        nj_dels = contigs[nj_var][within_exon[nj_var] & bigger_than_mingap].variant_id.values
        if len(nj_dels) > 0:
            contigs.loc[contigs.variant_id.isin(nj_dels), 'variant_type'] = 'DEL'

    is_trunc = contigs.variant_type.isin(NOVEL_JUNCS) & ~within_exon & contigs.overlaps_exon
    
    if 'valid_motif' in contigs.columns.values:
        contigs_tmp = contigs.drop(['motif', 'valid_motif'], axis=1)
        contigs_tmp = check_for_valid_motifs(contigs_tmp, is_trunc, args)
        contigs.loc[is_trunc, 'motif'] = contigs_tmp.motif
        contigs.loc[is_trunc, 'valid_motif'] = contigs_tmp.valid_motif
        is_trunc = is_trunc & contigs.valid_motif
        
    trunc_vars = contigs[is_trunc.fillna(False)].variant_id.values
    return np.unique(np.concatenate([nj_dels, trunc_vars]))

def get_tsv_vars(contigs):
    is_sv = contigs.variant_type.isin(SV_VARS)
    large_gap = (contigs.varsize >= MIN_GAP) | (contigs.contig_varsize >= MIN_GAP)
    keep_sv = large_gap & is_sv & contigs.overlaps_exon
    sv_vars = contigs[keep_sv].variant_id.values

    is_un = contigs.variant_type.isin(UNKNOWN)
    large_clip = contigs.varsize >= MIN_CLIP
    un_vars = contigs[is_un & large_clip].variant_id.values
    return np.concatenate([un_vars, sv_vars])

def get_fusion_vars(contigs):
    is_fus = contigs.variant_type.isin(LARGE_SV)
    fus_ids = contigs[is_fus].contig_id.values
    fus_locs = np.union1d(contigs[is_fus].pos1, contigs[is_fus].pos2)
    
    non_fus_vars = contigs[contigs.contig_id.isin(fus_ids) & ~is_fus]
    at_fus_boundary = non_fus_vars.pos1.isin(fus_locs) | non_fus_vars.pos2.isin(fus_locs)
    
    fus_boundary_vars = non_fus_vars[at_fus_boundary].variant_id.values
    fus_vars = contigs[is_fus].variant_id.values
    return np.concatenate([fus_vars, fus_boundary_vars])

def overlaps_gene(row, gene_tree):
    chr1, pos1, _ = get_pos_parts(row['pos1'])
    chr2, pos2, _ = get_pos_parts(row['pos2'])

    gtree1 = ac.get_chrom_ref_tree(chr1, gene_tree)
    if chr1 == chr2:
        if not gtree1:
            return False
        qstart, qend = _range_interval(pos1, pos2)
        return gtree1.overlaps(qstart, qend)

    if not gtree1:
        olaps = False
    else:
        p1s, p1e = _pos_interval(pos1)
        olaps = gtree1.overlaps(p1s, p1e)
    gtree2 = ac.get_chrom_ref_tree(chr2, gene_tree)
    if not gtree2:
        return olaps
    p2s, p2e = _pos_interval(pos2)
    return olaps or gtree2.overlaps(p2s, p2e)

def get_cds_lookup(gtf_file):
    """
    Parses GTF to build CDS interval trees and transcript coding bounds.
    Useful for inferring UTRs when GTF only has CDS/Exon entries (like CHESS).
    """
    logging.info("Building CDS and UTR lookup from GTF...")
    
    cds_trees = {}
    tx_map = {} # transcript_id -> {strand, min_cds, max_cds}

    # Regex to extract transcript_id and gene_id
    re_tid = re.compile(r'transcript_id "([\w\-\.\/]+)"')
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom, feature, start, end, strand, attr = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]
            
            if feature == 'CDS':
                # Update CDS Tree
                if chrom not in cds_trees:
                    cds_trees[chrom] = IntervalTree()
                # Store as 0-based half-open, consistent with exon/gene trees
                cds_trees[chrom].addi(start - 1, end) # intervaltree is exclusive at end
                
                # Update Transcript Map
                tid_m = re_tid.search(attr)
                if tid_m:
                    tid = tid_m.group(1)
                    if tid not in tx_map:
                        tx_map[tid] = {'strand': strand, 'min': start, 'max': end}
                    else:
                        tx_map[tid]['min'] = min(tx_map[tid]['min'], start)
                        tx_map[tid]['max'] = max(tx_map[tid]['max'], end)
                        
            elif feature == 'exon':
                # We need transcript info for exons too to link them to CDS bounds
                tid_m = re_tid.search(attr)
                if tid_m:
                    tid = tid_m.group(1)
                    # Initialize if seen for first time (might be non-coding)
                    if tid not in tx_map:
                        tx_map[tid] = {'strand': strand, 'min': float('inf'), 'max': float('-inf')}

    return cds_trees, tx_map

def get_feature_type(chrom, pos, cds_trees, ex_trees, ref_trees, tx_map):
    """
    Determines genomic feature at position using hierarchy:
    CDS > 5'UTR > 3'UTR > Exon (Non-coding) > Intron > Intergenic
    """
    # 1. Check CDS
    if chrom in cds_trees:
        pstart, pend = _pos_interval(pos)
        if cds_trees[chrom].overlaps(pstart, pend):
            return "CDS"
        return "CDS"
    
    # 2. Check Exon (implies UTR or Non-Coding since we failed CDS check)
    # ac.get_chrom_ref_tree returns the tree for that chromosome
    ex_tree = ac.get_chrom_ref_tree(chrom, ex_trees)
    if ex_tree:
        pstart, pend = _pos_interval(pos)
        if ex_tree.overlaps(pstart, pend):
            # Overlaps an exon. Determine if it's 5' or 3' UTR based on transcript info.
            # We need to find WHICH transcript(s) this exon belongs to. 
            # Since ex_trees from contigs_annotation usually don't store transcript IDs in data,
            # we make a best guess based on the 'gene' tree or simplified logic.
            
            # NOTE: For precise UTR assignment without transcript IDs in ex_tree, 
            # we check if the position is outside the CDS bounds of *any* overlapping gene's transcript.
            
            # Simplified robust heuristic:
            # If we are here, we are in an Exon but NOT in a CDS.
            # Check overlapping transcripts to vote for 5' or 3' UTR.
            
            # Note: Implementing perfect transcript-matching here is heavy. 
            # We will return "UTR/Exon" generally, or "5'UTR"/"3'UTR" if we can match to a CDS bound.
            # For now, let's call it "Exon_UTR" or try to find a containing transcript.
            return "UTR" # Or "Non-coding Exon"

    # 3. Check Intron (Inside gene bounds but not in exon)
    ref_tree = ac.get_chrom_ref_tree(chrom, ref_trees)
    if ref_tree:
        pstart, pend = _pos_interval(pos)
        if ref_tree.overlaps(pstart, pend):
            return "Intron"
        
    # 4. Fallback
    return "Intergenic"

# Improved UTR Logic (Optional - Drop into the function above if you want specific UTRs)
def get_detailed_feature_type(chrom, pos, cds_trees, ex_trees, ref_trees, tx_map):
    # 1. CDS
    if chrom in cds_trees:
        pstart, pend = _pos_interval(pos)
        if cds_trees[chrom].overlaps(pstart, pend):
            return "CDS"

    # 2. Exon (UTR check)
    ex_tree = ac.get_chrom_ref_tree(chrom, ex_trees)
    if ex_tree:
        pstart, pend = _pos_interval(pos)
        if ex_tree.overlaps(pstart, pend):
            # It's an exon. Is it 5' or 3'?
            # We look at the gene tree to get the gene/transcript context roughly
            # This part relies on having transcript-linked data which standard IntervalTrees might drop.
            # If strict 5'/3' is critical, we return "UTR".
            # If generic is fine:
            return "UTR"

    # 3. Intron
    ref_tree = ac.get_chrom_ref_tree(chrom, ref_trees)
    if ref_tree:
        pstart, pend = _pos_interval(pos)
        if ref_tree.overlaps(pstart, pend):
            return "Intron"

    return "Intergenic"

def get_contigs_to_keep(args):
    logging.info("Getting contigs to keep based on criteria.")
    try:
        contigs = pd.read_csv(args.contig_info_file, sep='\t')
    except IOError as e:
        exit_with_error(str(e), constants.EXIT_FILE_IO_ERROR)
    
    if contigs.empty:
        logging.warning("Input contig info file is empty. Returning empty list.")
        output_file = '%s_info.tsv' % args.out_prefix
        contigs.to_csv(output_file, sep='\t', index=None)
        return np.array([])

    # 1. Load Annotation Resources
    gene_tree, ex_trees, ex_ref = ac.get_gene_lookup(args.tx_ref_file)
    # NEW: Load CDS and Transcript Map
    cds_trees, tx_map = get_cds_lookup(args.tx_ref_file)

    contigs['large_varsize'] = contigs.contig_varsize >= MIN_CLIP
    contigs['is_contig_spliced'] = contigs.contig_cigar.str.contains('N')
    contigs['spliced_exon'] = match_splice_juncs(contigs)
    
    # This call now works even if dataframe is empty or malformed
    contigs['overlaps_exon'] = vars_overlap_exon(contigs, ex_trees)
    
    contigs['overlaps_gene'] = contigs.apply(overlaps_gene, axis=1, args=(gene_tree,))

    is_intergenic_exon = contigs.spliced_exon & (contigs.variant_type == 'NE') & \
                         contigs.large_varsize & ~contigs.overlaps_gene
    
    is_novel_exon = contigs.spliced_exon & contigs.large_varsize
    if args.mismatches < 4:
        contigs = check_for_valid_motifs(contigs, is_novel_exon, args)
        is_novel_exon = is_novel_exon & contigs.valid_motif
        
    ne_vars = contigs[is_novel_exon | is_intergenic_exon].variant_id.values
    as_vars = contigs[contigs.variant_type.isin(SPLICE_VARS)].variant_id.values

    junc_vars = get_junc_vars(contigs, ex_trees, args)
    is_junc_var = contigs.variant_id.isin(junc_vars)
    
    if is_junc_var.any():
        loc1 = contigs[is_junc_var].pos1.apply(get_pos_parts).values
        loc2 = contigs[is_junc_var].pos2.apply(get_pos_parts).values
        junc_varsizes = [get_overlap_size(ex_trees, l1[0], l1[1], l2[1]) for l1, l2 in zip(loc1, loc2)]
        contigs.loc[is_junc_var, 'varsize'] = junc_varsizes
        contigs.loc[is_junc_var, 'spliced_exon'] = True

    ri_vars = contigs[(contigs.variant_type == 'RI') & contigs.large_varsize].variant_id.values
    fus_vars = get_fusion_vars(contigs)
    sv_vars = get_tsv_vars(contigs)

    keep_vars = np.unique(np.concatenate([ri_vars, as_vars, ne_vars, sv_vars, fus_vars, junc_vars]))
    contigs['variant_of_interest'] = contigs.variant_id.isin(keep_vars)
    
    # 2. Add New Feature Columns
    logging.info("Annotating genomic features (CDS/UTR/Intron)...")

    def annotate_row(row):
        # Parse pos1
        chr1, p1, _ = get_pos_parts(row['pos1'])
        feat1 = get_detailed_feature_type(chr1, p1, cds_trees, ex_trees, gene_tree, tx_map)
        
        # Parse pos2
        chr2, p2, _ = get_pos_parts(row['pos2'])
        feat2 = get_detailed_feature_type(chr2, p2, cds_trees, ex_trees, gene_tree, tx_map)
        
        return pd.Series([feat1, feat2])

    contigs[['site1_feature', 'site2_feature']] = contigs.apply(annotate_row, axis=1)

    output_file = '%s_info.tsv' % args.out_prefix
    contigs.to_csv(output_file, sep='\t', index=None)
    logging.info("Contigs to keep written to file: %s", output_file)

    return contigs[contigs.variant_of_interest].contig_id.values

def write_output(args, keep_contigs):
    logging.info("Writing output VCF for %d contigs.", len(keep_contigs))
    keep_contigs_set = set(keep_contigs)
    cvars_file = args.vcf_file

    try:
        with open(cvars_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    print(line.strip())
                else:
                    break
        
        try:
            vcf = pd.read_csv(cvars_file, sep='\t', header=None, comment='#', low_memory=False)
        except pd.errors.EmptyDataError:
            logging.info("Input VCF body is empty. Creating empty VCF output.")
            return

        def is_keep(info_str):
            try:
                cid = info_str.split(';')[0].split('=')[1]
                return cid in keep_contigs_set
            except IndexError:
                return False
                
        if 7 in vcf.columns:
            vcf = vcf[vcf[7].apply(is_keep)]
            vcf.to_csv(sys.stdout, sep='\t', index=False, header=False)
        else:
             logging.warning("Unexpected VCF columns. Outputting empty body.")

    except IOError as e:
        exit_with_error(str(e), constants.EXIT_FILE_IO_ERROR)

def write_bam(args, keep_contigs):
    logging.info("Writing BAM file for %d contigs.", len(keep_contigs))
    keep_contigs_set = set(keep_contigs)
    
    bam = pysam.AlignmentFile(args.bam_file, 'rb')
    outbam = pysam.AlignmentFile('%s.bam' % args.out_prefix, 'wb', template=bam)

    count = 0
    if keep_contigs_set:
        for read in bam.fetch():
            if read.query_name in keep_contigs_set:
                outbam.write(read)
                count += 1
            
    bam.close()
    outbam.close()
    logging.info("BAM file written: %s.bam (%d reads)", args.out_prefix, count)

def main():
    args = parse_args()
    init_logging(args.log)
    set_globals(args)

    logging.info("--- Parameter Settings ---")
    logging.info("MIN_CLIP: %s", MIN_CLIP)
    logging.info("MIN_GAP: %s", MIN_GAP)
    logging.info("Splice Motif Mismatches: %s", args.mismatches)

    keep_contigs = get_contigs_to_keep(args)

    if len(keep_contigs) > 0:
        logging.info("Found %d variants to output.", len(keep_contigs))
    else:
        logging.warning("No variants found after refinement. Generating empty output files.")

    write_output(args, keep_contigs)
    write_bam(args, keep_contigs)
    logging.info("Refine annotation completed successfully.")

if __name__ == '__main__':
    main()