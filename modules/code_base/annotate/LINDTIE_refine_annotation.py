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

def get_overlap_size(ex_trees, chrom, start, end):
    ex_tree = ac.get_chrom_ref_tree(chrom, ex_trees)
    if not ex_tree: return float('nan')

    olap_left = ex_tree.overlaps(start - 1, start)
    olap_right = ex_tree.overlaps(end, end + 1)
    
    if not olap_left and not olap_right:
        return float('nan')
    
    match_left = ex_tree.overlap(start - 1, start)
    match_right = ex_tree.overlap(end, end + 1)
    
    if match_left == match_right or (not olap_left or not olap_right):
        olap_se = match_left if olap_left else match_right
        es, ee = [(x[0], x[1]) for x in olap_se][0]
        return min(ee, end) - start if start >= es else end - max(es, start)
    else:
        es1, ee1 = [(x[0], x[1]) for x in match_left][0]
        size1 = min(ee1, end) - start if start >= es1 else end - max(es1, start)
        
        es2, ee2 = [(x[0], x[1]) for x in match_right][0]
        size2 = min(ee2, end) - start if start >= es2 else end - max(es2, start)
        return max(size1, size2)

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
    if get_pos_parts(sv['pos2'])[0] != chr1: return False
    chr2, end, s2 = get_pos_parts(sv['pos2'])

    ex_tree = ac.get_chrom_ref_tree(chr1, ex_trees)
    if ex_tree:
        olap1 = ex_tree.overlap(start, start+1)
        olap2 = ex_tree.overlap(end, end+1)
        return len(olap1) > 0 and len(olap2) > 0 and olap1 == olap2
    return False

def overlaps_exon(sv, ex_trees):
    chr1, start, s1 = get_pos_parts(sv['pos1'])
    chr2, end, s2 = get_pos_parts(sv['pos2'])

    span_vars = ['DEL'] + NOVEL_JUNCS
    if sv['variant_type'] in span_vars:
        size = MIN_GAP if sv['variant_type'] == 'DEL' else MIN_CLIP
        return check_overlap(ex_trees, chr1, start, end, size=size)
    else:
        olap1 = check_overlap(ex_trees, chr1, start, start + 1)
        olap2 = check_overlap(ex_trees, chr2, end, end + 1)
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
    vars_overlaps_exon = np.array([False] * len(contigs))
    exonic_var = contigs.variant_type.isin(['EE', 'AS'])
    
    # FIXED: Update array directly instead of trying to modify DF with .loc
    vars_overlaps_exon[exonic_var] = True
    
    check_idx = ~exonic_var
    if check_idx.any():
        vars_overlaps_exon[check_idx] = contigs[check_idx].apply(overlaps_exon, axis=1, args=(ex_trees,))
    return vars_overlaps_exon

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
        return gtree1.overlaps(pos1, pos2) if gtree1 else False
    
    olaps = gtree1.overlaps(pos1, pos1 + 1) if gtree1 else False
    gtree2 = ac.get_chrom_ref_tree(chr2, gene_tree)
    return olaps or (gtree2.overlaps(pos2, pos2 + 1) if gtree2 else False)

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

    gene_tree, ex_trees, ex_ref = ac.get_gene_lookup(args.tx_ref_file)
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