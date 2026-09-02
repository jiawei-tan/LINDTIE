'''
Module      : LINDTIE_post_process
Description : Filters and collates novel variant information
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
Portability : POSIX
Performs post-processing of novel variant calls by collating annotation, 
expression, and VAF information; filtering variants by type, gene list, 
sequence quality, and biological criteria; selecting the most plausible 
variant type per contig using a scoring system; and producing fully ranked, 
formatted variant tables. It integrates differential expression, contig 
sequence context, motif validation, and read support to generate a refined, 
biologically informed set of variant annotations.

Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
'''

import numpy as np
import pandas as pd
import re
import sys
import logging
import os
import pysam
from typing import Dict, List, Tuple, Optional, Union
from LINDTIE_refine_annotation import get_pos_parts, get_block_seqs
from argparse import ArgumentParser
from utils import init_logging, exit_with_error
from pybedtools import BedTool
import constants

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1
BED_COLS = ['contig', 'start', 'end', 'name', 'score', 'strand', 'tStart', 'tEnd', 'itemRgb']
SPLIT_LEN = 10 # split variants longer than this many base-pairs into two separate junctions to count reads for
VAR_SEQ_LEN = 40 # extract this many bp for each variant (N / 2 in each direction)

# Subset of vaf_estimates.txt columns to merge (align with LINDTIE_estimate_VAF.R).
# Fisher / chi-square and read-count columns were previously dropped here.
VAF_ESTIMATE_MERGE_COLUMNS = [
    'contig_id',
    'TPM', 'sum_WT_TPM',
    'VAF_sum_WT_TPM',
    'control_gene_count_total', 'control_gene_count_wt',
    'control_contig_count',
    'num_reads_case', 'case_gene_count_wt',
    'fisher_p_val', 'odds_ratio',
]

# Global variables for minimum thresholds
MIN_CLIP = None
MIN_GAP = None

def set_globals(args):
    """
    Set global variables based on command line arguments
    """
    global MIN_CLIP
    global MIN_GAP

    if args.minClip:
        MIN_CLIP = args.minClip
    else:
        MIN_CLIP = constants.DEFAULT_MIN_CLIP

    if args.minGap:
        MIN_GAP = args.minGap
    else:
        MIN_GAP = constants.DEFAULT_MIN_GAP

def get_num_criteria(vtype):
    """
    Return the number of criteria for each variant type based on the criteria table
    
    Args:
        vtype: Variant type string
        
    Returns:
        Integer number of criteria for the variant type
    """
    criteria_counts = {
        'FUS': 3,   # has_clipping, size>min_clip, overlaps_gene
        'IGR': 3,   # has_clipping, size>min_clip, overlaps_gene
        'DEL': 4,   # is_contig_spliced, size>min_gap, overlaps_gene, overlaps_exon
        'INS': 4,   # is_contig_spliced, size>min_gap, overlaps_gene, overlaps_exon
        'AS': 3,    # is_contig_spliced, overlaps_gene, overlaps_exon
        'NEJ': 6,   # is_contig_spliced, size>min_gap, overlaps_gene, overlaps_exon, spliced_exon, valid_motif
        'PNJ': 6,   # is_contig_spliced, size>min_gap, overlaps_gene, overlaps_exon, spliced_exon, valid_motif
        'RI': 4,    # is_contig_spliced, size>min_clip, overlaps_gene, overlaps_exon
        'NE': 5,    # is_contig_spliced, size>min_clip, overlaps_gene, overlaps_exon, spliced_exon
        'EE': 6,    # is_contig_spliced, size>min_clip, overlaps_gene, overlaps_exon, spliced_exon, valid_motif
        'UN': 3     # has_soft_clip, size>min_clip, overlaps_gene
    }
    return criteria_counts.get(vtype, 0)

def calculate_variant_score(variant_row):
    """
    Integrated Scoring System: Structural Confidence + Functional Impact
    """
    vtype = variant_row.get('variant_type', '')
    
    # --- 1. SETUP VARIABLES ---
    try:
        size = abs(int(variant_row.get('varsize', 0)))
        reads = int(variant_row.get('num_reads_case', 0))
        vaf = float(variant_row.get('VAF_sum_WT_TPM', 0))
        contig_varsize = abs(int(variant_row.get('contig_varsize', 0)))
    except (ValueError, TypeError):
        size, reads, vaf, contig_varsize = 0, 0, 0.0, 0

    # Genomic Context fields
    feat1 = variant_row.get('site1_feature', 'Intergenic')
    feat2 = variant_row.get('site2_feature', 'Intergenic')
    is_coding = (feat1 == 'CDS') or (feat2 == 'CDS')

    # --- 2. BASE SCORES (User Optimized) ---
    base_scores = {
        'FUS': 100,  'IGR': 90,   # Major rearrangements
        'DEL': 80,   'INS': 80,   # Indels (potential frameshifts)
        'NEJ': 60,   'AS': 60,    # Splicing anomalies
        'PNJ': 60,   'EE': 40,    # Partial/Extended
        'RI': 40,    'NE': 40,    # Retention/Novel Exons
        'UN': 10                  # Unknown
    }
    score = base_scores.get(vtype, 20)

    # --- 3. STRUCTURAL SCORING (Confidence Checks) ---
    # Apply +/- 15 points based on structural expectations
    
    # helper for splicing checks
    def check_splicing(row, check_exon_overlap=True, expect_exon_overlap=True):
        s = 0
        s += 15 if bool(row.get('is_contig_spliced', False)) else -15
        s += 15 if bool(row.get('overlaps_gene', False)) else -15
        if check_exon_overlap:
            overlaps = bool(row.get('overlaps_exon', False))
            s += 15 if (overlaps == expect_exon_overlap) else -15
        return s

    if vtype in ['INS', 'DEL']:
        # Indels should generally overlap genes/exons to be relevant here
        score += 15 if size > MIN_GAP else -15
        score += check_splicing(variant_row, check_exon_overlap=True, expect_exon_overlap=True)

    elif vtype in ['FUS', 'IGR']:
        # Fusion Logic
        cigar_tuples = variant_row.get('cigar', [])
        has_clipping = False
        if isinstance(cigar_tuples, list):
            # Check constants for hard/soft clip codes
            has_clipping = any([op in [constants.CIGAR.get('hard-clip', 5), 
                                       constants.CIGAR.get('soft-clip', 4)] 
                                and val >= MIN_CLIP for op, val in cigar_tuples])
        
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        
        score += 15 if size > MIN_CLIP else -15
        score += 15 if has_clipping else -15
        score += 15 if overlaps_gene else -15

        # Fusion Refinement (Penalty for Read-Throughs)
        chr1 = str(variant_row.get('chr1', ''))
        chr2 = str(variant_row.get('chr2', ''))
        if chr1 != chr2 and chr1 and chr2:
            score += 40  # Bonus: True Translocation
        elif chr1 == chr2:
            score -= 20  # Penalty: Likely Read-Through (False Positive)

    elif vtype == 'UN':
        cigar_tuples = variant_row.get('cigar', [])
        has_soft_clip = False
        if isinstance(cigar_tuples, list):
             has_soft_clip = any([op == constants.CIGAR.get('soft-clip', 4) 
                                  and val >= MIN_CLIP for op, val in cigar_tuples])

        score += 15 if size > MIN_CLIP else -15
        score += 15 if has_soft_clip else -15
        score += 15 if bool(variant_row.get('overlaps_gene', False)) else -15

    elif vtype in ['EE', 'NE', 'AS', 'NEJ', 'PNJ', 'RI']:
        # Splicing Logic
        # NE and RI should NOT overlap existing exons
        expect_overlap = False if vtype in ['NE', 'RI'] else True
        score += check_splicing(variant_row, check_exon_overlap=True, expect_exon_overlap=expect_overlap)

        # Check for novel spliced exons
        if vtype in ['NE', 'EE', 'NEJ', 'PNJ']:
            score += 15 if bool(variant_row.get('spliced_exon', False)) else -15
        
        # Motif Check (Critical)
        valid_motif = variant_row.get('valid_motif', None)
        if valid_motif is True:
            score += 15
        elif valid_motif is False:
            score -= 30 # Harsh penalty for invalid motifs
            
        # Size Check
        threshold = MIN_GAP if vtype in ['NEJ', 'PNJ'] else MIN_CLIP
        score += 15 if size > threshold else -15

    # --- 4. FUNCTIONAL IMPACT SCORING (The "Killer" Features) ---
    
    # A. Frameshift Bonus (Only if Coding + Indel + Not divisible by 3)
    if is_coding and vtype in ['DEL', 'INS']:
        if size % 3 != 0:
            score += 50  # Massive Bonus: Frameshift
        else:
            score += 20  # Bonus: In-frame Indel

    # B. Retained Intron in CDS (Stop Codon Risk)
    # This rescues your low base score (40) -> 80 if it's dangerous
    if vtype == 'RI' and is_coding:
        score += 40

    # C. Location Context
    # Where does this variant land?
    features_touched = [feat1, feat2]
    if 'CDS' in features_touched:
        score += 20
    elif 'UTR' in features_touched:
        score += 10
    elif 'Intron' in features_touched:
        score -= 10  # Penalty: Intronic noise
    elif 'Intergenic' in features_touched:
        score -= 20  # Penalty: Mapping noise

    # --- 5. EVIDENCE WEIGHTING ---
    
    # Read Support
    if reads < 3:
        score -= 50
    elif reads >= 20:
        score += 20
        
    # VAF (VAF_sum_WT_TPM)
    if vaf < 0.05:
        score -= 30
    elif vaf > 0.2:
        score += 15

    # Contig Support
    if contig_varsize >= 50:
        score += 15
    elif contig_varsize >= 20:
        score += 10
    
    return max(0, score)

def select_best_variant_per_contig(contigs_df):
    """
    Select the best variant type for each contig using improved scoring
    
    Args:
        contigs_df: DataFrame with all variants (before filtering)
        
    Returns:
        DataFrame with one row per contig (best variant selected)
    """
    result_variants = []
    
    # Group by contig_id to handle multiple annotations per contig
    for contig_id, group in contigs_df.groupby('contig_id'):
        if len(group) == 1:
            # Single annotation - just add other_variant_type column
            variant = group.iloc[0].copy()
            variant['other_variant_type'] = ''
            result_variants.append(variant)
        else:
            # Multiple annotations - select best one using scoring
            group_with_scores = []
            
            for idx, row in group.iterrows():
                score = calculate_variant_score(row)
                variant_of_interest = bool(row.get('variant_of_interest', False))
                group_with_scores.append((score, row, variant_of_interest))
                
                # Log scoring details for debugging
                logging.debug(f"Contig {contig_id}, variant {row['variant_type']}: "
                            f"score={score}, variant_of_interest={variant_of_interest}")
            
            # Sort by score (highest first)
            group_with_scores.sort(key=lambda x: x[0], reverse=True)
            
            # Select primary variant (highest score)
            primary_variant = group_with_scores[0][1].copy()
            primary_score = group_with_scores[0][0]
            primary_interest = group_with_scores[0][2]
            
            # Collect other variant types
            other_types = [row['variant_type'] for score, row, interest in group_with_scores[1:]]
            primary_variant['other_variant_type'] = '|'.join(other_types)
            
            # Log selection details
            logging.info(f"Contig {contig_id}: Selected {primary_variant['variant_type']} "
                        f"(score={primary_score}, variant_of_interest={primary_interest}) "
                        f"over alternatives: {other_types}")
            
            result_variants.append(primary_variant)
    
    # Convert back to DataFrame
    return pd.DataFrame(result_variants)

def rank_variants_per_contig(contigs_df):
    """
    Keep all variants and add ranking information per contig using the same
    scoring system as selection. Highest score gets rank 1.

    Args:
        contigs_df: DataFrame with all variants (before filtering/collapse)

    Returns:
        DataFrame with added columns:
            - variant_score (float)
            - rank_within_contig (int)
            - is_primary (bool)
    """
    ranked = contigs_df.copy()
    # Calculate scores for every row
    ranked['variant_score'] = ranked.apply(calculate_variant_score, axis=1)
    # Rank within each contig_id (1 = best)
    ranked['rank_within_contig'] = (
        ranked.groupby('contig_id')['variant_score']
        .rank(method='first', ascending=False)
        .astype(int)
    )
    ranked['is_primary'] = ranked['rank_within_contig'] == 1
    return ranked

def validate_scoring_system(contigs_df):
    """
    Validate that the scoring system is working correctly and provide summary statistics
    
    Args:
        contigs_df: DataFrame after variant selection
    """
    logging.info("=== SCORING SYSTEM VALIDATION ===")
    
    if contigs_df.empty:
        logging.info("No contigs to validate (dataframe is empty).")
        logging.info("=== SCORING SYSTEM VALIDATION COMPLETE ===")
        return

    # Check if column exists
    if 'variant_of_interest' not in contigs_df.columns:
        logging.warning("'variant_of_interest' column missing from dataframe. Skipping detailed validation.")
        logging.info("=== SCORING SYSTEM VALIDATION COMPLETE ===")
        return
    
    # Check how many variants_of_interest were selected
    variants_of_interest = contigs_df[contigs_df['variant_of_interest'] == True]
    logging.info(f"Variants marked as 'variant_of_interest': {len(variants_of_interest)}")
    
    if len(variants_of_interest) > 0:
        logging.info("Variant types for variants_of_interest:")
        for vtype, count in variants_of_interest['variant_type'].value_counts().items():
            logging.info(f"  {vtype}: {count}")
    
    # Check variant type distribution
    logging.info("Final variant type distribution:")
    for vtype, count in contigs_df['variant_type'].value_counts().items():
        logging.info(f"  {vtype}: {count}")
    
    # Check contigs with multiple variant types
    multi_type_contigs = contigs_df[contigs_df['other_variant_type'] != '']
    logging.info(f"Contigs with multiple variant types: {len(multi_type_contigs)}")
    
    if len(multi_type_contigs) > 0:
        logging.info("Examples of variant type combinations:")
        for _, row in multi_type_contigs.head(3).iterrows():
            logging.info(f"  Contig {row['contig_id']}: Primary={row['variant_type']}, "
                        f"Others={row['other_variant_type']}")
            
            # Show scoring breakdown for the first few examples
            if _ < 2:  # Only show first 2 examples to avoid log spam
                breakdown = get_detailed_scoring_breakdown(row)
                logging.info(f"    Scoring breakdown: {breakdown}")
    
    # Check for any potential issues
    logging.info("=== SCORING SYSTEM VALIDATION COMPLETE ===")

def get_detailed_scoring_breakdown(variant_row):
    """
    Get detailed breakdown of scoring for a variant (matches new integrated logic)
    """
    breakdown = {}
    
    # --- 1. SETUP ---
    vtype = variant_row.get('variant_type', '')
    try:
        size = abs(int(variant_row.get('varsize', 0)))
        reads = int(variant_row.get('num_reads_case', 0))
        vaf = float(variant_row.get('VAF_sum_WT_TPM', 0))
        contig_varsize = abs(int(variant_row.get('contig_varsize', 0)))
    except:
        size, reads, vaf, contig_varsize = 0, 0, 0.0, 0

    feat1 = variant_row.get('site1_feature', 'Intergenic')
    feat2 = variant_row.get('site2_feature', 'Intergenic')
    is_coding = (feat1 == 'CDS') or (feat2 == 'CDS')

    # --- 2. BASE SCORE ---
    base_scores = {
        'FUS': 100,  'IGR': 90, 'DEL': 80, 'INS': 80,
        'NEJ': 60,   'AS': 60,  'PNJ': 60, 'EE': 40,
        'RI': 40,    'NE': 40,  'UN': 10
    }
    breakdown['base_score'] = base_scores.get(vtype, 20)
    
    # --- 3. STRUCTURAL SCORE ---
    # We calculate the net points gained/lost from structural checks
    struct_score = 0
    
    # Helper to mirror main function logic
    def check_splicing_score(row, check_exon=True, expect_exon=True):
        s = 0
        s += 15 if bool(row.get('is_contig_spliced', False)) else -15
        s += 15 if bool(row.get('overlaps_gene', False)) else -15
        if check_exon:
            overlaps = bool(row.get('overlaps_exon', False))
            s += 15 if (overlaps == expect_exon) else -15
        return s

    if vtype in ['INS', 'DEL']:
        struct_score += 15 if size > constants.DEFAULT_MIN_GAP else -15 # Approx constants if global not available
        struct_score += check_splicing_score(variant_row, True, True)
        
    elif vtype in ['FUS', 'IGR']:
        cigar = variant_row.get('cigar', [])
        has_clipping = False
        if isinstance(cigar, list):
             has_clipping = any([op in [4, 5] and val >= constants.DEFAULT_MIN_CLIP for op, val in cigar])
        
        struct_score += 15 if size > constants.DEFAULT_MIN_CLIP else -15
        struct_score += 15 if has_clipping else -15
        struct_score += 15 if bool(variant_row.get('overlaps_gene', False)) else -15
        
        # Translocation Bonus/Penalty
        chr1 = str(variant_row.get('chr1', ''))
        chr2 = str(variant_row.get('chr2', ''))
        if chr1 != chr2 and chr1 and chr2:
            breakdown['translocation_bonus'] = 40
            struct_score += 40
        elif chr1 == chr2:
            breakdown['read_through_penalty'] = -20
            struct_score -= 20

    elif vtype in ['EE', 'NE', 'AS', 'NEJ', 'PNJ', 'RI']:
        expect_overlap = False if vtype in ['NE', 'RI'] else True
        struct_score += check_splicing_score(variant_row, True, expect_overlap)
        
        if vtype in ['NE', 'EE', 'NEJ', 'PNJ']:
            struct_score += 15 if bool(variant_row.get('spliced_exon', False)) else -15
            
        valid_motif = variant_row.get('valid_motif', None)
        if valid_motif is True: struct_score += 15
        elif valid_motif is False: struct_score -= 30
        
        struct_score += 15 if size > constants.DEFAULT_MIN_CLIP else -15 # Simplified check

    elif vtype == 'UN':
         struct_score += 15 if size > constants.DEFAULT_MIN_CLIP else -15
         struct_score += 15 if bool(variant_row.get('overlaps_gene', False)) else -15
         # Soft clip check skipped for brevity in breakdown, assumming +/- 15 balance

    breakdown['structural_score'] = struct_score

    # --- 4. FUNCTIONAL SCORE ---
    func_score = 0
    
    # Frameshift
    if is_coding and vtype in ['DEL', 'INS']:
        if size % 3 != 0:
            func_score += 50
            breakdown['frameshift'] = 'Yes (+50)'
        else:
            func_score += 20
            breakdown['frameshift'] = 'In-frame (+20)'
            
    # RI in CDS
    if vtype == 'RI' and is_coding:
        func_score += 40
        breakdown['coding_retention'] = 'Yes (+40)'

    # Location
    feats = [feat1, feat2]
    if 'CDS' in feats: func_score += 20
    elif 'UTR' in feats: func_score += 10
    elif 'Intron' in feats: func_score -= 10
    elif 'Intergenic' in feats: func_score -= 20
    
    breakdown['functional_score'] = func_score
    breakdown['location_context'] = f"{feat1}/{feat2}"

    # --- 5. EVIDENCE SCORE ---
    ev_score = 0
    if reads < 3: ev_score -= 50
    elif reads >= 20: ev_score += 20
    
    if vaf < 0.05: ev_score -= 30
    elif vaf > 0.2: ev_score += 15
    
    if contig_varsize >= 50: ev_score += 15
    elif contig_varsize >= 20: ev_score += 10
    
    breakdown['evidence_score'] = ev_score

    # --- TOTAL ---
    breakdown['total_score'] = max(0, breakdown['base_score'] + struct_score + func_score + ev_score)
    
    return breakdown

def check_scoring_consistency(contigs_df):
    """
    Check if the scoring system is working consistently
    
    Args:
        contigs_df: DataFrame after variant selection
        
    Returns:
        Boolean indicating if scoring is consistent
    """
    logging.info("=== SCORING CONSISTENCY CHECK ===")
    
    # Note: Since variant_of_interest priority has been removed,
    # all variants are now scored purely based on biological characteristics
    # This function now serves as a placeholder for future consistency checks
    
    logging.info("Scoring system now uses pure biological criteria scoring")
    logging.info("All variants scored equally based on their characteristics")
    
    return True

def parse_args(args):
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Make supertranscript reference'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='sample',
                        metavar='SAMPLE',
                        type=str,
                        help='''Sample name.''')
    parser.add_argument(dest='contig_info',
                        metavar='CONTIG_INFO',
                        type=str,
                        help='''Contig information for novel contigs.''')
    parser.add_argument(dest='contig_fasta',
                        metavar='CONTIG_FASTA',
                        type=str,
                        help='''Fasta file containing contig sequences.''')
    parser.add_argument('--cosmic_tier_data',
                        metavar='COSMIC_TIER',
                        type=str,
                        default=None,
                        help='''Path to COSMIC Cancer Gene Census TSV (GENE_SYMBOL & TIER).''')
    parser.add_argument('--run_de',
                        type=str,
                        default='true',
                        help='''Whether differential expression was run (true/false).''')
    parser.add_argument('--single_sample_min_vaf',
                        type=float,
                        default=0.1,
                        help='''Minimum VAF_sum_WT_TPM to keep a variant when RUN_DE is false '''
                             '''(default: 0.1).''')
    parser.add_argument(dest='de_results',
                        metavar='DE_RESULTS',
                        type=str,
                        help='''Differential expression results.''')
    parser.add_argument(dest='vaf_estimates',
                        metavar='VAF_ESTIMATES',
                        type=str,
                        help='''VAF estimates file.''')
    parser.add_argument('--gene_filter',
                        metavar='GENE_FILTER',
                        type=str,
                        default='',
                        help='''File containing list of genes (one per line) to keep (filter out others).''')
    parser.add_argument('--var_filter',
                        metavar='VAR_FILTER',
                        type=str,
                        nargs='+',
                        help='''Types of variant to keep.''')
    parser.add_argument('--minClip',
                        type=int,
                        help='''Minimum length for a hard/soft clip (default: %(default)s).''')
    parser.add_argument('--minGap',
                        type=int,
                        help='''Minimum gap size for size-dependent scoring (default: %(default)s).''')
    parser.add_argument('--all_variants_out',
                        metavar='ALL_VARIANTS_TSV',
                        type=str,
                        default='',
                        help='''Optional path to write all variants with per-contig ranking (TSV).''')
    parser.add_argument('--discard_out',
                        metavar='DISCARD_TSV',
                        type=str,
                        default='',
                        help='''Optional path to write contigs discarded due to polyA/T or dinucleotide repeats in seq1/seq2.''')
    parser.add_argument('--detect_viral_integration',
                        type=str,
                        default='false',
                        help='''Whether to keep viral integration variants (true/false). Default: false.''')
    return parser.parse_args(args)

def get_all_genes(overlapping_genes):
    if isinstance(overlapping_genes, str):
        genes = overlapping_genes.split(':')
        genes = [gene.split('|') for gene in genes]
        genes = [g for gene in genes for g in gene if g != '']
        return genes
    else:
        return []

def filter_by_gene(contigs, gene_filter):
    genelist = gene_filter[0].values
    overlapping_genes = contigs.overlapping_genes.apply([lambda og: get_all_genes(og)])
    overlapping_genes = overlapping_genes.apply([lambda og: len(np.intersect1d(np.array(og), genelist)) > 0])
    contigs = contigs[overlapping_genes.values]
    return contigs

def add_de_info(contigs, de_results):
    de_results = de_results.rename(columns={'transcript_id': 'contig_id'})
    contigs = pd.merge(contigs, de_results, on='contig_id')
    # Trim leading/trailing spaces off every column name:
    contigs.columns = contigs.columns.str.strip()
    
    logging.info(
        'DE analysis results: %d contigs with DE info, %d total contigs',
        len(contigs.dropna(subset=['logFC'])), len(contigs)
    )
    
    return contigs

def load_cosmic_tiers(filepath):
    """
    Load COSMIC tier and fusion information.
    Expected columns: GENE_SYMBOL, TIER, COSMIC_fusion
    Returns: (tier_map, fusion_map)
    """
    if not filepath or not os.path.exists(filepath):
        return {}, {}
    
    try:
        df = pd.read_csv(filepath, sep='\t')
        # Normalize columns just in case
        df.columns = df.columns.str.upper().str.strip()
        
        if 'GENE_SYMBOL' not in df.columns or 'TIER' not in df.columns:
            logging.warning("COSMIC file missing 'GENE_SYMBOL' or 'TIER' columns.")
            return {}, {}

        # Normalize TIER values to canonical strings (e.g. '1', '2') so the
        # downstream `'1' in tiers_found` / `'2' in tiers_found` comparisons
        # in get_contig_tier() work. If TIER has any NaN values pandas infers
        # the column as float64, which would otherwise yield '1.0' / '2.0'.
        def _normalize_tier(v):
            if pd.isna(v):
                return None
            try:
                return str(int(float(v)))
            except (ValueError, TypeError):
                return str(v).strip()

        tier_series = df['TIER'].apply(_normalize_tier)
        tier_mask = tier_series.notna()
        tier_map = dict(zip(df.loc[tier_mask, 'GENE_SYMBOL'], tier_series[tier_mask]))

        fusion_map = {}
        if 'COSMIC_FUSION' in df.columns:
            fusion_map = dict(zip(df['GENE_SYMBOL'], df['COSMIC_FUSION']))

        return tier_map, fusion_map
    except Exception as e:
        logging.warning(f"Failed to load COSMIC file: {e}")
        return {}, {}

def get_contig_tier(overlapping_genes_str, tier_map):
    """
    Determines the highest COSMIC tier for a set of overlapping genes.
    Tier 1 > Tier 2 > None
    """
    if not tier_map:
        return ""
        
    genes = get_all_genes(overlapping_genes_str)
    tiers_found = []
    
    for gene in genes:
        # Check exact match or clean gene name
        if gene in tier_map:
            tiers_found.append(str(tier_map[gene]))
            
    # Prioritize Tier 1 over Tier 2
    if '1' in tiers_found:
        return '1'
    elif '2' in tiers_found:
        return '2'
    
    return ""

def add_cosmic_info(contigs, cosmic_filepath):
    """
    Adds COSMIC_tier and COSMIC_fusion columns to contigs DataFrame.
    """
    tier_map, fusion_map = load_cosmic_tiers(cosmic_filepath)
    
    if not tier_map:
        logging.info("No COSMIC tier information loaded.")
        contigs['COSMIC_tier'] = ''
    else:
        logging.info(f"Annotating with {len(tier_map)} COSMIC genes.")
        contigs['COSMIC_tier'] = contigs['overlapping_genes'].apply(lambda x: get_contig_tier(x, tier_map))
    
    # Default fusion annotation if column missing or no data
    if not fusion_map:
        contigs['COSMIC_fusion'] = ''
    else:
        contigs['COSMIC_fusion'] = contigs['overlapping_genes'].apply(
            lambda x: 'Yes' if any(g in fusion_map and str(fusion_map[g]).strip().upper() == 'YES' for g in get_all_genes(x)) else 'No'
        )
    return contigs

def get_short_gene_name(gene_string):
    '''
    Extract gene names from the overlapping_genes string
    '''
    if pd.isna(gene_string) or gene_string == '':
        return 'NA'
    
    # Split by ':' and '|' to get individual gene names
    genes = []
    for gene_group in str(gene_string).split(':'):
        genes.extend(gene_group.split('|'))
    
    # Filter out empty strings and take first few genes
    valid_genes = [g for g in genes if g and g != '']
    if not valid_genes:
        return 'NA'
    
    # Return first 3 genes joined by '|'
    return '|'.join(valid_genes[:3])

def add_other_variant_types_by_contig(contigs):
    """
    Consolidate multiple variant annotations per contig using improved scoring system
    
    This function replaces the old priority-based selection with a score-based approach
    that considers biological significance, variant size, and evidence strength.
    """
    logging.info("Starting improved variant type consolidation using scoring system...")
    
    # FIX: Handle empty input gracefully to avoid KeyError later
    if contigs.empty:
        logging.info("Input contigs dataframe is empty. Skipping consolidation.")
        contigs = contigs.copy()
        contigs['other_variant_type'] = ''
        return contigs

    original_count = len(contigs)
    
    # Apply the improved variant selection per contig
    consolidated_contigs = select_best_variant_per_contig(contigs)
    
    final_count = len(consolidated_contigs)
    
    # Count contigs with multiple variant types
    multiple_types_count = sum(consolidated_contigs['other_variant_type'] != '')
    
    logging.info(f"Improved variant consolidation summary:")
    logging.info(f"  Original variants: {original_count}")
    logging.info(f"  Contigs with multiple variant types: {multiple_types_count}")
    logging.info(f"  Final consolidated variants: {final_count}")
    logging.info(f"  Reduced {original_count - final_count} redundant variants")
    
    # Show examples of contigs with multiple variant types
    multi_type_contigs = consolidated_contigs[consolidated_contigs['other_variant_type'] != '']
    if len(multi_type_contigs) > 0:
        type_combinations = multi_type_contigs.groupby(['variant_type', 'other_variant_type']).size().reset_index(name='count')
        logging.info(f"  Examples of consolidated variant type combinations:")
        for _, row in type_combinations.head(5).iterrows():
            logging.info(f"    Primary: {row['variant_type']}, Other: {row['other_variant_type']} ({row['count']} contigs)")
    
    # Show size distribution of selected variants
    size_dist = consolidated_contigs.groupby('variant_type')['varsize'].agg(['count', 'mean', 'min', 'max'])
    logging.info("  Size distribution of selected variants:")
    for vtype, stats in size_dist.iterrows():
        if stats['count'] > 0:
            logging.info(f"    {vtype}: {int(stats['count'])} variants, mean size: {stats['mean']:.1f}bp")
    
    return consolidated_contigs

def _build_repeat_regexes():
    '''
    Build compiled regexes for homopolymer polyA/T and dinucleotide repeats.
    '''
    # PolyA/T of length >= 10
    poly_regex = re.compile(r'(?:A{10,}|T{10,})', re.IGNORECASE)
    # All dinucleotides with two different bases (e.g., AT, GC, TC, etc.), repeated >= 10 times
    dinucs = [x + y for x in 'ACGT' for y in 'ACGT' if x != y]
    dinuc_patterns = [f'(?:{d})' + '{10,}' for d in dinucs]
    dinuc_regex = re.compile('(?:' + '|'.join(dinuc_patterns) + ')', re.IGNORECASE)
    return poly_regex, dinuc_regex

def _is_repetitive_sequence(seq, poly_regex, dinuc_regex):
    '''
    Return True if sequence contains polyA/T or dinucleotide repeats.
    '''
    if not isinstance(seq, str) or seq == '':
        return False
    return bool(poly_regex.search(seq) or dinuc_regex.search(seq))

def ensure_seq_columns(df):
    '''
    Ensure sequence columns exist so final outputs always include them.
    '''
    for col in ['seq_loc1', 'seq_loc2', 'seq1', 'seq2']:
        if col not in df.columns:
            df[col] = ''
    return df


def prepare_vaf_merge_table(vafs):
    '''
    Dedupe VAF estimates and keep all merge columns present in the file.
    '''
    if vafs is None or vafs.empty:
        return vafs
    if 'contig_id' not in vafs.columns:
        logging.warning('VAF file has no contig_id column; using full table as-is.')
        return vafs.drop_duplicates()
    use_cols = ['contig_id'] + [
        c for c in VAF_ESTIMATE_MERGE_COLUMNS
        if c in vafs.columns and c != 'contig_id'
    ]
    return vafs[use_cols].drop_duplicates()


def merge_vaf_into_contigs(contigs, vafs):
    '''
    Left-merge VAF columns on contig_id. For columns present in both tables
    (e.g. num_reads_case, which exists on contigs from the DE/Oarfish path and in the
    VAF file), prefer the VAF value but fall back to the existing contigs value where the
    VAF table has no row for that contig. This keeps Oarfish read counts for variants that
    are not variant_of_interest and are therefore absent from vaf_estimates.txt (which
    would otherwise leave num_reads_case as NaN/NA in the final results).
    '''
    if vafs is None or vafs.empty or 'contig_id' not in vafs.columns:
        return contigs
    overlap = [c for c in vafs.columns if c != 'contig_id' and c in contigs.columns]
    fallback = None
    if overlap:
        # These columns are constant within a contig; dedupe so the backfill index is unique.
        fallback = contigs[['contig_id'] + overlap].drop_duplicates('contig_id').set_index('contig_id')
        contigs = contigs.drop(columns=overlap, errors='ignore')
    merged = contigs.merge(vafs, on='contig_id', how='left')
    if overlap:
        for col in overlap:
            backfill = merged['contig_id'].map(fallback[col])
            merged[col] = merged[col].where(merged[col].notna(), backfill)
    return merged

def filter_variants_with_chr(df, label, detect_viral_integration=False):
    '''
    Keep variants where chr1 or chr2 contains lowercase "chr".
    If detect_viral_integration is False, drop cases where only one side has "chr".
    '''
    if df.empty:
        return df
    chr1_has_chr = df['chr1'].astype(str).str.contains('chr', case=True, na=False)
    chr2_has_chr = df['chr2'].astype(str).str.contains('chr', case=True, na=False)
    any_chr = chr1_has_chr | chr2_has_chr
    filtered_out = (~any_chr).sum()
    if filtered_out > 0:
        logging.info("Filtered %d %s variants without 'chr' in chr1/chr2.", filtered_out, label)
    df = df[any_chr]
    if not detect_viral_integration:
        viral_mask = chr1_has_chr ^ chr2_has_chr
        viral_filtered = viral_mask.sum()
        if viral_filtered > 0:
            logging.info("Filtered %d %s viral-integration variants.", viral_filtered, label)
        df = df[~viral_mask]
    return df

def filter_contigs_by_sequence_repeats(contigs, contig_fasta):
    '''
    Extract seq1/seq2 using get_variant_seq, then split contigs into keep/discard
    based on presence of polyA/T or dinucleotide repeats in either sequence.
    Returns (kept_df, discarded_df). Writing is handled later after full formatting
    to ensure identical columns as final output.
    '''
    logging.info('Running sequence-based filtering for polyA/T and dinucleotide repeats...')
    # Ensure sequences are available
    contigs_with_seq = get_variant_seq(contigs.copy(), contig_fasta)

    poly_regex, dinuc_regex = _build_repeat_regexes()
    has_repeat_seq1 = contigs_with_seq['seq1'].apply(lambda s: _is_repetitive_sequence(s, poly_regex, dinuc_regex)) if 'seq1' in contigs_with_seq.columns else pd.Series(False, index=contigs_with_seq.index)
    has_repeat_seq2 = contigs_with_seq['seq2'].apply(lambda s: _is_repetitive_sequence(s, poly_regex, dinuc_regex)) if 'seq2' in contigs_with_seq.columns else pd.Series(False, index=contigs_with_seq.index)
    to_discard_mask = (has_repeat_seq1 | has_repeat_seq2).fillna(False)

    discarded = contigs_with_seq[to_discard_mask]
    kept = contigs_with_seq[~to_discard_mask]

    logging.info('Sequence filtering: %d discarded, %d kept (total %d).', len(discarded), len(kept), len(contigs))

    return kept, discarded

def get_variant_seq(contigs, contig_fasta):
    '''
    Extract variant sequences from contig fasta
    '''
    logging.info("Extracting variant sequences...")
    
    start = round(VAR_SEQ_LEN / 2)
    end = VAR_SEQ_LEN - start
    contig_info = []
    var_bed = []

    logging.info(f"Number of contigs before extraction: {len(contigs)}")

    for idx, loc in contigs.iterrows():
        # Window around cpos
        pos1 = max(0, int(loc.cpos) - start)
        pos2 = min(int(loc.contig_len), int(loc.cpos) + end)
        # Expand/adjust window to desired length without exceeding bounds
        if pos2 - pos1 < VAR_SEQ_LEN:
            # Prefer extend to the right if at start, else to the left if at end
            if pos1 == 0:
                pos2 = min(int(loc.contig_len), VAR_SEQ_LEN)
            elif pos2 == int(loc.contig_len):
                pos1 = max(0, int(loc.contig_len) - VAR_SEQ_LEN)
            else:
                # Centered window fallback
                mid = int(loc.cpos)
                pos1 = max(0, mid - start)
                pos2 = min(int(loc.contig_len), mid + end)
        var_bed.append([loc.contig_id, pos1, pos2])

        loc2_present = False
        if loc.contig_varsize > 0 and loc.variant_type != 'UN':
            var_end = int(loc.cpos) + int(loc.contig_varsize)
            pos3 = max(0, var_end - start)
            pos4 = min(int(loc.contig_len), var_end + end)
            if pos4 - pos3 < VAR_SEQ_LEN:
                if pos3 == 0:
                    pos4 = min(int(loc.contig_len), VAR_SEQ_LEN)
                elif pos4 == int(loc.contig_len):
                    pos3 = max(0, int(loc.contig_len) - VAR_SEQ_LEN)
                else:
                    mid2 = var_end
                    pos3 = max(0, mid2 - start)
                    pos4 = min(int(loc.contig_len), mid2 + end)
            var_bed.append([loc.contig_id, pos3, pos4])
            loc2_present = True

        seq1 = '%s:%d-%d' % (loc.contig_id, pos1, pos2)
        seq2 = '%s:%d-%d' % (loc.contig_id, pos3, pos4) if loc2_present else ''
        contig_info.append([loc.variant_id, seq1, seq2])

    var_bed = pd.DataFrame(var_bed, columns=['contig', 'start', 'end']).drop_duplicates()
    g = BedTool.from_dataframe(var_bed).remove_invalid()
    g = g.sequence(fi=contig_fasta)
    vs = get_block_seqs(g)

    logging.info(f"Retrieved {len(vs)} sequences from FASTA")

    var_seqs = pd.DataFrame.from_dict(vs, columns=['seq'], orient='index')

    contig_info = pd.DataFrame(contig_info, columns=['variant_id', 'seq_loc1', 'seq_loc2'])
    logging.info(f"Number of unique variant_id in contig_info: {contig_info['variant_id'].nunique()}")

    logging.info(f"Merging contig_info ({len(contig_info)}) with sequences")

    # Merge sequences with contig_info
    contig_info = contig_info.merge(var_seqs, left_on='seq_loc1', right_index=True, how='left')
    contig_info = contig_info.merge(var_seqs, left_on='seq_loc2', right_index=True, how='left', suffixes=['1', '2'])

    logging.info(f"Final contig_info size: {contig_info.shape}")

    # Sanity logging: how many sequences matched
    matched_seq1 = contig_info['seq1'].notna().sum() if 'seq1' in contig_info.columns else 0
    matched_seq2 = contig_info['seq2'].notna().sum() if 'seq2' in contig_info.columns else 0
    logging.info(f"Sequence merge matches: seq1={matched_seq1}, seq2={matched_seq2}")

    # Check for duplicates and aggregate if necessary
    dup_check = contig_info.variant_id.value_counts()
    if dup_check.max() > 10:
        logging.warning(f"High duplication detected in variant_id! Max duplicates: {dup_check.max()}")
        contig_info = contig_info.groupby('variant_id').agg({
            'seq_loc1': lambda x: ','.join(x.dropna().unique()),
            'seq_loc2': lambda x: ','.join(x.dropna().unique()),
            'seq1': lambda x: ','.join(x.dropna().unique()),
            'seq2': lambda x: ','.join(x.dropna().unique())
        }).reset_index()
        logging.info(f"Aggregated contig_info size: {contig_info.shape}")

    logging.info("Merging sequence info with original contigs...")
    contigs = contigs.merge(contig_info, on='variant_id', how='left')
    logging.info(f"Contigs shape after merge: {contigs.shape}")
    return contigs

def reformat_fields(contigs):
    '''
    Extract chrom, pos and strand fields.
    Reorder fields for clarity.
    Sort by p value.
    '''
    pos1 = contigs.pos1.apply(get_pos_parts).values
    pos2 = contigs.pos2.apply(get_pos_parts).values
    chr1, pos1, str1 = zip(*pos1)
    chr2, pos2, str2 = zip(*pos2)
    contigs['chr1'], contigs['pos1'], contigs['strand1'] = chr1, pos1, str1
    contigs['chr2'], contigs['pos2'], contigs['strand2'] = chr2, pos2, str2
    ran_de = 'logFC' in contigs.columns.values

    # Final column order when DE columns are available
    final_cols_de = [
        'chr1', 'pos1', 'strand1', 'site1_feature',
        'chr2', 'pos2', 'strand2', 'site2_feature',
        'variant_type', 'other_variant_type', 'overlapping_genes', 'sample', 'varsize',
        'supporting_read_count', 'supporting_read_count_reliability', 'other_supporting_read_count',
        'spanning_reads', 'junction_VAF', 'num_reads_case', 'total_num_reads_controls',
        'VAF_sum_WT_TPM',
        'logFC', 'FDR', 'PValue',
        'case_gene_count_wt', 'control_gene_count_wt', 'control_contig_count', 'control_gene_count_total',
        'TPM', 'sum_WT_TPM', 'fisher_p_val', 'odds_ratio',
        'large_varsize', 'is_contig_spliced', 'spliced_exon', 'overlaps_exon', 'overlaps_gene',
        'motif', 'valid_motif', 'COSMIC_tier', 'COSMIC_fusion',
        'vars_in_contig', 'contig_id', 'variant_id', 'partner_id', 'contig_varsize', 'cpos',
        'unique_contig_ID', 'contig_len', 'contig_cigar', 'seq_loc1', 'seq_loc2', 'seq1', 'seq2', 'variant_score'
    ]

    # Final column order when DE columns are not available (same as final_cols_de minus logFC/FDR/PValue)
    final_cols_node = [
        'chr1', 'pos1', 'strand1', 'site1_feature',
        'chr2', 'pos2', 'strand2', 'site2_feature',
        'variant_type', 'other_variant_type', 'overlapping_genes', 'sample', 'varsize',
        'supporting_read_count', 'supporting_read_count_reliability', 'other_supporting_read_count',
        'spanning_reads', 'junction_VAF', 'num_reads_case', 'total_num_reads_controls',
        'VAF_sum_WT_TPM',
        'case_gene_count_wt', 'control_gene_count_wt', 'control_contig_count', 'control_gene_count_total',
        'TPM', 'sum_WT_TPM', 'fisher_p_val', 'odds_ratio',
        'large_varsize', 'is_contig_spliced', 'spliced_exon', 'overlaps_exon', 'overlaps_gene',
        'motif', 'valid_motif', 'COSMIC_tier', 'COSMIC_fusion',
        'vars_in_contig', 'contig_id', 'variant_id', 'partner_id', 'contig_varsize', 'cpos',
        'unique_contig_ID', 'contig_len', 'contig_cigar', 'seq_loc1', 'seq_loc2', 'seq1', 'seq2', 'variant_score'
    ]

    # Sort by PValue (ascending) if available, otherwise by VAF_sum_WT_TPM (descending)
    if 'PValue' in contigs.columns:
        contigs['PValue'] = pd.to_numeric(contigs['PValue'], errors='coerce')
        contigs = contigs.sort_values(by='PValue', ascending=True, na_position='last')
    elif 'VAF_sum_WT_TPM' in contigs.columns:
        contigs['VAF_sum_WT_TPM'] = pd.to_numeric(contigs['VAF_sum_WT_TPM'], errors='coerce')
        contigs = contigs.sort_values('VAF_sum_WT_TPM', ascending=False, na_position='last')

    # Enforce full column set: missing columns (e.g. no DE / older VAF file) become NA
    final_cols = final_cols_de if ran_de else final_cols_node
    for c in final_cols:
        if c not in contigs.columns:
            contigs[c] = pd.NA

    # Read-count columns: emit as clean integers (no decimals). The supporting-read
    # columns come from a left merge, so they are float-with-NaN until cast to Int64.
    for _c in ('num_reads_case', 'case_gene_count_wt', 'supporting_read_count',
               'spanning_reads'):
        if _c in contigs.columns:
            contigs[_c] = pd.to_numeric(contigs[_c], errors='coerce').round().astype('Int64')

    contigs = contigs[final_cols]

    return contigs

def main():
    args = parse_args(sys.argv[1:])
    init_logging(args.log)
    
    # Set global variables for minimum thresholds
    set_globals(args)
    
    logging.info("Starting post-processing script with improved variant selection...")
    logging.info(f"Using MIN_CLIP={MIN_CLIP}, MIN_GAP={MIN_GAP}")

    try:
        logging.info("Loading input files...")
        contigs = pd.read_csv(args.contig_info, sep='\t', low_memory=False).fillna('')
        logging.info(f"Loaded contig information from {args.contig_info}.")
        
        # Log warning if input file itself is empty
        if contigs.empty:
            logging.warning("Input contig information file is empty.")
        
        de_results = pd.read_csv(args.de_results, sep='\t', low_memory=False)
        logging.info(f"Loaded differential expression results from {args.de_results}.")

        vafs = pd.read_csv(args.vaf_estimates, sep='\t', low_memory=False)
        logging.info(f"Loaded VAF estimates from {args.vaf_estimates}.")
        vafs = prepare_vaf_merge_table(vafs)
        logging.info(f"Processed VAF data: {len(vafs)} unique records.")

        # Load gene filter if provided
        gene_filter = []
        if args.gene_filter and os.path.exists(args.gene_filter):
            gene_filter = pd.read_csv(args.gene_filter, header=None)
            logging.info(f"Loaded gene filter with {len(gene_filter)} genes.")
        else:
            logging.info("No gene filter provided or file not found.")

    except Exception as e:
        logging.error(f"Error loading input files: {str(e)}")
        sys.exit(EXIT_FILE_IO_ERROR)

    logging.info(f"Initial contigs loaded: {len(contigs)}")
    logging.info(f"Unique contigs: {contigs['contig_id'].nunique()}")
    
    # STEP 1: Sequence-based filtering BEFORE scoring
    # Extract sequences, discard contigs with polyA/T or dinucleotide repeats in seq1/seq2
    contigs['sample'] = args.sample
    kept_contigs, discarded_contigs = filter_contigs_by_sequence_repeats(contigs, args.contig_fasta)
    # Compute variants-per-contig from the original (pre-filter) contigs
    vars_per_contig = contigs.groupby('contig_id', as_index=False).agg(
        {'variant_id': lambda x: len(np.unique(x))}
    ).rename({'variant_id': 'vars_in_contig'}, axis=1)
    # Attach vars_in_contig to filtered subsets
    kept_contigs = kept_contigs.merge(vars_per_contig, on='contig_id', how='left')
    discarded_contigs = discarded_contigs.merge(vars_per_contig, on='contig_id', how='left')
    # Ensure sequence columns are present and preserved going forward
    kept_contigs = ensure_seq_columns(kept_contigs)
    discarded_contigs = ensure_seq_columns(discarded_contigs)

    # Keep a copy of kept variants (pre-consolidation) for optional ranked export later
    all_variants = kept_contigs.copy()

    # STEP 2: Apply improved variant consolidation using scoring system on remaining contigs
    logging.info("Processing %d filtered variants using improved scoring system", len(kept_contigs))
    contigs = add_other_variant_types_by_contig(kept_contigs)
    logging.info(f"After variant consolidation: {len(contigs)} variants")

    logging.info("Added variants per contig information.")

    # STEP 3: Apply variant type filtering if specified
    if args.var_filter:
        contigs = contigs[contigs.variant_type.apply(lambda v: v in args.var_filter).values]
        logging.info(f"Filtered contigs by variant type: {args.var_filter}. Remaining: {len(contigs)}")

    # STEP 4: Apply gene filtering if specified
    if len(gene_filter) > 0:
        contigs = filter_by_gene(contigs, gene_filter)
        logging.info(f"Filtered contigs by gene list. Remaining: {len(contigs)}")

    if len(contigs) == 0:
        logging.warning('No variants found after filtering. Generating empty output and exiting.')
        contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')
        logging.info("Post-processing completed successfully (no variants found).")
        sys.exit()

    # STEP 5: Add DE and VAF information
    logging.info('Adding DE and VAF info...')
    contigs = add_de_info(contigs, de_results)
    contigs = merge_vaf_into_contigs(contigs, vafs)

    # Score only once num_reads_case (from add_de_info) and VAF_sum_WT_TPM (from the VAF merge)
    # are present. Scored any earlier, both .get() calls fall back to their 0 default and the
    # evidence block applies a flat -80 to every row (reads < 3 and vaf < 0.05), carrying no
    # information. varsize/contig_varsize come from the contig info file and are fine either way.
    #
    # select_best_variant_per_contig still scores pre-merge, because it has to pick a surviving
    # variant before the per-contig collapse. That is harmless: both columns are per-contig, so
    # the -80 lands equally on every candidate and cancels out of a within-contig ranking.
    contigs['variant_score'] = contigs.apply(calculate_variant_score, axis=1)
    validate_scoring_system(contigs)

    # Add COSMIC Annotation
    contigs = add_cosmic_info(contigs, args.cosmic_tier_data)

    # Optional: Log how many hits we found
    hits = contigs[contigs['COSMIC_tier'] != ''].shape[0]
    logging.info(f"Found COSMIC annotations for {hits} contigs.")

    logging.info("DE, VAF and COSMIC information added.")

    # ---------------------------------------------------------
    # NEW: Conditional Filtering for Single Sample Mode
    # ---------------------------------------------------------
    if str(args.run_de).lower() == 'false':
        logging.info("Single Sample Mode detected (RUN_DE=false). Applying strict COSMIC & VAF filtering.")
        
        initial_count = len(contigs)
        
        # DEBUG: Log what values we are seeing before filtering
        logging.info(f"DEBUG: Unique COSMIC Tiers found: {contigs['COSMIC_tier'].unique()}")
        logging.info(f"DEBUG: VAF_sum_WT_TPM stats: Min={contigs['VAF_sum_WT_TPM'].min()}, "
                     f"Max={contigs['VAF_sum_WT_TPM'].max()}")

        # 1. Clean VAF
        contigs['VAF_sum_WT_TPM'] = pd.to_numeric(contigs['VAF_sum_WT_TPM'], errors='coerce').fillna(0)

        # 2. Robust Tier Check
        # Allows matching '1', 'Tier 1', '1.0', etc.
        def is_tier_valid(val):
            val_str = str(val).upper()
            return '1' in val_str or '2' in val_str

        has_cosmic = contigs['COSMIC_tier'].apply(is_tier_valid)
        high_vaf = contigs['VAF_sum_WT_TPM'] >= args.single_sample_min_vaf

        # 3. Apply Filter
        contigs = contigs[has_cosmic & high_vaf]
        
        filtered_count = len(contigs)
        logging.info(f"Filtered {initial_count - filtered_count} variants. Remaining: {filtered_count}")
        
        if filtered_count == 0:
            logging.warning("No variants passed the Single Sample filters. Check if COSMIC file matched any genes.")
    # ---------------------------------------------------------

    logging.info("DE, VAF and COSMIC information added and filtered.")

    # STEP 6: Generate unique contig IDs; sequences already present from filtering step
    short_gnames = contigs.overlapping_genes.map(str).apply(get_short_gene_name)
    contig_ids, samples = contigs.contig_id, contigs['sample']
    con_names = ['|'.join([s, cid, sg]) for cid, s, sg in zip(contig_ids, samples, short_gnames)]
    contigs['unique_contig_ID'] = con_names
    logging.info("Generated unique contig IDs.")
    # Do not re-extract sequences here to avoid overwriting; keep existing seq columns

    # STEP 7: Final output formatting
    logging.info('Outputting to CSV')
    contigs = ensure_seq_columns(contigs)
    contigs = reformat_fields(contigs)

    # Filter out variants where neither chr1 nor chr2 contains "chr"
    detect_viral = str(args.detect_viral_integration).lower() == 'true'
    contigs = filter_variants_with_chr(contigs, "final", detect_viral_integration=detect_viral)

    # If requested, produce a discard TSV with identical columns to final output
    if getattr(args, 'discard_out', '') and not discarded_contigs.empty:
        try:
            logging.info('Preparing discard output with final columns...')
            # Add DE/VAF info
            discarded_full = add_de_info(discarded_contigs.copy(), de_results)
            discarded_full = merge_vaf_into_contigs(discarded_full, vafs)
            discarded_full = add_cosmic_info(discarded_full, args.cosmic_tier_data)
            # other_variant_type may not exist; ensure it exists for reformat
            if 'other_variant_type' not in discarded_full.columns:
                discarded_full['other_variant_type'] = ''
            # Score after the merges above, unconditionally: a variant_score carried over from
            # the pre-merge frame would have been computed without VAF_sum_WT_TPM.
            discarded_full['variant_score'] = discarded_full.apply(calculate_variant_score, axis=1)
            # Generate unique contig IDs and sequences
            short_gnames_disc = discarded_full.overlapping_genes.map(str).apply(get_short_gene_name)
            disc_ids, disc_samples = discarded_full.contig_id, discarded_full['sample']
            disc_names = ['|'.join([s, cid, sg]) for cid, s, sg in zip(disc_ids, disc_samples, short_gnames_disc)]
            discarded_full['unique_contig_ID'] = disc_names
            # Sequences already present; ensure columns exist
            discarded_full = ensure_seq_columns(discarded_full)
            # Reformat to final columns and write
            discarded_final = reformat_fields(discarded_full)
            discarded_final = filter_variants_with_chr(
                discarded_final, "discarded", detect_viral_integration=detect_viral
            )
            discarded_final.to_csv(args.discard_out, index=False, sep='\t', na_rep='NA')
            logging.info(f"Wrote discard file to {args.discard_out}")
        except Exception as e:
            logging.warning(f"Failed to write discard output: {str(e)}")

    # If requested, produce an all-variants ranked TSV with the same column order
    if getattr(args, 'all_variants_out', '') and not all_variants.empty:
        try:
            # Add DE/VAF info to all_variants first: rank_variants_per_contig scores each row,
            # and the VAF term of that score needs VAF_sum_WT_TPM to already be merged in.
            ranked_all = add_de_info(all_variants, de_results)
            ranked_all = merge_vaf_into_contigs(ranked_all, vafs)
            # Compute ranking on the filtered set (pre-consolidation copy)
            ranked_all = rank_variants_per_contig(ranked_all)
            ranked_all = add_cosmic_info(ranked_all, args.cosmic_tier_data)
            # Generate unique contig IDs to align with main output style
            short_gnames_all = ranked_all.overlapping_genes.map(str).apply(get_short_gene_name)
            contig_ids_all, samples_all = ranked_all.contig_id, ranked_all['sample']
            con_names_all = ['|'.join([s, cid, sg]) for cid, s, sg in zip(contig_ids_all, samples_all, short_gnames_all)]
            ranked_all['unique_contig_ID'] = con_names_all
            # Sequences already present from filtering step; ensure columns exist
            ranked_all = ensure_seq_columns(ranked_all)
            # Ensure 'other_variant_type' exists even though we didn't consolidate
            if 'other_variant_type' not in ranked_all.columns:
                ranked_all['other_variant_type'] = ''
            # Recompute positional fields for ordering compatibility
            ranked_all = ranked_all.copy()
            pos1_all = ranked_all.pos1.apply(get_pos_parts).values
            pos2_all = ranked_all.pos2.apply(get_pos_parts).values
            chr1_all, pos1_vals_all, str1_all = zip(*pos1_all)
            chr2_all, pos2_vals_all, str2_all = zip(*pos2_all)
            ranked_all['chr1'], ranked_all['pos1'], ranked_all['strand1'] = chr1_all, pos1_vals_all, str1_all
            ranked_all['chr2'], ranked_all['pos2'], ranked_all['strand2'] = chr2_all, pos2_vals_all, str2_all
            ranked_all = filter_variants_with_chr(
                ranked_all, "all-variants ranked", detect_viral_integration=detect_viral
            )
            # Sort before selecting columns, so the fallback key does not have to survive
            # the column selection below (by PValue if available else VAF_sum_WT_TPM desc).
            if 'PValue' in ranked_all.columns:
                ranked_all['PValue'] = pd.to_numeric(ranked_all['PValue'], errors='coerce')
                ranked_all = ranked_all.sort_values(by='PValue', ascending=True, na_position='last')
            elif 'VAF_sum_WT_TPM' in ranked_all.columns:
                ranked_all['VAF_sum_WT_TPM'] = pd.to_numeric(ranked_all['VAF_sum_WT_TPM'], errors='coerce')
                ranked_all = ranked_all.sort_values('VAF_sum_WT_TPM', ascending=False, na_position='last')
            # Compose final ordering to strictly follow the main output's column order
            base_cols = list(contigs.columns)
            # Define ranking columns that should be preserved (avoid duplicates)
            ranking_cols = [
                c for c in ['variant_score', 'rank_within_contig', 'is_primary']
                if c in ranked_all.columns and c not in base_cols
            ]
            # Start with the exact column order from the main output
            final_cols = [c for c in base_cols if c in ranked_all.columns]
            # Add ranking columns at the end
            final_cols.extend(ranking_cols)
            # Select only these columns
            ranked_all = ranked_all[final_cols]
            # Oarfish read-count columns: emit as clean integers (no decimals).
            for _c in ('num_reads_case', 'case_gene_count_wt'):
                if _c in ranked_all.columns:
                    ranked_all[_c] = pd.to_numeric(ranked_all[_c], errors='coerce').round().astype('Int64')
            # Write file
            ranked_all.to_csv(args.all_variants_out, index=False, sep='\t', na_rep='NA')
            logging.info(f"Wrote all-variants ranked table to {args.all_variants_out}")
        except Exception as e:
            logging.warning(f"Failed to write all-variants ranked output: {str(e)}")
    
    # Final summary
    final_variant_counts = contigs['variant_type'].value_counts()
    logging.info("FINAL RESULTS SUMMARY:")
    logging.info(f"Total variants output: {len(contigs)}")
    for vtype, count in final_variant_counts.items():
        logging.info(f"  {vtype}: {count}")
    
    contigs.to_csv(sys.stdout, index=False, sep='\t', na_rep='NA')
    logging.info("Post-processing completed successfully.")

if __name__ == '__main__':
    main()