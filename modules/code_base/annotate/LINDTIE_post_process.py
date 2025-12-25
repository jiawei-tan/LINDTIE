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
    Calculate a score for variant type selection based on multiple criteria
    Higher score = more likely to be the correct primary variant type
    
    Args:
        variant_row: Dictionary or pandas Series with variant information
    """
    vtype = variant_row.get('variant_type', '')
    size = abs(int(variant_row.get('varsize', 0)))
    contig_varsize = abs(int(variant_row.get('contig_varsize', 0)))
    
    # Base scores for different variant types based on biological significance
    base_scores = {
        # High-impact structural variants
        'FUS': 140,   # Fusions - highest priority for cancer/disease
        'IGR': 140,   # Intra-genic rearrangements
        
        # Size-dependent structural variants (score depends on size)
        'DEL': 80,    # Deletions - important, size-dependent
        'INS': 80,    # Insertions - important, size-dependent  
        
        # Splicing variants (moderate impact)
        'AS': 60,     # Alternative splicing
        'NEJ': 60,    # Novel exon junctions
        'PNJ': 60,    # Partial novel junctions
        
        # Novel sequences (lower priority unless large)
        'RI': 40,     # Retained introns
        'NE': 40,     # Novel exons
        'EE': 40,     # Extended exons
        
        # Unknown/unclear (lowest priority)
        'UN': 10      # Unknown
    }
    
    score = base_scores.get(vtype, 0)
    
    # FUSION PRIORITY: Give fusion variants absolute priority
    if vtype in ['FUS', 'IGR']:
        return 1000  # Maximum score to ensure fusion variants always win
    
    # Variant-specific scoring based on expected characteristics from the criteria table
    if vtype == 'INS':  # Insertion
        # Expected: Spliced Contig=Y, Varsize>min_gap, Overlaps Gene=Y, Overlaps Exon=Y
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        
        # Size check: should be > min_gap
        if size > MIN_GAP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
        
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
                
        if overlaps_exon:
            score += 15  # Bonus for expected exon overlap
        else:
            score -= 15  # Penalty for missing exon overlap
    
    elif vtype == 'DEL':  # Deletion
        # Expected: Spliced Contig=Y, Varsize>min_gap, Overlaps Gene=Y, Overlaps Exon=Y
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        
        # Size check: should be > min_gap
        if size > MIN_GAP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
        
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
                
        if overlaps_exon:
            score += 15  # Bonus for expected exon overlap
        else:
            score -= 15  # Penalty for missing exon overlap
    
    elif vtype == 'FUS':  # Fusion
        # Expected: Clipping=hard/soft, Varsize>min_clip, Overlaps Gene=Y
        has_clipping = any([op in [constants.CIGAR['hard-clip'], constants.CIGAR['soft-clip']] 
                           and val >= MIN_CLIP for op, val in variant_row.get('cigar', [])])
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        
        # Size check: should be > min_clip
        if size > MIN_CLIP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if has_clipping:
            score += 15  # Bonus for expected clipping
        else:
            score -= 15  # Penalty for missing clipping
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
    
    elif vtype == 'IGR':  # Intergenic Rearrangement
        # Expected: Clipping=hard/soft, Varsize>min_clip, Overlaps Gene=Y, within same gene
        has_clipping = any([op in [constants.CIGAR['hard-clip'], constants.CIGAR['soft-clip']] 
                           and val >= MIN_CLIP for op, val in variant_row.get('cigar', [])])
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        
        # Size check: should be > min_clip
        if size > MIN_CLIP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if has_clipping:
            score += 15  # Bonus for expected clipping
        else:
            score -= 15  # Penalty for missing clipping
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
    
    elif vtype == 'UN':  # Unknown
        # Expected: Clipping=soft, Varsize>min_clip, Overlaps Gene=Y
        has_soft_clip = any([op == constants.CIGAR['soft-clip'] 
                            and val >= MIN_CLIP for op, val in variant_row.get('cigar', [])])
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        
        # Size check: should be > min_clip
        if size > MIN_CLIP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if has_soft_clip:
            score += 15  # Bonus for expected soft clipping
        else:
            score -= 15  # Penalty for missing soft clipping
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
    
    elif vtype == 'EE':  # Exon-Exon
        # Expected: Spliced Contig=Y, Varsize>min_clip, Overlaps Gene=Y, Overlaps Exon=N, Spliced Exon=Y, Valid Motif=Y
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        spliced_exon = bool(variant_row.get('spliced_exon', False))
        valid_motif = variant_row.get('valid_motif', None)
        
        # Size check: should be > min_clip
        if size > MIN_CLIP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
                
        if not overlaps_exon:  # Should NOT overlap exon
            score += 15  # Bonus for correct exon non-overlap
        else:
            score -= 15  # Penalty for unexpected exon overlap
                
        if spliced_exon:
            score += 15  # Bonus for expected spliced exon
        else:
            score -= 15  # Penalty for missing spliced exon
                
        if valid_motif is True:
            score += 15  # Bonus for expected valid motif
        elif valid_motif is False:
            score -= 15  # Penalty for invalid motif
    
    elif vtype == 'NE':  # Novel Exon
        # Expected: Spliced Contig=Y, Varsize>min_clip, Overlaps Gene=N, Overlaps Exon=N, Spliced Exon=Y
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        spliced_exon = bool(variant_row.get('spliced_exon', False))
        
        # Size check: should be > min_clip
        if size > MIN_CLIP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
                
        if not overlaps_gene:  # Should NOT overlap gene
            score += 15  # Bonus for correct gene non-overlap
        else:
            score -= 15  # Penalty for unexpected gene overlap
                
        if not overlaps_exon:  # Should NOT overlap exon
            score += 15  # Bonus for correct exon non-overlap
        else:
            score -= 15  # Penalty for unexpected exon overlap
                
        if spliced_exon:
            score += 15  # Bonus for expected spliced exon
        else:
            score -= 15  # Penalty for missing spliced exon
    
    elif vtype == 'RI':  # Retained Intron
        # Expected: Spliced Contig=Y, Varsize>min_clip, Overlaps Gene=Y, Overlaps Exon=N
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        
        # Size check: should be > min_clip
        if size > MIN_CLIP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
                
        if not overlaps_exon:  # Should NOT overlap exon
            score += 15  # Bonus for correct exon non-overlap
        else:
            score -= 15  # Penalty for unexpected exon overlap
    
    elif vtype == 'AS':  # Alternative Splicing
        # Expected: Spliced Contig=Y, Overlaps Gene=Y, Overlaps Exon=Y, novel splice between 2 known sites
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
                
        if overlaps_exon:
            score += 15  # Bonus for expected exon overlap
        else:
            score -= 15  # Penalty for missing exon overlap
    
    elif vtype == 'NEJ':  # Novel Exon Junction
        # Expected: Spliced Contig=Y, Varsize>min_gap, Overlaps Gene=Y, Overlaps Exon=Y, Spliced Exon=Y, Valid Motif=Y
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        spliced_exon = bool(variant_row.get('spliced_exon', False))
        valid_motif = variant_row.get('valid_motif', None)
        
        # Size check: should be > min_gap
        if size > MIN_GAP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
                
        if overlaps_exon:
            score += 15  # Bonus for expected exon overlap
        else:
            score -= 15  # Penalty for missing exon overlap
                
        if spliced_exon:
            score += 15  # Bonus for expected spliced exon
        else:
            score -= 15  # Penalty for missing spliced exon
                
        if valid_motif is True:
            score += 15  # Bonus for expected valid motif
        elif valid_motif is False:
            score -= 15  # Penalty for invalid motif
    
    elif vtype == 'PNJ':  # Partial Novel Junction
        # Expected: Spliced Contig=Y, Varsize>min_gap, Overlaps Gene=Y, Overlaps Exon=Y, Spliced Exon=Y, Valid Motif=Y
        is_contig_spliced = bool(variant_row.get('is_contig_spliced', False))
        overlaps_gene = bool(variant_row.get('overlaps_gene', False))
        overlaps_exon = bool(variant_row.get('overlaps_exon', False))
        spliced_exon = bool(variant_row.get('spliced_exon', False))
        valid_motif = variant_row.get('valid_motif', None)
        
        # Size check: should be > min_gap
        if size > MIN_GAP:
            score += 15  # Bonus for meeting size criteria
        else:
            score -= 15  # Penalty for not meeting size criteria
        
        if is_contig_spliced:
            score += 15  # Bonus for expected splicing
        else:
            score -= 15  # Penalty for missing splicing
                
        if overlaps_gene:
            score += 15  # Bonus for expected gene overlap
        else:
            score -= 15  # Penalty for missing gene overlap
                
        if overlaps_exon:
            score += 15  # Bonus for expected exon overlap
        else:
            score -= 15  # Penalty for missing exon overlap
                
        if spliced_exon:
            score += 15  # Bonus for expected spliced exon
        else:
            score -= 15  # Penalty for missing spliced exon
                
        if valid_motif is True:
            score += 15  # Bonus for expected valid motif
        elif valid_motif is False:
            score -= 15  # Penalty for invalid motif
    
    # Contig evidence strength bonus (applies to all variants)
    if contig_varsize > 0:
        # Variants that affect contig sequence get bonus
        if contig_varsize >= 50:
            score += 15
        elif contig_varsize >= 20:
            score += 10
        elif contig_varsize >= 10:
            score += 5
    
    # Additional scoring factors
    # Read support bonus
    num_reads_case = variant_row.get('num_reads_case', 0)
    if num_reads_case > 0:
        if num_reads_case >= 10:
            score += 10  # High read support
        elif num_reads_case >= 5:
            score += 5   # Medium read support
    
    # VAF bonus (if available)
    vaf = variant_row.get('VAF', 0)
    if vaf > 0:
        if vaf >= 0.3:  # High VAF
            score += 10
        elif vaf >= 0.1:  # Medium VAF
            score += 5
    
    # NORMALIZED SCORING: Convert to 0-100 scale for fair competition
    # This ensures variants with different numbers of criteria compete fairly
    num_criteria = get_num_criteria(vtype)
    if num_criteria > 0:
        # Calculate maximum possible score for this variant type
        # Base score + (all criteria met: +15 each) + (contig bonus: +15) + (read bonus: +10) + (VAF bonus: +10)
        max_possible_score = base_scores.get(vtype, 0) + (num_criteria * 15) + 15 + 10 + 10
        
        # Calculate minimum possible score for this variant type  
        # Base score + (all criteria failed: -15 each) + (no contig bonus: +0) + (no read bonus: +0) + (no VAF bonus: +0)
        min_possible_score = base_scores.get(vtype, 0) - (num_criteria * 15)
        
        # Normalize score to 0-100 scale
        if max_possible_score > min_possible_score:
            normalized_score = ((score - min_possible_score) / (max_possible_score - min_possible_score)) * 100
            # Ensure score is within 0-100 bounds
            normalized_score = max(0, min(100, normalized_score))
            logging.debug(f"Variant {vtype}: Raw score={score}, Min={min_possible_score}, Max={max_possible_score}, Normalized={normalized_score:.1f}")
            return normalized_score
        else:
            # Fallback if max and min are the same
            return 50.0  # Middle score
    
    # Fallback for unknown variant types
    return score

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
    Get detailed breakdown of scoring for a variant (for debugging)
    
    Args:
        variant_row: Dictionary or pandas Series with variant information
        
    Returns:
        Dictionary with scoring breakdown
    """
    breakdown = {}
    
    # Base score
    vtype = variant_row.get('variant_type', '')
    base_scores = {
        'FUS': 140, 'IGR': 140, 'DEL': 80, 'INS': 80,
        'AS': 60, 'NEJ': 60, 'PNJ': 60, 'RI': 40, 'NE': 40, 'EE': 40, 'UN': 10
    }
    breakdown['base_score'] = base_scores.get(vtype, 0)
    
    # FUSION PRIORITY BONUS
    fusion_priority_bonus = 1000 if vtype in ['FUS', 'IGR'] else 0
    breakdown['fusion_priority_bonus'] = fusion_priority_bonus
    
    # Size bonuses (now handled in main scoring, but showing for debugging)
    size = abs(int(variant_row.get('varsize', 0)))
    size_bonus = 0
    # Note: Size scoring is now handled in the main variant-specific criteria
    # This is just for debugging display
    breakdown['size_bonus'] = size_bonus
    
    # Contig evidence bonus
    contig_varsize = abs(int(variant_row.get('contig_varsize', 0)))
    contig_bonus = 0
    if contig_varsize >= 50:
        contig_bonus = 15
    elif contig_varsize >= 20:
        contig_bonus = 10
    elif contig_varsize >= 10:
        contig_bonus = 5
    breakdown['contig_bonus'] = contig_bonus
    
    # Read support bonus
    num_reads_case = variant_row.get('num_reads_case', 0)
    read_bonus = 0
    if num_reads_case >= 10:
        read_bonus = 10
    elif num_reads_case >= 5:
        read_bonus = 5
    breakdown['read_support_bonus'] = read_bonus
    
    # VAF bonus
    vaf = variant_row.get('VAF', 0)
    vaf_bonus = 0
    if vaf >= 0.3:
        vaf_bonus = 10
    elif vaf >= 0.1:
        vaf_bonus = 5
    breakdown['vaf_bonus'] = vaf_bonus
    
    # Total score
    breakdown['total_score'] = sum([
        breakdown['base_score'],
        breakdown['fusion_priority_bonus'], # Added fusion priority bonus
        breakdown['size_bonus'],
        breakdown['contig_bonus'],
        breakdown['read_support_bonus'],
        breakdown['vaf_bonus']
    ])
    
    # NORMALIZED SCORING INFORMATION
    num_criteria = get_num_criteria(vtype)
    if num_criteria > 0:
        # Calculate normalization parameters
        max_possible_score = breakdown['base_score'] + (num_criteria * 15) + 15 + 10 + 10
        min_possible_score = breakdown['base_score'] - (num_criteria * 15)
        
        breakdown['num_criteria'] = num_criteria
        breakdown['max_possible_score'] = max_possible_score
        breakdown['min_possible_score'] = min_possible_score
        
        # Calculate normalized score
        if max_possible_score > min_possible_score:
            normalized_score = ((breakdown['total_score'] - min_possible_score) / (max_possible_score - min_possible_score)) * 100
            normalized_score = max(0, min(100, normalized_score))
            breakdown['normalized_score'] = round(normalized_score, 1)
        else:
            breakdown['normalized_score'] = 50.0
    else:
        breakdown['num_criteria'] = 0
        breakdown['max_possible_score'] = breakdown['total_score']
        breakdown['min_possible_score'] = breakdown['total_score']
        breakdown['normalized_score'] = breakdown['total_score']
    
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
        'chr1', 'pos1', 'strand1',
        'chr2', 'pos2', 'strand2',
        'variant_type', 'other_variant_type', 'overlapping_genes', 'sample',
        'variant_id', 'partner_id', 'vars_in_contig',
        'varsize', 'contig_varsize', 'cpos',
        'TPM', 'mean_WT_TPM','VAF', 'logFC', 'FDR', 'PValue', 'num_reads_case', 'total_num_reads_controls',
        'large_varsize', 'is_contig_spliced', 'spliced_exon', 'overlaps_exon', 'overlaps_gene',
        'motif', 'valid_motif',
        'contig_id', 'unique_contig_ID', 'contig_len', 'contig_cigar',
        'seq_loc1', 'seq_loc2', 'seq1', 'seq2'
    ]

    # Final column order when DE columns are not available
    final_cols_node = [
        'chr1', 'pos1', 'strand1',
        'chr2', 'pos2', 'strand2',
        'variant_type', 'other_variant_type', 'overlapping_genes', 'sample',
        'variant_id', 'partner_id', 'vars_in_contig',
        'varsize', 'contig_varsize', 'cpos',  
        'TPM', 'mean_WT_TPM', 'VAF',
        'large_varsize', 'is_contig_spliced', 'spliced_exon', 'overlaps_exon', 'overlaps_gene',
        'motif', 'valid_motif',
        'contig_id', 'unique_contig_ID', 'contig_len', 'contig_cigar',
        'seq_loc1', 'seq_loc2', 'seq1', 'seq2'
    ]

    # Sort by PValue (ascending) if available, otherwise by VAF (descending)
    if 'PValue' in contigs.columns:
        contigs['PValue'] = pd.to_numeric(contigs['PValue'], errors='coerce')
        contigs = contigs.sort_values(by='PValue', ascending=True, na_position='last')
    elif 'VAF' in contigs.columns:
        contigs['VAF'] = pd.to_numeric(contigs['VAF'], errors='coerce')
        contigs = contigs.sort_values('VAF', ascending=False, na_position='last')

    # Select columns according to availability
    final_cols = final_cols_de if ran_de else final_cols_node
    cols_present = [c for c in final_cols if c in contigs.columns]
    contigs = contigs[cols_present]

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
        
        # Process VAF data to keep only necessary columns
        if 'TPM' in vafs.columns and 'mean_WT_TPM' in vafs.columns:
            vafs = vafs[['contig_id', 'TPM', 'mean_WT_TPM', 'VAF']].drop_duplicates()
        else:
            # Keep all columns if the expected ones are not present
            vafs = vafs.drop_duplicates()
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
    # Ensure sequence columns are present and preserved going forward
    kept_contigs = ensure_seq_columns(kept_contigs)
    discarded_contigs = ensure_seq_columns(discarded_contigs)

    # Keep a copy of kept variants (pre-consolidation) for optional ranked export later
    all_variants = kept_contigs.copy()

    # STEP 2: Apply improved variant consolidation using scoring system on remaining contigs
    logging.info("Processing %d filtered variants using improved scoring system", len(kept_contigs))
    contigs = add_other_variant_types_by_contig(kept_contigs)
    logging.info(f"After variant consolidation: {len(contigs)} variants")

    # Validate the scoring system worked correctly
    validate_scoring_system(contigs)

    # Calculate variants per contig (needed for downstream processing)
    vars_per_contig = contigs.groupby('contig_id', as_index=False)
    vars_per_contig = vars_per_contig.agg({'variant_id': lambda x: len(np.unique(x))})
    vars_per_contig = vars_per_contig.rename({'variant_id': 'vars_in_contig'}, axis=1)
    contigs = contigs.merge(vars_per_contig, on='contig_id')
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
    contigs = pd.merge(contigs, vafs, on='contig_id', how='left')
    logging.info("DE and VAF information added.")

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

    # If requested, produce a discard TSV with identical columns to final output
    if getattr(args, 'discard_out', '') and not discarded_contigs.empty:
        try:
            logging.info('Preparing discard output with final columns...')
            # Add DE/VAF info
            discarded_full = add_de_info(discarded_contigs.copy(), de_results)
            discarded_full = pd.merge(discarded_full, vafs, on='contig_id', how='left')
            # other_variant_type may not exist; ensure it exists for reformat
            if 'other_variant_type' not in discarded_full.columns:
                discarded_full['other_variant_type'] = ''
            # Generate unique contig IDs and sequences
            short_gnames_disc = discarded_full.overlapping_genes.map(str).apply(get_short_gene_name)
            disc_ids, disc_samples = discarded_full.contig_id, discarded_full['sample']
            disc_names = ['|'.join([s, cid, sg]) for cid, s, sg in zip(disc_ids, disc_samples, short_gnames_disc)]
            discarded_full['unique_contig_ID'] = disc_names
            # Sequences already present; ensure columns exist
            discarded_full = ensure_seq_columns(discarded_full)
            # Reformat to final columns and write
            discarded_final = reformat_fields(discarded_full)
            discarded_final.to_csv(args.discard_out, index=False, sep='\t', na_rep='NA')
            logging.info(f"Wrote discard file to {args.discard_out}")
        except Exception as e:
            logging.warning(f"Failed to write discard output: {str(e)}")

    # If requested, produce an all-variants ranked TSV with the same column order
    if getattr(args, 'all_variants_out', '') and not all_variants.empty:
        try:
            # Compute ranking on the filtered set (pre-consolidation copy)
            ranked_all = rank_variants_per_contig(all_variants)
            # Add DE/VAF info to all_variants
            ranked_all = add_de_info(ranked_all, de_results)
            ranked_all = pd.merge(ranked_all, vafs, on='contig_id', how='left')
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
            # Compose final ordering to strictly follow the main output's column order
            base_cols = list(contigs.columns)
            # Define ranking columns that should be preserved
            ranking_cols = [c for c in ['variant_score', 'rank_within_contig', 'is_primary'] if c in ranked_all.columns]
            # Filter out extra DE columns that aren't in the original output
            extra_de_cols = ['variant_of_interest', 'logCPM', 'F', 'len', 'Overdispersion']
            # Start with the exact column order from the main output
            final_cols = [c for c in base_cols if c in ranked_all.columns]
            # Add ranking columns at the end
            final_cols.extend(ranking_cols)
            # Select only these columns
            ranked_all = ranked_all[final_cols]
            # Sort to match original output behavior (by PValue if available else VAF desc)
            if 'PValue' in ranked_all.columns:
                ranked_all['PValue'] = pd.to_numeric(ranked_all['PValue'], errors='coerce')
                ranked_all = ranked_all.sort_values(by='PValue', ascending=True, na_position='last')
            elif 'VAF' in ranked_all.columns:
                ranked_all['VAF'] = pd.to_numeric(ranked_all['VAF'], errors='coerce')
                ranked_all = ranked_all.sort_values('VAF', ascending=False, na_position='last')
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