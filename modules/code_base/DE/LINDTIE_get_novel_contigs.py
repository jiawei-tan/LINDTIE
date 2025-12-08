'''
Module      : LINDTIE_get_novel_contigs
Description : Gets novel transcripts from transcript count matrix
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
Portability : POSIX
Reads the transcript count matrix, identifies novel transcripts,
and writes the novel transcripts to a tab-separated value (TSV) file.

Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
'''

import pandas as pd
import numpy as np
import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser

def parse_args(args):
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Get novel transcripts'
    parser = ArgumentParser(description = description)
    parser.add_argument(dest='tcm_file',
                        metavar='TCM_FILE',
                        type=str,
                        help='''transcript_counts_matrix.tsv file.''')
    parser.add_argument(dest='ref_tx_fasta',
                        metavar='REF_FASTA',
                        type=str,
                        help='''Transcriptome reference fasta.''')
    parser.add_argument(dest='denovo_fasta',
                        metavar='DENOVO_FASTA',
                        type=str,
                        help='''Sample de novo filtered fasta assembly file.''')
    return parser.parse_args(args)

def get_ref_txs(ref_tx_fasta):
    '''
    Get all reference transcript IDs
    from reference transcriptome fasta
    '''
    handle = open(ref_tx_fasta, 'r')
    ref_txs = []
    for record in SeqIO.parse(handle, 'fasta'):
        ref_txs.append(record.id)
    handle.close()
    return ref_txs

def get_novel_transcripts(tcm, ref_txs, outdir):
    '''
    Process transcript count matrix to identify novel transcripts
    and prepare output table
    '''
    print('Processing transcript count matrix...')
    
    # Create a copy of the transcript count matrix
    novel_tx = tcm.copy()
    
    # Identify novel transcripts (not in reference)
    novel_tx['is_novel'] = ~novel_tx['transcript_id'].isin(ref_txs)
    
    # Filter to only novel transcripts
    print('Filtering for novel transcripts...')
    novel_tx = novel_tx[novel_tx['is_novel']]
    
    if len(novel_tx) == 0:
        print('Warning: No novel transcripts found!')
        # Create empty dataframe with expected columns
        novel_tx = pd.DataFrame(columns=['contig', 'case_reads'])
    else:
        # Get the first sample column (assuming it's the case sample)
        sample_columns = [col for col in tcm.columns if col != 'transcript_id']
        case_sample = sample_columns[0]  # Use first sample as case
        
        # Prepare output table
        novel_tx = novel_tx[['transcript_id', case_sample]].copy()
        novel_tx = novel_tx.rename(columns={'transcript_id': 'contig', case_sample: 'case_reads'})
    
    # Write full transcript table (including reference transcripts)
    print('Writing full transcript table...')
    full_tx_table = tcm.copy()
    full_tx_table['is_novel'] = ~full_tx_table['transcript_id'].isin(ref_txs)
    full_tx_table.to_csv('%s/transcript_table.txt' % outdir, sep='\t', index=False)
    
    return novel_tx

def main():
    args = parse_args(sys.argv[1:])
    try:
        print('Reading in transcript count matrix file...')
        tcm = pd.read_csv(args.tcm_file, sep = '\t')
        print('Fetching reference transcripts...')
        ref_txs = get_ref_txs(args.ref_tx_fasta)
    except IOError as message:
        print("{} ERROR: {}, exiting".format("get_novel_contigs", message), file=sys.stderr)
        sys.exit(1)

    # get novel transcripts table
    outdir = os.path.dirname(args.tcm_file)
    novel_tx = get_novel_transcripts(tcm, ref_txs, outdir)
    
    # write dummy DE file
    novel_tx.to_csv('%s/transcript_de.txt' % outdir, sep = '\t', index = False)

if __name__ == '__main__':
    main()
