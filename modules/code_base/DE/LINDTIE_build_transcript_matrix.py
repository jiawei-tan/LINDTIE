'''
Module      : LINDTIE_build_transcript_matrix
Description : Builds transcript count matrix from Oarfish .quant files
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
Portability : POSIX
Combines all quant files into a transcript count matrix, 
and writes the matrix to a tab-separated value (TSV) file.

Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
'''

import pandas as pd
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description='Build transcript count matrix from Oarfish .quant files')
    parser.add_argument('inputs', nargs='*', help="quant files, sample names (comma-separated), and output filename")
    args = parser.parse_args()

    if len(args.inputs) < 3:
        parser.error("At least 3 arguments required: quant files, sample names (comma-separated), and output filename.")

    return args

def load_quant_file(filepath):
    '''
    Load Oarfish .quant file with columns: tname, len, num_reads
    '''
    df = pd.read_csv(filepath, sep='\t')
    if 'tname' not in df.columns or 'num_reads' not in df.columns:
        raise ValueError(f"Missing required columns in file: {filepath}")
    return df[['tname', 'num_reads']].set_index('tname')

def build_transcript_matrix(quant_files, sample_names):
    '''
    Combine all quant files into a transcript count matrix
    '''
    count_dfs = []
    for file, sample in zip(quant_files, sample_names):
        df = load_quant_file(file)
        df.columns = [sample]
        count_dfs.append(df)

    count_matrix = pd.concat(count_dfs, axis=1).fillna(0).astype(int)
    count_matrix.index.name = 'transcript_id'
    return count_matrix.reset_index()

def main():
    args = parse_args()
    quant_files = args.inputs[:-2]
    sample_names = args.inputs[-2].split(',')
    outfile = args.inputs[-1]

    if len(quant_files) != len(sample_names):
        raise ValueError("The number of quant files must match the number of sample names.")

    print('Loading Oarfish transcript counts...')
    count_matrix = build_transcript_matrix(quant_files, sample_names)

    print(f'Writing transcript count matrix to {outfile}...')
    count_matrix.to_csv(outfile, sep='\t', index=False)
    print('Done.')

if __name__ == '__main__':
    main()
