#!/usr/bin/env python3

import itertools, sys
from Bio import SeqIO
import pandas as pd
from pathlib import Path

"""Script simply returns the codon usage frequency (columns) for each transcript
(rows) of a given FASTA file of ORFs.

Can be run as a standalone tool."""

# Returns codon usage for a single ORF.
def get_codon_usage(seq):
    codon_dict = {''.join(c):0 for c in list(itertools.product(['G','A','C','T'], repeat=3))}
    total_codons = 0
    mlen = len(seq) - len(seq)%3
    for n in range(0, mlen, 3):
        if seq[n:n+3] in codon_dict.keys():
            total_codons += 1
            codon_dict[seq[n:n+3].upper()] += 1
    for k, v in codon_dict.items():
        codon_dict[k] = round(v*100/total_codons,3)
    return codon_dict

# Generates the codon usage table for all ORFs in a fasta file.
def codon_counts_fasta(fasta_file, out_dir, taxon_code, training = True, RF = True):
    Path(f'{out_dir}Codon_Counts').mkdir(exist_ok = True, parents = True)

    if training and RF:
        cc_tsv = f'{out_dir}Codon_Counts/{taxon_code}.CodonCounts.RandomRF.tsv'
    elif training and not RF:
        cc_tsv = f'{out_dir}Codon_Counts/{taxon_code}.CodonCounts.RefContamEval.tsv'
    else:
        cc_tsv = f'{out_dir}Codon_Counts/{taxon_code}.CodonCounts.tsv'

    seq_db = {}
    for seq_rec in SeqIO.parse(fasta_file,'fasta'):
        seq_db[seq_rec.id] = get_codon_usage(f'{seq_rec.seq}')

    df = pd.DataFrame(seq_db)
    # Transform dataframe, genes as rows rather than columns
    df = df.T
    df.index.name = 'Gene'
    df = df.reset_index(level=0)

    if training:
        if RF:
            # Label genes in the +1 reading frame as 1, all others as 0.
            rf = df['Gene'].map(lambda x: 1 if x.endswith("RF1") else 0)
            df.insert(loc=1, column='RF', value=rf)
        else:
            stype = df['Gene'].map(lambda x: 1 if x.endswith("_Ref") else 0)
            df.insert(loc=1, column='REF', value=stype)
        # Randomly shuffles the rows, to make subsampling a bit easier later.
        # Important as there is no need to sample the ENTIRE training set...
        df = df.sample(frac=1).reset_index(drop=True)

    df.to_csv(cc_tsv, sep = '\t', index = False)

    return cc_tsv

if __name__ == '__main__':
    if len(sys.argv[1:]) == 3:
        fasta_file = sys.argv[1]
        taxon_code = sys.argv[2]
        out_dir = sys.argv[3]
    else:
        print('Usage:\n    python3 codon_counts.py [FASTA-FILE] [TAXON-CODE] [OUTPUT-NAME]\n')
        sys.exit(1)

    codon_counts_fasta(fasta_file, out_dir, taxon_code, training = True)
