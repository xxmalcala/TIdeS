#!/usr/bin/env python3

"""
Provides a quick summary of the basic composition of a given FASTA file of
nucleotide sequences.

Returns a table of GC1, GC2, GC3, GC12, and GC3 at four-fold degenerate sites.
"""

import sys
from pathlib import Path

import pandas as pd


def eval_gcode(gencode: str) -> list:
    supp_codes = {
        '1':['TGA','TAG','TAA'], '6':['TGA'], '10':['TAG','TAA']
    }
    if gencode in supp_codes:
        return supp_codes[gencode]
    else:
        print('Error: Unsupported Genetic Code Provided.')
        print('For stop-codon recognition, the only supported codes are:')
        print('--- 1  [TGA, TAG, TAA]')
        print('--- 6  [TGA]')
        print('--- 10 [TAG, TAA]')
        sys.exit()


def calc_gc_content(gc_val: str) -> str:
    return f'{100*(gc_val.count("G") + gc_val.count("C"))/len(gc_val):.2f}'


def calc_gc_stats(seq: str, stop_codons: list) -> list:
    gc_stats = {'GC1':[''], 'GC2':[''], 'GC3':[''], 'GC12':[''], 'GC3_DGN':['']}
    gc1 = ''
    gc2 = ''
    gc3 = ''
    gc12 = ''
    gc3_dgn = ''
    degenerate_codons = ['TC','CT','CC','CG','AC','GT','GC','GG']
    for n in range(0, len(seq), 3):
        if seq[n:n+3] not in stop_codons:
            gc_stats['GC1'][0] += seq[n]
            gc_stats['GC2'][0] += seq[n+1]
            gc_stats['GC3'][0] += seq[n+2]
            gc_stats['GC12'][0] += seq[n:n+2]
            if seq[n:n+2] in degenerate_codons:
                gc_stats['GC3_DGN'][0] += seq[n+2]
    for k, v in gc_stats.items():
        gc_stats[k] = calc_gc_content(v[0])
    gc_stats['Seq_Length'] = len(seq)
    return gc_stats


def parse_fasta_file(fasta_file: str, gencode: str) -> dict:
    seq_comp_dict = {}
    out_tsv = f'ORF_Composition/{fasta_file.rpartition(".fa")[0]}.CompSummary.tsv'

    stop_codons = eval_gcode(gencode)

    for i in open(fasta_file).read().split('>')[1:]:
        seq_id = i.partition('\n')[0]
        seq = i.partition('\n')[2].replace('\n','').upper()
        seq_comp_dict[seq_id] = calc_gc_stats(seq, stop_codons)

    df = pd.DataFrame.from_dict(seq_comp_dict, orient = 'index')
    df.index.name = 'ORF'
    df.to_csv(out_tsv, sep = '\t')

if __name__ == '__main__':
    if len(sys.argv[1:]) == 2:
        fasta_file = sys.argv[1]
        gcode = sys.argv[2]
    elif len(sys.argv[1:]) == 1:
        fasta_file = sys.argv[1]
        gcode = 1
    else:
        print('Usage:\n    python3 orf_composition.py [FASTA-FILE] [GENETIC-CODE [1, 6, 10]]\n')
        sys.exit()

    Path('ORF_Composition').mkdir(exist_ok = True)

    parse_fasta_file(fasta_file, f'{gcode}')
