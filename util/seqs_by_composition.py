#!/usr/bin/env python3

import argparse
import random
import sys

import pandas as pd

from Bio import SeqIO

def collect_args():
    parser = argparse.ArgumentParser(description = f'Prepares txt files for classifying predicted ORFs.',
            usage=argparse.SUPPRESS, add_help = False,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    g = parser.add_argument_group('General Options', description = (
    '''--fin (-f)      input file in FASTA format\n'''
    '''--tsv (-t)      TSV-output from orf_composition.py\n'''
    '''--out (-o)      output-name [e.g., taxon-name, cluster-1]\n'''
    '''--min_gc12      minimum %GC12 (default = 0)\n'''
    '''--max_gc12      maximum %GC12 (default = 100)\n'''
    '''--min_gc3       minimum %GC3 (default = 0)\n'''
    '''--max_gc3       maximum %GC3 (default = 100)\n'''
    '''--min_gc3d      minimum %GC3 at 4-fold degenerate sites (default = 0)\n'''
    '''--max_gc3d      minimum %GC3 at 4-fold degenerate sites (default = 100)\n'''
    '''--nseqs         number of sequences to randomly select (default = 100)\n'''
    '''--reps          number of random replicates to generate (default = 1)\n'''
    '''--help (-h)     show this help message and exit\n'''))

    g.add_argument('--help', '-h', action="help", help=argparse.SUPPRESS)

    g.add_argument('--fin', '-f', action = 'store',
        metavar = ('[FASTA File]'), type = str, required = True,
        help = argparse.SUPPRESS)

    g.add_argument('--tsv', '-t', action = 'store',
        metavar = ('[TSV File]'), type = str, required = True,
        help = argparse.SUPPRESS)

    g.add_argument('--out','-o', action = 'store',
        type = str, required = True,
        help = argparse.SUPPRESS)

    g.add_argument('--min_gc12', action = 'store',
        default = 0, type = float,
        help = argparse.SUPPRESS)

    g.add_argument('--max_gc12', action = 'store',
        default = 100, type = float,
        help = argparse.SUPPRESS)

    g.add_argument('--min_gc3', action = 'store',
        default = 0, type = float,
        help = argparse.SUPPRESS)

    g.add_argument('--max_gc3', action = 'store',
        default = 100, type = float,
        help = argparse.SUPPRESS)

    g.add_argument('--min_gc3d', action = 'store',
        default = 0, type = float,
        help = argparse.SUPPRESS)

    g.add_argument('--max_gc3d', action = 'store',
        default = 100, type = float,
        help = argparse.SUPPRESS)

    g.add_argument('--nseqs', action = 'store',
        default = 100, type = int,
        help = argparse.SUPPRESS)

    g.add_argument('--reps', action = 'store',
        default = 1, type = int,
        help = argparse.SUPPRESS)

    args = parser.parse_args()

    return args

def grab_seqs(tsv_file: str, min_3: float = 0, max_3: float = 100, min_12:float = 0, max_12: float = 100, min_3d: float = 0, max_3d: float = 100):

    df = pd.read_table(tsv_file)
    df2 = df.loc[(df['GC12'] >= min_12) & (df['GC12'] <= max_12)
                & (df['GC3'] >= min_3) & (df['GC3'] <= max_3)
                & (df['GC3_DGN'] >= min_3d) & (df['GC3_DGN'] <= max_3d)]

    return list(df2.ORF)


def store_info(filtered_seqs: list, fasta_file: str, out_name: str, nseqs: int = 100, replicate: int = 1, all_seqs: bool = False):
    if all_seqs:
        out_fas = f'{fasta_file.rpartition(".fas")[0]}.{out_name}.FilteredSeqs.fasta'
        tmp_seqs = [i for i in SeqIO.parse(fasta_file,'fasta') if i.description in filtered_seqs]
        SeqIO.write(tmp_seqs, out_fas, 'fasta')
        return None

    filt_seqs_annot = [f'{i.split(" ")[0]}\t{out_name}\n' for i in filtered_seqs]
    out_txt = f'{fasta_file.rpartition(".fas")[0]}.{out_name}.FilteredSeqs'

    if replicate > 1:
        for n in range(1, replicate+1):
            random.shuffle(filt_seqs_annot)
            if nseqs > 0:
                with open(f'{out_txt}.Rep{n}.txt','w+') as w:
                    w.write(''.join(filt_seqs_annot[:nseqs]))
    else:
        if nseqs > 0:
            with open(f'{out_txt}.txt','w+') as w:
                w.write(''.join(filt_seqs_annot[:nseqs]))

    return None

if __name__ == '__main__':
    args = collect_args()

    filt_seqs = grab_seqs(
                    args.tsv,
                    args.min_gc3,
                    args.max_gc3,
                    args.min_gc12,
                    args.max_gc12,
                    args.min_gc3d,
                    args.max_gc3d)

    store_info(
        filt_seqs,
        args.fin,
        args.out,
        args.nseqs,
        args.reps)
