#!/usr/bin/env python3


import sys

from collections import defaultdict

from Bio import SeqIO


def eval_cds_fas(fasta_file):
    comp_seqs = defaultdict(list)

    outf = fasta_file.rpartition("/")[-1].split(".fas")[0]

    for i in SeqIO.parse(fasta_file, 'fasta'):
        if (
            len(i.seq)%3 == 0 and
            f'{i.seq[:3]}'.upper() == 'ATG' and
            f'{i.seq[-3:]}'.upper() in ['TGA','TAA','TAG']
            ):
            comp_seqs[i.id.rpartition('.')[0]].append((len(i.seq), i))

    else:
        outf += '.LongestSeq.fas'
        top_seqs = [max(s, key = lambda x: x[0])[1] for s in comp_seqs.values()]

    SeqIO.write(top_seqs, outf, 'fasta')



if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
    except:
        print('\nUsage:\n    python3 get_longest [FASTA-FILE]\n')
        sys.exit()

    eval_cds_fas(fasta_file)
