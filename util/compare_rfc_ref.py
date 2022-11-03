#!/usr/bin/env python3

import os, sys
import subprocess
from Bio import SeqIO


def dmnd_blast(fasta_file, ref_genome):
    dmnd_cmd = 'diamond blastx -k 1 ' \
        f'-q {fasta_file} ' \
        f'-d {ref_genome} ' \
        f'-o {fasta_file.split("/")[-1].split(".fas")[0]}.RefGenome.tsv ' \
        '-f 6 qseqid qlen sseqid slen pident length qstart qend sstart send qframe'

    os.system(dmnd_cmd)
    return f'{fasta_file.split("/")[-1].split(".fas")[0]}.RefGenome.tsv'

def eval_hits(fasta_file, ref_genome):
    tsv = dmnd_blast(fasta_file, ref_genome)
    hits = [i for i in open(tsv).readlines()]
    unique_hits = set([i.split('\t')[0] for i in hits])
    frame_hits = set([i.split('\t')[0] for i in hits if i.endswith('\t1\n')])
    num_seqs = len([i for i in SeqIO.parse(fasta_file,'fasta')])
    print(f'Number of seqs: {num_seqs}')
    print(f'Number of hits: {len(unique_hits)}')
    print(f'Number of hits in-frame: {len(frame_hits)}')

if __name__ == '__main__':
    if len(sys.argv[1:]) == 2:
        fas = sys.argv[1]
        rg = sys.argv[2]
    else:
        print('Usage:\n    python3 compare_rfc_ref.py [RFC-PRED-FASTA-FILE] '
            '[REF-GENOME-AA-FASTA]\n')
        sys.exit(1)

    eval_hits(fas, rg)
