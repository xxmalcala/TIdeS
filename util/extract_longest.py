#!/usr/bin/env python3
import sys
from collections import defaultdict
from Bio import SeqIO

def get_longest(fasta_file):
    txpt_dict = defaultdict(list)

    for i in SeqIO.parse(fasta_file,'fasta'):
        txpt_dict[i.id.split('.pORF')[0]].append(i)

    with open(f'{fasta_file.split("/")[-1].split(".fas")[0]}.Longest.fas','w+') as w:
        for k, v in txpt_dict.items():
            lngst = max(v, key=lambda x: len(x.seq))
            w.write(f'>{lngst.id}\n{lngst.seq}\n')

if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
    except:
        print('Usage:\n    python3 extract_longest.py [FASTA-FILE]')
        sys.exit()

    get_longest(fasta_file)
