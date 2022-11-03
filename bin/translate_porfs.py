#!/usr/bin/env python3

import shutil
from Bio import SeqIO

"""Translates the pORFs using a given (and supported) translation table."""

def prep_final_pORFs(pORF_fas, ttable, taxon_dir):
    print(f'Translating pORFs using translation table: {ttable}')
    shutil.copy2(pORF_fas, taxon_dir)

    pORF_aa_fas = pORF_fas.rpartition("/")[-1].replace(".fas",".AA.fas")
    with open(f'{taxon_dir}{pORF_aa_fas}','w+') as w:
        for i in SeqIO.parse(pORF_fas, 'fasta'):
            w.write(f'>{i.description}\n{i.seq[:-3].translate(table=int(ttable))}\n')

    return f'{taxon_dir}{pORF_aa_fas.replace("AA.","")}', f'{taxon_dir}{pORF_aa_fas}'
