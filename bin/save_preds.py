#!/usr/bin/env python3

import shutil
from Bio import SeqIO

from bin import orf_call as oc

def finalize_outputs(txn_code: str,  ref_orf_fas: str, imp_orfs: list, minor_orfs: list, gcode: str = '1', min_len: int = 300, contam = False):
    if contam:
        imp_orf_fas = f'{txn_code}_TIdeS/Classified/{txn_code}.TIdeS.fas'
        mnr_orf_fas = f'{txn_code}_TIdeS/Classified/{txn_code}.NonTarget.TIdeS.fas'
    else:
        imp_orf_fas = f'{txn_code}_TIdeS/Classified/{txn_code}.TIdeS.fas'
        mnr_orf_fas = f'{txn_code}_TIdeS/Classified/{txn_code}.TIdeS.AllPreds.fas'

    tds_ps = {}

    ref_dict = {i.id:i for i in SeqIO.parse(ref_orf_fas, 'fasta')}

    with open(imp_orf_fas, 'w+') as w:
        for i in imp_orfs:
            w.write(f'>{ref_dict[i].description}\n{ref_dict[i].seq}\n')
            tds_ps[ref_dict[i].description] = f'{ref_dict[i].seq}'

    with open(mnr_orf_fas, 'w+') as w:
        for i in minor_orfs:
            w.write(f'>{ref_dict[i].description}\n{ref_dict[i].seq}\n')

    with open(imp_orf_fas.replace("fas","AA.fas"), 'w+') as w:
        tds_psa = oc.translate_seqs(tds_ps, gcode)
        for k, v in tds_psa.items():
                w.write(f'>{k}\n{v}\n')

    shutil.copy2(imp_orf_fas, f'{txn_code}_TIdeS/')
    shutil.copy2(imp_orf_fas.replace("fas","AA.fas"), f'{txn_code}_TIdeS/')
