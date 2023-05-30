#!/usr/bin/env python3

"""Methods to generate output FASTA-files based on the classification steps from
other modules."""

import shutil, sys
from collections import defaultdict

from Bio import SeqIO

from bin import orf_call as oc

def finalize_outputs(
                    txn_code: str,
                    ref_orf_fas: str,
                    imp_orfs: list,
                    minor_orfs: list,
                    gcode: str = '1',
                    min_len: int = 300,
                    contam = False) -> None:
    tds_ps = {}
    ref_dict = {i.id:i for i in SeqIO.parse(ref_orf_fas, 'fasta')}
    if contam:
        imp_orf_fas = f'{txn_code}_TIdeS/Classified/{txn_code}.TIdeS.fas'
        mnr_orf_fas = f'{txn_code}_TIdeS/Classified/{txn_code}.TIdeS.NonTarget.fas'
    else:
        imp_orf_fas = f'{txn_code}_TIdeS/Classified/{txn_code}.TIdeS.fas'

    with open(imp_orf_fas, 'w+') as w:
        if contam:
            for i in imp_orfs:
                w.write(f'>{ref_dict[i].description}\n{ref_dict[i].seq}\n')
                tds_ps[ref_dict[i].description] = f'{ref_dict[i].seq}'

            with open(mnr_orf_fas, 'w+') as w:
                for i in minor_orfs:
                    w.write(f'>{ref_dict[i].description}\n{ref_dict[i].seq}\n')

        else:
            top_orfs = filter_tides_preds(ref_dict, imp_orfs)

            for i in top_orfs:
                w.write(f'>{i.description}\n{i.seq}\n')
                tds_ps[i.description] = f'{i.seq}'

    with open(imp_orf_fas.replace("fas","AA.fas"), 'w+') as w:
        tds_psa = oc.translate_seqs(tds_ps, gcode)

        for k, v in tds_psa.items():
                w.write(f'>{k}\n{v}\n')

    shutil.copy2(imp_orf_fas, f'{txn_code}_TIdeS/')
    shutil.copy2(imp_orf_fas.replace("fas","AA.fas"), f'{txn_code}_TIdeS/')


def filter_tides_preds(ref_dict: dict, imp_orfs: list) -> list:
    tmp_dict = defaultdict(list)
    top_orfs = []

    for k, v in ref_dict.items():
        if k in imp_orfs:
            tmp_dict[k.split('.pORF')[0]].append(v)

    for k, v in tmp_dict.items():
        if len(v) > 1:
            maxlen = max([len(i) for i in v])
            top_orfs += [i for i in v if len(i) == maxlen]

        else:
            top_orfs += v

    return top_orfs
