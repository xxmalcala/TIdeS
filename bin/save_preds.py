#!/usr/bin/env python3

"""
Something something...


Dependencies include: Optuna, Scikit-learn, and XGBoost.
"""

import pickle, sys

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


"""Need to add a log file with params, so that if given those data, then all options
are set to be the same...

This can be stored in the TIdeS.pkl file too..."""


def save_model(taxon_code, overlap, kmer, step, cvec, clf):
    mdl_pkl = f'{taxon_code}_TIdeS/{taxon_code}.TIdeS.pkl'

    with open(mdl_pkl, 'wb') as fout:
        pickle.dump((overlap, kmer, step, cvec, clf), fout)

def save_seqs(taxon_code, all_orfs, clf_summary, gcode, single_best = True, contam = False):
    if not contam:
        ntd_seqs = []
        aa_seqs = []
        if single_best:
            out_ntd = f'{taxon_code}_TIdeS/{taxon_code}.TIdeS.fasta'
            for k, v in clf_summary[1].items():
                s = v[0][1]

                ntd_seqs.append(SeqRecord(Seq(all_orfs[s]), s.split()[0], '', s))

                aa_seqs.append(SeqRecord(Seq(all_orfs[s]).translate(table = gcode), s.split()[0], '', s))

        else:
            out_ntd = f'{taxon_code}_TIdeS/{taxon_code}.TIdeS.All_CRF.fasta'

        out_aa = out_ntd.replace("TIdeS.","TIdeS.AA.")

        SeqIO.write(ntd_seqs, out_ntd, 'fasta')
        SeqIO.write(aa_seqs, out_aa, 'fasta')

    else:
        class_seqs = defaultdict(list)

        for k, v in all_orfs.items():
            class_seqs[clf_summary[k][0]].append(SeqRecord(Seq(v), k, '', ''))

        for k, v in class_seqs.items():
            out_ntd =  f'{taxon_code}_TIdeS/{taxon_code}.TIdeS.{k}.fasta'
            SeqIO.write(v, out_ntd, 'fasta')
