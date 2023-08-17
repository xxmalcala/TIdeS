#!/usr/bin/env python3


import sys

from collections import defaultdict

from Bio import SeqIO


def eval_cds_fas(fasta_file):
    comp_seqs = defaultdict(list)
    td_predict = False

    outf = fasta_file.rpartition("/")[-1].split(".fas")[0]

    for i in SeqIO.parse(fasta_file, 'fasta'):
        scr = None
        if not (
            len(i.seq)%3 == 0 and
            f'{i.seq[:3]}'.upper() == 'ATG' and
            f'{i.seq[-3:]}'.upper() in ['TGA','TAA','TAG']
            ):
            continue

        if 'score=' in i.description:
            td_predict = True
            scr = float(i.description.rpartition('score=')[-1].split()[0])

        comp_seqs[i.id.rpartition('.')[0]].append((len(i.seq), i, scr))

    if td_predict:
        outf += '.MaxScore.fas'
        top_seqs = [max(s, key = lambda x: x[-1])[1] for s in comp_seqs.values()]

    else:
        outf += '.LongestSeq.fas'
        top_seqs = [max(s, key = lambda x: x[0])[1] for s in comp_seqs.values()]

    SeqIO.write(top_seqs, outf, 'fasta')

# def eval_cds_fas(fasta_file):
#     out_name = f'{fasta_file}.MaxScore.fas'
#     scores = defaultdict(list)
#     comp_scores = defaultdict(list)
#     for i in SeqIO.parse(fasta_file,'fasta'):
#         scr = float(i.description.split("score=")[1].split()[0])
#         scores[i.id.split('.p')[0]].append((scr, i))
#         if f'{i.seq}'[:3].upper() == 'ATG' and f'{i.seq}'[-3:].upper() in ['TGA','TAG','TAA'] and len(i.seq)%3 == 0:
#             comp_scores[i.id.split('.p')[0]].append((scr, i))
#
#
#     max_scr_seqs = []
#     for k, v in scores.items():
#         v.sort(key=lambda x: -float(x[0]))
#         max_scr_seqs.append(v[0][1])
#     SeqIO.write(max_scr_seqs, out_name, 'fasta')
#
#
#     max_scr_seqs = []
#     for k, v in comp_scores.items():
#         v.sort(key=lambda x: -float(x[0]))
#         max_scr_seqs.append(v[0][1])
#     SeqIO.write(max_scr_seqs, out_name.replace("MaxScore.fas","MaxScore.Comp.fas"), 'fasta')


if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
    except:
        print('\nUsage:\n    python3 best_td_predict.py [TransDecoder.cds-FILE]\n')
        sys.exit()

    eval_cds_fas(fasta_file)
