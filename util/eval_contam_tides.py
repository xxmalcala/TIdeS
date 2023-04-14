#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np


import seaborn as sns

def parse_tables(preds, ref):
    tp, fp, tn, fn = 0,0,0,0
    tides_preds = {i.split('\t')[0]:i.split('\t')[1].rstrip().replace("-","") for i in open(preds).readlines()[1:]}
    ref_curated = {i.split('\t')[0]:i.split('\t')[1].rstrip().replace("-","") for i in open(ref).readlines()[1:]}
    for k, v in ref_curated.items():
        try:
            pred = tides_preds[k]
        except:
            continue
        if v == pred == 'Target':
            tp += 1
        elif v == pred == 'NonTarget':
            tn += 1
        elif v != pred and pred == 'Target':
            fp += 1
        elif v != pred and pred == 'NonTarget':
            fn += 1
    return np.array([[tp, fp],[fn, tn]])




def save_confusion_matrix(preds, ref):
    matrix_orig = parse_tables(preds, ref)
    matrix = matrix_orig.astype('float') / matrix_orig.sum(axis=1)[:, np.newaxis]*100
    plt.figure(figsize=(7,4))
    sns.set(font_scale=0.8)
    sns.heatmap(matrix, annot=True, annot_kws={'size':8}, fmt=".2f",
            cmap=plt.cm.Blues, linewidths=0.2, cbar=False)

    class_names = ['Target-Seqs', 'Non-Target-Seqs']
    tick_marks = np.arange(len(class_names))+0.5
    tick_marks2 = tick_marks
    plt.xticks(tick_marks, class_names, rotation=0)
    plt.yticks(tick_marks2, class_names, rotation=0)

    plt.title('Confusion Matrix for Random Forest Model (Percentages)')
    plt.ylabel('True Labels')
    plt.xlabel('TIdeS Predicted Labels')
    plt.tight_layout()
    contingency_png = preds.split("/")[-1].replace('.tsv','.ContingencyTable.png')
    plt.savefig(contingency_png, dpi = 300)

    contingency_tsv = preds.split("/")[-1].replace('.tsv','.ContingencyTable.tsv')

    with open(contingency_tsv,'w+') as w:
        w.write('\tTarget\tNon-Target\n')
        w.write(f'Target\t{matrix_orig[0][0]}\t{matrix_orig[0][1]}\n')
        w.write(f'Non-Target\t{matrix_orig[1][0]}\t{matrix_orig[1][1]}\n')

if __name__ == '__main__':
    try:
        tides_preds = sys.argv[1]
        scored_hits = sys.argv[2]
    except:
        print('Usage:\n    python3 eval_contam_tides.py [TIdeS-Predicted-Classification] [Reference-Sister-Summary]\n')
        sys.exit(1)

    save_confusion_matrix(tides_preds, scored_hits)
