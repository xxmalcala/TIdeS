#!/usr/bin/env python3

"""
Controls the training and prediction steps for TIdeS.

Dependencies include: Optuna, Scikit-learn.
"""

import pickle
import random
import sys
import time
from collections import defaultdict
from datetime import timedelta
from pathlib import Path

import optuna
import numpy as np

from sklearn.dummy import DummyClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import classification_report

def objective_svm(trial, X_features, X_labels, threads):
    train_x, test_x, train_y, test_y = train_test_split(
                                        X_features,
                                        X_labels,
                                        shuffle = True,
                                        test_size = 0.25,
                                        random_state = random.randint(0,10000),
                                        stratify = X_labels
                                        )
    params = {
        'C': trial.suggest_float('C', 1e-8, 10),
        'kernel': 'rbf',
        'random_state': random.randint(0,10000)
        }

    clf = SVC(**params)
    clf_score = cross_val_score(
                    clf,
                    train_x,
                    train_y,
                    n_jobs = threads,
                    cv = 5
                    ).mean()

    return clf_score


def train_model(train_data, threads, contam):
    X_features, X_labels = train_data[1:]
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    # optuna.logging.set_verbosity(optuna.logging.INFO)

    study = optuna.create_study(direction='maximize')

    study.optimize(lambda trial: objective_svm(trial, X_features, X_labels, threads), n_trials = 100)

    opt_trial = study.best_trial

    opt_clf = SVC(probability = True, random_state = random.randint(0,10000), **opt_trial.params)

    dm_clf = DummyClassifier(strategy = "stratified", random_state = random.randint(0,10000))

    train_x, test_x, train_y, test_y = train_test_split(
                                        X_features,
                                        X_labels,
                                        shuffle = True,
                                        test_size = 0.25,
                                        random_state = random.randint(0,10000),
                                        stratify = X_labels
                                        )

    opt_clf.fit(train_x, train_y)
    dm_clf.fit(train_x, train_y)

    opt_preds = opt_clf.predict(test_x)
    dm_preds = dm_clf.predict(test_x)

    opt_clf_report = classification_report(test_y, opt_preds)
    dm_clf_report = classification_report(test_y, dm_preds)

    return opt_clf, study, (opt_clf_report, dm_clf_report)


def parse_qpreds(q_names, q_preds, eval_names, contam = False):
    summary = defaultdict(list)

    for n in range(0, len(q_preds)):
        if contam:
            summary[q_names[n]] += [f'{eval_names[np.argmax(q_preds[n])]}']
            for i in q_preds[n]:
                summary[q_names[n]] += [f'{i:.3f}']

        else:
            summary[q_names[n].split('.pORF')[0]].append([
                                                    f'{q_names[n].split()[0]}',
                                                    f'{q_names[n]}',
                                                    f'{eval_names[np.argmax(q_preds[n])]}',
                                                    f'{q_preds[n][0]:.4f}',
                                                    f'{q_preds[n][1]:.4f}'
                                                    ])
    return summary


def distance(co1, co2):
    # Capture distance between two sets of points
    return ((abs(co1[0] - co2[0])**2) + (abs(co1[1] - co2[1])**2))**0.5


def best_preds(pred_list):
    preds_to_keep = []
    crf_evals = sorted([i for i in pred_list if 'CRF' in i], key = lambda x: (-float(x[-1]), -int(x[1].split()[2].split(':')[-1])))

    if crf_evals:

        preds_to_keep = [crf_evals[0]]
        coords_eval = [tuple(sorted(int(i) for i in crf_evals[0][1].split(':')[-1].split('(')[0].split('-')))]

        for pred in crf_evals[1:]:
            crds = sorted(int(i) for i in pred[1].split(':')[-1].split('(')[0].split('-'))
            nearest = min(coords_eval, key = lambda x: distance(x, crds))

            if crds[0] <= nearest[1] and crds[1] >= nearest[0]:
                continue

            else:
                preds_to_keep.append(pred)
                coords_eval.append(tuple(crds))

    return preds_to_keep


def classify_orfs(taxon_code: str, start_time, train_data, query_data, threads: int = -1, model = None, contam = False, verb = True):

    clf_dir = f'{taxon_code}_TIdeS/Classification/'
    clf_tsv = f'{clf_dir}{taxon_code}.TIdeS_Classification.tsv'
    clf_tp_sb_tsv = f'{clf_dir}{taxon_code}.TIdeS_Classification.SingleBest.tsv'
    clf_scr_tsv = f'{clf_dir}{taxon_code}.TIdeS_Classification_Report.txt'
    clf_stdy = f'{clf_dir}{taxon_code}.TIdeS.OptunaStudy.pkl'

    non_contam_header = 'Transcript\tpORF\tFull-pORF-Name\tEvaluation\tProb-OOF\tProb-CRF\n'

    Path(clf_dir).mkdir(exist_ok = True, parents = True)

    query_names, query_features = query_data

    if model:
        trained_clf = model

    else:
        if verb:
            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Training TIdeS')

        trained_clf, opt_study, clf_scrs = train_model(train_data, threads, contam)

        with open(clf_scr_tsv, 'w+') as w:
            w.write('TIdeS-Classification-Report\n')
            w.write(f'{clf_scrs[0]}\n')
            w.write('"Dummy"-Classification-Report\n')
            w.write(f'{clf_scrs[1]}')

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Classifying ORFs')

    trained_clf_qpreds = trained_clf.predict_proba(query_data[1])


    if not contam:
        eval_names = ['OOF', 'CRF']
    else:

        eval_names = trained_clf.classes_
        tmp = '\t'.join([f'Prob-{name}' for name in eval_names])
        contam_header = f'ORF\tEvaluation\t{tmp}\n'

    query_summary = parse_qpreds(query_data[0], trained_clf_qpreds, eval_names, contam)

    with open(clf_tsv, 'w+') as w:
        if contam:
            w.write(contam_header)
            for k, v in query_summary.items():
                x = '\t'.join(v)
                w.write(f'{k}\t{x}\n')

        else:
            w.write(non_contam_header)
            for k, v in query_summary.items():
                for i in v:
                    x = '\t'.join(i)
                    w.write(f'{k}\t{x}\n')
    if not model:
        with open(clf_stdy, 'wb') as fout:
            pickle.dump(opt_study, fout)

    if not contam:
        query_best_preds = {}
        for k, v in query_summary.items():
            bps = best_preds(v)
            if bps:
                query_best_preds[k] = bps

        with open(clf_tp_sb_tsv, 'w+') as w:
            w.write(non_contam_header)
            for k, v in query_best_preds.items():
                if v:
                    x = '\t'.join(v[0])
                    w.write(f'{k}\t{x}\n')

        return (query_summary, query_best_preds), trained_clf

    return query_summary, trained_clf
