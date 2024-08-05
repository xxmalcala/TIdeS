#!/usr/bin/env python3

"""
Controls the training and prediction steps for TIdeS.

Dependencies include: Optuna, NumPy, Scikit-learn.
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
from sklearn.model_selection import train_test_split, cross_val_score, TunedThresholdClassifierCV
from sklearn.metrics import classification_report, make_scorer, f1_score


def objective_svm(trial, X_features, X_labels, threads):
    """
    Runs and evaluates Optuna trials

    Parameters
    ----------
    trial:       Optuna trial to run/evaluate
    X_features:  training features for hyperparameter tuning
    X_labels:    feature labels for evaluating hyperparameter tunings
    threads:     number of threads to use for training

    Returns
    ----------
    clf_score:  classifier tuning score
    """

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
        'random_state': random.randint(0,10000),
        'tol': trial.suggest_float('tol', 1e-6, 1e-2),
        }

    clf = SVC(**params)

    clf_score = cross_val_score(
                    clf,
                    train_x,
                    train_y,
                    n_jobs = int(threads),
                    cv = 5
                    ).mean()
    return clf_score


def train_model(train_data, threads, contam):
    """
    Manages SVC model training, calling upon Optuna

    Parameters
    ----------
    train_data:  dataframe with training data features and associated labels
    threads:     number of threads to use for training
    contam:      True if evaluating contamination

    Returns
    ----------
    opt_clf:       hyperparameter optimized SVC model
    study:         outputs of the Optuna tuning
    *_clf_report:  classifier performance summaries
    """

    X_features, X_labels = train_data[1:]

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    # optuna.logging.set_verbosity(optuna.logging.INFO)

    study = optuna.create_study(direction='maximize')

    study.optimize(lambda trial: objective_svm(trial, X_features, X_labels, threads), n_jobs = int(threads), n_trials = 100)

    opt_trial = study.best_trial
    # print(study.best_params)

    opt_clf = SVC(probability = True, cache_size = 500, random_state = random.randint(0,10000), **opt_trial.params)

    dm_clf = DummyClassifier(strategy = "stratified", random_state = random.randint(0,10000))

    train_x, test_x, train_y, test_y = train_test_split(
                                        X_features,
                                        X_labels,
                                        shuffle = True,
                                        test_size = 0.25,
                                        random_state = random.randint(0,10000),
                                        stratify = X_labels
                                        )


    if not contam:
        scorer = make_scorer(f1_score, pos_label = 1)
        opt_clf = TunedThresholdClassifierCV(opt_clf, scoring = scorer, n_jobs = int(threads))

        scorer(opt_clf.fit(train_x, train_y), train_x, train_y)

    else:
        opt_clf.fit(train_x, train_y)

    opt_preds = opt_clf.predict(test_x)

    dm_clf.fit(train_x, train_y)

    dm_preds = dm_clf.predict(test_x)

    opt_clf_report = classification_report(test_y, opt_preds)
    dm_clf_report = classification_report(test_y, dm_preds)

    return opt_clf, study, (opt_clf_report, dm_clf_report)


def parse_qpreds(q_names, q_preds, eval_names, contam = False):
    """
    Parses the query ORFs and returns the prediction summary from the SVC

    Parameters
    ----------
    q_names:     Query ORF names
    q_preds:     Probability of the query ORF being in or out of frame
    eval_names:  List of classification categories (binary for ORFs, potentially multi-class for contamination)
    contam:      Whether the evaluations are for contamination

    Returns
    ----------
    summary:  summary of query ORF classifications
    """

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
    """
    Capture distance between two sets of points

    Parameters
    ----------
    co1:  position 1
    co2:  position 2

    Returns
    ----------
    distance between positions/coordinates
    """

    return ((abs(co1[0] - co2[0])**2) + (abs(co1[1] - co2[1])**2))**0.5


def best_preds(pred_list):
    """
    Keeps the best predictions per transcript. For ORF-calling, prioritized by probability of being in-frame,
    then by length (longer chosen over shorter). For contam, simply by classified group.

    Parameters
    ----------
    pred_list:  summary of predictions

    Returns
    ----------
    preds_to_keep:  list of best predictions
    """

    preds_to_keep = []

    crf_evals = sorted([i for i in pred_list if 'CRF' in i], key = lambda x: (-float(x[-1]), -int(x[1].split()[2].split(':')[-1])))
    # crf_evals = sorted([i for i in pred_list if 'CRF' in i], key = lambda x: (-float(x[-1]), -int(x[1].split('_')[4])))

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


def classify_orfs(taxon_code: str, start_time, train_data, query_data, threads: int = -1, model = None, training: bool = False, contam: bool = False, verb: bool = True):
    """
    Manages the ORF classification

    Parameters
    ----------
    taxon_code:  species/taxon name or abbreviated code
    start_time:  start time for current step
    train_data:  dataframe with training data features and associated labels
    query_data:  dataframe with query ORF features and ORF names
    threads:     number of threads to use
    model:       user-provided pre-trained model
    contam:      True if evaluating contamination
    verb:        print verbose messages

    Returns
    ----------
    query_summary:  prediction summary for the query datasets
    trained_clf:    optuna optimized SVC model
    """

    clf_dir = f'{taxon_code}_TIdeS/Classification/'
    clf_tsv = f'{clf_dir}{taxon_code}.TIdeS_Classification.tsv'
    clf_tp_sb_tsv = f'{clf_dir}{taxon_code}.TIdeS_Classification.SingleBest.tsv'
    clf_scr_tsv = f'{clf_dir}{taxon_code}.TIdeS_Classification_Report.txt'
    clf_stdy = f'{clf_dir}{taxon_code}.TIdeS.OptunaStudy.pkl'

    non_contam_header = 'Transcript\tpORF\tFull-pORF-Name\tEvaluation\tProb-OOF\tProb-CRF\n'

    Path(clf_dir).mkdir(exist_ok = True, parents = True)

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

        with open(clf_stdy, 'wb') as fout:
            pickle.dump(opt_study, fout)


    if not training:

        if verb:
            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Classifying ORFs')

        query_names, query_features = query_data

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

    else:
        return None, trained_clf
