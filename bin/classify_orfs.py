#!/usr/bin/env python3

"""Performs binary classification of given ORFs using supervised training of random forests.


Dependencies include: Sci-Kit Learn."""

import pickle, time
from datetime import timedelta
from pathlib import Path
from numpy import argmax
import numpy as np
import pandas as pd

import optuna

import sklearn.ensemble
import sklearn.metrics
import sklearn.model_selection
import sklearn.naive_bayes

def extract_features(feature_dict: dict, contam: bool, training: bool) -> pd.core.frame.DataFrame:
    feature_df = pd.DataFrame.from_dict(feature_dict, orient = 'index')

    if not training:
        return feature_df

    elif contam:
        y_ref_df = feature_df.index.map(lambda x: 1 if x.split('_')[-1] == 'Target' else 0)

        return feature_df, y_ref_df

    else:
        y_ref_df = feature_df.index.map(lambda x: 1 if x[-1] == '1' else 0)

        return feature_df, y_ref_df


def objective(trial, X_train, y_train, threads, contam):
    if contam:
        mnb_params = {
            'alpha': trial.suggest_float('alpha', 1e-5, 1.0, log = True)
        }
        clf = sklearn.naive_bayes.MultinomialNB(**mnb_params)
    else:
        rf_params = {
            'n_estimators': trial.suggest_int('n_estimators', 100, 500, 100),
            'max_depth': trial.suggest_int('max_depth', 4, 16, 2),
            'max_features': trial.suggest_categorical('max_features', choices=['sqrt', 'log2']),
            'min_samples_split': trial.suggest_int('min_samples_split', 2, 5),
            'criterion': trial.suggest_categorical('criterion', ['gini', 'entropy']),
            'n_jobs': int(threads),
            'oob_score': True,
            'random_state': 121219
            }

        clf = sklearn.ensemble.RandomForestClassifier(**rf_params)

    clf_score = sklearn.model_selection.cross_val_score(
                                            clf,
                                            X_train,
                                            y_train,
                                            n_jobs = threads,
                                            cv = 5
                                            ).mean()

    return clf_score


def train_clf(X_ref_df, y_ref_df, threads: int = -1, contam = False):
    X_train, X_test, y_train, y_test = \
        sklearn.model_selection.train_test_split(
                                    X_ref_df,
                                    y_ref_df,
                                    test_size = 0.3,
                                    random_state = 42,
                                    stratify = y_ref_df)

    optuna.logging.set_verbosity(optuna.logging.INFO)
    study = optuna.create_study(direction = 'maximize')
    study.optimize(lambda trial: objective(trial, X_train, y_train, threads, contam), n_trials = 50)

    trial = study.best_trial
    if contam:
        opt_clf = sklearn.naive_bayes.MultinomialNB(**trial.params)
    else:
        opt_clf = sklearn.ensemble.RandomForestClassifier(random_state = 42, **trial.params)

    opt_clf.fit(X_train, y_train)

    return opt_clf, trial


def save_clf(clf_model, txn_code: str) -> None:
    clf_pkl = f'{txn_code}_TIdeS/{txn_code}.TIdeS.pkl'
    pickle.dump(clf_model, open(clf_pkl, 'wb+'))


def process_predictions(
        txn_code: str,
        query_preds,
        query_df: pd.core.frame.DataFrame,
        contam: bool = False) -> dict:

    tds_dir = f'{txn_code}_TIdeS/Classified/'
    Path(tds_dir).mkdir(parents = True, exist_ok = True)

    pos, neg = 'CRF', 'OOF'
    if contam:
        pos, neg = 'Target', 'Non-Target'

    pos_pred_seqs = []
    contam_seqs = []

    tds_tsv = f'{tds_dir}{txn_code}.TIdeS_Summary.tsv'

    with open(tds_tsv, 'w+') as w:
        w.write(f'Gene\tTIdes_Prediction\tProb_{neg}\tProb_{pos}\n')
        for n in range(len(query_preds)):
            if argmax(query_preds[n]) == 0:
                w.write(f'{query_df.index[n]}\t{neg}\t'
                        f'{query_preds[n][0]*100:.2f}\t'
                        f'{query_preds[n][1]*100:.2f}\n')
                if contam:
                    contam_seqs.append(query_df.index[n])
            else:
                pos_pred_seqs.append(query_df.index[n])
                w.write(f'{query_df.index[n]}\t{pos}\t'
                        f'{query_preds[n][0]*100:.2f}\t'
                        f'{query_preds[n][1]*100:.2f}\n')

    if contam:
        return pos_pred_seqs, contam_seqs

    return best_porf(pos_pred_seqs, tds_tsv)


def best_porf(pos_pred_seqs, tds_tsv):
    final_pORFs = []

    from collections import defaultdict
    top_porfs = defaultdict(list)

    for line in open(tds_tsv).readlines():
        if '\tCRF\t' in line:
            txpt = line.split('\t')[0].split('.pORF')[0]
            top_porfs[txpt].append((line.split('\t')[0],float(line.split('\t')[-1])))

    for k, v in top_porfs.items():
        top_conf = max(v,key=lambda item:item[1])[-1]
        final_pORFs += [i[0] for i in v if i[1] == top_conf]

    return pos_pred_seqs, final_pORFs


def classify_orfs(
                txn_code: str,
                train_orfs: dict,
                query_orfs: dict,
                start_time,
                pretrained: str = None,
                threads: int = 4,
                verb: bool = True,
                contam: bool = False) -> dict:

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Extracting sequence features')

    query_df = extract_features(query_orfs, contam, False)

    if not pretrained:
        X_ref_df, y_ref_df = extract_features(train_orfs, contam, True)

        if verb:
            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Training TIdeS')

        # trained_rfc, finished_trial = train_rfc(X_ref_df, y_ref_df, threads)
        trained_clf, finished_trial = train_clf(X_ref_df, y_ref_df, threads, contam)
    else:
        trained_clf = pickle.load(open(pretrained, 'rb'))

    save_clf(trained_clf, txn_code)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Classifying ORFs')

    query_preds = trained_clf.predict_proba(query_df)

    if contam:
        target_seqs, contam_seqs =  process_predictions(txn_code, query_preds, query_df, contam)

        return target_seqs, contam_seqs

    else:
        return process_predictions(txn_code, query_preds, query_df)
