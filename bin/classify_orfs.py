#!/usr/bin/env python3

import pickle, time
import pandas as pd

from datetime import timedelta
from numpy import argmax
from pathlib import Path

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import (train_test_split, GridSearchCV)


def extract_features(feature_dict: dict, contam: bool, training: bool):
    feature_df = pd.DataFrame.from_dict(feature_dict, orient = 'index')
    if training:
        if contam:
            y_ref_df = feature_df.index.map(lambda x: 1 if x.split('_')[-1] == 'Target' else 0)

        else:
            y_ref_df = feature_df.index.map(lambda x: 1 if x[-1] == '1' else 0)

        return feature_df, y_ref_df

    return feature_df



def train_rfc(X_ref_df, y_ref_df, threads: int):
    X_train, X_test, y_train, y_test = train_test_split(X_ref_df,
                                                        y_ref_df,
                                                        test_size = 0.2,
                                                        random_state = 42,
                                                        stratify = y_ref_df)
    param_grid = {
        'n_estimators': [300, 500],
        'max_features': ['sqrt', 'log2'],
        'min_samples_split': [2, 5],
        'max_depth' : [4, 8, 12, 16],
        'criterion' :['gini', 'entropy'],
        'oob_score': [True],
        'n_jobs': [int(threads)]}

    CV_rfc = GridSearchCV(
        estimator = RandomForestClassifier(random_state = 42),
        param_grid = param_grid,
        cv = 5,
        n_jobs = int(threads))

    CV_rfc.fit(X_train, y_train)

    rfc = RandomForestClassifier(random_state = 42, **CV_rfc.best_params_)

    return rfc.fit(X_train, y_train)


def save_rfc(rfc_model, txn_code):
    rfc_pkl = f'{txn_code}_TIdeS/{txn_code}.TIdeS.pkl'
    pickle.dump(rfc_model, open(rfc_pkl, 'wb+'))


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

        trained_rfc = train_rfc(X_ref_df, y_ref_df, threads)

    else:
        trained_rfc = pickle.load(open(pretrained, 'rb'))

    save_rfc(trained_rfc, txn_code)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Classifying ORFs')

    query_preds = trained_rfc.predict_proba(query_df)

    if contam:
        target_seqs, contam_seqs =  process_predictions(txn_code, query_preds, query_df, contam)

        return target_seqs, contam_seqs

    else:
        return process_predictions(txn_code, query_preds, query_df)
