import sys
import pandas as pd
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import (train_test_split,
                                    GridSearchCV)

"""Uses Sci-Kit Learn's Random Forest Classifier to identify the best putative
ORFs (currently from transcriptomes) based on similarity to ORFs with
reasonable hits to a diverse eukaryotic protein database.

This is intended to provide a more robust means to analyze young lineage-specific
genes derived from traditionally "noisy" transcriptomic (RNA-seq) data."""


def parse_known_ORF_table(train_orf_tsv, contam = False):
    orf_df = pd.read_table(train_orf_tsv, header = 0)

    orf_df = orf_df.sample(frac=1)
    # orf_df = orf_df.sort_values(orf_df.columns[1])

    # Grab the appropriate data from dataframe (col-1 is the categorization)
    orf_rfs = orf_df[orf_df.columns[1]]
    orf_cdns = orf_df.loc[:,orf_df.columns[2]:]

    # Put a check in for smallest category for RFC_Model
    min_y_cnts = min(orf_rfs.value_counts())
    if min_y_cnts < 100:
        print('WARNING: Too few categorized data to train confidently.'
            ' Results likely to be biased.')

    # return orf_cdns[:500], orf_rfs[:500]
    return orf_cdns, orf_rfs


def train_rfc(train_orf_tsv):
    x, y = parse_known_ORF_table(train_orf_tsv)

    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.2, random_state = 42)

    # Parameters to try with GridSearchCV.
    param_grid = {
        'n_estimators': [200, 500, 1000],
        'max_features': ['sqrt', 0.2, 0.3, 0.4],
        'min_samples_split': [5],
        'max_depth' : [4,6,8,10],
        'criterion' :['gini', 'entropy'],
        'n_jobs': [-1]}

    # Might change to nested-gridsearch (usure)
    CV_rfc = GridSearchCV(
        estimator = RandomForestClassifier(random_state = 42),
        param_grid = param_grid,
        cv = 5,
        n_jobs = -1)

    CV_rfc.fit(x_train, y_train)
    print(CV_rfc.best_params_)

    rfc = RandomForestClassifier(random_state = 42, **CV_rfc.best_params_)
    rfc.fit(x_train, y_train)

    return rfc


def rfc_classify_query(query_orf_tsv, taxon_code, taxon_dir, rfc, RF = True):
    """
    Classifies the query data into categories using the trained RFC.

    :param q_df: pORF pandas codon count dataframe.
    :param taxon_code: Taxon name (for naming files).
    :param taxon_dir: Current TIdeS working directory.
    :param rfc: Trained random forest classifier.
    :param rep_num: Current TIdeS round.
    :param RF: Classification mode for pORFs (True) or contamination (False).
    :return tides_out_dir: Directory for output TIdeS predictions.
    :return rfc_call_tsv: TSV file of TIdeS's predictions.
    """

    tides_out_dir = f'{taxon_dir}TIdeS_Predicted/'
    Path(tides_out_dir).mkdir(parents = True, exist_ok = True)

    best_q_preds = defaultdict(list)

    q_df = pd.read_table(query_orf_tsv)
    q_genes = q_df['Gene']
    q_cdns = q_df.loc[:,'GGG':]

    q_preds = rfc.predict_proba(q_cdns)
    q_preds_brief = rfc.predict(q_cdns)

    ### Need to link all pORFs with their probs, not just "positives"
    positive_pred_pos = [p for p, c in enumerate(q_preds_brief) if c == 1]

    q_pred_eval = q_preds[positive_pred_pos]
    q_genes_eval = q_genes[positive_pred_pos]

    preds_to_eval = [(i,j) for i,j in zip(q_genes_eval, q_pred_eval)]

    for i in preds_to_eval:
        best_q_preds[i[0].split("_Len")[0]].append(i)

        if RF:
            rfc_call_tsv = f'{tides_out_dir}{taxon_code}.ORF_RF_Call_ALL.RFC.tsv'
        else:
            rfc_call_tsv = f'{tides_out_dir}{taxon_code}.NonContam.RFC.tsv'

    with open(rfc_call_tsv,'w+') as w:
        w.write('Gene\tRFC_Call\tORF_Prob_Cat_0\tORF_Prob_Cat_1\n')
        for k, v in best_q_preds.items():
            g = max(v, key=lambda x: x[-1][-1])
            if g[-1][-1] > g[-1][0]:
                ev = '1'
            else:
                ev = '0'
            w.write(f'{g[0]}\t{ev}\t{g[1][0]:.3f}\t{g[1][-1]:.3f}\n')

    return tides_out_dir, rfc_call_tsv


def rfc_filter_orfs(fasta_file, taxon_code, tides_out_dir, rfc_call_tsv, rfc = None, RF = True):
    """
    Filters pORFs/transcripts from a single FASTA file based on TIdeS predictions.

    :param fasta_file: FASTA file of pORFs.
    :param taxon_code: Taxon name (for naming files).
    :param tides_out_dir: TIdeS directory for the taxon.
    :param rep_num: Current TIdeS round.
    :param rfc_call_tsv: TSV-file of TIdeS predictions.
    :param RF: Classification mode for pORFs (True) or contamination (False).
    :return: TIdeS predicted ORFs (FASTA format).
    """

    rfc_seqs = [i.split('\t')[0] for i in open(rfc_call_tsv).readlines()[1:]]

    rfc_filt_fas = f'{tides_out_dir}{taxon_code}.TIdeS_Pred.fas'

    with open(rfc_filt_fas, 'w+') as w:
        for i in SeqIO.parse(fasta_file,'fasta'):
            if i.id in rfc_seqs:
                w.write(f'>{i.description}\n{i.seq}\n')

    if not RF:
        contam_fas = f'{tides_out_dir}{taxon_code}.TIdeS_Pred.Contam.fas'
        with open(contam_fas, 'w+') as w:
            for i in SeqIO.parse(fasta_file,'fasta'):
                if i.id not in rfc_seqs:
                    w.write(f'>{i.description}\n{i.seq}\n')

    return rfc_filt_fas


def save_rfc_model(taxon_dir, rfc):
    import pickle
    model_file = f'{taxon_dir}{taxon_dir.rstrip("/")}.RFC_Model.pkl'
    with open(model_file, 'wb') as rfcm:
        pickle.dump(rfc, rfcm)


def classify_seqs(train_orf_tsv, query_orf_tsv, fasta_file, taxon_dir, taxon_code, rfc, train = False, RF = True):
    """
    Calls on fucntions for classifying pORFs using Sci-Kit Learn's Random Forest Classifier.

    :param train_orf_tsv: training codon count table.
    :param query_orf_tsv: pORF codon count table.
    :param fasta_file: FASTA file of pORFs.
    :param taxon_dir: Current TIdeS working directory.
    :param taxon_code: Taxon name (for naming files).
    :param RF: Classification mode for pORFs (True) or contamination (False).
    :return: TIdeS predicted ORFs (FASTA format).
    """

    if not rfc:
        print('Training TIdeS.')
        trained_rfc = train_rfc(train_orf_tsv)
    else:
        import pickle
        trained_rfc = pickle.load(open(rfc, 'rb'))

    save_rfc_model(taxon_dir, trained_rfc)

    if not train:
        print('Classifying pORFs.')
        tds_dir, rfc_pred_tsv = rfc_classify_query(
                                    query_orf_tsv,
                                    taxon_code,
                                    taxon_dir,
                                    trained_rfc,
                                    RF)

        rfc_fas = rfc_filter_orfs(
                    fasta_file,
                    taxon_code,
                    tds_dir,
                    rfc_pred_tsv,
                    RF)

        return rfc_fas

    else:
        return None
