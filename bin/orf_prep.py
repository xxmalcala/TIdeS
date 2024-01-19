#!/usr/bin/env python3

"""
Prepares query and training ORFs for subsequent classification.

Dependencies include: BioPython, Scikit-learn.
"""

import os, random, sys

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sklearn.feature_extraction.text import CountVectorizer


def trim_seq(seq: str, rf: int) -> str:
    """
    Trims sequences to ensure multiple of 3.

    Parameters
    ----------
    seq:  nucleotide sequence to trim
    rf:   sequence's reading frame

    Returns sequence with length divisble by 3.
    """

    temp = seq[rf-1:]

    return temp[:len(temp)-len(temp)%3]


def randomize_orientation(ref_orf_db: dict, taxon_code: str) -> dict:
    """
    Randomly changes the orientation in a subset of the training ORFs.

    Parameters
    ----------
    ref_orf_db:  dictionary of all the training ORFs
    taxon_code:  species/taxon name or abbreviated code

    Returns
    ----------
    rnd_ornt:  dictionary of all training ORFs and their classification (0/1)

    """

    rnd_ornt = {}
    rnd_ornt_seqs = []

    # Shuffle the dictionary keys
    seq_subsamp = random.sample(list(ref_orf_db.keys()), int(len(ref_orf_db)))


    for n in range(0, len(seq_subsamp)):
        ref_seq_name = f'{seq_subsamp[n]}_RF1'
        ref_seq = ref_orf_db[seq_subsamp[n]]


        rf = random.sample([-3,-2,-1, 2, 3], 1)[0]
        if rf > 0:
            rnd_seq_name = f'{seq_subsamp[n]}_RF{rf}'
            rnd_seq = trim_seq(ref_orf_db[seq_subsamp[n]], rf)

        else:
            rnd_seq_name = f'{seq_subsamp[n]}_RF{abs(rf) + 3}'
            rnd_seq = trim_seq(f'{Seq(ref_orf_db[seq_subsamp[n]]).reverse_complement()}', abs(rf))

        rnd_ornt[ref_seq_name] = [1, ref_seq]
        rnd_ornt[rnd_seq_name] = [0, rnd_seq]

        rnd_ornt_seqs.append(SeqRecord(Seq(ref_seq), ref_seq_name, '',''))
        rnd_ornt_seqs.append(SeqRecord(Seq(rnd_seq), rnd_seq_name, '',''))

        # Keep 50% of the training ORFs in the correct orientation
        # if n < len(seq_subsamp) * 0.5:
        #     rnd_seq_name = f'{seq_subsamp[n]}_RF1'
        #     rnd_seq = ref_orf_db[seq_subsamp[n]]
        #     rnd_ornt[rnd_seq_name] = [1, rnd_seq]
        #
        # else:
        #     # Randomly select "incorrect" reading frames
        #     rf = random.sample([-3,-2,-1, 2, 3], 1)[0]
        #
        #     if rf > 0:
        #         rnd_seq_name = f'{seq_subsamp[n]}_RF{rf}'
        #         rnd_seq = trim_seq(ref_orf_db[seq_subsamp[n]], rf)
        #
        #     else:
        #         rnd_seq_name = f'{seq_subsamp[n]}_RF{abs(rf) + 3}'
        #         rnd_seq = trim_seq(f'{Seq(ref_orf_db[seq_subsamp[n]]).reverse_complement()}', abs(rf))
        #
        #     rnd_ornt[rnd_seq_name] = [0, rnd_seq]
        #
        #     rnd_ornt_seqs.append(SeqRecord(Seq(rnd_seq), rnd_seq_name, '',''))

    train_fas = f'{taxon_code}_TIdeS/ORF_Calling/{taxon_code}.TrainingORFs.fas'

    SeqIO.write(rnd_ornt_seqs, train_fas, 'fasta')

    return rnd_ornt


def chunk_seq(seq: str, contam: bool, overlap: bool, kmer: int, step: int):
    """
    Chunks the sequences into appropriate "words" by kmer-length and with/out
    overlapping kmers with a given step size.

    Parameters
    ----------
    seq:      nucleotide sequence to trim
    contam:   training ORFs are for contamination calling
    overlap:  overlapping kmers permitted if True
    kmer:     kmer size (number of nt) to convert sequences into
    step:     nucleotide distance to consider for overlapping kmers

    Returns kmers of the given sequence with/out overlapping (and appropriate steps)
    """

    if overlap:
        if not step:
            step = int(kmer/2)
        kmer_list = [seq[n:n+kmer] for n in range(0, len(seq)-(kmer-1), step)]

    else:
        kmer_list = [seq[n:n+kmer] for n in range(0, len(seq), kmer)]

    return ' '.join(kmer_list)


def kmer_ngram_counts(train_orfs_dict: dict, query_orfs_dict: dict, taxon_code: str, partial: bool = False, cvec = None, contam: bool = False, overlap: bool = False, kmer: int = 3, step = None):
    """
    Generates the training and query arrays for subsequent classification.

    Parameters
    ----------
    train_orfs_dict:  dictionary of training ORFs
    query_orfs_dict:  dictionary of all putative ORFs (partials included)
    taxon_code:       species/taxon name or abbreviated code
    partial:          evaluate partial ORFs
    cvec:             pretrained CountVectorizer from prior run
    contam:           training ORFs are for contamination calling
    overlap:          overlapping kmers permitted if True
    kmer:             kmer size (number of nt) to convert sequences into
    step:             nucleotide distance to consider for overlapping kmers

    Returns
    ----------
    train_labels:  sequence names for the training dataset
    X_train_orfs:  classification-ready data for prepared training ORFs
    train_class:   labels (0/1) for the training dataset (X_train_orfs)
    query_labels:  sequence names for the query pORF dataset
    X_query_orfs:  classification-ready data for prepared query pORFs
    cvec:          trained CountVectorizer
    """

    train_labels, train_seqs, train_class = [],[],[]
    query_labels, query_seqs = [],[]

    for k, v in query_orfs_dict.items():
        if not contam:
            if not partial and 'orf_type:complete' in k:
                    query_seqs.append(chunk_seq(v, contam, overlap, kmer, step))
                    query_labels.append(k)
            else:
                query_seqs.append(chunk_seq(v, contam, overlap, kmer, step))
                query_labels.append(k)
        else:
            query_seqs.append(chunk_seq(v, contam, overlap, kmer, step))
            query_labels.append(k)

    if cvec:
        X_query_orfs = np.array([i/sum(i) for i in cvec.transform(query_seqs).toarray()])
        return None, (query_labels, X_query_orfs), cvec

    if not contam:
        train_orfs_dict = randomize_orientation(train_orfs_dict, taxon_code)

    for k, v in train_orfs_dict.items():
        train_seqs.append(chunk_seq(v[1], contam, overlap, kmer, step))
        train_labels.append(k)
        train_class.append(v[0])

    cvec = CountVectorizer(ngram_range = (1,1))

    X_query_orfs = np.array([i/sum(i) for i in cvec.fit_transform(query_seqs).toarray()])
    X_train_orfs = np.array([i/sum(i) for i in cvec.transform(train_seqs).toarray()])

    return (train_labels, X_train_orfs, train_class), (query_labels, X_query_orfs), cvec
