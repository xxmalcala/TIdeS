#!/usr/bin/env python3

"""
Prepares folder(s) and filters the transcriptome by length and putative
rRNA sequences.

Note that rRNA filtering is performed by Barrnap. Taxa with unusual rRNAs, such
as Foraminifera, will need to have their rRNAs removed independently.

Dependencies include: Barrnap, BioPython, CD-HIT.
"""

import shutil, subprocess, sys, time
from datetime import timedelta
from pathlib import Path

from Bio import SeqIO


def filt_len(fasta_file: str, taxon_code: str, out_dir: str, min_len: int, max_len: int) -> str:
    """
    Remove sequences shorter than a minimum length

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store filtered data
    min_len:     minimum ORF length to consider [default is 300nt]

    Returns
    ----------
    filt_len_fasta:  FASTA formatted file with sequences above a given minimum length
    """

    filt_len_fasta = f'{out_dir}/{taxon_code}.{min_len}bp.fas'

    t_cnt = 0
    updated_seqs = []

    for seq_rec in SeqIO.parse(fasta_file, 'fasta'):
        if min_len <= len(seq_rec.seq) <= max_len:
            t_cnt += 1

            updated_seq_name = f'{taxon_code}_XX_Transcript_{t_cnt}_Len_{len(seq_rec.seq)}'

            if 'cov' in seq_rec.id:
                cov_val = float(seq_rec.id.split("cov_")[1].split("_")[0])
                updated_seq_name += f'_Cov_{cov_val:.2f}'

            seq_rec.id = updated_seq_name
            seq_rec.description = ''
            seq_rec.seq = seq_rec.seq.upper()

            updated_seqs.append(seq_rec)

    SeqIO.write(updated_seqs, filt_len_fasta, 'fasta')

    return filt_len_fasta


def run_barrnap(fasta_file: str, threads: int) -> list:
    """
    Run Barrnap to remove easily identifiable rRNAs.

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    threads:     number of cpu threads to use

    Returns
    ----------
    rRNA_seqs:  list of putative rRNA sequence names
    """

    rRNA_seqs = []
    kdm = ['bac','arc','mito','euk']

    for k in kdm:
        bnp_cmd = ['barrnap', '--kingdom', k, '--threads', f'{threads}', fasta_file]

        bnp_rslt = subprocess.run(
                        bnp_cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True
                        )

        rRNA_seqs += [i.split('\t')[0] for i in bnp_rslt.stdout.split('\n') if '_XX_' in i]

    return rRNA_seqs


def remove_rRNA(fasta_file: str, taxon_code: str, out_dir: str, min_len: int, threads: int) -> str:
    """
    Remove putative rRNA sequences from FASTA file

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store filtered data
    min_len:     minimum ORF length to consider [default is 300nt]
    threads:     number of cpu threads to use

    Returns
    ----------
    rRNA_clean_fas:  FASTA formatted file without putative rRNA sequences
    """

    rRNA_seqs = run_barrnap(fasta_file, threads)

    rRNA_clean_fas = f'{out_dir}{taxon_code}.{min_len}bp.Filt_rRNA.fas'
    rRNA_fas = f'{out_dir}{taxon_code}.{min_len}bp.rRNA_Seqs.fas'

    rRNA_contam, clean_seqs = [],[]

    # bin sequences as either putative rRNAs or "clean"
    for i in SeqIO.parse(fasta_file, 'fasta'):
        if i.id in rRNA_seqs:
            rRNA_contam.append(i)
        else:
            clean_seqs.append(i)

    SeqIO.write(rRNA_contam, rRNA_fas, 'fasta')
    SeqIO.write(clean_seqs, rRNA_clean_fas, 'fasta')

    return rRNA_clean_fas


def clust_txps(fasta_file: str, taxon_code: str, out_dir: str, pid: float, threads: int) -> str:
    """
    Remove redundant sequences using a given minium percent identity

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store clustered transcripts
    pid:         minimum percent identity to cluster redundant sequences
    threads:     number of cpu threads to use

    Returns
    ----------
    clust_fas:  FASTA formatted file with redundant transcripts removed
    """

    clust_fas = f'{out_dir}{fasta_file.split("/")[-1].split("Filt")[0]}{int(pid*100)}pid.Cluster.fas'

    cd_hit_cmd = f'cd-hit-est -T {threads} '\
                    '-G 0 ' \
                    f'-c {pid} ' \
                    '-aS 1.0 ' \
                    '-aL .005 ' \
                    f'-i {fasta_file} ' \
                    f'-o {clust_fas}'

    cdht_rslt = subprocess.run(
                    cd_hit_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                    )

    return clust_fas


def prep_contam(fasta_file: str, taxon_code: str, sis_smry: str, model: str, start_time, verb: bool = True) -> dict:

    """
    Prepares target/non-target sequences from a FASTA formatted file using a user
    provided summary table.

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    sis_smry:    tab-delimited file with user-defined "target" and "non-target" sequences
    pretrained:  random forest model from previous TIdeS run
    start_time:  initial timestamp to track runtime
    verb:        verbose print statements

    Returns
    ----------
    targ_seqs:      dictionary of exemplar target sequences considered to capture
    non_targ_seqs:  dictionary of exemplar non-target sequences to avoid
    """

    backup_dir = f'{taxon_code}_TIdeS/Original/'
    eval_seq_fas = f'{backup_dir}{taxon_code}.EvalSeqs.fas'
    query_orfs = {}

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Backing up data for {taxon_code}')

    prep_dir(backup_dir)

    shutil.copy2(fasta_file, backup_dir)

    if not model:
        shutil.copy2(sis_smry, backup_dir)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Parsing ORF classifications')

    if not model:
        train_orfs = {i.split('\t')[0]:i.split('\t')[1].rstrip() for i in open(sis_smry).readlines()}
        contam_seqs = {i:0 for i in list(set(train_orfs.values()))}

    # Group sequences by user classification
    for i in SeqIO.parse(fasta_file, 'fasta'):
        if model:
            query_orfs[i.id] = f'{i.seq}'
        else:
            if i.id in train_orfs.keys():
                contam_seqs[train_orfs[i.id]] += 1
                train_orfs[i.id] = [train_orfs[i.id], f'{i.seq}']
            else:
                query_orfs[i.id] = f'{i.seq}'

    if model:
        return None, query_orfs

    for k, v in contam_seqs.items():
        if v < 50:
            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Error: ' \
                f'Fewer than 50 sequences were labeled as: {k}')

            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Exiting TIdeS: Too few training seqs.\n')
            # sys.exit()

        elif 50 <= v < 100:
                print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Warning: ' \
                    f'Fewer than 100 sequences were labeled as: {k}')

        else:
            pass

    print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Classification summary: ')
    for k, v in contam_seqs.items():
            print(f'{" "*11}----- {v} sequences were labeled as: {k}')

    return train_orfs, query_orfs


def prep_dir(new_dir: str) -> None:
    """
    Creates a new directory.

    Parameters
    ----------
    new_dir:  directory to create
    """

    Path(new_dir).mkdir(parents = True, exist_ok = True)


def filter_transcripts(fasta_file: str, taxon_code: str, start_time, min_len: int = 300, max_len: int = 10000, threads: int = 4, pid: float = 0.97, verb: bool = True) -> str:
    """
    Performs all the initial filtering steps to ease classification

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    start_time:  initial timestamp to track runtime
    min_len:     minimum ORF length to consider [default is 300nt]
    threads:     number of cpu threads to use
    pid:         minimum percent identity to cluster redundat sequences
    verb:        verbose print statements

    Returns
    ----------
    post_clst_fas:  FASTA formatted file of filtered transcripts for ORF-calling
    """


    backup_dir = f'{taxon_code}_TIdeS/Original/'
    filt_len_dir = f'{taxon_code}_TIdeS/Filter_Steps/Length_Filter/'
    clust_dir = f'{taxon_code}_TIdeS/Filter_Steps/Clustering/'
    filt_rRNA_dir = f'{taxon_code}_TIdeS/Filter_Steps/rRNA_Filter/'

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Backing up data for {taxon_code}')

    prep_dir(backup_dir)
    shutil.copy2(fasta_file, backup_dir)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Removing short transcripts')

    prep_dir(filt_len_dir)
    mlen_fas = filt_len(fasta_file, taxon_code, filt_len_dir, min_len, max_len)


    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Removing rRNA contamination')

    prep_dir(filt_rRNA_dir)
    rRNA_free = remove_rRNA(mlen_fas, taxon_code, filt_rRNA_dir, min_len, threads)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Filtering redundant transcripts')

    prep_dir(clust_dir)
    post_clst_fas = clust_txps(rRNA_free, taxon_code, clust_dir, pid, threads)

    return post_clst_fas
