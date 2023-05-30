#!/usr/bin/env python3

"""Prepares folder(s) and filters the transcriptome by length and putative
rRNA sequences.

Note that rRNA filtering is performed by Barrnap. Taxa with unusual rRNAs, such
as Foraminifera, will need to have their rRNAs removed independently.

Dependencies include: Barrnap, CD-HIT."""


import shutil, subprocess, sys, time
from datetime import timedelta
from pathlib import Path

from Bio import SeqIO


def prep_dir(new_dir: str) -> None:
    Path(new_dir).mkdir(parents = True, exist_ok = True)


def filt_len(fasta_file: str,
            taxon_code: str,
            out_dir: str,
            min_len: int) -> str:
    """
    Remove sequences shorter than a minimum length

    Parameters
    ----------
    fasta_file: FASTA formatted file
    taxon_code: species/taxon name or abbreviated code
    out_dir: output directory to store filtered data
    min_len: minimum ORF/transcript length to consider

    Returns length-filtered FASTA formatted file
    """

    len_filt_fasta = f'{out_dir}/{taxon_code}.{min_len}bp.fas'

    with open(len_filt_fasta, 'w+') as w:
        t_cnt = 1
        for seq_rec in SeqIO.parse(fasta_file, 'fasta'):
            if len(seq_rec.seq) >= min_len:
                seq_name = f'>{taxon_code}_XX_Transcript_{t_cnt}_Len_{len(seq_rec.seq)}'

                if 'cov' in seq_rec.id:
                    cov_val = float(seq_rec.id.split("cov_")[1].split("_")[0])
                    seq_name += f'_Cov_{cov_val:.2f}'

                w.write(f'{seq_name}\n{seq_rec.seq}\n')
                t_cnt += 1

    return len_filt_fasta


def clust_txps(fasta_file: str,
                taxon_code: str,
                out_dir: str,
                pid: float,
                threads: int) -> str:
    """
    Remove redundant sequences using a given minium percent identity

    Parameters
    ----------
    fasta_file: FASTA formatted file
    taxon_code: species/taxon name or abbreviated code
    out_dir: output directory to store filtered data
    pid: minimum percent identity to cluster redundat sequences
    threads: number of cpu threads to use

    Returns a redundant sequence filtered FASTA formatted file
    """

    clust_fas = f'{out_dir}{fasta_file.split("/")[-1].split("Filt")[0]}{int(pid*100)}pid.Cluster.fas'

    cd_hit_cmd = f'cd-hit-est -T {threads} '\
                    '-G 0 ' \
                    f'-c {pid} ' \
                    '-aS 1.0 ' \
                    '-aL .005 ' \
                    f'-i {fasta_file} ' \
                    f'-o {clust_fas}'

    cdht_rslt = subprocess.run(cd_hit_cmd.split(),
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)

    return clust_fas


def run_barrnap(fasta_file: str,
                threads: int) -> list:
    """
    Run Barrnap to remove easily identifiable rRNAs

    Parameters
    ----------
    fasta_file: FASTA formatted file
    threads: number of cpu threads to use

    Returns a list of putative rRNA sequences
    """

    rRNA_seqs = []
    kdm = ['bac','arc','mito','euk']

    for k in kdm:
        bnp_cmd = ['barrnap', '--kingdom', k, '--threads', f'{threads}', fasta_file]

        bnp_rslt = subprocess.run(bnp_cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)

        rRNA_seqs += [i.split('\t')[0] for i in bnp_rslt.stdout.split('\n') if '_XX_' in i]

    return rRNA_seqs


def remove_rRNA(fasta_file: str,
                taxon_code: str,
                out_dir: str,
                min_len: int,
                threads: int) -> str:
    """
    Remove putative rRNA sequences from FASTA file

    Parameters
    ----------
    fasta_file: FASTA formatted file
    taxon_code: species/taxon name or abbreviated code
    out_dir: output directory to store filtered data
    min_len: minimum ORF/transcript length to consider
    threads: number of cpu threads to use

    Returns an rRNA-free FASTA formatted file
    """

    rRNA_seqs = run_barrnap(fasta_file, threads)

    rRNA_clean_fas = f'{out_dir}{taxon_code}.{min_len}bp.Filt_rRNA.fas'
    rRNA_fas = f'{out_dir}{taxon_code}.{min_len}bp.rRNA_Seqs.fas'

    rRNA_contam, clean_seqs = [],[]

    # bin sequences as either putative rRNAs or "clean"
    for i in SeqIO.parse(fasta_file, 'fasta'):
        if i.id in rRNA_seqs:
            rRNA_contam.append(f'>{i.id}\n{i.seq}')
        else:
            clean_seqs.append(f'>{i.id}\n{i.seq}')

    with open(rRNA_fas,'w+') as w:
        w.write('\n'.join(rRNA_contam))

    with open(rRNA_clean_fas,'w+') as w:
        w.write('\n'.join(clean_seqs))

    return rRNA_clean_fas


def prep_contam(fasta_file: str,
                taxon_code: str,
                sis_smry: str,
                pretrained: str,
                start_time,
                verb: bool = True) -> dict:

    """
    Prepares target/non-target sequences from a FASTA formatted file using a user
    provided summary table.

    Parameters
    ----------
    fasta_file: FASTA formatted file
    taxon_code: species/taxon name or abbreviated code
    sis_smry: tab-delimited file with user-defined "target" and "non-target" sequences
    pretrained: random forest model from previous TIdeS run
    start_time: initial timestamp to track runtime
    verb: verbose print statements

    Returns two lists based on the 'sis_smry' text file: "target" sequences and "non-target" sequences
    """

    backup_dir = f'{taxon_code}_TIdeS/Original/'
    eval_seq_fas = f'{backup_dir}{taxon_code}.EvalSeqs.fas'

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Backing up data for {taxon_code}')

    prep_dir(backup_dir)

    shutil.copy2(fasta_file, backup_dir)

    if pretrained:
        return None, None

    shutil.copy2(sis_smry, backup_dir)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Parsing target/non-target assignments')

    seq_summary = {i.split('\t')[0]:i.rstrip().split('\t')[1].lower() for i in open(sis_smry).readlines()}

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Saving target/non-target sequences')

    targ_seqs, non_targ_seqs = {}, {}

    # Create FASTA files for annotated 'target' and 'non-target' sequences for training
    with open(eval_seq_fas, 'w+') as w:
        for i in SeqIO.parse(fasta_file, 'fasta'):
            try:
                if seq_summary[i.id] in ['target', 'host', 'positive', '1']:
                    w.write(f'>{i.id}_Target\n{i.seq.upper()}\n')
                    targ_seqs[f'{i.id}_Target'] = f'{i.seq.upper()}'

                else:
                    w.write(f'>{i.id}_NonTarget\n{i.seq.upper()}\n')
                    non_targ_seqs[f'{i.id}_NonTarget'] = f'{i.seq.upper()}'

            except KeyError:
                continue

    if min([len(targ_seqs), len(non_targ_seqs)]) < 50:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Error: ' \
            'Fewer than 50 examples each of target and non-target sequences were assigned.')

        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Exiting TIdeS: Too Few Training Seqs.\n')
        sys.exit()

    elif 50 <= min([len(targ_seqs), len(non_targ_seqs)]) < 100:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Warning: ' \
            'Fewer than 100 target and 100 non-target sequences were assigned.')

    else:
        pass

    return targ_seqs, non_targ_seqs


def filter_transcripts(fasta_file: str,
                        taxon_code: str,
                        start_time,
                        min_len: int = 300,
                        threads: int = 4,
                        pid: float = 0.97,
                        verb: bool = True) -> str:
    """
    Performs all the initial filtering steps to ease classification

    Parameters
    ----------
    fasta_file: FASTA formatted file
    txn_code: species/taxon name or abbreviated code
    start_time: initial timestamp to track runtime
    min_len: minimum ORF/transcript length to consider
    threads: number of cpu threads to use
    pid: minimum percent identity to cluster redundat sequences
    verb: verbose print statements

    Returns the 'final' filtered FASTA formatted file
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
    mlen_fas = filt_len(fasta_file, taxon_code, filt_len_dir, min_len)


    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Removing rRNA contamination')

    prep_dir(filt_rRNA_dir)
    rRNA_free = remove_rRNA(mlen_fas, taxon_code, filt_rRNA_dir, min_len, threads)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Filtering redundant transcripts')

    prep_dir(clust_dir)
    post_clst_fas = clust_txps(rRNA_free, taxon_code, clust_dir, pid, threads)

    return post_clst_fas
