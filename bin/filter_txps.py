#!/usr/bin/env python3

import shutil, subprocess, time
from datetime import timedelta
from pathlib import Path
from Bio import SeqIO

"""Script prepares folder(s) and filters the transcriptome by length and putative
rRNA sequences.

Note that rRNA filtering is performed by Barrnap. Taxa with unusual rRNAs, such
as Foraminifera, will need to have their rRNAs removed independently.

Dependencies include: Barrnap, CD-HIT."""


def prep_dir(new_dir: str) -> None:
    Path(new_dir).mkdir(parents = True, exist_ok = True)


def filt_len(fasta_file: str, taxon_code: str, out_dir: str, min_len: int) -> str:
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


def clust_txps(fasta_file: str, taxon_code: str, out_dir: str, pid: float, threads: int) -> str:
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


def run_barrnap(fasta_file: str, threads: int) -> list:
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


def remove_rRNA(fasta_file: str, taxon_code: str, out_dir: str, min_len: int, threads: int) -> str:
    rRNA_seqs = run_barrnap(fasta_file, threads)

    rRNA_clean_fas = f'{out_dir}{taxon_code}.{min_len}bp.Filt_rRNA.fas'
    rRNA_fas = f'{out_dir}{taxon_code}.{min_len}bp.rRNA_Seqs.fas'

    rRNA_contam, clean_seqs = [],[]

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


def filter_transcripts(txm_fas: str, txn_code: str, start_time, min_len: int = 300, threads: int = 4, pid: float = 0.97, verb: bool = True):
    backup_dir = f'{txn_code}_TIdeS/Original/'
    filt_len_dir = f'{txn_code}_TIdeS/Filter_Steps/Length_Filter/'
    clust_dir = f'{txn_code}_TIdeS/Filter_Steps/Clustering/'
    filt_rRNA_dir = f'{txn_code}_TIdeS/Filter_Steps/rRNA_Filter/'

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Backing up data for {txn_code}')

    prep_dir(backup_dir)
    shutil.copy2(txm_fas, backup_dir)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Removing short transcripts')

    prep_dir(filt_len_dir)
    mlen_fas = filt_len(txm_fas, txn_code, filt_len_dir, min_len)


    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Removing rRNA contamination')

    prep_dir(filt_rRNA_dir)
    rRNA_free = remove_rRNA(mlen_fas, txn_code, filt_rRNA_dir, min_len, threads)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Filtering redundant transcripts')

    prep_dir(clust_dir)
    post_clst_fas = clust_txps(rRNA_free, txn_code, clust_dir, pid, threads)

    return post_clst_fas
