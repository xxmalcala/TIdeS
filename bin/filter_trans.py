#!/usr/bin/env python3

import glob, shutil, subprocess, sys
from pathlib import Path
from Bio import SeqIO

"""Script prepares folder(s) and filters the transcriptome by length and putative
rRNA sequences.

Note that rRNA filtering is performed by Barrnap. Taxa with unusual rRNAs, such
as Foraminifera, will need to have their rRNAs removed independently.

Dependencies include: Barrnap, CD-HIT."""


def filter_length(transcriptome, taxon_code, out_dir, min_length = 300):
    filt_dir = f'{out_dir}/Filter_Steps/Length_Filter/'
    Path(filt_dir).mkdir(exist_ok = True, parents = True)
    len_filt_fasta = f'{filt_dir}/{taxon_code}.{min_length}bp.fas'

    with open(len_filt_fasta, 'w+') as w:
        t_count = 1

        for seq_rec in SeqIO.parse(transcriptome, 'fasta'):
            if len(seq_rec.seq) >= min_length:
                w.write(f'>{taxon_code}_XX_Transcript_{t_count}_Len_{len(seq_rec.seq)}' \
                        f'\n{seq_rec.seq}\n')
                t_count += 1

    return len_filt_fasta


def clust_transcripts(len_filt_fasta, taxon_code, out_dir, min_length, threads = 4):
    clust_dir = f'{out_dir}/Filter_Steps/Clustering/'
    Path(clust_dir).mkdir(exist_ok = True, parents = True)
    clust_fas = f'{clust_dir}{taxon_code}.{min_length}bp.Cluster.fas'

    cd_hit_cmd = f'cd-hit-est -T {threads} '\
                    '-G 0 ' \
                    '-c 0.97 ' \
                    '-aS 1.0 ' \
                    '-aL .005 ' \
                    f'-i {len_filt_fasta} ' \
                    f'-o {clust_fas}'

    cd_hit_call = subprocess.call(cd_hit_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return clust_fas

def identify_rRNA(clust_fas, taxon_code, out_dir, threads = 4):
    rRNA_dir = f'{out_dir}Filter_Steps/rRNA_Filter/'
    Path(rRNA_dir).mkdir(exist_ok = True, parents = True)

    rRNA_db = ['bac','arc','euk','mito']

    for ribo_type in rRNA_db:
        rRNA_fas = f'{rRNA_dir}/{taxon_code}.rRNA_{ribo_type}.fas'
        barrnap_cmd = f'barrnap --threads {threads} '\
                        f'--kingdom {ribo_type} ' \
                        f'--quiet ' \
                        f'--outseq {rRNA_fas} ' \
                        f'{clust_fas}'

        barrnap_call = subprocess.call(barrnap_cmd, shell=True,
            stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return rRNA_dir


def remove_rRNA(len_filt_fasta, taxon_code, rRNA_dir, min_length):
    filt_dir = rRNA_dir.rstrip('rRNA_Filter/')
    rRNA_filt_fasta = f'{filt_dir}/{taxon_code}.{min_length}bp.Filt_rRNA.fas'
    rRNA_seqs = []

    for fas in glob.glob(f'{rRNA_dir}*rRNA*fas'):
        for seq_rec in SeqIO.parse(f'{fas}','fasta'):
            rRNA_seqs.append(seq_rec.id.split(':')[2])

    with open(rRNA_filt_fasta, 'w+') as w:
        trxpt_num = 1
        for seq_rec in SeqIO.parse(len_filt_fasta,'fasta'):
            if seq_rec.id not in set(rRNA_seqs):

                if 'cov' in seq_rec.id.lower():
                    cov = seq_rec.id.lower().split('cov_')[1].split('.')[0]
                    new_name = f'>{taxon_code}_XX_Transcript_{trxpt_num}_Len_' \
                            f'{len(seq_rec.seq)}_Cov_{cov}'

                else:
                    new_name = f'>{taxon_code}_XX_Transcript_{trxpt_num}_Len_' \
                            f'{len(seq_rec.seq)}'

                w.write(f'{new_name}\n{seq_rec.seq}\n')
                trxpt_num += 1

    return rRNA_filt_fasta


def filter_txpts(fasta_file, taxon_code, min_length, threads):
    out_dir = f'{taxon_code}_TIdeS/'
    Path(f'{out_dir}Original/').mkdir(exist_ok = True, parents = True)
    shutil.copy2(fasta_file, f'{out_dir}Original/')

    mlen_fas = filter_length(
                fasta_file,
                taxon_code,
                out_dir,
                min_length)

    cl_fas = clust_transcripts(
                mlen_fas,
                taxon_code,
                out_dir,
                min_length,
                threads)

    rRNA_dir = identify_rRNA(
                cl_fas,
                taxon_code,
                out_dir,
                threads)

    filt_trxme = remove_rRNA(
                    cl_fas,
                    taxon_code,
                    rRNA_dir,
                    min_length)
                    
    return filt_trxme, out_dir

# Update with optional ARGPARSE!
if __name__ == '__main__':
    if len(sys.argv[1:]) >= 2:
        fasta_file = sys.argv[1]
        taxon_code = sys.argv[2]
        try:
            min_length = int(sys.argv[3])
        except IndexError:
            min_length = 300
        try:
            threads = int(sys.argv[4])
        except IndexError:
            threads = 4
    else:
        print('Usage:\n    python3 filter_trans.py [FASTA-FILE] [TAXON-CODE] '
                '[MIN-TRANSCRIPT-LENGTH (default = 300bp)]\n')
        sys.exit(1)

    filter_txpts(fasta_file, taxon_code, min_length = 300, threads = 4)
    # print(filt_trxme)
