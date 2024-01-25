#!/usr/bin/env python3

"""
Prepares folder(s) and filters the transcriptome by length and putative
rRNA sequences.

Note that rRNA filtering is performed by Barrnap. Taxa with unusual rRNAs, such
as Foraminifera, will need to have their rRNAs removed independently.

Dependencies include: Barrnap, BioPython, CD-HIT.
"""

import numpy, random, shutil, subprocess, sys, time
from datetime import timedelta
from pathlib import Path

from Bio import SeqIO
from ete3 import NCBITaxa


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
    filt_len_conv_tbl = f'{out_dir}/{taxon_code}.{min_len}bp.NameConversion.tsv'

    t_cnt = 0
    updated_seqs = []
    name_convert = []

    for i in SeqIO.parse(fasta_file, 'fasta'):
        if min_len <= len(i.seq) <= max_len:
            t_cnt += 1

            updated_seq_name = f'{taxon_code}|Transcript_{t_cnt}_Len_{len(i.seq)}'

            if 'cov' in i.id:
                cov_val = float(i.id.split("cov_")[1].split("_")[0])
                updated_seq_name += f'_Cov_{cov_val:.2f}'

            name_convert.append(f'{updated_seq_name}\t{i.description}\n')
            i.id = updated_seq_name
            i.description = ''
            i.seq = i.seq.upper()

            updated_seqs.append(i)

    SeqIO.write(updated_seqs, filt_len_fasta, 'fasta')

    with open(filt_len_conv_tbl,'w+') as w:
        w.write('Updated-Name\tOriginal-Name\n')
        w.write(''.join(name_convert))

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

        rRNA_seqs += [i.split('\t')[0] for i in bnp_rslt.stdout.split('\n') if '|' in i]

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


def run_kraken2(fasta_file: str, taxon_code: str, out_dir: str, kraken_db: str, threads: int) -> str:
    """
    Run Kraken2 to identify putative non-eukaryotic contamination

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store Kraken2's classification reports
    kraken_db:   path to kraken2-database directory
    threads:     number of cpu threads to use

    Returns
    ----------
    krak_out:   kraken2 taxonomic assignment report
    """

    krak_out = f'{out_dir}{taxon_code}.Kraken2_Txnmy.txt'

    krak_cmd = 'kraken2 --confidence 0.1 ' \
                f'--db {kraken_db} ' \
                f'--threads {threads} ' \
                f'--output {krak_out} ' \
                f'{fasta_file}'

    krak_rslt = subprocess.run(
                    krak_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                    )

    return krak_out


def assign_taxonomy(fasta_file: str, taxon_code: str, out_dir: str, kraken_db: str, threads: int) -> dict:
    """
    Runs Kraken2 and assigns taxonomy to the classified sequences.

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store Kraken2's classification reports
    kraken_db:   path to kraken2-database directory
    threads:     number of cpu threads to use

    Returns
    ----------
    krak_seq_txnmy: dictionary of sequences and their broad classification
    """

    kraken_output = run_kraken2(fasta_file, taxon_code, out_dir, kraken_db, threads)

    krak_seq_txnmy = {}

    ncbi = NCBITaxa()

    if not Path(f'{Path.home()}/.etetoolkit/taxa.sqlite').is_file():
        ncbi.update_taxonomy_database()

    with open(kraken_output, 'r') as fi:
        for ln in fi:
            if ln[0] == 'C':
                taxid = ln.split('\t')[2]
                txnmy = ncbi.get_taxid_translator(ncbi.get_lineage(taxid))
                if 'Bacteria' in txnmy.values():
                    krak_seq_txnmy[ln.split('\t')[1]] = 'Bacteria'
                elif 'Archaea' in txnmy.values():
                    krak_seq_txnmy[ln.split('\t')[1]] = 'Archaea'
                elif 'Viruses' in txnmy.values():
                    krak_seq_txnmy[ln.split('\t')[1]] = 'Virus'
                else:
                    krak_seq_txnmy[ln.split('\t')[1]] = 'PutativeEuk'

            else:
                krak_seq_txnmy[ln.split('\t')[1]] = 'PutativeEuk'

    return krak_seq_txnmy


def kraken_contam_eval(fasta_file: str, taxon_code: str, out_dir: str, kraken_db: str, threads: int) -> str:
    """
    Takes the taxonomic assignments to identify sequences for training TIdeS

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store Kraken2's classification reports and
                 sequence classification table
    kraken_db:   path to kraken2-database directory
    threads:     number of cpu threads to use

    Returns
    ----------
    seq_eval_summary: table of annotated sequences for training
    """


    seq_eval_summary = f'{out_dir}{taxon_code}.KrakenEval.TIdeS.txt'

    krak_seq_txnmy = assign_taxonomy(fasta_file, taxon_code, out_dir, kraken_db, threads)

    euk_seqs = []
    non_euk_seqs = []
    for k, v in krak_seq_txnmy.items():
        if v != 'PutativeEuk':
            non_euk_seqs.append(f'{k}\tNon-Eukaryotic\n')
        else:
            euk_seqs.append(f'{k}\t{taxon_code}\n')

    random.shuffle(non_euk_seqs)
    random.shuffle(euk_seqs)

    fin_data = non_euk_seqs[:min([len(non_euk_seqs), len(euk_seqs), 200])]
    fin_data += euk_seqs[:min([len(non_euk_seqs), len(euk_seqs), 200])]

    if not fin_data:
        print(f' ----- No sequences were annotated as "non-eukaryotic"')
        print(f' ----- Quitting TIdeS')
        sys.exit(1)

    with open(seq_eval_summary, 'w+') as w:
        w.write(''.join(fin_data))

    return seq_eval_summary


def remove_non_euk(fasta_file: str, taxon_code: str, out_dir: str, kraken_db: str, threads: int) -> str:
    """
    Remove putative non-eukaryotic contamination from FASTA file

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store annotated sequences (euk and non-euk)
    kraken_db:   path to kraken2-database directory
    threads:     number of cpu threads to use

    Returns
    ----------
    peuk_fasta': FASTA formatted file of putative eukaryotic transcripts
    """

    peuk_fasta = f'{out_dir}{taxon_code}.Kraken_Txnmy.PutativeEuk.fasta'

    krak_seq_txnmy = assign_taxonomy(fasta_file, taxon_code, out_dir, kraken_db, threads)

    non_euk_seqs = {'Bacteria':[],'Viruses':[], 'Archaea':[]}
    peuk_seqs = []

    for i in SeqIO.parse(fasta_file,'fasta'):
        if i.id not in krak_seq_txnmy:
            peuk_seqs.append(i)
        else:
            if krak_seq_txnmy[i.id] == 'PutativeEuk':
                peuk_seqs.append(i)
            else:
                non_euk_seqs[krak_seq_txnmy[i.id]].append(i)

    for k, v in non_euk_seqs.items():
        if v:
            outf = f'{out_dir}{taxon_code}.Kraken_Txnmy.{k}.fasta'
            SeqIO.write(v, outf, 'fasta')

    if not peuk_seqs:
        print(f'ERROR: No putative eukaryotic sequences are present for {taxon_code} from {fasta_file}')
        sys.exit(1)

    SeqIO.write(peuk_seqs, peuk_fasta, 'fasta')

    return peuk_fasta


def clust_txps(fasta_file: str, taxon_code: str, out_dir: str, pid: float, mem: int, threads: int) -> str:
    """
    Remove redundant sequences using a given minium percent identity

    Parameters
    ----------
    fasta_file:  FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    out_dir:     output directory to store clustered transcripts
    pid:         minimum percent identity to cluster redundant sequences
    mem:         memory limit (MB) for cd-hit (default: 800, unlimited: 0)
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
                    f'-M {int(mem)} ' \
                    f'-i {fasta_file} ' \
                    f'-o {clust_fas}'

    cdht_rslt = subprocess.run(
                    cd_hit_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                    )

    return clust_fas


def prep_contam(fasta_file: str, taxon_code: str, sis_smry: str, kraken_db: str, model: str, start_time, threads: int, verb: bool) -> dict:

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
        if kraken_db:
            if verb:
                print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Classifying ORFs with Kraken2')

            sis_smry = kraken_contam_eval(fasta_file, taxon_code, backup_dir, kraken_db, threads)
        else:
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

    if len(contam_seqs) < 2:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Error: ' \
            f'Missing at least TWO categories to run TIdeS')
        sys.exit(1)

    for k, v in contam_seqs.items():
        if v < 25:
            print(contam_seqs)
            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Error: ' \
                f'Fewer than 25 sequences were labeled as: {k}')

            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Exiting TIdeS: Too few training seqs.\n')
            sys.exit(1)

        elif 25 <= v < 100:
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


def filter_transcripts(fasta_file: str, taxon_code: str, start_time, kraken_db: str, skip_filter: bool = False, min_len: int = 300, max_len: int = 10000, threads: int = 4, pid: float = 0.97, mem: int = 2000, verb: bool = True) -> str:
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
    krak_dir = f'{taxon_code}_TIdeS/Filter_Steps/Kraken_Filter/'
    filt_rRNA_dir = f'{taxon_code}_TIdeS/Filter_Steps/rRNA_Filter/'

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Backing up data for {taxon_code}')

    prep_dir(backup_dir)
    shutil.copy2(fasta_file, backup_dir)

    if skip_filter:
        return f'{backup_dir}{fasta_file.rpartition("/")[-1]}'

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
    post_clst = clust_txps(rRNA_free, taxon_code, clust_dir, pid, mem, threads)

    if verb and kraken_db:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Removing non-eukaryotic transcripts')

        prep_dir(krak_dir)
        post_krak = remove_non_euk(post_clst, taxon_code, krak_dir, kraken_db, threads)

        return post_krak

    return post_clst
