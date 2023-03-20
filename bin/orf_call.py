#!/usr/bin/env python3

import re, subprocess, sys, time
import numpy as np

from datetime import timedelta
from itertools import product
from pathlib import Path
from random import sample

from Bio import SeqIO
from Bio.Seq import Seq


def eval_gcode_ttable(translation_table: str) -> list:
    ncbi_link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi'

    supported_tables = {'1': ['TGA','TAA','TAG','TCA','TTA','CTA',1],
                        'universal': ['TGA','TAA','TAG','TCA','TTA','CTA',1],
                        '4': ['TAA','TAG','TTA','CTA',4],
                        'blepharisma': ['TAA','TAG','TTA','CTA',4],
                        '6': ['TGA','TCA',6],
                        'ciliate': ['TGA','TCA',6],
                        '10': ['TAA','TAG','TTA','CTA',10],
                        'euplotes': ['TAA','TAG','TTA','CTA',10],
                        '12': ['TGA','TAA','TAG','TCA','TTA','CTA',12],
                        'alt-yeast': ['TGA','TAA','TAG','TCA','TTA','CTA',12],
                        '26': ['TGA','TAA','TAG','TCA','TTA','CTA',26],
                        'pacysolen': ['TGA','TAA','TAG','TCA','TTA','CTA',26],
                        '29': ['TGA','TCA',29],
                        'mesodinium': ['TGA','TCA',29],
                        '30': ['TGA','TCA',30],
                        'peritrch': ['TGA','TCA',30]}



    if translation_table not in supported_tables.keys():
        tmp_tbl = list(supported_tables.keys())
        sppt_tbls = [', '.join(i) for i in list(zip(tmp_tbl[0::2], tmp_tbl[1::2]))]
        print('\nCheck the chosen genetic code/translation table.\n\n' \
            'Supported tables include:\n    ' +
            '\n    '.join(sppt_tbls))

        print(f'\nFor more information, visit: {ncbi_link}\n')
        return None

    else:
        return supported_tables[f'{translation_table}']


def bin_pos(cdn_pos: list, rf: dict, std: str) -> dict:
    rf_adj = 0
    if std == '-':
        rf_adj += 3
    for i in cdn_pos:
        if i % 3 == 0:
            rf[1 + rf_adj].append(i)
        elif i % 3 == 1:
            rf[2 + rf_adj].append(i)
        else:
            rf[3 + rf_adj].append(i)
    return rf


def eval_coords(strt_coords: np.ndarray, stp_coords:  np.ndarray, rframe: int, min_len: int) -> list:
    porf_coords = []

    if len(stp_coords) == 2:
        try:
            if rframe > 3:
                new_start = strt_coords[strt_coords > stp_coords[1]].max()
            else:
                new_start = strt_coords[strt_coords < stp_coords[1]].min()
            if abs(new_start - stp_coords[1]) >= min_len:
                porf_coords.append((new_start, stp_coords[1]))

        except ValueError:
            pass

    else:
        valid_stops = [(stp_coords[i], stp_coords[i+1]) for i, j
                        in enumerate(stp_coords[1:] - stp_coords[:-1]) if j >= min_len]

        for pos in valid_stops:
            try:
                if rframe > 3:
                    tmp_start = [i for i in strt_coords[strt_coords <= pos[1]] if (i - pos[0]) >= min_len]
                    porf_coords.append((max(tmp_start), pos[0]))
                    porf_coords += [(tmp_start[n+1], pos[0]) for n in range(len(tmp_start)-1)
                                    if (tmp_start[n+1] - tmp_start[n]) > 14
                                    and (pos[1] - tmp_start[n+1]) >= min_len]

                else:
                    tmp_start = [i for i in strt_coords[strt_coords >= pos[0]] if (pos[1] - i) >= min_len]
                    porf_coords.append((min(tmp_start), pos[1]))
                    porf_coords += [(tmp_start[n+1], pos[1]) for n in range(len(tmp_start)-1)
                                    if (tmp_start[n+1] - tmp_start[n]) > 14
                                    and (pos[1] - tmp_start[n+1]) >= min_len]

            except ValueError:
                pass

    return porf_coords


def get_start_stop_codons(seq: str, stp_cdns: list, min_len: int) -> dict:
    strt_rf = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[]}
    stp_rf = {1:[0], 2:[0], 3:[0], 4:[len(seq)], 5:[len(seq)], 6:[len(seq)]}

    stp_cdns_p = []
    stp_cdns_n = []

    for n in range(len(stp_cdns)):
        if n < len(stp_cdns)/2:
            stp_cdns_p += [m.end() for m in re.finditer(stp_cdns[n], seq)]
        else:
            stp_cdns_n += [m.start() for m in re.finditer(stp_cdns[n], seq)]

    atg_cdns_p = [m.start() for m in re.finditer('ATG', seq)]
    atg_cdns_n = [m.end() for m in re.finditer('CAT', seq)]

    stp_rf = bin_pos(stp_cdns_p, stp_rf, '+')
    stp_rf = bin_pos(stp_cdns_n, stp_rf, '-')

    strt_rf = bin_pos(atg_cdns_p, strt_rf, '+')
    strt_rf = bin_pos(atg_cdns_n, strt_rf, '-')

    return strt_rf, stp_rf


def grab_inframe_orfs(seq: str, porf_coords: list, strd: str) -> list:
    inframe_orfs = []

    for pos in porf_coords:
        if strd == '+':
            orf = seq[pos[0]:pos[1]]
            inframe_orfs.append(f':{pos[0]+1}-{pos[1]+1}({strd})\n{orf}\n')

        else:
            orf = f'{Seq(seq[pos[1]:pos[0]]).reverse_complement()}'
            inframe_orfs.append(f':{pos[0]+1}-{pos[1]+1}({strd})\n{orf}\n')

    return inframe_orfs


def extract_orfs(seq: str, stp_cdns: list, min_len: int, strand: str) -> list:
    strt_rf, stp_rf = get_start_stop_codons(seq, stp_cdns, min_len)

    final_pORFs = []

    if strand == 'plus':
        del stp_rf[4]
        del stp_rf[5]
        del stp_rf[6]

    elif strand == 'minus':
        del stp_rf[1]
        del stp_rf[2]
        del stp_rf[3]
    else:
        pass

    for k in stp_rf.keys():
        stp_coords = list(set(stp_rf[k]))
        strt_coords = strt_rf[k]

        if stp_coords and strt_coords:
            stp_coords.sort()
            stp_coords = np.array(stp_coords)
            strt_coords = np.array(strt_coords)

            porf_coords = eval_coords(strt_coords, stp_coords, k, min_len)

            if porf_coords:
                strd = '+' if k < 4 else '-'
                final_pORFs += grab_inframe_orfs(seq, porf_coords, strd)

    return final_pORFs


def translate_seqs(seq_dict: dict, gcode: str) -> dict:
    ttable = eval_gcode_ttable(f'{gcode}'.lower())[-1]
    return {k:f'{Seq(v[:-3]).translate(ttable)}' for k, v in seq_dict.items()}


def call_all_orfs(fasta_file: str, txn_code: str, gcode: str, min_len: int, strand: str) -> str:
    orf_fas = f'{txn_code}_TIdeS/ORF_Calling/{txn_code}.{min_len}bp.CompORFs.fas'

    stp_cdns = eval_gcode_ttable(str(gcode).lower())
    orf_calls = {}

    with open(orf_fas, 'w+') as w:
        for seq_rec in SeqIO.parse(fasta_file, 'fasta'):
            orf_num = 1
            seq_porfs = extract_orfs(f'{seq_rec.seq}', stp_cdns[:-1], min_len, strand)

            if seq_porfs:
                for porf in seq_porfs:
                    orf_calls[f'{seq_rec.id}.pORF{orf_num}'] = porf.split('\n')[1]
                    w.write(f'>{seq_rec.id}.pORF{orf_num} {seq_rec.id}{porf}')
                    orf_num += 1

    return orf_calls, orf_fas


def trim_seq(seq: str, rf: int) -> str:
    temp = seq[rf-1:]
    return temp[:len(temp)-len(temp)%3]


def ref_orf_dmnd(
                query: str,
                txn_code: str,
                dmnd_db: str,
                min_len: int,
                evalue: float,
                threads: int,
                intermed: bool) -> list:

    outfmt = '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe'

    if intermed:
        tsv_out = f'{txn_code}_TIdeS/ORF_Calling/{txn_code}.{min_len}bp.DMND_REF_ORFs.tsv'
        outfmt = f'{outfmt} -o {tsv_out}'

    dmnd_cmd = f'diamond blastx ' \
                f'-p {threads} ' \
                f'-e {evalue} ' \
                f'-q {query} ' \
                f'-d {dmnd_db} ' \
                f'-k 1 ' \
                f'-f {outfmt}'

    dmnd_rslt = subprocess.run(dmnd_cmd.split(),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)
    if intermed:
        return [i.rstrip() for i in open(outfmt.split()[-1]).readlines()]

    else:
        return dmnd_rslt.stdout.split('\n')[:-1]


def generate_ref_orfs(
                    fasta_file: str,
                    txn_code: str,
                    start_time,
                    dmnd_db: str,
                    min_len: int,
                    evalue: float,
                    threads: int,
                    verb: bool ) -> dict:

    ref_orf_fas = f'{txn_code}_TIdeS/ORF_Calling/{txn_code}.{min_len}bp.DMND_REF_ORFs.fas'

    ref_orfs = {}

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Running DIAMOND BLASTX against protein database')

    dmnd_hits = {i.split('\t')[0]:i.split('\t') for i in
                    ref_orf_dmnd(fasta_file, txn_code, dmnd_db, min_len, evalue, threads, True)}

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Generating training ORFs for {txn_code}')

    for seq_rec in SeqIO.parse(fasta_file,'fasta'):
        if seq_rec.id in dmnd_hits:
            orf_pos = [int(j) for j in dmnd_hits[seq_rec.id][4:6]]
            orf = seq_rec.seq[min(orf_pos)-1:max(orf_pos)]

            if orf_pos[0] > orf_pos[1]:
                ref_orfs[seq_rec.id] = f'{Seq(orf).reverse_complement()}'

            else:
                ref_orfs[seq_rec.id] = f'{orf}'

    rnd_orf_orient = randomize_orientation(ref_orfs)

    with open(ref_orf_fas, 'w+') as w:
        for k, v in ref_orfs.items():
            w.write(f'>{k}\n{v}\n')

    with open(ref_orf_fas.replace("fas","RandomRF.fas"), 'w+') as w:
        for k, v in rnd_orf_orient.items():
            w.write(f'>{k}\n{v}\n')

    return ref_orfs, rnd_orf_orient


def randomize_orientation(ref_orf_dict: dict) -> dict:
    seq_subsamp = sample(list(ref_orf_dict.keys()), int(len(ref_orf_dict)))
    rnd_ornt = {}

    for n in range(0, len(seq_subsamp)):
        if n < len(seq_subsamp) * 0.5:
            rnd_ornt[f'{seq_subsamp[n]}_RF1'] = ref_orf_dict[seq_subsamp[n]]

        else:
            rf = sample([-3,-2,-1, 2, 3], 1)[0]
            if rf > 0:
                rnd_ornt[f'{seq_subsamp[n]}_RF{rf}'] = trim_seq(ref_orf_dict[seq_subsamp[n]], rf)

            else:
                temp_seq = f'{Seq(ref_orf_dict[seq_subsamp[n]]).reverse_complement()}'
                rnd_ornt[f'{seq_subsamp[n]}_RF{abs(rf) + 3}'] = trim_seq(temp_seq, abs(rf))

    return rnd_ornt


def freq_counts(orf_dict: dict, kmer: int = 3) -> dict:
    orf_kmer_dict = {}
    kmer_set = [''.join(i) for i in product(['A', 'T', 'G', 'C'], repeat = kmer)]

    for k, v in orf_dict.items():
        kmer_list = [v[n:n+kmer] for n in range(0, len(v), kmer)]
        orf_kmer_dict[k] = {kmer: (kmer_list.count(kmer) / len(kmer_list)) for kmer in kmer_set}

    return orf_kmer_dict


def generate_contam_calls(
                        targ_seqs: dict,
                        non_targ_seqs: dict,
                        fasta_file: str,
                        txn_code: str,
                        pretrained: str,
                        start_time,
                        gcode: str = '1',
                        kmer: int = 3,
                        verb: bool = True):

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Preparing training and query ORFs for {txn_code}')

    tmp_query = {}
    for i in SeqIO.parse(fasta_file,'fasta'):
        tmp_query[i.id] = f'{i.seq}'

        if len(i)%3 > 0:
            tmp_query[i.id] = f'{i.seq}'[:-(len(i)%3)]

    query_orfs = freq_counts(tmp_query)
    ref_orfs = None

    if not pretrained:
        ref_orfs = freq_counts(targ_seqs)
        ref_orfs.update(freq_counts(non_targ_seqs))


    return query_orfs, ref_orfs


def generate_orf_calls(
                    fasta_file: str,
                    txn_code: str,
                    start_time,
                    dmnd_db: str,
                    gcode: str = '1',
                    pretrained = None,
                    min_len: int = 300,
                    evalue: float = 1e-30,
                    threads: int = 1,
                    strand: str = 'both',
                    kmer: int = 3,
                    verb: bool = True):

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Extracting complete putative ORFs for {txn_code}')

    Path(f'{txn_code}_TIdeS/ORF_Calling/').mkdir(exist_ok = True, parents = True)

    stp_cdns = eval_gcode_ttable(str(gcode))[:-1]
    init_query_orfs, comp_orf_fas = call_all_orfs(fasta_file, txn_code, gcode, min_len, strand)

    if pretrained:
        rnd_orfs = None

        if verb:
            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Preparing query ORFs for {txn_code}')

    else:


        tmp_ref_orfs, tmp_rnd_rfs = generate_ref_orfs(fasta_file,
                                                        txn_code,
                                                        start_time,
                                                        dmnd_db,
                                                        min_len,
                                                        evalue,
                                                        threads,
                                                        verb)
        if verb:
            print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Preparing training and query ORFs for {txn_code}')

        rnd_orfs = freq_counts(tmp_rnd_rfs)

    query_orfs = freq_counts(init_query_orfs)

    return query_orfs, rnd_orfs, comp_orf_fas
