#!/usr/bin/env python3

"""
Extracts putative ORFs from a given transcriptome's FASTA-formatted file,
using a user-supplied genetic code and which meet a minimum length criterion.

Returns the putative ORFs and their sequences from all transcripts.


Dependencies include: BioPython, NumPy.
"""

import os, re, subprocess, sys, time

from datetime import timedelta
from pathlib import Path

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def gcode_start_stops(gcode: str) -> list:
    """
    Provides a list of stop codons from a given genetic code (i.e. translation table).

    Parameters
    ----------
    gcode:  NCBI numerical genetic code

    Returns list of valid stop codons on plus and minus strands as well as translation table
    """

    ncbi_link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi'

    stop_codon_by_gcodes = {
        '1': [['TAA','TAG','TGA'],['TTA','CTA','TCA'], 1],
        '4': [['TAA','TAG'],['TTA','CTA'], 4],
        '6': [['TGA'],['TCA'], 6],
        '10': [['TAA','TAG'],['TTA','CTA'], 10],
        '12': [['TAA','TAG','TGA'],['TTA','CTA','TCA'], 12],
        '26': [['TAA','TAG','TGA'],['TTA','CTA','TCA'], 26],
        '29': [['TGA'],['TCA'], 29],
        '30': [['TGA'],['TCA'], 30],
        'universal': [['TAA','TAG','TGA'],['TTA','CTA','TCA'], 1],
        'blepharisma': [['TAA','TAG'],['TTA','CTA'], 4],
        'ciliate': [['TGA'],['TCA'], 6],
        'euplotes': [['TAA','TAG'],['TTA','CTA'], 10],
        'alt-yeast': [['TAA','TAG','TGA'],['TTA','CTA','TCA'], 12],
        'pachysolen': [['TAA','TAG','TGA'],['TTA','CTA','TCA'], 26],
        'mesodinium': [['TGA'],['TCA'], 29],
        'peritrich': [['TGA'],['TCA'], 30]
        }

    if str(gcode).lower() not in stop_codon_by_gcodes.keys():
        print('\nCheck the chosen genetic code/translation table.\n\n' \
            'Supported tables include:\n    ' +
            ',\n    '.join(list(stop_codon_by_gcodes.keys())))

        print(f'\nFor more information, visit: {ncbi_link}\n')

        return None

    return stop_codon_by_gcodes[str(gcode)]


def bin_pos(cdn_pos: list, rf: dict, strand: str = '+') -> dict:
    """
    Organizes all codon positions of interest by reading frame and strand.

    Parameters
    ----------
    cdn_pos:  list of codon positions
    rf:       reading-frame dictionary to populate
    strand:   plus or minus strand

    Returns
    ----------
    rf: dictionary of all codon positions by reading frame
    """

    rf_adj = 0
    if strand != '+':
        rf_adj += 3

    for i in cdn_pos:
        if i % 3 == 0:
            rf[1 + rf_adj].append(i)
        elif i % 3 == 1:
            rf[2 + rf_adj].append(i)
        else:
            rf[3 + rf_adj].append(i)

    return rf


def capture_start_stop_pos(seq: str, fstop: list, rstop: list) -> dict:
    """
    Identifies all start and stop codon position for a given sequence in all
    reading frames.

    Parameters
    ----------
    seq:    single transcript/sequence
    fstop:  list of stop codons in the plus strand orientation
    rstop:  list of stop codons in the minus strand orientation

    Returns
    ----------
    start_rf:  dictionary of start codon positions by reading frame
    stop_rf:   dictionary of stop codon positions by reading frame
    """

    stop_cdns_p, stop_cdns_m = [], []

    # Capture indices of start codons from both strands
    start_cdns_p = [m.start() for m in re.finditer('ATG', seq)]
    start_cdns_m = [m.end() for m in re.finditer('CAT', seq)]

    # Capture indices of all stop codons from both strands
    for n in range(len(fstop)):
        stop_cdns_p += [m.end() for m in re.finditer(fstop[n], seq)]
        stop_cdns_m += [m.start() for m in re.finditer(rstop[n], seq)]

    start_rf = bin_pos(start_cdns_p, {1:[], 2:[], 3:[], 4:[], 5:[], 6:[]})
    start_rf = bin_pos(start_cdns_m, start_rf, '-')

    stop_rf = bin_pos(stop_cdns_p, {1:[], 2:[], 3:[], 4:[], 5:[], 6:[]})
    stop_rf = bin_pos(stop_cdns_m, stop_rf, '-')

    return start_rf, stop_rf


def check_start_dist(tmp_start: list, stop_coord: int, rf: int) -> list:
    """
    Assesses putative start codon positions, keeping a set of start codons at
    least 15nt apart from one another. Note that this starts with the first
    codon (by position in array) in the "tmp_start" array of start codons positions.

    Parameters
    ----------
    tmp_start:   numpy array of start codon positions
    stop_coord:  stop codon coordinate for this putative ORF
    rf:          current reading frame

    Returns list of pORF coordinates for a single transcript fragment between stop codons
    """

    stp_pos = stop_coord[1]
    if rf > 3:
        stp_pos = stop_coord[0]

    # Grab the initial start codon coordinate and stop codon for this pORF
    porf_coords = [(tmp_start[0], stp_pos)]

    # Start tracking the last valid start codon
    last_valid = tmp_start[0]
    tmp_start = tmp_start[1:]

    # Track whether to skip the next start codon coordinate in the array
    skip_pos = False
    for n in range(len(tmp_start)):
        if skip_pos:
            skip_pos = False
            continue

        # Check that current start coordinate is at least 15nt from previous valid start codon
        elif abs(tmp_start[n] - last_valid) > 14:
            porf_coords.append((tmp_start[n], stp_pos))
            last_valid = tmp_start[n]

            # If final coordinate, no further comparison needed
            if n != len(tmp_start) - 1:

                # Evaluate skipping next codon
                if abs(tmp_start[n] - tmp_start[n+1]) < 14:
                    skip_pos = True
                else:
                    skip_pos = False
        else:
            pass

    return list(set(porf_coords))


def check_partial_pORFs(start_coords: np.ndarray, stop_coords: np.ndarray, rf: int, min_len: int, seq_len: int) -> list:
    """
    Checks for partial pORFs if no start codons could create a complete pORF
    (5' partial ORF), or no stop codon could create a complete pORF (3' partial ORF).

    Parameters
    ----------
    start_coords:  numpy array of start codon positions
    stop_coords:   numpy array of stop codon positions
    rf:            current reading frame
    min_len:       minimum ORF length to consider [default is 300nt]

    Returns
    ----------
     partial_coords:  list of partial pORF coordinates in the given reading frame
    """

    partial_coords = []
    if rf > 3:
        seq_start = rf - 4

    else:
        seq_start = rf - 1

    seq_end = seq_len - (seq_len - seq_start)%3

    if stop_coords.size == start_coords.size == 0:
        partial_coords.append((seq_start, seq_end))

    # Captures partial pORFs lacking a start codon in the current reading frame
    elif stop_coords.size != start_coords.size == 0:
        if rf > 3:
            if (seq_end - stop_coords[-1]) >= min_len:
                partial_coords.append((stop_coords[-1], seq_end))
        else:
            if (stop_coords[0] - seq_start) >= min_len:
                partial_coords.append((seq_start, stop_coords[0]))

    # Captures partial 3' pORFs lacking a stop codon in the current reading frame
    elif start_coords.size != stop_coords.size == 0:
        if rf > 3:
            tmp_start = sorted(start_coords[start_coords >= min_len])
        else:
            tmp_start = sorted(start_coords[abs(seq_len - start_coords) >= min_len])
        if tmp_start:
            partial_coords += check_start_dist(tmp_start, [seq_start, seq_end], rf)

    else:
        # Keep first and last stop codons to check for partial ORF ends
        eval_stops = sorted([min(stop_coords), max(stop_coords)])

        if rf > 3:
            # Capture partial 5' ORFs
            if (seq_end - eval_stops[1]) >= min_len:
                if eval_stops[1] > max(start_coords):
                    partial_coords.append((eval_stops[1], seq_end))
                elif abs(eval_stops[1] - max(start_coords)) < min_len:
                    partial_coords.append((eval_stops[1], seq_end))

            # Capture partial 3' ORFs
            eval_start = [i for i in start_coords if i <= eval_stops[0] and abs(i - seq_start) >= min_len]
            if eval_start:
                partial_coords += check_start_dist(eval_start, [seq_start, seq_end], rf)

        else:
            # Capture partial 5' ORFs
            if (eval_stops[0] - seq_start) >= min_len:
                if eval_stops[0] < min(start_coords):
                    partial_coords.append((seq_start, eval_stops[0]))
                elif abs(eval_stops[0] - min(start_coords)) < min_len:
                    partial_coords.append((seq_start, eval_stops[0]))

            # Capture partial 3' ORFs
            eval_start = [i for i in start_coords if i >= eval_stops[-1] and abs(i - seq_end) >= min_len]
            if eval_start:
                partial_coords += check_start_dist(eval_start, [seq_start, seq_end], rf)

    return partial_coords


def eval_pORF_coords(start_coords: np.ndarray, stop_coords: np.ndarray, rf: int, min_len: int, seq_len: int):
    """
    Evaluates combinations of start and stop coordinates to return a list of
    putative ORF (pORF) coordinates from a given transcript.

    Parameters
    ----------
    start_coords:  numpy array of start codon positions
    stop_coords:   numpy array of stop codon positions
    rf:            current reading frame
    min_len:       minimum ORF length to consider [default is 300nt]

    Returns
    ----------
    porf_coords:  list of pORF coordinates in the given reading frame
    """

    porf_coords = []
    valid_stops = []

    if rf > 3:
        seq_start = rf - 4
        seq_end = seq_len - (seq_len - seq_start)%3
        if stop_coords.size > 0 and abs(stop_coords[-1] - seq_end) >= min_len:
            valid_stops.append((stop_coords[-1], seq_end))

    else:
        seq_start = rf - 1
        if stop_coords.size  > 0 and abs(stop_coords[0] - seq_start) >= min_len:
            valid_stops.append((seq_start, stop_coords[0]))

    # keep windows between stop codons that are at least "min_len" nucleotides apart
    valid_stops += [(j, stop_coords[i+1]) for i, j in enumerate(stop_coords[:-1]) if abs(j - stop_coords[i+1]) >= min_len]

    if valid_stops:
        for stop_coord in valid_stops:
            if rf > 3:
                # Captures all start codon positions between stop codons that would
                # create ORFs larger than the minimum ORF size
                tmp_start = sorted(
                    [i for i in start_coords[start_coords <= stop_coord[1]] if (i - stop_coord[0]) >= min_len],
                    reverse = True
                    )
            else:
                # captures all start codon positions between stop codons that would
                # create ORFs larger than the minimum ORF size
                tmp_start = sorted(
                    [i for i in start_coords[start_coords >= stop_coord[0]] if (stop_coord[1] - i) >= min_len]
                    )
            if tmp_start:
                porf_coords += check_start_dist(tmp_start, stop_coord, rf)

    porf_coords += check_partial_pORFs(start_coords, stop_coords, rf, min_len, seq_len)

    return porf_coords


def prepare_pORF_seqs(seq_num: int, seq_name: str, seq: str, coord: list, fstop: list, rf: int):
    """
    Prepares the pORFs from a set of transcript coordinates.

    Parameters
    ----------
    seq_num:   the current pORF number to append to name
    seq_name:  original transcript name
    seq:       the transcript nucleotide sequence
    coord:     ORF start-stop coordinates
    fstop:     all stop codons on the plus strand/orienatation
    rf:        current reading frame

    Returns
    ----------
    porf:            ORF nucleotide sequence
    porf_name:       name of the putative ORF
    porf_full_name:  the detailed pORF name (include completeness, length, coordinates)
    """

    porf = seq[min(coord):max(coord)]
    porf_len = f'orf_length:{len(porf)}'
    porf_name = f'{seq_name}.pORF{seq_num}'

    if rf > 3:
        porf = f'{Seq(porf).reverse_complement()}'

        if porf[:3] == 'ATG' and porf[-3:] in fstop:
            comp = 'orf_type:complete'

        elif porf[:3] != 'ATG' and porf[-3:] in fstop:
            comp = 'orf_type:5prime_partial'

        else:
            comp = 'orf_type:3prime_partial'

        porf_full_name = f'{porf_name} {comp} {porf_len} {seq_name}:{max(coord)}-{min(coord)+1}(-)'

    else:
        if porf[:3] == 'ATG' and porf[-3:] in fstop:
            comp = 'orf_type:complete'

        elif porf[:3] != 'ATG' and porf[-3:] in fstop:
            comp = 'orf_type:5prime_partial'

        elif porf[:3] == 'ATG' and porf[-3:] not in fstop:
            comp = 'orf_type:3prime_partial'

        else:
            comp = 'orf_type:incomplete'

        porf_full_name = f'{porf_name} {comp} {porf_len} {seq_name}:{coord[0]+1}-{coord[1]}(+)'

    return porf_name, porf_full_name, porf


def extract_orfs_from_seq(seq_name: str, seq: str, fstop, rstop, min_len: int, stranded: str = None):
    """
    Extracts all the putative ORFs greater than some minimum length for all reading
    frames (unless stranded RNA-seq/assembly) for a given transcript.

    Parameters
    ----------
    query:     name of the "filtered" transcriptome's FASTA formatted file
    seq:       the transcript nucleotide sequence
    fstop:     all stop codons on the plus strand
    rstop:     all stop codons on the minus strand (rev. complement of fstop)
    min_len:   minimum ORF length to consider [default is 300nt]
    stranded:  designates from which strand (if not both) to call putative ORFs

    Returns dictionary of all putative ORFs (partials included)
    """

    all_pORFs = {}

    start_rf, stop_rf = capture_start_stop_pos(seq, fstop, rstop)

    # Remove reading frames that shouldn't be considered due to
    # stranded-sequencing and assembly protocols
    if stranded and stranded.lower() == 'plus':
        del stop_rf[4]
        del stop_rf[5]
        del stop_rf[6]

    elif stranded and stranded.lower() == 'minus':
        del stop_rf[1]
        del stop_rf[2]
        del stop_rf[3]

    else:
        pass

    porf_num = 1

    for k in stop_rf.keys():
        stop_coords = np.array(sorted(stop_rf[k]))
        start_coords = np.array(sorted(start_rf[k]))

        porf_coords = eval_pORF_coords(
                        start_coords,
                        stop_coords,
                        k,
                        min_len,
                        len(seq))
        # If a given reading frame has pORFs, then prepare them...
        if porf_coords:
            for coord in porf_coords:
                prepped_porf = prepare_pORF_seqs(porf_num, seq_name, seq, coord, fstop, k)

                all_pORFs[prepped_porf[0]] = prepped_porf[1:]

                porf_num += 1

    return all_pORFs


def ref_orf_dmnd(query: str, taxon_code: str, dmnd_db: str, min_len: int, evalue: float, threads: int, intermed: bool) -> list:
    """
    Runs DIAMOND BLASTX against a reference protein database.

    Parameters
    ----------
    query:       name of the "filtered" transcriptome's FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    dmnd_db:     path to a protein database for DIAMOND
    min_len:     minimum ORF length to consider
    evalue:      maximum e-value to consider hits from DIAMOND
    threads:     number of threads to use
    intermed:    keep intermediate TSV outputs from DIAMOND

    Returns list of "BLASTX" hits for parsing
    """

    outfmt = '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe'

    if intermed:
        tsv_out = f'{taxon_code}_TIdeS/ORF_Calling/{taxon_code}.{min_len}bp.DMND_REF_ORFs.tsv'
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


def generate_ref_orfs(fasta_file: str, taxon_code: str, start_time, dmnd_db: str, min_len: int, evalue: float, threads: int, verb: bool) -> dict:
    """
    Generate reference ORFs from outputs of DIAMOND BLASTX against a reference
    protein database.

    Parameters
    ----------
    query:       name of the "filtered" transcriptome's FASTA formatted file
    fasta_file:  name of the "filtered" transcriptome's FASTA formatted file
    gcode:       NCBI numerical genetic code
    rstop:       all stop codons on the minus strand (rev. complement of fstop)
    min_len:     minimum ORF length to consider [default is 300nt]
    stranded:    designates from which strand (if not both) to call putative ORFs

    Returns
    ----------
    ref_orf_db: dictionary of reference ORFs
    """

    ref_orf_db = {}
    ref_orfs = []

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Running DIAMOND BLASTX against protein database')

    dmnd_hits = {i.split('\t')[0]:i.split('\t') for i in
                    ref_orf_dmnd(fasta_file, taxon_code, dmnd_db, min_len, evalue, threads, True)}

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Generating training ORFs for {taxon_code}')

    for seq_rec in SeqIO.parse(fasta_file,'fasta'):
        if seq_rec.id in dmnd_hits:
            orf_pos = [int(j) for j in dmnd_hits[seq_rec.id][4:6]]
            orf = f'{seq_rec.seq[min(orf_pos)-1:max(orf_pos)]}'

            if orf_pos[0] > orf_pos[1]:
                ref_orf_db[seq_rec.id] = f'{Seq(orf).reverse_complement()}'
                ref_orfs.append(SeqRecord(Seq(orf).reverse_complement(),seq_rec.id,'',''))

            else:
                ref_orf_db[seq_rec.id] = orf
                ref_orfs.append(SeqRecord(Seq(orf),seq_rec.id,'',''))

    return ref_orf_db


def capture_pORFs(fasta_file: str, taxon_code: str, start_time, gcode: int = 1, min_len: int = 300, partial: bool = False, stranded: str = None, verb: bool = True):
    """
    Manages the capture of all putative ORFs from all transcripts in a given
    FASTA file.

    Parameters
    ----------
    fasta_file:  name of the "filtered" transcriptome's FASTA formatted file
    taxon_code:  species/taxon name or abbreviated code
    start_time:  initial timestamp to track runtime
    gcode:       NCBI numerical genetic code
    min_len:     minimum ORF length to consider [default is 300nt]
    stranded:    designates from which strand (if not both) to call putative ORFs
    verb:        verbose print statements

    Returns
    ----------
    pORF_db:     dictionary of all putative ORFs (partials included)
    porf_fasta:  FASTA formatted file of all putative ORFs
    """

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-start_time))}]  Extracting putative ORFs for {taxon_code}')

    porf_dir = f'{taxon_code}_TIdeS/ORF_Calling/'
    porf_fasta = f'{porf_dir}{taxon_code}.{min_len}bp.pORFs.fas'

    Path(porf_dir).mkdir(exist_ok = True, parents = True)

    gcode_info = gcode_start_stops(gcode)

    pORF_db = {}
    all_pORFs = []

    for i in SeqIO.parse(fasta_file, 'fasta'):
        seq_porfs = extract_orfs_from_seq(
                        i.id,
                        f'{i.seq}',
                        gcode_info[0],
                        gcode_info[1],
                        min_len,
                        stranded)
        if seq_porfs:
            for i in seq_porfs.values():
                if not partial and 'orf_type:complete' not in i[0]:
                    continue
                all_pORFs.append(SeqRecord(Seq(i[1]), i[0].split()[0], '', i[0]))

                pORF_db[i[0]] = i[1]

    SeqIO.write(all_pORFs, porf_fasta, 'fasta')

    return pORF_db, porf_fasta
