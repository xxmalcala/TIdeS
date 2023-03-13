#!/usr/bin/env python3

import argparse, glob, os, shutil, sys, time

from argparse import RawTextHelpFormatter, SUPPRESS
from datetime import timedelta
from Bio import SeqIO

from bin import filter_txps as ft
from bin import orf_call as oc
from bin import classify_orfs as cl
from bin import save_preds as sp


def collect_args():
    parser = argparse.ArgumentParser(description = f'{ascii_msg()}\n     TIdeS ' \
            'leverages supervised machine-learning to hunt\n     ' \
            'and select "optimal" putative ORFs (pORFs) from a given\n     ' \
            'transcriptome and can infer ORFs from a target taxon\n     ' \
            'from "noisy" (e.g. scRNA-seq) datasets.',
            formatter_class=RawTextHelpFormatter)

    # Set of required options for the script
    required_arg_group = parser.add_argument_group('Input-Output Options')

    required_arg_group.add_argument('--fin', '-i', action = 'store',
        metavar = '[FASTA File]', type = str, required = True,
        help = ' FASTA file of transcripts or ORFs.\n\n')

    required_arg_group.add_argument('--taxon-name','-n', nargs='+',
        action = 'store', metavar = '[Taxon]', type = str, required = True,
        help=" Taxon name or PhyloToL taxon-code.\n")

    porf_arg_group = parser.add_argument_group('TIdeS ORF-Calling Options')

    porf_arg_group.add_argument('--db','-d', action = 'store',
        metavar = '[Protein Database]', type = str,
        help = ' Path to protein database.\n\n')\

    porf_arg_group.add_argument('--genetic-code','-g', action = 'store',
        default = 1, metavar = '[Genetic-Code]', type = int,
        help = ' NCBI supported translation table (default = 1).\n\n')

    porf_arg_group.add_argument('--evalue','-e', action = 'store',
        default = 1e-30, metavar = '[e-value]', type = float,
        help = ' E-value for reference ORF calling (default = 1e-30)\n\n')

    # porf_arg_group.add_argument('--stranded','-st', action = 'store_true',
    #     help = ' Transcripts are stranded (consider RFs +1 to +3).\n\n')

    # porf_arg_group.add_argument('--orfs','-orfs', action = 'store_true',
    #     help = ' FASTA file input contains user-defined "true" ORFs to be used for training.\n\n')


    contam_arg_group = parser.add_argument_group('TIdeS Contamination-Calling Options')

    contam_arg_group.add_argument('--contam','-c', #action = 'store_true',
        nargs='?', const = 'All',
        help = 'Default accounts for all non-same-minor clade sequences,\n' \
        'otherwise provide a comma-separated list of contaminant\ntaxon-codes. ' \
        '(e.g. EE_cr_Gthe,Sr_st).\n\n')

    contam_arg_group.add_argument('--sister-table','-s', action = 'store',
        metavar = '[Sister-Relationship Table]', type = str,
        help = ' Run contamination pipeline... (TO UPDATE)\n')

    optional_arg_group = parser.add_argument_group('General Options')
    optional_arg_group.add_argument('--min-len','-l', action = 'store',
        default = 300, metavar = '[bp]', type = int,
        help = ' Minimum ORF size (bp) to consider (default = 300)\n\n')

    optional_arg_group.add_argument('--threads','-p', action = 'store',
        default = 4, metavar = '[Threads]', type = int,
        help = ' Number of threads (default = 4)\n')

    optional_arg_group.add_argument('--rfc','-r', action = 'store',
        metavar = '[Trained-RFC]', type = str, default = None,
        help = ' Previous TIdeS trained RFC model.\n')

    # optional_arg_group.add_argument('--train-rfc','-train-rfc', action = 'store_true',
    #     help = ' Only runs the RFC training.\n\n')

    optional_arg_group.add_argument('--gzip','-gz', action = 'store_true',
        help = ' Compress output folder and data.\n\n')


    # Ensures that just script description provided if no arguments provided
    if len(sys.argv[1:]) == 0:
        print (parser.description+'\n')
        sys.exit()

    args = parser.parse_args()
    # print(args)
    # sys.exit()
    args.taxon_name = '_'.join(args.taxon_name)

    return args



def ascii_msg():
    msg = """      _____ ___    _     ___
     |_   _|_ _|__| |___/ __|
       | |  | |/ _` / -_)__ \\
       |_| |___\__,_\___|___/

     Version 1.2.0
    """
    return msg

### Change default threads to 2!!!
def predict_orfs(fasta_file: str,
                taxon_code:str,
                dmnd_db: str,
                gcode: str = '1',
                pretrained = None,
                min_len:int = 300,
                pid: float = 0.97,
                threads: int = 24,
                intermed: bool = True,
                verb: bool = True) -> None:

    sttime = time.time()

    if verb:
        print(ascii_msg())

    if verb:
        print('#------ Preparing Transcriptome Data -------#')

    filt_fas = ft.filter_transcripts(fasta_file,
                                        taxon_code,
                                        sttime,
                                        min_len,
                                        threads,
                                        pid,
                                        verb)

    if verb:
        print('\n#------- Calling Training and pORFs --------#')

    query_orfs, rnd_orfs = oc.generate_orf_calls(filt_fas,
                                                    taxon_code,
                                                    sttime,
                                                    dmnd_db,
                                                    gcode,
                                                    min_len,
                                                    threads,
                                                    intermed,
                                                    3,
                                                    verb)

    if verb:
        print('\n#----------- ORF Classification ------------#')

    pos_tides_preds, final_preds = cl.classify_orfs(taxon_code,
                                                    rnd_orfs,
                                                    query_orfs,
                                                    sttime,
                                                    pretrained,
                                                    threads,
                                                    verb,
                                                    False)

    if verb:
        print('\n#---------- Saving TIdeS Outputs -----------#')

    comp_orf_fas = f'{taxon_code}_TIdeS/ORF_Calling/{taxon_code}.{min_len}bp.CompORFs.fas'


    sp.finalize_outputs(taxon_code,
                        comp_orf_fas,
                        final_preds,
                        pos_tides_preds,
                        gcode,
                        min_len)

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Finished running TIdeS')
        print('\n#--------- Thanks for using TIdeS! ---------#\n')



def eval_contam(fasta_file: str,
                taxon_code: str,
                sister_summary: str,
                pretrained: str = None,
                threads: int = 24,
                intermed: bool = True,
                verb: bool = True) -> None:

    sttime = time.time()

    if verb:
        print(ascii_msg())

    sttime = time.time()

    if verb:
        print('#---- Preparing User-Assessed ORF Data -----#')

    tg_seqs, ntg_seqs = ft.prep_contam(fasta_file,
                                taxon_code,
                                sister_summary,
                                sttime,
                                verb)

    query_orfs, ref_orfs = oc.generate_contam_calls(tg_seqs,
                                                    ntg_seqs,
                                                    fasta_file,
                                                    taxon_code,
                                                    sttime,
                                                    gcode,
                                                    3,
                                                    verb)

    if verb:
        print('\n#----------- ORF Classification ------------#')

    tg_preds, ntg_preds = cl.classify_orfs(taxon_code,
                                    ref_orfs,
                                    query_orfs,
                                    sttime,
                                    pretrained,
                                    threads,
                                    verb,
                                    True)

    if verb:
        print('\n#---------- Saving TIdeS Outputs -----------#')

    tds_tg_pf = f'{taxon_code}_TIdeS/{taxon_code}.TIdeS.fas'
    tds_ntg_pf = f'{taxon_code}_TIdeS/{taxon_code}.NonTarget.TIdeS.fas'
    tds_tg_ps, tds_ntg_ps = {}, {}

    with open(tds_tg_pf, 'w+') as w:
        for i in SeqIO.parse(fasta_file,'fasta'):
            if i.id in tg_preds:
                w.write(f'>{i.description}\n{i.seq}\n')
                tds_tg_ps[i.description] = f'{i.seq}'

    with open(tds_ntg_pf, 'w+') as w:
        for i in SeqIO.parse(fasta_file,'fasta'):
            if i.id in ntg_preds:
                w.write(f'>{i.description}\n{i.seq}\n')
                tds_ntg_ps[i.description] = f'{i.seq}'

    with open(tds_tg_pf.replace("fas","AA.fas"), 'w+') as w:
        tds_psa = oc.translate_seqs(tds_tg_ps, gcode)
        for k, v in tds_psa.items():
                w.write(f'{k}\n{v}\n')

    with open(tds_ntg_pf.replace("fas","AA.fas"), 'w+') as w:
        tds_psa = oc.translate_seqs(tds_ntg_ps, gcode)
        for k, v in tds_psa.items():
                w.write(f'{k}\n{v}\n')

    for fas in glob.glob(f'{taxon_code}_TIdeS/{taxon_code}.*fas'):
        shutil.copy2(fas, f'{taxon_code}_TIdeS/Classified/')

    if verb:
        print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Finished running TIdeS')
        print('\n#--------- Thanks for using TIdeS! ---------#\n')


if __name__ == '__main__':

    """
    Notes to self...

    FASTA file of all predictions (ORF-calling) should be viewed as intermediate
    files. -- FLAG this

    Add 'single-best' option for ORFs to break ties when ORF calls have identical
    proba values (i.e. return longest when prob for pORF1 == pORF2)... All preds is intermed file.

    Add stranded-seq option. Need to find dataset for this to be certain...

    Add argparse related stuff...

    Also put save-preds functionality into eval_contam pipeline.

    Reset threads to 1 or 2 as defaults...
    """


    # args = check_args()
    # print('Under Construction!')
    # sys.exit()

    args = collect_args()

    # if 1 < len(sys.argv[1:]) < 5:
    #     fasta_file = sys.argv[1]
    #     taxon_code = sys.argv[2]
    #     try:
    #         gcode = sys.argv[3]
    #         dmnd_db = sys.argv[4]
    #     except:
    #         gcode = '1'
    #         dmnd_db = 'tides_aa_db.dmnd'
    #
    # else:
    #     print("\nTemporary Usage:\n    python3 tides.py [TRANSCRIPTOME] [TAXON-NAME] [TRANSLATION-TABLE]\n")
    #     sys.exit()

    contam = False

    if not oc.eval_gcode_ttable(gcode):
        sys.exit(1)

    if not contam:
        predict_orfs(fasta_file, taxon_code, dmnd_db)
    else:
        eval_contam(fasta_file, taxon_code, sister_summary)
