#!/usr/bin/env python3

import argparse, glob, os, shutil, sys, time

from datetime import timedelta
from Bio import SeqIO

from bin import filter_txps as ft
from bin import orf_call as oc
from bin import classify_orfs as cl
from bin import save_preds as sp


def collect_args():

    parser = argparse.ArgumentParser(description = f'{ascii_msg()}\n     TIdeS '
            'leverages supervised machine-learning to hunt\n     '
            'and select "optimal" putative ORFs (pORFs) from a given\n     '
            'transcriptome and can infer ORFs from a target taxon\n     '
            'from "noisy" (e.g. scRNA-seq) datasets.\n\nUsage:\n  tides.py '
            '[options] --in [FASTA file] --taxon [taxon name]',
            usage=argparse.SUPPRESS, add_help = False,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    g = parser.add_argument_group('General Options', description = (
    '''--fin (-f)             input file in FASTA format\n'''
    '''--taxon (-n)          taxon-name or PhyloToL taxon-code\n'''
    '''--threads (-p)        number of CPU threads (default = 1)\n'''
    '''--model (-m)          previously trained TIdeS model (".pkl" file)\n'''
    '''--quiet (-q)          no console output\n'''
    '''--gzip (-gz)          gzip TIdeS output\n'''
    '''--help (-h)           show this help message and exit'''))

    g.add_argument('--help', '-h', action="help", help=argparse.SUPPRESS)

    g.add_argument('--fin', '-f', action = 'store',
        metavar = ('[FASTA File]'), type = str, required = True,
        help = argparse.SUPPRESS)

    g.add_argument('--taxon','-n', nargs='+',
        action = 'store', metavar = '[Taxon]', type = str, required = True,
        help = argparse.SUPPRESS)

    g.add_argument('--threads','-p', action = 'store',
        default = 1, metavar = '[Threads]', type = int,
        help = argparse.SUPPRESS)

    g.add_argument('--model','-m', action = 'store',
        metavar = '[Trained-RFC]', type = str, default = None,
        help = argparse.SUPPRESS)

    g.add_argument('--quiet','-q', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--gzip','-gz', action = 'store_true',
        help = argparse.SUPPRESS)

    porf = parser.add_argument_group('ORF-Calling Options', description = (
    '''--db (-d)             protien database (FASTA or DIAMOND format)\n'''
    '''--id (-id)            minimum % identity to remove redundant transcripts (default = 97)\n'''
    '''--min-orf (-l)        minimum ORF length (bp) to evaluate (default = 300)\n'''
    '''--evalue (-e)         maximum e-value to infer reference ORFs (default = 1e-30)\n'''
    '''--gencode (-gc)       genetic code to use to translate ORFs\n'''
    '''--strand              query strands to call ORFs (both/minus/plus, default = both)'''))

    porf.add_argument('--db','-d', action = 'store',
        metavar = '[Protein Database]', type = str,
        help = argparse.SUPPRESS)

    porf.add_argument('--pid', '--id', '-id', action = 'store',
        default = 97, metavar = '[perc-identity]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--min-orf','-l', action = 'store',
        default = 300, metavar = '[bp]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--evalue','-e', action = 'store',
        default = 1e-30, metavar = '[e-value]', type = float,
        help = argparse.SUPPRESS)

    porf.add_argument('--gencode','-gc', action = 'store',
        default = '1', metavar = '[Genetic-Code]', type = str,
        help = argparse.SUPPRESS)

    porf.add_argument('--strand', choices = ['both', 'plus', 'minus'], action = 'store',
        default = 'both', type = str, help = argparse.SUPPRESS)

    cntm = parser.add_argument_group('Contamination-Calling Options', description = (
    '''--contam (-c)         evalute ORFs for contamination (may be discontinued)\n'''
    '''--sister-table (-s)   table of sequence annotated as target or contamination\n\n'''))

    cntm.add_argument('--contam','-c', action = 'store_true',
        help = argparse.SUPPRESS)

    cntm.add_argument('--sister-table','-s', action = 'store',
        metavar = '[Sister-Relationship Table]', type = str,
        help = argparse.SUPPRESS)


    # Ensures that just script description provided if no arguments provided
    if len(sys.argv[1:]) == 0:
        print (parser.description+'\n')
        sys.exit()

    args = parser.parse_args()

    args.taxon = '_'.join(args.taxon)

    return args


def ascii_msg():
    msg = """      _____ ___    _     ___
     |_   _|_ _|__| |___/ __|
       | |  | |/ _` / -_)__ \\
       |_| |___\__,_\___|___/

     Version 1.2.0
    """
    return msg


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

    args = collect_args()

    if not oc.eval_gcode_ttable(args.gencode):
        sys.exit(1)

    args.pid = args.pid / 100

    if not args.contam:
        if not args.db:
            print('Missing diamond "blast" database -- use the prep_tides_db.sh '
                    'script in the Tides/util/ directory.')
            sys.exit(1)

        predict_orfs(args.fin,
                    args.taxon,
                    args.db,
                    args.gencode,
                    args.model,
                    args.min_orf,
                    args.pid,
                    args.threads,
                    True,
                    not args.quiet)
    else:
        eval_contam(args.fin,
                    args.taxon,
                    args.sister_table,
                    args.model,
                    args.threads,
                    True,
                    not args.quiet)
