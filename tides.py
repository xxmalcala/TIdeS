#!/usr/bin/env python3

import argparse, glob, os, shutil, sys, time

from datetime import timedelta
from Bio import SeqIO

from bin import filter_txps as ft
from bin import orf_call as oc
from bin import classify_orfs as cl
from bin import save_preds as sp


def collect_args():

    parser = argparse.ArgumentParser(description = f'{ascii_logo_vsn()}\n{usage_msg()}',
            usage=argparse.SUPPRESS, add_help = False,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    g = parser.add_argument_group('General Options', description = (
    '''--fin (-f)            input file in FASTA format\n'''
    '''--taxon (-n)          taxon-name or PhyloToL taxon-code\n'''
    '''--threads (-p)        number of CPU threads (default = 1)\n'''
    '''--model (-m)          previously trained TIdeS model (".pkl" file)\n'''
    '''--quiet (-q)          no console output\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n'''
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
    '''--strand (-s)         query strands to call ORFs (both/minus/plus, default = both)'''))

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

    porf.add_argument('--strand', '-s', choices = ['both', 'plus', 'minus'], action = 'store',
        default = 'both', type = str, help = argparse.SUPPRESS)

    cntm = parser.add_argument_group('Contamination-Calling Options', description = (
    '''--contam (-c)         evalute ORFs for contamination (may be discontinued)\n'''
    '''--sister-table        table of sequences annotated as target or contamination\n\n'''))

    cntm.add_argument('--contam','-c', action = 'store_true',
        help = argparse.SUPPRESS)

    cntm.add_argument('--sister-table', action = 'store',
        metavar = '[Sister-Relationship Table]', type = str,
        help = argparse.SUPPRESS)


    # Ensures that just script description provided if no arguments provided
    if len(sys.argv[1:]) == 0:
        print(ascii_logo_vsn())
        print(detail_msg())
        print(f'{usage_msg()}\n')
        sys.exit()

    args = parser.parse_args()

    args.taxon = '_'.join(args.taxon)

    return args


def ascii_logo_vsn():
    alv_msg = """      _____ ___    _     ___
     |_   _|_ _|__| |___/ __|
       | |  | |/ _` / -_)__ \\
       |_| |___\__,_\___|___/

     Version 1.2.0
    """
    return alv_msg


def detail_msg():
    dtl_msg = """     TIdeS leverages supervised machine-learning to hunt
     and select "optimal" putative ORFs (pORFs) from a
     transcriptome and can infer ORFs from a target taxon
     from "noisy" (e.g. scRNA-seq) datasets.\n"""
    return dtl_msg


def usage_msg():
    return """Usage:\n    tides.py [options] --fin [FASTA file] --taxon [taxon name]"""


def predict_orfs(
                fasta_file: str,
                taxon_code:str,
                dmnd_db: str,
                gcode: str = '1',
                pretrained = None,
                min_len:int = 300,
                pid: float = 0.97,
                evalue: float = 1e-30,
                threads: int = 1,
                strand:str = 'both',
                verb: bool = True) -> None:

    sttime = time.time()

    if verb:
        print('#------ Preparing Transcriptome Data -------#')

    filt_fas = ft.filter_transcripts(
                                    fasta_file,
                                    taxon_code,
                                    sttime,
                                    min_len,
                                    threads,
                                    pid,
                                    verb)

    if verb:
        if not pretrained:
            print('\n#------- Calling Training and pORFs --------#')
        else:
            print('\n#------------- Calling pORFs ---------------#')

    query_orfs, rnd_orfs, comp_orf_fas = oc.generate_orf_calls(
                                                            filt_fas,
                                                            taxon_code,
                                                            sttime,
                                                            dmnd_db,
                                                            gcode,
                                                            pretrained,
                                                            min_len,
                                                            evalue,
                                                            threads,
                                                            strand,
                                                            3,
                                                            verb)

    if verb:
        print('\n#----------- ORF Classification ------------#')

    pos_tides_preds, final_preds = cl.classify_orfs(
                                                taxon_code,
                                                rnd_orfs,
                                                query_orfs,
                                                sttime,
                                                pretrained,
                                                threads,
                                                verb,
                                                False)

    if verb:
        print('\n#---------- Saving TIdeS Outputs -----------#')

    sp.finalize_outputs(
                    taxon_code,
                    comp_orf_fas,
                    final_preds,
                    pos_tides_preds,
                    gcode,
                    min_len)

    return sttime


def eval_contam(
            fasta_file: str,
            taxon_code: str,
            sister_summary: str,
            gcode: str = '1',
            pretrained: str = None,
            min_len: int = 300,
            threads: int = 1,
            verb: bool = True) -> None:

    sttime = time.time()

    if verb:
        print('#---- Preparing User-Assessed ORF Data -----#')

    tg_seqs, ntg_seqs = ft.prep_contam(
                                fasta_file,
                                taxon_code,
                                sister_summary,
                                pretrained,
                                sttime,
                                verb)

    query_orfs, ref_orfs = oc.generate_contam_calls(
                                                tg_seqs,
                                                ntg_seqs,
                                                fasta_file,
                                                taxon_code,
                                                pretrained,
                                                sttime,
                                                gcode,
                                                3,
                                                verb)

    if verb:
        print('\n#----------- ORF Classification ------------#')

    tg_preds, ntg_preds = cl.classify_orfs(
                                        taxon_code,
                                        ref_orfs,
                                        query_orfs,
                                        sttime,
                                        pretrained,
                                        threads,
                                        verb,
                                        True)

    if verb:
        print('\n#---------- Saving TIdeS Outputs -----------#')


    sp.finalize_outputs(
                    taxon_code,
                    fasta_file,
                    ntg_preds,
                    tg_preds,
                    gcode,
                    min_len,
                    True)

    return sttime


if __name__ == '__main__':

    """Note to self:
    Add method to check the sister-summary table format...

    Requirements:
    -- Tab-delimited
    -- Target/Non-Target OR 1/0
    """

    args = collect_args()

    if not oc.eval_gcode_ttable(args.gencode):
        sys.exit(1)

    args.pid = args.pid / 100

    if not args.contam:
        if not args.db and not args.model:
            print('\nERROR: Missing diamond "blast" database -- use the "prep_tides_db.sh" '
                    'script in the Tides/util/ directory.\n')
            sys.exit(1)

        if not args.quiet:
            print(ascii_logo_vsn())

        sttime = predict_orfs(args.fin,
                                args.taxon,
                                args.db,
                                args.gencode,
                                args.model,
                                args.min_orf,
                                args.pid,
                                args.evalue,
                                args.threads,
                                args.strand,
                                not args.quiet)
    else:
        if not args.quiet:
            print(ascii_logo_vsn())

        sttime = eval_contam(args.fin,
                            args.taxon,
                            args.sister_table,
                            args.gencode,
                            args.model,
                            args.min_orf,
                            args.threads,
                            not args.quiet)


    if args.gzip:
        if not args.quiet:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Compressing TIdeS Outputs')

        os.system(f'tar -zcf {args.taxon}_TIdeS.tar.gz {args.taxon}_TIdeS/')

    if not args.quiet:
        print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Finished running TIdeS')
        print('\n#--------- Thanks for using TIdeS! ---------#\n')
