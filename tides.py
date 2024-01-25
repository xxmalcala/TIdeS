#!/usr/bin/env python3

"""
Transcript Identification and Selection (TIdeS) is a too for bespoke ORF calling
and sequence classification based on limited information and classic machine learning
approaches.

Overall, the intent was to aid in ORF-calling and decontamination of single-cell
transcriptomes from largely uncultivable microeukaryotes, as the popular tool
TransDecoder underperforms in these taxa.

TIdeS does learn best from quality data defined by the user, and will generate
trained models to speed up ORF-calling and ORF-classification from similar taxa
or replicates of the same taxon.
"""

import argparse, glob, os, pickle, shutil, sys, time
from datetime import timedelta
from pathlib import Path

from tides.bin import filt_seqs as ft
from tides.bin import orf_call as oc
from tides.bin import orf_prep as op
from tides.bin import classify_orfs as co
from tides.bin import save_preds as sp


def collect_args():

    parser = argparse.ArgumentParser(description = f'{ascii_logo_vsn()}\n{usage_msg()}',
            usage=argparse.SUPPRESS, add_help = False,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    g = parser.add_argument_group('General Options', description = (
    '''--fin (-i)            input file in FASTA format\n'''
    '''--taxon (-o)          taxon or output name\n'''
    '''--threads (-t)        number of CPU threads (default = 1)\n'''
    '''--kraken (-k)         kraken2 database to identify and filter\n'''
    '''                      non-eukaryotic sequences\n'''
    '''--no-filter           skip the rRNA and transcript clustering steps\n'''
    '''--model (-m)          previously trained TIdeS model (".pkl" file)\n'''
    '''--kmer                kmer size for generating sequence features (default = 3)\n'''
    '''--overlap             permit overlapping kmers (see --kmer)\n'''
    '''--step                step-size for overlapping kmers\n'''
    '''                      (default is kmer-length/2)\n'''
    '''--clean               remove intermediate filter-step files\n'''
    '''--quiet (-q)          no console output\n'''
    '''--gzip (-gz)          tar and gzip TIdeS output\n'''
    '''--help (-h)           show this help message and exit'''))

    g.add_argument('--help', '-h', action="help", help = argparse.SUPPRESS)

    g.add_argument('--fin', '-i', action = 'store',
        metavar = ('[FASTA File]'), type = str,
        help = argparse.SUPPRESS)

    g.add_argument('--taxon','-o', nargs='+',
        action = 'store', metavar = '[Taxon]', type = str,
        help = argparse.SUPPRESS)

    g.add_argument('--threads','-t', action = 'store',
        default = 1, metavar = '[Threads]', type = int,
        help = argparse.SUPPRESS)

    g.add_argument('--kraken','-k', action = 'store',
        metavar = '[Kraken2 Database]', type = str, default = None,
        help = argparse.SUPPRESS)

    g.add_argument('--no-filter', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--model','-m', action = 'store',
        metavar = '[Trained-RFC]', type = str, default = None,
        help = argparse.SUPPRESS)

    g.add_argument('--kmer', action = 'store',
        default = 3, metavar = '[kmer]', type = int,
        help = argparse.SUPPRESS)

    g.add_argument('--overlap', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--step', action = 'store',
        type = int, default = None,
        help = argparse.SUPPRESS)

    g.add_argument('--clean', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--quiet','-q', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--gzip','-gz', action = 'store_true',
        help = argparse.SUPPRESS)

    g.add_argument('--version', action = 'store_true',
        help = argparse.SUPPRESS)

    porf = parser.add_argument_group('ORF-Calling Options', description = (
    '''--db (-d)             protein database (FASTA or DIAMOND format)\n'''
    '''--partial (-p)        evaluate partial ORFs as well\n'''
    '''--id (-id)            minimum % identity to remove redundant transcripts\n'''
    '''                      (default = 97)\n'''
    '''--memory              memory limit (MB) for CD-HIT (default = 2000, unlimited = 0)\n'''
    '''--min-orf (-l)        minimum ORF length (bp) to evaluate (default = 300)\n'''
    '''--max-orf (-ml)       maximum ORF length (bp) to evaluate (default = 10000)\n'''
    '''--evalue (-e)         maximum e-value to infer reference ORFs\n'''
    '''                      (default = 1e-30)\n'''
    '''--gencode (-g)        genetic code to use to translate ORFs\n'''
    '''--strand (-s)         query strands to call ORFs\n'''
    '''                      (both/minus/plus, default = both)'''))

    porf.add_argument('--db', '-d', action = 'store',
        metavar = '[Protein Database]', type = str,
        help = argparse.SUPPRESS)

    porf.add_argument('--partial','-p', action = 'store_true',
        help = argparse.SUPPRESS)

    porf.add_argument('--pid', '--id', '-id', action = 'store',
        default = 97, metavar = '[perc-identity]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--memory', action = 'store',
        default = 2000, metavar = '[memory-limit]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--min-orf', '-l', action = 'store',
        default = 300, metavar = '[bp]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--max-orf', '-ml', action = 'store',
        default = 10000, metavar = '[bp]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--evalue', '-e', action = 'store',
        default = 1e-30, metavar = '[e-value]', type = float,

        help = argparse.SUPPRESS)

    porf.add_argument('--gencode', '-g', action = 'store',
        default = '1', metavar = '[Genetic-Code]', type = str,
        help = argparse.SUPPRESS)

    porf.add_argument('--strand', '-s', choices = ['both', 'plus', 'minus'], action = 'store',
        default = 'both', type = str, help = argparse.SUPPRESS)

    cntm = parser.add_argument_group('Contamination-Calling Options', description = (
    '''--contam (-c)         table of annotated sequences\n'''))

    cntm.add_argument('--contam','-c', nargs = '?', const = True,
        help = argparse.SUPPRESS)

    # Ensures that just script description provided if no arguments provided
    if len(sys.argv[1:]) == 0:
        print(ascii_logo_vsn())
        print(detail_msg())
        print(f'{usage_msg()}\n')
        sys.exit()

    args = parser.parse_args()

    if args.version:
        print(ascii_logo_vsn())
        sys.exit()

    elif not args.taxon or not args.fin:
        print('\nMissing required inputs: FASTA-File and Output Name')
        print(f'\n{usage_msg()}\n')
        sys.exit()

    args.taxon = '_'.join(args.taxon)

    return args


def ascii_logo_vsn():
    alv_msg = """      _____ ___    _     ___
     |_   _|_ _|__| |___/ __|
       | |  | |/ _` / -_)__ \\
       |_| |___\__,_\___|___/

     Version 1.1.4
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


def check_dependencies(args):
    dpnds = [
        ('BioPython', 'Bio.SeqIO'), ('Optuna', 'optuna'), ('Pandas', 'pandas'),
        ('Scikit-Learn', 'sklearn.svm'), ('ete3', 'ete3'), ('Barrnap', 'barrnap'),
        ('CD-HIT', 'cd-hit-est'), ('DIAMOND', 'diamond'), ('Kraken2', 'kraken2')]

    dpnd_status = {}
    for n in range(9):
        if n < 5:
            try:
                __import__(dpnds[n][1])
                dpnd_status[dpnds[n][0]] = 'Check'
            except:
                dpnd_status[dpnds[n][0]] = 'Error'
        else:
            if shutil.which(dpnds[n][1]):
                dpnd_status[dpnds[n][0]] = 'Check'
            else:
                if not args.contam:
                    if not args.kraken:
                        if n == 8:
                            dpnd_status[dpnds[n][0]] = 'Warning'
                        else:
                            dpnd_status[dpnds[n][0]] = 'Error'
                else:
                    if n == 8 and args.kraken:
                        dpnd_status[dpnds[n][0]] = 'Error'
                    else:
                        dpnd_status[dpnds[n][0]] = 'Warning'

    if 'Error' in dpnd_status.values():
        print('Error: Status of all dependencies:')
        for k, v in dpnd_status.items():
            print(f'    {k}:  {v}')
        sys.exit()


def predict_orfs(
        fasta_file: str,
        taxon_code:str,
        dmnd_db: str,
        kraken_db: str,
        partial: bool,
        filter: bool,
        gcode: str,
        kmer: int,
        overlap: bool,
        step: int,
        model: str,
        min_len: int,
        max_len: int,
        pid: float,
        memory: int,
        evalue: float,
        threads: int,
        strand: str,
        verb: bool) -> None:

    """
    Predicts in-frame Open Reading Frames (ORFs) from a given transcriptome.

    Parameters
    ----------
    fasta_file:  FASTA formatted transcriptome file
    taxon_code:  species/taxon name or abbreviated code
    dmnd_db:     path to a protein database for DIAMOND
    partial:     evaluate partial ORFs
    gcode:       genetic code used for ORF-calling and translation steps
    kmer:        kmer size (number of nt) to convert sequences into
    overlap:     overlapping kmers permitted if True
    step:        step-size for overlapping kmers
    model:       pickle file from previous TIdeS run
    min_len:     minimum ORF length to consider
    max_len:     maximum ORF length to consider
    pid:         percent identity (0-1.0) for removing redundant sequences
    memory:      memory limit (MB) for CD-HIT
    evalue:      maximum e-value to keep hits from DIAMOND
    threads:     number of threads to use
    strand:      designate strand(s) for ORF calling
    verb:        verbose print statements

    Returns current time to track overall runtime.
    """

    sttime = time.time()
    ref_orfs = ref_orf_fas = cvec = clf = None
    stop_codons, rstop_codons, ttable = oc.gcode_start_stops(gcode)

    if overlap and not step:
        step = int(kmer/2)

    """NEED to compress log of options into pickle file too ... then can pass all
    relevant info to the appropriate functions, like the kmer, overlap, step, blah blah..."""

    if verb:
        print('#------ Preparing Transcriptome Data -------#')

    filt_fas = ft.filter_transcripts(
                    fasta_file,
                    taxon_code,
                    sttime,
                    kraken_db,
                    filter,
                    min_len,
                    max_len,
                    threads,
                    pid,
                    memory,
                    verb
                    )


    if verb:
        if not model:
            print('\n#------- Calling Training and pORFs --------#')
        else:
            print('\n#------------- Calling pORFs ---------------#')

    putative_orfs, porf_fas = oc.capture_pORFs(
                                    filt_fas,
                                    taxon_code,
                                    sttime,
                                    gcode,
                                    min_len,
                                    partial,
                                    strand,
                                    verb
                                    )

    if not model:
        ref_orfs = oc.generate_ref_orfs(
                                    filt_fas,
                                    taxon_code,
                                    sttime,
                                    dmnd_db,
                                    min_len,
                                    evalue,
                                    threads,
                                    verb
                                    )


    if model:
        with open(model, 'rb') as f:
            overlap, kmer, step, cvec, clf = pickle.load(f)

    if verb:
        if not model:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Preparing training and query ORFs for {taxon_code}')
        else:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Preparing query ORFs for {taxon_code}')

    train_data, query_data, cvec = op.kmer_ngram_counts(
                                    ref_orfs,
                                    putative_orfs,
                                    taxon_code,
                                    partial,
                                    cvec,
                                    False,
                                    overlap,
                                    kmer,
                                    None
                                    )

    if verb:
        print('\n#----------- ORF Classification ------------#')

    clf_summary, clf = co.classify_orfs(
                        taxon_code,
                        sttime,
                        train_data,
                        query_data,
                        threads,
                        clf,
                        False,
                        verb
                        )

    if verb:
        print('\n#---------- Saving TIdeS Outputs -----------#')
        print(f'[{timedelta(seconds=round(time.time()-sttime))}] '
                ' Making FASTA files and storing TIdeS model')

    sp.save_model(
        taxon_code,
        overlap,
        kmer,
        step,
        cvec,
        clf
        )

    sp.save_seqs(
        taxon_code,
        putative_orfs,
        clf_summary,
        ttable,
        True
        )

    return sttime


def eval_contam(fasta_file: str,
                taxon_code: str,
                contam_list: str,
                kraken_db: list,
                gcode: str,
                kmer: int,
                overlap: bool,
                step: int,
                model: str,
                threads: int,
                verb: bool) -> None:

    """
    Binary classification of target versus non-target sequences

    Parameters
    ----------
    fasta_file: FASTA formatted transcriptome file
    taxon_code: species/taxon name or abbreviated code
    sister_summary: tab-delimited file with user-defined "target" and "non-target" sequences
    gcode: genetic code used for ORF-calling and translation steps
    pretrained: random forest model from previous TIdeS run
    min_len: minimum ORF length to consider
    threads: number of cpu threads to use
    verb: verbose print statements

    Returns current time to track overall runtime
    """


    sttime = time.time()
    cvec = clf = None
    stop_codons, rstop_codons, ttable = oc.gcode_start_stops(gcode)

    if overlap and not step:
        step = int(kmer/2)

    if verb:
        print('#---- Preparing User-Assessed ORF Data -----#')

    train_orfs, query_orfs = ft.prep_contam(
                                fasta_file,
                                taxon_code,
                                contam_list,
                                kraken_db,
                                model,
                                sttime,
                                threads,
                                verb
                                )

    if model:
        with open(model, 'rb') as f:
            overlap, kmer, step, cvec, clf = pickle.load(f)

    if verb:
        if not model:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Preparing training and query ORFs for {taxon_code}')
        else:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Preparing query ORFs for {taxon_code}')

    train_data, query_data, cvec = op.kmer_ngram_counts(
                                    train_orfs,
                                    query_orfs,
                                    taxon_code,
                                    False,
                                    cvec,
                                    True,
                                    overlap,
                                    kmer,
                                    step
                                    )


    if verb:
        print('\n#----------- ORF Classification ------------#')

    clf_summary, clf = co.classify_orfs(
                        taxon_code,
                        sttime,
                        train_data,
                        query_data,
                        threads,
                        clf,
                        True,
                        verb
                        )

    if verb:
        print('\n#---------- Saving TIdeS Outputs -----------#')
        print(f'[{timedelta(seconds=round(time.time()-sttime))}] '
                ' Making FASTA files and storing TIdeS model')

    sp.save_model(
        taxon_code,
        overlap,
        kmer,
        step,
        cvec,
        clf
        )

    sp.save_seqs(
        taxon_code,
        query_orfs,
        clf_summary,
        ttable,
        True,
        True
        )

    return sttime

def main():

    """Note to self:
    Add method to check the sister-summary table format...

    Requirements:
    -- Tab-delimited
    -- Target/Non-Target OR 1/0

    Also provide options for multi-classes (e.g., not solely current binary approach)...
    """

    args = collect_args()

    if not oc.gcode_start_stops(args.gencode):
        sys.exit(1)

    check_dependencies(args)

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
                                args.kraken,
                                args.partial,
                                args.no_filter,
                                args.gencode,
                                args.kmer,
                                args.overlap,
                                args.step,
                                args.model,
                                args.min_orf,
                                args.max_orf,
                                args.pid,
                                args.evalue,
                                args.threads,
                                args.strand,
                                not args.quiet)
    else:
        if not args.quiet:
            print(ascii_logo_vsn())

        if not isinstance(args.contam, str) and not args.kraken and not args.model:
            print(args)
            print('ERROR: Please include the path to a file of annotated sequences to classify OR\n' \
                'the path to a suitable Kraken2 database for automatic annotation of\n' \
                'non-eukaryotic sequences.\n')
            sys.exit(1)

        elif isinstance(args.contam, str) and not args.kraken and not Path(args.contam).is_file():
            print(f'Could not find the file {args.contam}. Please check that the correct\n' \
                'path is provided.\n')
            sys.exit(1)

        sttime = eval_contam(args.fin,
                            args.taxon,
                            args.contam,
                            args.kraken,
                            args.gencode,
                            args.kmer,
                            args.overlap,
                            args.step,
                            args.model,
                            args.threads,
                            not args.quiet)

    if args.clean:
        if not args.quiet:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Removing intermediate files')

        if Path(f'{args.taxon}_TIdeS/Filter_Steps/').exists():
            os.system(f'rm -rf {args.taxon}_TIdeS/Filter_Steps')
        if Path(f'{args.taxon}_TIdeS/ORF_Calling/').exists():
            os.system(f'rm -rf {args.taxon}_TIdeS/ORF_Calling/')

    if args.gzip:
        if not args.quiet:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Compressing TIdeS Outputs')

        os.system(f'tar -zcf {args.taxon}_TIdeS.tar.gz {args.taxon}_TIdeS/')

    if not args.quiet:
        print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Finished running TIdeS')
        print('\n#--------- Thanks for using TIdeS! ---------#\n')
