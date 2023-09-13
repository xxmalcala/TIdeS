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

from Bio import SeqIO

from bin import filt_seqs as ft
from bin import orf_call as oc
from bin import orf_prep as op
from bin import classify_orfs as co
from bin import save_preds as sp


def collect_args():

    parser = argparse.ArgumentParser(description = f'{ascii_logo_vsn()}\n{usage_msg()}',
            usage=argparse.SUPPRESS, add_help = False,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    g = parser.add_argument_group('General Options', description = (
    '''--fin (-f)            input file in FASTA format\n'''
    '''--taxon (-n)          taxon-name\n'''
    '''--threads (-t)        number of CPU threads (default = 1)\n'''
    '''--model (-m)          previously trained TIdeS model (".pkl" file)\n'''
    '''--kmer (-k)           kmer size for generating sequence features (default = 3)\n'''
    '''--overlap (-ov)       permit overlapping kmers (see --kmer)\n'''
    '''--step                step-size for overlapping kmers (default is kmer-length/2)\n'''
    '''--clean               remove intermediate filter-step files\n'''
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

    g.add_argument('--threads','-t', action = 'store',
        default = 1, metavar = '[Threads]', type = int,
        help = argparse.SUPPRESS)

    g.add_argument('--model','-m', action = 'store',
        metavar = '[Trained-RFC]', type = str, default = None,
        help = argparse.SUPPRESS)

    g.add_argument('--kmer','-k', action = 'store',
        default = 3, metavar = '[kmer]', type = int,
        help = argparse.SUPPRESS)

    g.add_argument('--overlap','-ov', action = 'store_true',
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

    porf = parser.add_argument_group('ORF-Calling Options', description = (
    '''--db (-d)             protien database (FASTA or DIAMOND format)\n'''
    '''--partial (-p)        evaluate partial ORFs as well\n'''
    '''--id (-id)            minimum % identity to remove redundant transcripts (default = 97)\n'''
    '''--min-orf (-l)        minimum ORF length (bp) to evaluate (default = 300)\n'''
    '''--max-orf (-ml)       maximum ORF length (bp) to evaluate (default = 10000)\n'''
    '''--evalue (-e)         maximum e-value to infer reference ORFs (default = 1e-30)\n'''
    '''--gencode (-gc)       genetic code to use to translate ORFs\n'''
    '''--strand (-s)         query strands to call ORFs (both/minus/plus, default = both)'''))

    porf.add_argument('--db','-d', action = 'store',
        metavar = '[Protein Database]', type = str,
        help = argparse.SUPPRESS)

    g.add_argument('--partial','-p', action = 'store_true',
        help = argparse.SUPPRESS)

    porf.add_argument('--pid', '--id', '-id', action = 'store',
        default = 97, metavar = '[perc-identity]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--min-orf','-l', action = 'store',
        default = 300, metavar = '[bp]', type = int,
        help = argparse.SUPPRESS)

    porf.add_argument('--max-orf','-ml', action = 'store',
        default = 10000, metavar = '[bp]', type = int,
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

    args.taxon = '_'.join(args.taxon)

    return args


def ascii_logo_vsn():
    alv_msg = """      _____ ___    _     ___
     |_   _|_ _|__| |___/ __|
       | |  | |/ _` / -_)__ \\
       |_| |___\__,_\___|___/

     Version 1.0.0
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
        partial: bool = False,
        gcode: str = '1',
        kmer: int = 3,
        overlap: bool = False,
        step = None,
        model = None,
        min_len: int = 300,
        max_len: int = 10000,
        pid: float = 0.97,
        evalue: float = 1e-30,
        threads: int = 1,
        strand:str = 'both',
        verb: bool = True
        ) -> None:

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
                    min_len,
                    max_len,
                    threads,
                    pid,
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

    # print(list(putative_orfs.keys())[0])
    # print(list(clf_summary[1].values())[0][0])

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
                gcode: str = '1',
                kmer: int = 3,
                overlap = True,
                step = None,
                model: str = None,
                threads: int = 1,
                verb: bool = True) -> None:

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
                                model,
                                sttime,
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


if __name__ == '__main__':

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
                                args.partial,
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
        # print(args)
        # sys.exit()
        if not args.quiet:
            print(ascii_logo_vsn())

        sttime = eval_contam(args.fin,
                            args.taxon,
                            args.contam,
                            args.gencode,
                            args.kmer,
                            args.overlap,
                            args.step,
                            args.model,
                            args.threads,
                            not args.quiet)

    if args.clean:
        os.system(f'rm -rf {args.taxon}_TIdeS/Filter_Steps/')

    if args.gzip:
        if not args.quiet:
            print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Compressing TIdeS Outputs')

        os.system(f'tar -zcf {args.taxon}_TIdeS.tar.gz {args.taxon}_TIdeS/')

    if not args.quiet:
        print(f'[{timedelta(seconds=round(time.time()-sttime))}]  Finished running TIdeS')
        print('\n#--------- Thanks for using TIdeS! ---------#\n')
