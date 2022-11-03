#!/usr/bin/env python3

import argparse, pickle, shutil, sys
from argparse import RawTextHelpFormatter, SUPPRESS
from pathlib import Path



def collect_args():

    parser = argparse.ArgumentParser(description = '\nTIdeS uses a Random Forest ' \
            'Classifier, trained on the composition of alignable ORFs,\n' \
            'to infer complete putative ORFs (pORFs) from the transcriptome.',
            formatter_class=RawTextHelpFormatter)

    # Set of required options for the script
    required_arg_group = parser.add_argument_group('Input-Output Options')

    required_arg_group.add_argument('--fin', '-i', action = 'store',
        metavar = '[FASTA File]', type = str, required = True,
        help = ' Transcriptome FASTA file.\n\n')

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

    porf_arg_group.add_argument('--stranded','-st', action = 'store_true',
        help = ' Transcripts are stranded (consider RFs +1 to +3).\n\n')

    porf_arg_group.add_argument('--orfs','-orfs', action = 'store',
        metavar = '[Training-ORF FASTA File]', type = str,
        help = ' User-defined "true" ORFs to be used for training.\n\n')

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

    optional_arg_group.add_argument('--train-rfc','-train-rfc', action = 'store_true',
        help = ' Only runs the RFC training.\n\n')

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


def classify_pORFs(args):
    from bin import (filter_trans as ft,
            dmnd_search as dms,
            random_orientation as rando,
            codon_counts as ccnt,
            call_orfs as co,
            rfc_classify as rfcc,
            translate_porfs as tp)

    co.eval_gcode_ttable(f'{args.genetic_code}')

    print(f'\nIntial Filtering of FASTA-FILE: {args.fin.split("/")[-1]}')
    filt_fas, taxon_dir = ft.filter_txpts(
                            args.fin,
                            args.taxon_name,
                            args.min_len,
                            args.threads)

    if not args.rfc:
        if not args.orfs:
            print(f'Generating set of reference ORFs.')
            ref_orfs = dms.extract_orf_hits(
                            filt_fas,
                            args.taxon_name,
                            args.db,
                            taxon_dir,
                            args.threads,
                            args.evalue)
        else:
            Path(f'{taxon_dir}ORF_Calling').mkdir(parents = True, exist_ok = True)
            ref_orfs = f'{taxon_dir}ORF_Calling/{args.orfs}'
            shutil.copy2(args.orfs, ref_orfs)


        print(f'Generating mixed orientation data for: {ref_orfs.split("/")[-1]}')
        ref_rand_orfs = rando.gen_rand_orientation(ref_orfs)

    print('Calling "complete" ORFs.')
    query_orf_fas = co.call_orfs(
                        filt_fas,
                        taxon_dir,
                        args.genetic_code,
                        True,
                        False,
                        args.min_len)

    if not args.rfc:
        print('Generating codon counts for training data.')
        train_orf_tsv = ccnt.codon_counts_fasta(
                                ref_rand_orfs,
                                taxon_dir,
                                args.taxon_name,
                                True)
    else:
        train_orf_tsv = None

    if not args.train_rfc:
        print('Generating codon counts for query data.')
        query_orf_tsv = ccnt.codon_counts_fasta(
                                query_orf_fas,
                                taxon_dir,
                                args.taxon_name,
                                False)
    else:
        query_orf_tsv = None
    # reorganize inputs for consistency...

    rfc_fas = rfcc.classify_seqs(
                train_orf_tsv,
                query_orf_tsv,
                query_orf_fas,
                taxon_dir,
                args.taxon_name,
                args.rfc,
                args.train_rfc)

    if not args.train_rfc:
        rfc_fin_fas, rfc_fin_aa = tp.prep_final_pORFs(
                                    rfc_fas,
                                    args.genetic_code,
                                    taxon_dir)

    return taxon_dir


def classify_contam(args):
    from bin import (codon_counts as ccnt,
                    parse_sister as ps,
                    rfc_classify as rfcc)

    # Parse tree-walking summary and mark known contam and "clean" for training.
    print('\nGathering marked contaminant and "clean" sequences for training TIdeS.')
    mlen_fas, ref_cntm_fas, taxon_dir = ps.bin_seqs(
                                args.fin,
                                args.taxon_name,
                                args.sister_table,
                                args.min_len)

    if not args.rfc:
        print('Generating codon counts for training data.')
        train_orf_tsv = ccnt.codon_counts_fasta(
                                ref_cntm_fas,
                                taxon_dir,
                                args.taxon_name,
                                True, False)
    else:
        train_orf_tsv = None

    if not args.train_rfc:
        query_orf_tsv = ccnt.codon_counts_fasta(
                                mlen_fas,
                                taxon_dir,
                                args.taxon_name,
                                False)
    else:
        query_orf_tsv = None

    # Use Random Forest Classifier to predict contaminant and "clean" genes.
    rfc_fas = rfcc.classify_seqs(
                train_orf_tsv,
                query_orf_tsv,
                mlen_fas,
                taxon_dir,
                args.taxon_name,
                args.rfc,
                args.train_rfc,
                False)

    if not args.train_rfc:
        shutil.copy2(rfc_fas, taxon_dir)
        shutil.copy2(rfc_fas.replace("fas","Contam.fas"), taxon_dir)

    return taxon_dir


if __name__ == '__main__':

    args = collect_args()

    if args.contam:
        print(args.contam)
        # sys.exit()
        taxon_dir = classify_contam(args)

    else:

        if (args.db or args.orfs or args.rfc or args.train_rfc):
            taxon_dir = classify_pORFs(args)
        # if args.db == args.rfc == None and args.train == False:
        else:
            print('\nERROR: If you are trying to use TIdeS to hunt for putative' \
                ' ORFs, make sure to include a protein database or a trained RFC ' \
                'from a prior TIdeS run.\n')
            print('Use "tides.py -h" for more information.\n')
            sys.exit()


    if not args.train_rfc:
        print(f'\nFinal TIdeS predictions can be found: {taxon_dir}\n\nThanks ' \
            'for using TIdeS!\n')
    else:
        print(f'\nFinal TIdeS trained RFC model can be found: {taxon_dir}\n\nThanks ' \
            'for using TIdeS!\n')

    if args.gzip:
        import subprocess
        tar_cmd = f'tar -zcvf {taxon_dir.rstrip("/")}.tar.gz {taxon_dir}'

        dmnd_call = subprocess.call(tar_cmd, shell=True,
            stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
