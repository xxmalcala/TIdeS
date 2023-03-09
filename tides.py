#!/usr/bin/env python3

import glob, os, shutil, sys, time
from datetime import timedelta
from Bio import SeqIO

from bin import filter_txps as ft
from bin import orf_call as oc
from bin import classify_orfs as cl
from bin import save_preds as sp


# def check_args() -> class:
#     pass



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

    Definitely "require" a genetic code even for contamination!!! This may help
    infer target vs non-self in weird taxa (e.g. ciliates)...

    Pickle some supporting files as well? Maybe make a log with parameters to
    improve repeatability? e.g. what translation table did we use last week
    """


    # args = check_args()
    # print('Under Construction!')
    # sys.exit()

    if 1 < len(sys.argv[1:]) < 5:
        fasta_file = sys.argv[1]
        taxon_code = sys.argv[2]
        try:
            gcode = sys.argv[3]
            dmnd_db = sys.argv[4]
        except:
            gcode = '1'
            dmnd_db = 'Prot_DB/tides_db.dmnd'

    else:
        print("\nTemporary Usage:\n    python3 tides.py [TRANSCRIPTOME] [TAXON-NAME] [TRANSLATION-TABLE]\n")
        sys.exit()

    if not oc.eval_gcode_ttable(gcode):
        sys.exit(1)


    predict_orfs(fasta_file, taxon_code, dmnd_db)
    #
    # eval_contam(fasta_file, taxon_code, sister_summary)
