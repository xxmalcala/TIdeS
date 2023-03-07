import os, shutil, sys, time
from datetime import timedelta
from Bio import SeqIO


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

    import filter_txps as ft
    import orf_call as oc
    import classify_orfs as cl

    if not oc.eval_gcode_ttable(gcode):
        sys.exit(1)

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

    query_orfs, rnd_orfs, stp_cdns = oc.generate_orf_calls(filt_fas,
                                                            taxon_code,
                                                            sttime,
                                                            dmnd_db,
                                                            gcode,
                                                            min_len,
                                                            threads,
                                                            intermed,
                                                            False,
                                                            3,
                                                            verb)

    if verb:
        print('\n#----------- ORF Classification ------------#')
    tides_preds = cl.classify_orfs(taxon_code,
                                    rnd_orfs,
                                    query_orfs,
                                    stp_cdns,
                                    sttime,
                                    pretrained,
                                    threads,
                                    verb,
                                    False)

    if verb:
        print('\n#---------- Saving TIdeS Outputs -----------#')

    tds_pf = f'{taxon_code}_TIdeS/{taxon_code}.TIdeS.fas'
    tds_ps = {}

    with open(tds_pf, 'w+') as w:
        orf_fas = f'{taxon_code}_TIdeS/ORF_Calling/{taxon_code}.{min_len}bp.CompORFs.fas'
        for i in SeqIO.parse(orf_fas,'fasta'):
            if i.id in tides_preds:
                w.write(f'>{i.description}\n{i.seq}\n')
                tds_ps[i.description] = f'{i.seq}'

    with open(tds_pf.replace("fas","AA.fas"), 'w+') as w:
        tds_psa = oc.translate_seqs(tds_ps, gcode)
        for k, v in tds_psa.items():
                w.write(f'{k}\n{v}\n')

    shutil.copy2(tds_pf, f'{taxon_code}_TIdeS/Classified/')

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



    predict_orfs(fasta_file, taxon_code, dmnd_db)
