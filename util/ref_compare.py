#!/usr/bin/env python3

import os, sys


def diamond_blast(query_fas, ref_fas, outname) -> str:
    out_tsv = f'{outname}.BLASTX_RefDB.tsv'
    stringent_tsv = out_tsv.replace('.tsv','.Stringent.tsv')
    dmnd_cmd_stringent = f'diamond blastx -q {query_fas} -d {ref_fas} -k 1 -o {stringent_tsv} ' \
                '--very-sensitive -f 6 qseqid sseqid qlen slen ' \
                'length pident qframe --quiet --query-cover 70 --subject-cover 70'

    dmnd_cmd = f'diamond blastx -q {query_fas} -d {ref_fas} -k 1 -o {out_tsv} ' \
                '--very-sensitive -f 6 qseqid sseqid qlen slen ' \
                'length pident qframe --quiet'

    os.system(dmnd_cmd)
    os.system(dmnd_cmd_stringent)

    return out_tsv, stringent_tsv

def parse_hits(default_hits, strict_hits, qseqs, outname) -> None:

    d_unique_txp_hits = [(i.split('.pORF')[0], int(i.split('\t')[-1])) for i in open(default_hits).readlines()]
    s_unique_txp_hits = [(i.split('.pORF')[0], int(i.split('\t')[-1])) for i in open(strict_hits).readlines()]

    d_uth = len(set([i[0] for i in d_unique_txp_hits]))
    s_uth = len(set([i[0] for i in s_unique_txp_hits]))

    d_crf_hits = len(set([i[0] for i in d_unique_txp_hits if i[-1] == 1]))
    s_crf_hits = len(set([i[0] for i in s_unique_txp_hits if i[-1] == 1]))

    with open(f'{outname}.CRF.Summary.tsv','w+') as w:
        w.write('Total-Query-Seqs\tStringency\tUnique-Transcript-Hits\tUnique-Txp-CRFs\n')
        w.write(f'{qseqs}\tDefault\t{d_uth}\t{d_crf_hits}\n')
        w.write(f'{qseqs}\tStrict\t{s_uth}\t{s_crf_hits}\n')

    print(f'Total query seqs:                 {qseqs}')
    print(f'Unique transcript hits (default): {d_uth}')
    print(f'Unique txpt CRF hits (default):   {d_crf_hits}')
    print(f'Unique transcript hits (srtict):  {s_uth}')
    print(f'Unique txpt CRF hits (strict):    {s_crf_hits}')

if __name__ == '__main__':
    try:
        query_fas = sys.argv[1]
        ref_fas = sys.argv[2]
        outname = sys.argv[3]
    except:
        print('Usage:\n    python3 ref_compare.py [QUERY-FASTA] [REF-FASTA] [OUTPUT-NAME]\n')
        sys.exit()

    print(f'\n{outname}')
    default_hits, strict_hits = diamond_blast(query_fas, ref_fas, outname)
    qseqs = len([i for i in open(query_fas).readlines() if i[0] == '>'])

    parse_hits(default_hits, strict_hits, qseqs, outname)
