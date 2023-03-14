#!/usr/bin/env python3

import os, sys


def diamond_blast(query_fas, ref_fas, outname) -> str:
    out_tsv = f'{outname}.BLASTX_RefDB.tsv'
    dmnd_cmd = f'diamond blastx -q {query_fas} -d {ref_fas} -k 1 -o {out_tsv} ' \
                '-f 6 qseqid sseqid qlen slen length pident qframe --quiet'
    os.system(dmnd_cmd)
    return out_tsv

def parse_hits(query_fas: str, ref_fas: str, outname: str) -> None:
    dmnd_hits = diamond_blast(query_fas, ref_fas, outname)

    qseqs = len([i for i in open(query_fas).readlines() if i[0] == '>'])
    unique_txp_hits = [(i.split('.pORF')[0], int(i.split('\t')[-1])) for i in open(dmnd_hits).readlines()]
    uth = len(set([i[0] for i in unique_txp_hits]))
    crf_hits = len(set([i[0] for i in unique_txp_hits if i[-1] == 1]))

    with open(f'{outname}.CRF.Summary.tsv','w+') as w:
        w.write('Total-Query-Seqs\tUnique-Transcript-Hits\tUnique-Txp-CRFs\n')
        w.write(f'{qseqs}\t{uth}\t{crf_hits}\n')

    print(f'Total query seqs:       {qseqs}')
    print(f'Unique transcript hits: {uth}')
    print(f'Unique txpt CRF hits:   {crf_hits}')

if __name__ == '__main__':
    try:
        query_fas = sys.argv[1]
        ref_fas = sys.argv[2]
        outname = sys.argv[3]
    except:
        print('Usage:\n    python3 ref_compare.py [QUERY-FASTA] [REF-FASTA] [OUTPUT-NAME]\n')
        sys.exit()

    parse_hits(query_fas, ref_fas, outname)
