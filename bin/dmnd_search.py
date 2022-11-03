import subprocess, sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

"""Script identifies 'high-quality' (e-value <= 1e-30) ORFs from the "BLASTX" comparisons of a
transcriptome against a given protein database.

These identified ORFs will then be used to train the ML models to discern the
most appropriate ORFs from the rest of the transcriptome.

Dependencies include: DIAMOND."""

# Run "BLASTX" with the transcriptome as the query against a given protein database.
def diamond_blastx(rRNA_filt_fasta, taxon_code, dmnd_db, out_dir, threads = 4, eval = 1e-30):
    dmnd_dir = f'{out_dir}/ORF_Calling/'
    dmnd_tsv = f'{dmnd_dir}{taxon_code}.DMND_REF_ORFs.tsv'
    Path(dmnd_dir).mkdir(parents = True, exist_ok = True)

    dmnd_cmd = 'diamond blastx ' \
            f'-q {rRNA_filt_fasta} ' \
            f'-d {dmnd_db} ' \
            f'-o {dmnd_tsv} ' \
            f'-e {eval} ' \
            f'-p {threads} ' \
            '-k 1 ' \
            '--query-cover 0.6 ' \
            '--subject-cover 0.6 ' \
            '-f 6 qseqid sseqid qstart qend qframe'

    dmnd_call = subprocess.call(dmnd_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return dmnd_tsv

# Extracts the (nucleotide) alignments from the DIAMOND outputs, from the transcriptome.
# Will be reference for 'in-frame' ORFs.
def extract_orf_hits(rRNA_filt_fasta, taxon_code, dmnd_db, out_dir, threads, evalue):

    dmnd_tsv = diamond_blastx(
                rRNA_filt_fasta,
                taxon_code,
                dmnd_db,
                out_dir,
                threads,
                evalue)

    hit_info = {i.split('\t')[0]:[int(j) for j in i.split('\t')[2:]] for i in open(dmnd_tsv).readlines()}
    seq_db = {i.id:[f'{i.seq}'] for i in SeqIO.parse(rRNA_filt_fasta,'fasta') if i.id in hit_info.keys()}

    for k,v in hit_info.items():
        qs = min(v[:2])-1
        qe = max(v[:2])
        seq = seq_db[k][0][qs:qe]

        if v[-1] < 0:
            seq = f'{Seq(seq).reverse_complement()}'

        seq_db[k].append(seq)

    with open(dmnd_tsv.replace("tsv","fas"),'w+') as w:
        for k, v in seq_db.items():
            w.write(f'>{k}\n{v[-1]}\n')

    return dmnd_tsv.replace("tsv","fas")
