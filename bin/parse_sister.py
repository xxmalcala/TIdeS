import shutil, sys
from pathlib import Path

from Bio import SeqIO

"""Script parses a phylogenetic "sisters" summary table generated from walking
through a series of phylogenetic trees.

It is recommended to use the PhyloToL naming convention and associated tree-walking
scripts to assist in generating this "sisters" summary table."""

# Option for specific relationships?

# Requires a taxon name and the phylogenetic sisters table for (ideally) several
# hundred phylogenetic trees. Returns "clean" and putative contaminant sequence
# names.

def parse_sister(taxon_code, sister_txt, contam = 'All', short = True):
    q_tax = [line for line in open(sister_txt).readlines() if line.startswith(taxon_code)]
    q_tax_ref, q_tax_contam = [],[]

    for i in q_tax:
        if 'same-minor' in i:
            q_tax_ref.append((i.split('\t')[2],i))
        if i.endswith('short\n'):
            # Ignore non-monophyletic?
            if 'non-monophyletic' not in i:
                # Adds all sequences that are not found in the same minor clade.
                if contam == 'All':
                    q_tax_contam.append((i.split('\t')[2],i))
                # Only adds sequences if they are found to be "sister" to any
                # discrete contaminating taxa/clades (e.g. EE_cr_Gthe, and/or Op_fu)
                elif any(cntm_tax in i for cntm_tax in contam.split(',')):
                    q_tax_contam.append((i.split('\t')[2],i))
    if len(q_tax_ref) < 100 or len(q_tax_contam) < 20:
        print('\nWARNING: Check your phylogenetic sister-relationship summary table.')
        print('-------- There are too few data to train accurately.')
        print('-------- Check your taxon-code/name (common issue) before proceeding.')
        print('\n\nQuitting TIdeS.\n')
        sys.exit()
    return q_tax_ref, q_tax_contam

# Simply backs up the dataset. Needs the transcriptome, taxon name, and phylogenetic
# sisters table. Outputs the directory for TIdeS to work in.
def backup_data(fasta_file, taxon_code, sister_txt,  min_len):
    out_dir = f'{taxon_code}_TIdeS/'
    Path(f'{out_dir}Original').mkdir(parents = True, exist_ok = True)

    min_len_fas = f'{out_dir}Original/{taxon_code}.{min_len}bp.fas'

    shutil.copy2(fasta_file, f'{out_dir}Original/')
    shutil.copy2(sister_txt, f'{out_dir}Original/')

    with open(min_len_fas,'w+') as w:
        for i in SeqIO.parse(fasta_file, 'fasta'):
            if len(i.seq) >= min_len:
                w.write(f'>{i.id}\n{i.seq}\n')

    return min_len_fas, out_dir

# Requires the transcriptome, taxon name, and phylogenetic sisters table.
# Returns a single FASTA file, with clean and contaminant sequences marked, and
# the working directory for TIdeS.
def bin_seqs(fasta_file, taxon_code, sister_txt, min_len = 300):
    mlen_fas, out_dir = backup_data(fasta_file, taxon_code, sister_txt, min_len)

    ctg_dir = f'{out_dir}Categorized/'
    Path(f'{ctg_dir}').mkdir(parents = True, exist_ok = True)

    ref_cntm_fas = f'{ctg_dir}{taxon_code}.RefContam.fas'

    ref, contam = parse_sister(taxon_code, sister_txt)

    # put some check for min number of refs/contams ... ~100 each?
    ref_seq_names = [i[0] for i in ref]
    contam_seq_names = [i[0] for i in contam]

    with open(ref_cntm_fas,'w+') as w:
        for i in SeqIO.parse(fasta_file,'fasta'):
            if len(i.seq) >= min_len:
                if i.id in ref_seq_names:
                    w.write(f'>{i.id}_Ref\n{i.seq}\n')
                elif i.id in contam_seq_names:
                    w.write(f'>{i.id}_Contam\n{i.seq}\n')

    return mlen_fas, ref_cntm_fas, out_dir
