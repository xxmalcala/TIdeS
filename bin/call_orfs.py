from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import re, sys

"""Currently grabs just the complete open reading frames (start + stop codons).

Will update to be more "flexible" at a later date."""

# Requires a translation table (integer). Returns list of associated stop codons.
def eval_gcode_ttable(translation_table):
    ncbi_link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi'
    supported_tables = {'1': ['TGA','TAA','TAG','TCA','TTA','CTA'],
                        '4': ['TAA','TAG','TTA','CTA'],
                        '6': ['TGA','TCA'], '10': ['TAA','TAG','TTA','CTA'],
                        '12': ['TGA','TAA','TAG','TCA','TTA','CTA'],
                        '26': ['TGA','TAA','TAG','TCA','TTA','CTA'],
                        '29': ['TGA','TCA'],'30': ['TGA','TCA']}

    if translation_table not in supported_tables.keys():
        print('\nCheck the chosen genetic code/translation table.' \
            '\nSupported tables include:\n\n    '+', '.join(supported_tables.keys()))
        print(f'\nFor more information, visit: {ncbi_link}\n')
        sys.exit(1)

    else:
        return supported_tables[f'{translation_table}']

# Identifies putative Open Reading Frames (ORFs) from a given pair of coordinates.
def extract_pORF(seq, fp, tp, orf_num, strand, comp, min_len):

    if abs(tp - fp) >= min_len:
        if strand == '+':
            initial_pORF = f'{seq.seq[fp:tp]}'
        else:
            initial_pORF = f'{seq.seq[fp:tp].reverse_complement()}'

        # Hunts for a start codon, so complete ORFs (start + stop codons) can be returned.
        updated_pORF, st_pos = check_start_pos(initial_pORF, min_len)

        # Checks that ORF is complete. Coded as "partial" pORFs may be supported
        # later. Unsure if/when that will happen.
        if comp:
            if updated_pORF[:3] == 'ATG' and updated_pORF[-3:] in ['TAA','TAG','TGA']:
                orf_num += 1

                if strand == '+':
                    fin_pORF_name = f'>{seq.id}.pORF{orf_num} {seq.id}:{fp+st_pos+1}-{tp}({strand})\n'
                else:
                    fin_pORF_name = f'>{seq.id}.pORF{orf_num} {seq.id}:{tp-st_pos}-{fp+1}({strand})\n'

                return f'{fin_pORF_name}{updated_pORF}\n', orf_num

            else:
                return None, orf_num

        else:
            return None, orf_num

    else:
        return None, orf_num


# Takes a list of stop codons and returns a list of putative ORF sequences.
def eval_pORFs(seq, stop_pos, orf_num, strand, comp, min_len):
    size_filt_pORFs = []

    # For transcripts/reading frames with a singular in-frame stop codon.
    if len(stop_pos) == 1:
        if strand == '+':
            tp = stop_pos[-1]+3
            fp = tp%3
        else:
            fp = stop_pos[-1]
            tp = len(seq) - len(seq[fp:])%3
        fin_pORF, orf_num = extract_pORF(seq, fp, tp, orf_num, strand, comp, min_len)
        if fin_pORF:
            size_filt_pORFs.append(fin_pORF)

    else:
        # Extract all ORFs from a particular transcript's current reading frame.
        for n in range(len(stop_pos)-1):
            if strand == '+':
                fp, tp = stop_pos[n]+3, stop_pos[n+1]+3
            else:
                fp, tp = stop_pos[n], stop_pos[n+1]
            fin_pORF, orf_num = extract_pORF(seq, fp, tp, orf_num, strand, comp, min_len)
            if fin_pORF:
                size_filt_pORFs.append(fin_pORF)

            # Basis for partial transcripts here. Unsupported.
            # else:
            #     orf_num += 1
            #     if strand == '+':
            #         fin_pORF_name = f'>{seq.id}.pORF{orf_num} {seq.id}:{fp+st_pos+1}-{tp}({strand})\n'
            #     else:
            #         fin_pORF_name = f'>{seq.id}.pORF{orf_num} {seq.id}:{tp-st_pos}-{fp+1}({strand})\n'
            #     size_filt_pORFs.append(f'{fin_pORF_name}{updated_pORF}\n')

    return size_filt_pORFs, orf_num

# Takes a sequence (as a string) and minimum length. Returns the putative ORF.
def check_start_pos(pORF, min_len):
    all_start = [m.start() for m in re.finditer('ATG', pORF) if m.start()%3 == 0]
    if all_start:
        # Ensure new ORF still meets criterion.
        if len(pORF[all_start[0]:]) >= min_len and len(pORF[all_start[0]:])%3 == 0:
            return pORF[all_start[0]:], all_start[0]
        else:
            return pORF, 0
    else:
        return pORF, 0

# Takes a transcript and its stop codons, returns positions of all stop codons
# for each reading frame.
def collect_stop_codons(seq, stp_cdns):
    stop_codon_pos = {'1':[],'2':[],'3':[], '4':[], '5':[], '6':[], 'what':[]}
    rf_adj = 0

    for n in range(len(stp_cdns)):
        stp_pos = [m.start() for m in re.finditer(stp_cdns[n], f'{seq.seq}')]
        if stp_pos:
            # Adjust the reading-frame based on number of stop codons -- needed
            # for alternative genetic codes (e.g. translation table 6)
            if n > (len(stp_cdns)/2 - 1):
                rf_adj = 3
                # print(stp_cdns[n], rf_adj)
            for p in stp_pos:
                if p%3 == 0:
                    stop_codon_pos[f'{1+rf_adj}'].append(p)
                elif (p-1)%3 == 0:
                    stop_codon_pos[f'{2+rf_adj}'].append(p)
                elif (p-2)%3 == 0:
                    stop_codon_pos[f'{3+rf_adj}'].append(p)
                else:
                    stop_codon_pos['what'].append(p)

    stop_codon_pos = {k: sorted(stop_codon_pos[k]) for k, v in stop_codon_pos.items() if v}

    return stop_codon_pos

# Grabs all pORFs from a transcript.
def all_orfs_from_seq(seq, stp_cdns, comp, min_len):
    all_orfs_from_seq = []
    orf_num = 0
    stop_codons = collect_stop_codons(seq, stp_cdns)

    for k, v in stop_codons.items():
        if int(k) > 3:
            strand = '-'
        else:
            strand = '+'

        orfs_from_rf, orf_num = eval_pORFs(seq, v, orf_num, strand, comp, min_len)
        all_orfs_from_seq += orfs_from_rf

    return all_orfs_from_seq

# Longest is UNSUPPORTED
# Gets all pORFs for all transcripts in a given FASTA file.
def call_orfs(fasta_file, out_dir, ttable = 1, comp = True, longest = False, min_len = 300):

    Path(f'{out_dir}ORF_Calling/').mkdir(parents = True, exist_ok = True)

    all_orf_fas = f'{fasta_file.rpartition("/")[-1].rpartition(".")[0]}.CompORFs.fas'
    stp_cdns = eval_gcode_ttable(f'{ttable}')
    all_orfs_from_transcriptome = []

    for i in SeqIO.parse(fasta_file,'fasta'):

        if len(i.seq) >= min_len:
            all_orfs_from_transcriptome += all_orfs_from_seq(i, stp_cdns, comp, min_len)

    with open(f'{out_dir}ORF_Calling/{all_orf_fas}','w+') as w:
        w.write(''.join(all_orfs_from_transcriptome))

    return f'{out_dir}ORF_Calling/{all_orf_fas}'
