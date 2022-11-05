import random, sys
from Bio import SeqIO


def get_seq(seq, rf):
    temp = f'{seq.seq}'[rf-1:]
    return temp[:len(temp)-len(temp)%3]

def change_reading_frame(seq, shuff = False):
    rfs = [-3,-2,-1, 2, 3]
    rf = random.sample(rfs, 1)[0]
    if rf > 0:
        new_seq = get_seq(seq,rf)
        if shuff:
            new_seq = list(new_seq)
            random.shuffle(new_seq)
            new_seq = ''.join(new_seq)
        final_seq = f'>{seq.id}_RF{rf}\n{new_seq}\n'
    else:
        new_seq = get_seq(seq.reverse_complement(), abs(rf))
        if shuff:
            new_seq = list(new_seq)
            random.shuffle(new_seq)
            new_seq = ''.join(new_seq)
        final_seq = f'>{seq.id}_RF{abs(rf)+3}\n{new_seq}\n'
    return final_seq

def gen_rand_orientation(fasta_file):
    initial_seqs = [seq_rec for seq_rec in SeqIO.parse(fasta_file, 'fasta')]
    random.shuffle(initial_seqs)
    quarter_initial = int(len(initial_seqs)*.25)
    half_initial = int(len(initial_seqs)*.5)
    final_seqs = [f'>{seq_rec.id}_RF1\n{seq_rec.seq}\n' for seq_rec in initial_seqs[:half_initial]]
    # final_seqs += [change_reading_frame(seq_rec, True) for seq_rec in initial_seqs[quarter_initial:-quarter_initial]]
    final_seqs += [change_reading_frame(seq_rec) for seq_rec in initial_seqs[-half_initial:]]


    with open(f'{fasta_file.split(".fa")[0]}.RandomRF.fas','w+') as w:
        w.write(''.join(final_seqs))

    return f'{fasta_file.split(".fa")[0]}.RandomRF.fas'

if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1]
    except IndexError:
        print('Usage:\n\n    python random_orientaton.py [FASTA-FILE]\n')
        sys.exit(1)

    gen_rand_orientation(fasta_file)
