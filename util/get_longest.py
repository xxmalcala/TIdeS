
from Bio import SeqIO
from collections import defaultdict
import sys

def get_longest(fasta_file):
    name = fasta_file.split('/')[-1][:10]
    x = defaultdict(str)
    for i in SeqIO.parse(fasta_file, 'fasta'):
        if len(i.seq) > len(x[i.id.split(".p")[0]]):
            x[i.id.split(".p")[0]] = f'{i.seq}'
    with open(f'{name}.LONGEST.fas','w+') as w:
        for k, v in x.items():
            w.write(f'>{k}\n{v}\n')

if __name__ == '__main__':
    if len(sys.argv[1:]) == 1:
        fas = sys.argv[1]
    else:
        print('python3 get_longest.py [FASTA-FILE]')
        sys.exit(1)

    get_longest(fas)
