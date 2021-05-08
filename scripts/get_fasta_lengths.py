from Bio import SeqIO
import argparse


parser= argparse.ArgumentParser(description= 'get lengths of a fasta as a 2 column csv id,length')
parser.add_argument('--inFasta', help='input fasta')
parser.add_argument('--outFile', help='output file')
args=parser.parse_args()

with open(args.outFile , 'w+') as of:
    fa=SeqIO.parse(args.inFasta, 'fastq')
    [of.write(','.join([record.id, str(len(record) )] ) +'\n' ) for record in fa]





