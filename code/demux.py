from Bio import SeqIO
import os
import gzip
from pathlib import Path
from argparse import ArgumentParser
import sys

# see https://github.com/timoast/basc/blob/main/scripts/demux.py
# https://doi.org/10.1016/j.molcel.2021.09.019

# Format:
# cell barcode in read name (ran sinto barcode)
# Tn5 barcode in first 8 bp of read1 and read2
# ME sequence after Tn5 barcode needs to be removed (remove first 42 bp according to the paper)

# parse args
parser = ArgumentParser(description='Read demultiplexer')
parser.add_argument('--read1', help='Path to read1')
parser.add_argument('--read2', help='Path to read2')
parser.add_argument('--tn5_i5', help='Path to Tn5 i5 index FASTA file')
parser.add_argument('--tn5_i7', help='Path to Tn5 i7 index FASTA file')
parser.add_argument('--output', help='Path to output directory')
args = parser.parse_args()


def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def open_fastq(f):
    if os.path.isfile(f):
        if f.endswith(".gz"):
            f_open = gzip.open(f, "rt") # rt read text mode (decodes binary gzip)
        else:
            f_open = open(f, "r")
        return(f_open)
    else:
        raise Exception("File not found")


def read_sample_barcodes(inpath):
    all_bc = dict()
    for i in SeqIO.parse(open(inpath),'fasta'):
        all_bc[str(i.seq)] = i.id
    return(all_bc)


def get_entry(f):
    return([f.readline(), f.readline(), f.readline(), f.readline()])


def extract_sequences(sequence, bc_len=8, trim=42):
    tn5_bc = sequence[:bc_len]
    trimmed = sequence[trim:]
    return((tn5_bc, trimmed))


# Load barcode sequences from fasta file
tn5_barcodes_i5 = read_sample_barcodes(args.tn5_i5)
tn5_barcodes_i7 = read_sample_barcodes(args.tn5_i7)

r1 = open_fastq(f=args.read1)
r2 = open_fastq(f=args.read2)

# Create output files for each barcode combination
outf = dict()
outpath = Path(args.output)
if not outpath.exists():
    os.mkdir(outpath)

i5_marks = [tn5_barcodes_i5[x] for x in tn5_barcodes_i5.keys()]
i7_marks = [tn5_barcodes_i7[x] for x in tn5_barcodes_i7.keys()]
i5_marks.append("unknown")
i7_marks.append("unknown")

for i5 in i5_marks:
    for i7 in i7_marks:
        bc =  i5 + "-" + i7
        # outputting uncompressed fastq is ~10x faster
        fname_1 = open(outpath / (bc + ".R1.fastq"), "w+")
        fname_2 = open(outpath / (bc + ".R2.fastq"), "w+")
        outf[bc] = (fname_1, fname_2)

        
# iterate over reads
x = 0
while True:
    r1_entry = get_entry(r1)
    r2_entry = get_entry(r2)

    if r1_entry[0] == '':
        break

    # Extract Tn5 barcode sequences for each read
    # first element = tn5 barcode, second element = trimmed read
    i5_seqs = extract_sequences(sequence=r1_entry[1])
    i7_seqs = extract_sequences(sequence=r2_entry[1])
    
    if i5_seqs[0] in tn5_barcodes_i5.keys():
        i5_mark = tn5_barcodes_i5[i5_seqs[0]]
    else:
        # compute hamming distance
        hams = [hamming_distance(i5_seqs[0], x) < 2 for x in tn5_barcodes_i5.keys()]
        if sum(hams) == 1:
            i5_mark = tn5_barcodes_i5[list(tn5_barcodes_i5.keys())[hams.index(True)]]
        else:
            i5_mark = "unknown"
    
    if i7_seqs[0] in tn5_barcodes_i7.keys():
        i7_mark = tn5_barcodes_i7[i7_seqs[0]]
    else:
        # compute hamming distance
        hams = [hamming_distance(i7_seqs[0], x) < 2 for x in tn5_barcodes_i7.keys()]
        if sum(hams) == 1:
            i7_mark = tn5_barcodes_i7[list(tn5_barcodes_i7.keys())[hams.index(True)]]
        else:
            i7_mark = "unknown"
    
    # insert trimmed read
    r1_entry[1] = i5_seqs[1]
    r2_entry[1] = i7_seqs[1]

    # write to file according to barcodes
    outfile = i5_mark + "-" + i7_mark
    r1_outf = outf[outfile][0]
    r2_outf = outf[outfile][1]

    r1_outf.write("".join(r1_entry))
    r2_outf.write("".join(r2_entry))
    
    x += 1
    if x % 1e6 == 0:
        print("Processed " + str(int(x/1e6)) + " million reads", file=sys.stderr, end="\r")

# close all files
for i in outf.keys():
    outf[i][0].close()
    outf[i][1].close()