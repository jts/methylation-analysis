import sys
import argparse
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import MutableSeq

description = """
Modify a reference genome to methylate positions according to the input recognition set
"""

parser = argparse.ArgumentParser(description=description, epilog='')
parser.add_argument('input', action='store', help='the input reference file')
parser.add_argument('--recognition', action='store', help='the recognition set')
args = parser.parse_args()

recognition_sites = list()
recognition_sites_methylated = list()

if args.recognition == "CpG":
    recognition_sites = ["CG"]
    recognition_sites_methylated = ["MG"]
elif args.recognition == "dam":
    recognition_sites = ["GATC"]
    recognition_sites_methylated = ["GMTC"]
elif args.recognition == "dcm":
    recognition_sites = ["CCAGG", "CCTGG"]
    recognition_sites_methylated = ["CMAGG", "CMTGG"]
else:
    sys.stderr.write("unknown recognition: " + args.recognition)
    sys.exit(1)

recognition_length = len(recognition_sites[0])

# 
for rec in SeqIO.parse(args.input, "fasta"):
    outseq = rec.seq.tomutable()
    for bi in xrange(0, len(rec) - 1):

        for s,m in zip(recognition_sites, recognition_sites_methylated):
            if str(rec.seq[bi:bi + recognition_length]) == s:
                outseq[bi:bi + recognition_length] = m
    rec.seq = outseq
    SeqIO.write(rec, sys.stdout, "fasta")
