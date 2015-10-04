import sys
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import MutableSeq

# Change all CG dinucleotides to MG
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    outseq = rec.seq.tomutable()
    for bi in xrange(0, len(rec) - 1):
        if str(rec.seq[bi:bi+2]) == "CG":
            outseq[bi] = 'M'
    rec.seq = outseq
    SeqIO.write(rec, sys.stdout, "fasta")
