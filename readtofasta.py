#!/usr/bin/env python
from __future__ import print_function
import argparse
import h5py
import sys
import os

def readevents(filename):
    events = []
    kmers = []
    sds = []
    f = h5py.File(filename,"r")

    calls = "BaseCalled_2D"
    if not "Analyses" in f or not "Basecall_2D_000" in f["Analyses"] or not calls in f["Analyses"]["Basecall_2D_000"]:
        raise KeyError("Not present in HDF5 file")

    if not "Fastq" in f["Analyses"]["Basecall_2D_000"][calls]:
        raise KeyError("Fastq not found in HDF5 file")
    data = f["Analyses"]["Basecall_2D_000"][calls]["Fastq"]
    lines =  data[()].splitlines()
    seq = lines[1].decode('utf-8')
    label = lines[0][1:].strip()
    labelSuffix = "_twodirections" 
    label = label+labelSuffix
    label = label.decode('utf-8')

    return label, seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('infile', type=str)
    parser.add_argument('-o','--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)") 
    parser.add_argument('-l','--linelength', type=int, default=60)

    args = parser.parse_args()

    try:
        if os.path.isfile(args.infile):
            label, sequence = readevents(args.infile)
            label += " "+arg.infile
            seqs = [(label, seq)]
        elif os.path.isdir(args.infile):
            seqs = []
            for f in os.listdir(args.infile):
               pathname = os.path.join(args.infile, f)
               if os.path.isfile(pathname) and pathname.endswith('.fast5'):
                   label, sequence = readevents(pathname)
                   seqs.append((label+" "+pathname, sequence))
        else:
            raise ValueError("Path is not a directory or file: "+args.infile)
    except KeyError:
        sys.exit(1)

    outf = args.output
    infname = args.infile
    llen = args.linelength

    for label, sequence in seqs:
        print(">"+label, file=outf)
        for i in range( len(sequence)//llen ):
            lstart = i*llen
            lend = min(lstart + llen, len(sequence))
            print( sequence[lstart:lend], file=outf )

    sys.exit(0)
