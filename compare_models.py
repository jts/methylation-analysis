#! /usr/bin/env python
import sys

def read_model(fn):
    fh = open(fn)
    model = dict()

    for line in fh:
        line = line.rstrip()

        # copy then skip header lines
        if line[0] == '#' or line.find("kmer") == 0:
            continue

        fields = line.split()
        model[fields[0]] = tuple(fields[1:6])
    return model

model_a = read_model(sys.argv[1])
model_b = read_model(sys.argv[2])

for key in sorted(model_a):
    values_a = model_a[key]

    if key in model_b:
        values_b = model_b[key]
        print "\t".join([key, values_a[0], values_b[0], values_a[1], values_b[1], str(float(values_a[0]) - float(values_b[0]))])

