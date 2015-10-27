#! /usr/bin/env python
import sys

def read_model(fn, print_header):
    fh = open(fn)
    model = dict()

    for line in fh:
        line = line.rstrip()

        # copy then skip header lines
        if line[0] == '#' or line.find("kmer") == 0:
            if print_header:
                print(line)
            continue

        fields = line.split()
        model[fields[0]] = tuple(fields[1:6])
    return model

model_m_fn = sys.argv[1]
model_u_fn = sys.argv[2]

model_m = read_model(model_m_fn, True)
model_u = read_model(model_u_fn, False)

for key in sorted(model_m):
    values = model_m[key]

    if key.find('M') == -1:
        # replace with u model data
        values = model_u[key]
    print(key + "\t" + "\t".join([str(f) for f in values]))

