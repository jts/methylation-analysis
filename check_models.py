#! /usr/bin/env python
import sys
import os
import operator
import itertools
import argparse
import fast5_reader
import math

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

base_fn = [ 'r7.3_complement_median68pA_pop1.model',  'r7.3_complement_median68pA_pop2.model', 'r7.3_template_median68pA.model' ]
base_models = dict()
checked_models = dict()

for fn in base_fn:
    base_models[fn] = read_model(fn)
    checked_models[fn] = 0

# Iterate over reads checking models
for directory in sys.argv[1:]:
    print("Checking", directory)
    for read_fn in os.listdir(directory):
        path = directory + "/" + read_fn
        if not os.path.isfile(path) or read_fn.find(".fast5") == -1:
            continue 

        for strand in range(0, 2):
            f = fast5_reader.File(path)
            if not f.have_model(strand):
                print('file [' + read_fn + '] has no model for strand [' + str(strand) + ']')
                sys.exit(0)

            m, a = f.get_model(strand)
            model_file = f.get_model_file(strand)

            # Shorten the model name by stripping /opt/chimera/model
            model_file = model_file.replace("/opt/chimaera/model/", "")

            # Replace / with _
            model_file = model_file.replace("/", "_")
            for e in m:
                kmer = e[0]
                mean = e[1]
                expected_mean = base_models[model_file][kmer][0]
                d = float(mean) - float(expected_mean)
                if abs(d) > 0.0001:
                    print("Model mismatch %s -- %s -- %s\n" % (read_fn, model_file, kmer))
            checked_models[model_file] += 1

for (k,v) in checked_models.items():
    print(k,v)
