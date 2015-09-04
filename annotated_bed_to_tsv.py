# /usr/bin/env python

import sys

print_header = True
num_fields = -1

for line in sys.stdin:
    fields = line.rstrip().split()
    kv_str = fields[3]
    kv_pairs = kv_str.split(';')

    keys = [x.split("=")[0] for x in kv_pairs]
    values = [x.split("=")[1] for x in kv_pairs]

    if print_header:
        header = ["chr", "start", "end"] + keys
        print "\t".join(header)
        print_header = False
        num_fields = len(header)

    output = fields[0:3] + values
    if len(output) != num_fields:
        print "Error, field length mismatch"
    print "\t".join(output)
