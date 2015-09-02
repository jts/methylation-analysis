#! /bin/bash

cat $1 | grep SITE | awk '{ print "chr" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 }'
