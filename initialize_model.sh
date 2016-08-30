#! /bin/bash

IN=$1
STRAND=$2
NAME=$3
KIT=$4
# Add the default model type ("ONT") to the initial model
cat <(echo -e "#model_name\t$NAME\n#type\tONT\n#kit\t$KIT\n#strand\t$STRAND") $IN
