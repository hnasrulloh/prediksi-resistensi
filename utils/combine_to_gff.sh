#!/bin/bash
# Run this script for the project directory
# Arguments
#   1: a file listing of ID genomes
#   2: a directory containing the gff3 files
#   3: a direcorty containing the fna files
#   4: a output directory

ID_LIST=$1
GFF_DIR=$2
FNA_DIR=$3
OUT_DIR=$4

mkdir -p $OUT_DIR

# Combine a coresponding fna and gff file into a single gff
for id in `cat $ID_LIST`; do
    cat $GFF_DIR/$id.PATRIC.gff > $OUT_DIR/$id.gff
    echo "##FASTA" >> $OUT_DIR/$id.gff
    cat $FNA_DIR/$id.fna >> $OUT_DIR/$id.gff
done

# Fix the incompatible syntax
sed -i -e 's/accn|//g' $OUT_DIR/*
