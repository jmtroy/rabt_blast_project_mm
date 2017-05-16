#!/bin/bash

#
# This script takes one tab separated file as input
# and take the column number where the chromosome is in
# and will convert that column in the file
# from 1,2,3... to chr1, chr2, chr3 ... for chromosome ids of mouse genomes
# 

# input file is the first and only argument.
FILE=$1
COLUMN=$2

# create a unique name for the temp file
dt=`date +%Y%m%d`
tm=`date +%H%M%S`
TMP=TMP_"$dt"_"$tm"

awk -F $'\t' -v col=$COLUMN 'BEGIN {OFS=FS} $col == "1" {$col = "chr1"} \
$col == "2" {$col = "chr2"} \
$col == "3" {$col = "chr3"} \
$col == "4" {$col = "chr4"} \
$col == "5" {$col = "chr5"} \
$col == "6" {$col = "chr6"} \
$col == "7" {$col = "chr7"} \
$col == "8" {$col = "chr8"} \
$col == "9" {$col = "chr9"} \
$col == "10" {$col = "chr10"} \
$col == "11" {$col = "chr11"} \
$col == "12" {$col = "chr12"} \
$col == "13" {$col = "chr13"} \
$col == "14" {$col = "chr14"} \
$col == "15" {$col = "chr15"} \
$col == "16" {$col = "chr16"} \
$col == "17" {$col = "chr17"} \
$col == "18" {$col = "chr18"} \
$col == "19" {$col = "chr19"} \
$col == "20" {$col = "chr20"} \
$col == "21" {$col = "chr21"} \
$col == "22" {$col = "chr22"} \
$col == "X" {$col = "chrX"} \
$col == "Y" {$col = "chrY"} \
$col == "MT" {$col = "chrM"} ;1' $FILE > $TMP && mv $TMP $FILE