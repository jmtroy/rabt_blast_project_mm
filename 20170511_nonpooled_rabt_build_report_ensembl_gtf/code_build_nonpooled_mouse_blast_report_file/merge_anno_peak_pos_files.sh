#!/bin/bash
# script to merge homer annotation and peaks files
ANNO_INPUT_FILE=$1
PEAK_INPUT_FILE=$2
PEAK_OUTPUT_BED_FILE=$3
OUTPUT_DATA_FOLDER=$4

module load bedtools/2.25.0

#
# get data from anno file
#
# remove first column
#
# then use awk to build a bed file
# col 2 > 1		chrom
# col 3 > 2		start
# col 4 > 3		end
# col 1 > 4		name
# "." > 5		value
# col 5 > 6		strand
# col 8 > 7		annotation (from Homer)
# col 24 to 8	mgi symbol from homer
TMP_ANNO_FILE=$OUTPUT_DATA_FOLDER/tmp_homer_anno.bed
tail -n +2 $ANNO_INPUT_FILE | awk -F"\t" -v OFS="\t" '{print $2, $3, $4, $1, ".", $5, $8, $16}' > $TMP_ANNO_FILE

#
# get date from peaks file
#
# remove columns beginning with '#'
#
# then use awk to build a bed file
# col 2 > 1		chrom
# col 3 > 2		start
# col 4 > 3		end
# col 1 > 4		name
# "." > 5		value
# col 5 > 6		strand
# col 8 > 7 	findPeaks score (from Homer)
# col 11 > 8	Fold Change vs Control (from Homer)
# col 12 > 9 	p-value vs Control (from Homer)
TMP_PEAK_FILE=$OUTPUT_DATA_FOLDER/tmp_homer_peak.bed
grep -v "^#" $PEAK_INPUT_FILE | awk -F"\t" -v OFS="\t" '{print $2, $3, $4, $1, ".", $5, $8, $11, $12}' > $TMP_PEAK_FILE

# the two files should have the same genomic coordinates, so require 100% overlap for the intersect
TMP_COMBINDED_FILE=$OUTPUT_DATA_FOLDER/tmp_homer_combinded.bed
bedtools intersect -a $TMP_ANNO_FILE -b $TMP_PEAK_FILE -wo -f 1.0 -F 1.0 > $TMP_COMBINDED_FILE

# now keep only the needed data
# 1 > 1 chrom
# 2 > 2 start
# 3 > 3 end
# 4 > 4 homer id/name
# 5 > 5 empty
# 6 > 6 strand / meaningless
# 7 > 7 annotation from homer
# 8 > 8 mgi symbol from homer
# 15 > 9 findPeaks score from homer
# 16 > 10 Fold Change vs Control (from Homer)
# 17 > 11 p-value vs Control (from Homer)
awk -F"\t" -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $15, $16, $17}' $TMP_COMBINDED_FILE > $PEAK_OUTPUT_BED_FILE

# now convert the chr1, chr2 to 1, 2 etc ...
TMP1="$OUTPUT_DATA_FOLDER"/tmp1_merge_peaks.bed
awk -F $'\t'  'BEGIN {OFS=FS} $1 == "chr1" {$1 = "1"} \
$1 == "chr2" {$1 = "2"} \
$1 == "chr3" {$1 = "3"} \
$1 == "chr4" {$1 = "4"} \
$1 == "chr5" {$1 = "5"} \
$1 == "chr6" {$1 = "6"} \
$1 == "chr7" {$1 = "7"} \
$1 == "chr8" {$1 = "8"} \
$1 == "chr9" {$1 = "9"} \
$1 == "chr10" {$1 = "10"} \
$1 == "chr11" {$1 = "11"} \
$1 == "chr12" {$1 = "12"} \
$1 == "chr13" {$1 = "13"} \
$1 == "chr14" {$1 = "14"} \
$1 == "chr15" {$1 = "15"} \
$1 == "chr16" {$1 = "16"} \
$1 == "chr17" {$1 = "17"} \
$1 == "chr18" {$1 = "18"} \
$1 == "chr19" {$1 = "19"} \
$1 == "chr20" {$1 = "20"} \
$1 == "chr21" {$1 = "21"} \
$1 == "chr22" {$1 = "22"} \
$1 == "chrX" {$1 = "X"} \
$1 == "chrY" {$1 = "Y"} ;1' $PEAK_OUTPUT_BED_FILE > $TMP1
# now move the temp file to the output file
mv $TMP1 $PEAK_OUTPUT_BED_FILE


# clean up temp files
rm $TMP_ANNO_FILE
rm $TMP_PEAK_FILE
rm $TMP_COMBINDED_FILE
