#!/bin/bash
#

OUTPUT_DATA_FOLDER=$1
PROJECT_INPUT_DATA_FOLDER=$2
BLAST_RESULTS_FILE=$3
DATA_SET_NAME=$4
AVG_RPKM=$5  ## this parameter is actually not used in this version of the script.
H3K4ME3_PEAK_FILE=$6

# genscan data retrieval was not run on biocluster, but the script and the data
# are in the folder below
# note we use the file with chrom names like 1,2,3,...
GENSCAN_DATA_FOLDER="$PROJECT_INPUT_DATA_FOLDER"/genscan_genes
#GENSCAN_BED_FILE="$GENSCAN_DATA_FOLDER"/genscan_mouse_genes_w_bed_cols_w_short_chrom_names.bed
#GENSCAN_BED_HDR="$GENSCAN_DATA_FOLDER"/genscan_mouse_genes_w_bed_cols_header.txt
GENSCAN_BED_FILE="$GENSCAN_DATA_FOLDER"/genscan_mouse_exons_w_short_chrom_names.bed
GENSCAN_BED_HDR="$GENSCAN_DATA_FOLDER"/genscan_mouse_exons_bed_cols_header.txt

# N-SCAN data retrieval was not run on biocluster, but the script and the data
# are in the folder below
# note we use the file with chrom names like 1,2,3,...
NSCAN_DATA_FOLDER="$PROJECT_INPUT_DATA_FOLDER"/N-SCAN
#NSCAN_BED_FILE="$NSCAN_DATA_FOLDER"/N-SCAN_mouse_genes_w_bed_cols_w_short_chrom_names.bed
#NSCAN_BED_HDR="$NSCAN_DATA_FOLDER"/N-SCAN_mouse_genes_w_bed_cols_header.txt
NSCAN_BED_FILE="$NSCAN_DATA_FOLDER"/N-SCAN_mouse_exons_w_short_chrom_names.bed
NSCAN_BED_HDR="$NSCAN_DATA_FOLDER"/N-SCAN_mouse_exons_bed_cols_header.txt

# ensemble TSS data retrieval was not run on biocluster, but the script and the data
# are in the folder below
TSS_DATA_FOLDER="$PROJECT_INPUT_DATA_FOLDER"/ensembl_TSS_data
TSS_BED_FILE="$TSS_DATA_FOLDER"/ncbi37_biomart_transcript.bed
TSS_BED_HDR="$TSS_DATA_FOLDER"/ncbi37_biomart_transcript_BED_cols_header.txt

# We will need a tmp file for some file manipulation below
TMP="$OUTPUT_DATA_FOLDER"/tmp.txt

# now that we have all the data, use bedtools to decorate the blastx results with TSS, genscan and N-SCAN data
BED_FILE="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME".bed
echo "$BED_FILE"
# remove comment lines beginning with '#' with grep
# remove the header line with tail -n +2
# now use awk to put the 12 standard bed file columns if front 
# for now we have have data for the first 6 bed columns chrom, start, stop, name, score(fpkm), strand and the rest are "."
grep '^#' -v $BLAST_RESULTS_FILE | tail -n +2 | awk -F"\t" -v OFS="\t" -v SET="$DATA_SET_NAME" '{print $1, $2, $3, $4,$5,$6, SET, $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' > "$BED_FILE"

# remove comment lines beginning with '#' with grep
# get the header row with head -n 1
# now use awk to put the 12 standard bed file columns if front of original data
# for now we have have data for the first 4 bed columns chrom, start, stop, name, the rest are "."
HDR_FILE="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr.txt
grep '^#' -v $BLAST_RESULTS_FILE | head -n 1 | awk -F"\t" -v OFS="\t" '{print "chrom","chromStart","chromEnd","name","rpkm","strand","data_source",$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' > $HDR_FILE

module load bedtools/2.21.0
# sort the blastx results and the tss data
SORTED_BED1="$OUTPUT_DATA_FOLDER"/sorted_bed1.bed
SORTED_BED2="$OUTPUT_DATA_FOLDER"/sorted_bed2.bed
bedtools sort -i $BED_FILE > $SORTED_BED1
bedtools sort -i $TSS_BED_FILE > $SORTED_BED2
# combine the blastx results and tss data into a new file
BED_w_TSS="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_w_TSS.bed
# for the bedtools closest command use -s for strand (consider strand) and -t first to use the first closest TSS in case of a tie
bedtools closest -a $SORTED_BED1 -b $SORTED_BED2 -s -t first > $BED_w_TSS
# combine the blastx header and tss header
HDR_w_TSS="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr_w_TSS.txt
paste $HDR_FILE $TSS_BED_HDR > $HDR_w_TSS

# sort the genscan data and combine it with the blastx & tss data
BED_w_TSS_GENSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_w_TSS_genscan.bed
bedtools sort -i $GENSCAN_BED_FILE > $SORTED_BED2
bedtools intersect -wao -a $BED_w_TSS -b $SORTED_BED2 -s > $BED_w_TSS_GENSCAN 2>bedtools_error.txt
# combine the genscan header with the blastx & tss header
HDR_w_TSS_GENSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr_w_TSS_genscan.txt
paste $HDR_w_TSS $GENSCAN_BED_HDR > $HDR_w_TSS_GENSCAN
awk -F"\t" -v OFS="\t" '{print $0,"overlap_[genscan]"}' $HDR_w_TSS_GENSCAN > $TMP && mv $TMP $HDR_w_TSS_GENSCAN

# sort the N-SCAN data and combine it with the blastx,tss & genscan data
BED_w_TSS_GENSCAN_NSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_w_TSS_genscan_nscan.bed
bedtools sort -i $NSCAN_BED_FILE > $SORTED_BED2
bedtools intersect -wao -a $BED_w_TSS_GENSCAN -b $SORTED_BED2 -s > $BED_w_TSS_GENSCAN_NSCAN 2>>bedtools_error.txt
# combine the N-SCAN header with the blastx & tss & genscan header
HDR_w_TSS_GENSCAN_NSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr_w_TSS_genscan_nscan.txt
paste $HDR_w_TSS_GENSCAN $NSCAN_BED_HDR > $HDR_w_TSS_GENSCAN_NSCAN
awk -F"\t" -v OFS="\t" '{print $0,"overlap_[N-SCAN]"}' $HDR_w_TSS_GENSCAN_NSCAN > $TMP && mv $TMP $HDR_w_TSS_GENSCAN_NSCAN

# Added March 14, 2017 - call script to add closest peak (h3k4me3) to file.  Also update header.
FILE_TO_UPDATE=$BED_w_TSS_GENSCAN_NSCAN
HDR_TO_UPDATE=$HDR_w_TSS_GENSCAN_NSCAN
sh closest_peak.sh $OUTPUT_DATA_FOLDER $FILE_TO_UPDATE $HDR_TO_UPDATE $H3K4ME3_PEAK_FILE $DATA_SET_NAME


################
#  Joe Troy 1/14/17 - comment out the code that adds the RPKM, as we have added already
#  in the pipeline.
# 
# # sort the RPKM data and combine it with the blastz, txx, genscan & nscan data
# BED_w_TSS_GENSCAN_NSCAN_RPKM="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_w_TSS_genscan_nscan_rpkm.bed
# bedtools sort -i $AVG_RPKM > $SORTED_BED2
# # use the option '-f 0.99 -r' for 99% overlap as the rpkm data should have exact coordinates as the BLAST data
# bedtools intersect -wao -f 0.99 -r -a $BED_w_TSS_GENSCAN_NSCAN -b $SORTED_BED2 -s > $BED_w_TSS_GENSCAN_NSCAN_RPKM 2>>bedtools_error.txt
# # combine the N-SCAN header with the blastx & tss & genscan header
# HDR_w_TSS_GENSCAN_NSCAN_RPKM="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr_w_TSS_genscan_nscan_rpkm.txt
# RPKM_HDR="$OUTPUT_DATA_FOLDER"/rpkm_header.txt
# echo -e "chrom[rpkm]\tchromStart[rpkm]\tchromEnd[rpkm]\tname[rpkm]\tvalue[rpkm]\tstrand[rpkm]\tthickStart[rpkm]\tthickEnd[rpkm]\titemRbg[rpkm]\tblockCount[rpkm]\tblockSizes[rpkm]\tblockStarts[rpkm]\tAVG_RPKM" > $RPKM_HDR
# paste $HDR_w_TSS_GENSCAN_NSCAN $RPKM_HDR > $HDR_w_TSS_GENSCAN_NSCAN_RPKM
# awk -F"\t" -v OFS="\t" '{print $0,"overlap_[RPKM]"}' $HDR_w_TSS_GENSCAN_NSCAN_RPKM > $TMP && mv $TMP $HDR_w_TSS_GENSCAN_NSCAN_RPKM


# put the final header and final data together and write it out.
FINAL_FILE="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_results.txt
cat $HDR_w_TSS_GENSCAN_NSCAN $BED_w_TSS_GENSCAN_NSCAN > $FINAL_FILE

# Finally we will add the data for each set together
ALL_DATA_SET_REPORT="$OUTPUT_DATA_FOLDER"/all_data_set_report.txt
if [ ! -f $ALL_DATA_SET_REPORT ];
then
	# if the all data set report does not exist, first add then header
	cat $HDR_w_TSS_GENSCAN_NSCAN > $ALL_DATA_SET_REPORT
fi
# now add the data to the report
cat $BED_w_TSS_GENSCAN_NSCAN >> $ALL_DATA_SET_REPORT

# change chrom names to chr1, chr2, etc.  (columns 1, 26, 40, 54, 62 and 70)
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $BED_w_TSS_GENSCAN_NSCAN 1
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $BED_w_TSS_GENSCAN_NSCAN 20
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $BED_w_TSS_GENSCAN_NSCAN 34
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $BED_w_TSS_GENSCAN_NSCAN 48
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $BED_w_TSS_GENSCAN_NSCAN 56
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $BED_w_TSS_GENSCAN_NSCAN 64

sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $FINAL_FILE 1
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $FINAL_FILE 20
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $FINAL_FILE 34
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $FINAL_FILE 48
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $FINAL_FILE 56
sh mm_ensembl_chrom_ids_to_UCSC_by_column.sh $FINAL_FILE 64

# move the important file to the data folder
mkdir -p "$OUTPUT_DATA_FOLDER"/data
mv $FINAL_FILE "$OUTPUT_DATA_FOLDER"/data/

# create a folder for UCSC custom tracks
mkdir -p "$OUTPUT_DATA_FOLDER"/UCSC_custom_tracks
EXON_BED="$OUTPUT_DATA_FOLDER"/UCSC_custom_tracks/"$DATA_SET_NAME"_custom_track.bed
echo "track name='$DATA_SET_NAME' description='$DATA_SET_NAME'" > $EXON_BED
cut -f 1-4 $BED_w_TSS_GENSCAN_NSCAN | sort -u >> $EXON_BED

echo "end of individual data script"