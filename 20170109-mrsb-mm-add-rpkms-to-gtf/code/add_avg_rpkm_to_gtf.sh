#!/bin/bash

# Joe Troy 12/11/2016 - Created
# Joe Troy 12/13/2016 - add simple report (in bed format) with just RPKM for A30C, FC30C and H30C
# Joe Troy 01/24/2017 - fixed H gtf input GTF_IN=/home/jmtroy2/jmt/projects/20170109-mrsb-mm-add-rpkms-to-gtf/project_input_data/H30C-combined-gtf/contained.combined.gtf
# Joe Troy 04/27/2017 - added to github, created new code folder, and changed to use rabt cuffcompare output

#
# Create RPKM counts from cuffcompare .gtf files and add back to .gtf file
#

# refer to the link below for analyzeRepeats.pl documentation
# template from documentation   http://homer.salk.edu/homer/ngs/analyzeRNA.html
# analyzeRepeats.pl <rna|repeats|gtf file> <genome> [options] -count [genes|exons|introns|cds|3utr|5utr] -d <Tag Directory> [Tag Directory 2] ... > <output file>
#
# Although the HOMER program is called analyzeRepeats.pl it will also quantify expression
# for a custom set

# load HOMER
module load homer

# set variable to HOMER's bin folder - needed to run HOMER programs
HOMERBIN=/home/apps/homer/homer-4.7/bin

OUTPUT_DATA_FOLDER=/home/jmtroy2/jmt/projects/20170109-mrsb-mm-add-rpkms-to-gtf/rabt_output
# remove and recreate output data folder
if [ -d "$OUTPUT_DATA_FOLDER" ]
then
  rm -R "$OUTPUT_DATA_FOLDER"
fi

# create / recreate the output folder
mkdir "$OUTPUT_DATA_FOLDER"

# define some temp files
TMP1="$OUTPUT_DATA_FOLDER"/tmp1.txt
TMP2="$OUTPUT_DATA_FOLDER"/tmp2.txt
TMP_HDR="$OUTPUT_DATA_FOLDER"/tmp_header.txt
TMP_DATA="$OUTPUT_DATA_FOLDER"/tmp_data.txt

####################################################### A30C
#
# Process Amygdala 30 minute control un-annotated exons
# define .gtf input
#
GTF_IN=/home/jmtroy2/jmt/projects/20170109-mrsb-mm-add-rpkms-to-gtf/project_input_data/A30C-combined-gtf/contained.combined.gtf_from_cufflinks_rabt/contained.combined.gtf

# defile tag directories
T1=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmA30C1_tags
T2=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmA30C2_tags
T3=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmA30C3_tags
T4=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmA30C4_tags
T5=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmA30C5_tags
OUT_PREFIX=A30C
#

RPKM_FILE="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_w_rpkm_counts.gtf
RPKM_TMP_FILE="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_w_rpkm_counts_tmp.txt
SORTED_GTF="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_sorted_gtf.gtf
# use Homer's analyzeRepeats to get the raw counts
# "$GTF_IN" is the input gtf file
# "mm9" is the genone, but for this we don't care
# "-strand -" Is the strand setting that is appropriate for TruSeq illumina processing - but would recommend always testing this
# "-count exons" count if feature in col 3 of the gtf is 'exon' (the should all be exons'
# "-noCondensingParts" do not roll up counts per each transcript, provide counts exon by exon even when a transcript has multiple exons
# "-raw" provide the RPKM (Reads Per Kilobase of transcript per Million mapped reads) value
# "  (head -n1 && tail -n +2 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 ) " Sort everything except the first line
# by chromosome, start, end and strand
#perl "$HOMERBIN"/analyzeRepeats.pl $GTF_IN mm9 -strand - -count exons -noCondensingParts -rpkm -d $T1 $T2 $T3 $T4 $T5 | (head -n1 && tail -n +2 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 ) > $TMP
perl "$HOMERBIN"/analyzeRepeats.pl $GTF_IN mm9 -strand - -count exons -noCondensingParts -rpkm -d $T1 $T2 $T3 $T4 $T5 > $TMP1

# remove header 
# then sort by chromosome, start, end, and strand
# then only print the average of the rpkm columns 
tail -n +2 $TMP1 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 > $RPKM_TMP_FILE
cat $RPKM_TMP_FILE | awk -F"\t" -v OFS="\t" -v SUM=0 -v AVG=0 '{SUM=$9+$10+$11+$12+$13;AVG=SUM/5;print AVG}' > $TMP2

# now sort the GTF_IN files by by chromosome, start, end, and strand
sort -t $'\t' -k 1,1 -k 4,4n -k 5,5n -k 7,7 $GTF_IN > $SORTED_GTF

# now we PASTE the tmp2 and the sorted_gtf file together to get a new .gtf file
# this will be the file used as input to blast.
paste $SORTED_GTF $TMP2 > $RPKM_FILE

####################################################### FC30C
#
# Process Frontal Cortex 30 minute control un-annotated exons
# define .gtf input
#
GTF_IN=/home/jmtroy2/jmt/projects/20170109-mrsb-mm-add-rpkms-to-gtf/project_input_data/FC30C-combined-gtf/contained.combined.gtf

# defile tag directories
T1=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmFC30C1_tags
T2=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmFC30C2_tags
T3=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmFC30C3_tags
T4=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmFC30C4_tags
T5=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmFC30C5_tags
OUT_PREFIX=FC30C
#

RPKM_FILE="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_w_rpkm_counts.gtf
RPKM_TMP_FILE="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_w_rpkm_counts_tmp.txt
SORTED_GTF="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_sorted_gtf.gtf
# use Homer's analyzeRepeats to get the raw counts
# "$GTF_IN" is the input gtf file
# "mm9" is the genone, but for this we don't care
# "-strand -" Is the strand setting that is appropriate for TruSeq illumina processing - but would recommend always testing this
# "-count exons" count if feature in col 3 of the gtf is 'exon' (the should all be exons'
# "-noCondensingParts" do not roll up counts per each transcript, provide counts exon by exon even when a transcript has multiple exons
# "-raw" provide the RPKM (Reads Per Kilobase of transcript per Million mapped reads) value
# "  (head -n1 && tail -n +2 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 ) " Sort everything except the first line
# by chromosome, start, end and strand
#perl "$HOMERBIN"/analyzeRepeats.pl $GTF_IN mm9 -strand - -count exons -noCondensingParts -rpkm -d $T1 $T2 $T3 $T4 $T5 | (head -n1 && tail -n +2 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 ) > $TMP
perl "$HOMERBIN"/analyzeRepeats.pl $GTF_IN mm9 -strand - -count exons -noCondensingParts -rpkm -d $T1 $T2 $T3 $T4 $T5 > $TMP1

# remove header 
# then sort by chromosome, start, end, and strand
# then only print the average of the rpkm columns 
tail -n +2 $TMP1 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 > $RPKM_TMP_FILE
cat $RPKM_TMP_FILE | awk -F"\t" -v OFS="\t" -v SUM=0 -v AVG=0 '{SUM=$9+$10+$11+$12+$13;AVG=SUM/5;print AVG}' > $TMP2

# now sort the GTF_IN files by by chromosome, start, end, and strand
sort -t $'\t' -k 1,1 -k 4,4n -k 5,5n -k 7,7 $GTF_IN > $SORTED_GTF

# now we PASTE the tmp2 and the sorted_gtf file together to get a new .gtf file
# this will be the file used as input to blast.
paste $SORTED_GTF $TMP2 > $RPKM_FILE

####################################################### H30C
#
# Process Hypothalamus 30 minute control un-annotated exons
# define .gtf input
#
GTF_IN=/home/jmtroy2/jmt/projects/20170109-mrsb-mm-add-rpkms-to-gtf/project_input_data/H30C-combined-gtf/contained.combined.gtf

# defile tag directories
T2=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmH30C2_tags
T3=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmH30C3_tags
T5=/home/jmtroy2/jmt/projects/20161207-temp-mrsb-homer-tag-dirs/mmH30C5_tags
OUT_PREFIX=H30C
#

RPKM_FILE="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_w_rpkm_counts.gtf
RPKM_TMP_FILE="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_w_rpkm_counts_tmp.txt
SORTED_GTF="$OUTPUT_DATA_FOLDER"/"$OUT_PREFIX"_sorted_gtf.gtf
# use Homer's analyzeRepeats to get the raw counts
# "$GTF_IN" is the input gtf file
# "mm9" is the genone, but for this we don't care
# "-strand -" Is the strand setting that is appropriate for TruSeq illumina processing - but would recommend always testing this
# "-count exons" count if feature in col 3 of the gtf is 'exon' (the should all be exons'
# "-noCondensingParts" do not roll up counts per each transcript, provide counts exon by exon even when a transcript has multiple exons
# "-raw" provide the RPKM (Reads Per Kilobase of transcript per Million mapped reads) value
# "  (head -n1 && tail -n +2 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 ) " Sort everything except the first line
# by chromosome, start, end and strand
#perl "$HOMERBIN"/analyzeRepeats.pl $GTF_IN mm9 -strand - -count exons -noCondensingParts -rpkm -d $T1 $T2 $T3 $T4 $T5 | (head -n1 && tail -n +2 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 ) > $TMP
perl "$HOMERBIN"/analyzeRepeats.pl $GTF_IN mm9 -strand - -count exons -noCondensingParts -rpkm -d $T2 $T3 $T5 > $TMP1

# remove header 
# then sort by chromosome, start, end, and strand
# then only print the average of the rpkm columns 
tail -n +2 $TMP1 | sort -t $'\t' -k 2,2 -k 3,3n -k 4,4n -k 5,5 > $RPKM_TMP_FILE
cat $RPKM_TMP_FILE | awk -F"\t" -v OFS="\t" -v SUM=0 -v AVG=0 '{SUM=$9+$10+$11;AVG=SUM/3;print AVG}' > $TMP2

# now sort the GTF_IN files by by chromosome, start, end, and strand
sort -t $'\t' -k 1,1 -k 4,4n -k 5,5n -k 7,7 $GTF_IN > $SORTED_GTF

# now we PASTE the tmp2 and the sorted_gtf file together to get a new .gtf file
# this will be the file used as input to blast.
paste $SORTED_GTF $TMP2 > $RPKM_FILE
