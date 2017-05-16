#!/bin/bash
#PBS -A simons
#PBS -l walltime=168:00:00
#PBS -l nodes=1:ppn=8
#PBS -N report
#PBS -o report.out
#PBS -e report.err
#PBS -m abe
#PBS -M jmtroy2@igb.illinois.edu

# IMPORTANT: use the below qsub command to run from the projects code folder
# qsub -v my_script_name=main_pooled_cufflinks_rabt.sh -S /bin/bash main_pooled_cufflinks_rabt.sh

# pool mouse bam files together and then run cufflinks

# change directory to torque working directory (the directory the script was run from)
cd $PBS_O_WORKDIR

# echo out script name and time
echo begin script "$my_script_name" at `date`

# Get project name, the project name is assumed to be same as the parent folder
# of the current folder.
# the current folder is assumed to be the code folder.
CURRENT_FOLDER=`pwd`
PARENT_FOLDER="$(dirname "$CURRENT_FOLDER")"
PROJECT_NAME="$(basename "$PARENT_FOLDER")"

# TO DO, set input_data_folder if needed
INPUT_DATA_FOLDER="/home/groups/simons/Joe/mm_blast_project/input_data"

# the PROJECT_FOLDER is in the simoms project foldes
PROJECT_FOLDER="$PARENT_FOLDER"
CODE_FOLDER="$CURRENT_FOLDER"
PROJECT_INPUT_DATA_FOLDER="$PROJECT_FOLDER"/project_input_data

## run id variable - output and intermediate files will go
## in the run id directory (for example /RUN_20130708_162650)
dt=`date +%Y%m%d`
tm=`date +%H%M%S`
RUN_ID=RUN_"$dt"_"$tm"

# TO DO, to help identify the contents of the outfolder, add something to the RUN_ID
RUN_ID="$RUN_ID"

## set a variable with the name of the directory the output (and interim data) will be placed, and then create the folders

# TO DO set the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# INTERIM_DATA_FOLDER="$PROJECT_FOLDER"/"$RUN_ID"/interim_data
# TO DO create the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# mkdir -p "$INTERIM_DATA_FOLDER" 

OUTPUT_DATA_FOLDER="$PROJECT_FOLDER"/output_pooled_cufflinks_rabt"_"$RUN_ID
SAVED_CODE="$OUTPUT_DATA_FOLDER"/saved_code
## (use -p option below) mkdir "$PROJECT_FOLDER"/"$RUN_ID"/output_data
mkdir -p "$OUTPUT_DATA_FOLDER" # the -p option will create any leading directories that do not already exist.
INTERIM_DATA_FOLDER="$OUTPUT_DATA_FOLDER"/interim_data
# not needed for this script # mkdir -p "$INTERIM_DATA_FOLDER"

# create a run log file in the output folder
# This file can be used to capture any log information you need in your script.
RUN_LOG_FILE="$OUTPUT_DATA_FOLDER"/"$RUN_ID"_LOG.txt
echo begin script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 
echo "The qsub job name: $PBS_JOBNAME" at `date` >> "$RUN_LOG_FILE"
echo "The qsub job id: $PBS_JOBID" at `date` >> "$RUN_LOG_FILE"
echo "The project folder is: $PROJECT_FOLDER" >> "$RUN_LOG_FILE"
echo "The code folder is: $CODE_FOLDER" >> "$RUN_LOG_FILE"
echo "The output folder is: $OUTPUT_DATA_FOLDER" >> "$RUN_LOG_FILE"

#################################################################################################
### BEGIN THE REAL WORK NOW THAT THE INITIAL HOUSE KEEPING IS DONE ##############################
#################################################################################################

A30C1=/home/groups/simons/Joe/mouse_A_30/A_30_CK1_ATTACTCG-TAATCTTA_L00M_R1_001.tophat/accepted_hits.bam
A30C2=/home/groups/simons/Joe/mouse_A_30/A_30_CK2_ATTACTCG-CAGGACGT_L00M_R1_001.tophat/accepted_hits.bam
A30C3=/home/groups/simons/Joe/mouse_A_30/A_30_CK3_ATTACTCG-GTACTGAC_L00M_R1_001.tophat/accepted_hits.bam
A30C4=/home/groups/simons/Joe/mouse_A_30/A_30_CK4_TCCGGAGA-TATAGCCT_L00M_R1_001.tophat/accepted_hits.bam
A30C5=/home/groups/simons/Joe/mouse_A_30/A_30_CK5_TCCGGAGA-ATAGAGGC_L00M_R1_001.tophat/accepted_hits.bam

# use samtools to merge the file above into one big freaking BAM file.
module load samtools/1.3.1
module load cufflinks/2.1.1

# samtools merge
# -r Attach an RG tag to each alignment. The tag value is inferred from file names.
# -b List of input BAM files, one file per line.


samtools merge -r $OUTPUT_DATA_FOLDER/A30C_merged_accepted_hits.bam \
$A30C1 \
$A30C2 \
$A30C3 \
$A30C4 \
$A30C5

# OUT_FLD=/home/groups/simons/Joe/mm_blast_project/20170411_cufflinks_pooled/output_950_build_mouse_blast_report_file_RUN_20170411_114340
# use OUT_FLD instead of $OUTPUT_DATA_FOLDER because we have already merged the file

MERGED_BAM=$OUTPUT_DATA_FOLDER/A30C_merged_accepted_hits.bam
CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/A30C_cufflinks_pooled_rabt
gtf=/home/groups/simons/Joe/ensembl37_NCBIM37_genome/genes.gtf
cufflinks -p 8 -o $CUFFLINKS_OUT --GTF-guide $gtf $MERGED_BAM

# do cuffcompare (see notes below)
GTF1="$CUFFLINKS_OUT"/transcripts.gtf
CUFFCOMPARE_OUT="$OUTPUT_DATA_FOLDER"/A30_cuffcompare_contained
# point to folder with fastq for each chromosome
CHR_FA=/home/mirrors/igenome/Mus_musculus/NCBI/build37.2/Sequence/Chromosomes
cuffcompare -V -o "$CUFFCOMPARE_OUT" -C -s "$CHR_FA" -r "$gtf" -R "$GTF1" 2> "$OUTPUT_DATA_FOLDER"/A30c_contained.log.txt

#
# FC30C
#

FC30C1=/home/groups/simons/Joe/mouse_FC_30/FC_30_CK1_ATTCAGAA-GGCTCTGA_L00M_R1_001.tophat/accepted_hits.bam
FC30C2=/home/groups/simons/Joe/mouse_FC_30/FC_30_CK2_ATTCAGAA-AGGCGAAG_L00M_R1_001.tophat/accepted_hits.bam
FC30C3=/home/groups/simons/Joe/mouse_FC_30/FC_30_CK3_ATTCAGAA-TAATCTTA_L00M_R1_001.tophat/accepted_hits.bam
FC30C4=/home/groups/simons/Joe/mouse_FC_30/FC_30_CK4_ATTCAGAA-CAGGACGT_L00M_R1_001.tophat/accepted_hits.bam
FC30C5=/home/groups/simons/Joe/mouse_FC_30/FC_30_CK5_ATTCAGAA-GTACTGAC_L00M_R1_001.tophat/accepted_hits.bam

# samtools merge
# -r Attach an RG tag to each alignment. The tag value is inferred from file names.
# -b List of input BAM files, one file per line.

MERGED_BAM=$OUTPUT_DATA_FOLDER/FC30C_merged_accepted_hits.bam

samtools merge -r $MERGED_BAM \
$FC30C1 \
$FC30C2 \
$FC30C3 \
$FC30C4 \
$FC30C5

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/FC30C_cufflinks_pooled_rabt
gtf=/home/groups/simons/Joe/ensembl37_NCBIM37_genome/genes.gtf
cufflinks -p 8 -o $CUFFLINKS_OUT --GTF-guide $gtf $MERGED_BAM

# do cuffcompare (see notes below)
GTF1="$CUFFLINKS_OUT"/transcripts.gtf
CUFFCOMPARE_OUT="$OUTPUT_DATA_FOLDER"/FC30_cuffcompare_contained
# point to folder with fastq for each chromosome
CHR_FA=/home/mirrors/igenome/Mus_musculus/NCBI/build37.2/Sequence/Chromosomes
cuffcompare -V -o "$CUFFCOMPARE_OUT" -C -s "$CHR_FA" -r "$gtf" -R "$GTF1" 2> "$CUFFCOMPARE_OUT"/FC30c_contained.log.txt


#
# H30C
#

H30C2=/home/groups/simons/Joe/mouse_H_30/H_30_CK2_CGGCTATG-GGCTCTGA_L00M_R1_001.tophat/accepted_hits.bam
H30C3=/home/groups/simons/Joe/mouse_H_30/H_30_CK3_CGGCTATG-AGGCGAAG_L00M_R1_001.tophat/accepted_hits.bam
H30C5=/home/groups/simons/Joe/mouse_H_30/H_30_CK5_CGGCTATG-TAATCTTA_L00M_R1_001.tophat/accepted_hits.bam

# samtools merge
# -r Attach an RG tag to each alignment. The tag value is inferred from file names.
# -b List of input BAM files, one file per line.

MERGED_BAM=$OUTPUT_DATA_FOLDER/H30C_merged_accepted_hits.bam

samtools merge -r $MERGED_BAM \
$H30C2 \
$H30C3 \
$H30C5

$CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/H30C_cufflinks_pooled_rabt

gtf=/home/groups/simons/Joe/ensembl37_NCBIM37_genome/genes.gtf
cufflinks -p 8 -o $CUFFLINKS_OUT --GTF-guide $gtf $MERGED_BAM

# do cuffcompare (see notes below)
GTF1="$CUFFLINKS_OUT"/transcripts.gtf
CUFFCOMPARE_OUT="$OUTPUT_DATA_FOLDER"/H30_cuffcompare_contained
# point to folder with fastq for each chromosome
CHR_FA=/home/mirrors/igenome/Mus_musculus/NCBI/build37.2/Sequence/Chromosomes
cuffcompare -V -o "$CUFFCOMPARE_OUT" -C -s "$CHR_FA" -r "$gtf" -R "$GTF1" 2> "$CUFFCOMPARE_OUT"/contained.log.txt

# NOTES on cuffcompare
# from cuffcompare document at http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html
#
# -o <outprefix>
# 
# All output files created by Cuffcompare will have this prefix
# (e.g. .loci, .tracking, etc.). If this option is not provided
# the default output prefix being used is: "cuffcmp"
# 
# -r
# 
# An optional “reference” annotation GFF file. Each sample is matched
# against this file, and sample isoforms are tagged as overlapping, matching, or novel
# where appropriate. See the refmap and tmap output file descriptions below.
# 
# -R
# 
# If -r was specified, this option causes cuffcompare to ignore reference
# transcripts that are not overlapped by any transcript in one of cuff1.gtf,…,cuffN.gtf.
# Useful for ignoring annotated transcripts that are not present in your RNA-Seq samples
# and thus adjusting the “sensitivity” calculation in the accuracy report
# written in the file
# 
# -s <seq_dir>
# 
# Causes cuffcompare to look into for fasta files with the underlying genomic
# sequences (one file per contig) against which your reads were aligned for
# some optional classification functions. For example, Cufflinks transcripts consisting
# mostly of lower-case bases are classified as repeats. Note that must contain
# one fasta file per reference chromosome, and each file must be named after
# the chromosome, and have a .fa or .fasta extension.
# 
# -C
# 
# Enables the “contained” transcripts to be also written in the .combined.gtffile,
# with the attribute "contained_in" showing the first container transfrag found. By
# default, without this option, cuffcompare does not write in that file isoforms that
# were found to be fully contained/covered (with the same compatible intron structure)
# by other transfrags in the same locus.
#
# -V
# 
# Cuffcompare is a little more verbose about what it’s doing, printing messages 
# to stderr, and it will also show warning messages about any inconsistencies 
# or potential issues found while reading the given GFF file(s).



#####################################################################################
### END OF THE REAL WORK - DO the FINAL HOUSE KEEPING  ##############################
#####################################################################################

# copy the contents of the current folder (with this script and other code) to the saved code folder
cp -R "$CODE_FOLDER" "$SAVED_CODE"
echo end script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 
