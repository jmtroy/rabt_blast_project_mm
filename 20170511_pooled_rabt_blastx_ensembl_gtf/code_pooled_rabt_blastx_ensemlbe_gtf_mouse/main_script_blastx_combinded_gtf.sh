#!/bin/bash
#PBS -A simons
#PBS -l walltime=168:00:00
#PBS -l nodes=1:ppn=1
#PBS -N blast
#PBS -o blast.out
#PBS -e blast.err
#PBS -m abe
#PBS -M jmtroy2@igb.illinois.edu

# IMPORTANT: use the below qsub command to run from the code directory
# qsub -v my_script_name=main_script_blastx_combinded_gtf.sh -S /bin/bash main_script_blastx_combinded_gtf.sh

#
#  This script...
#
#	Does a grep on the combined gtf file from cufflinks to find intergenic transcripts
#   Note that the combinded gtf file has been updated with exon RPKMs using homer software
#		either grep on 'class_code "u"'
#		or grep on 'class_code "i"'
#
#	use bedops to convert the intergenic .gtf file to an intergenic .bed file
#
#	use repeats data to filter out exons with repeats
#
#	use pseudo gene data to filter out exons that overlap with pseudo genes
#
#	filter out exons with RKPM less than 1
#
#	blastx the remaining exons against human proteins
#
#

# change directory to torque working directory (the directory "qsub -S /bin/bash main_script.sh" was run from)
cd $PBS_O_WORKDIR

# set a variable with this scripts name
echo begin script "$my_script_name" at `date`

# Get project name, the project name is assumed to be same as the parent folder
# of the current folder.
# the current folder is assumed to be the code folder.
CURRENT_FOLDER=`pwd`
PARENT_FOLDER="$(dirname "$CURRENT_FOLDER")"
PROJECT_NAME="$(basename "$PARENT_FOLDER")"

# TO DO, set input_data_folder if needed
# input data folder not needed for this script 
# INPUT_DATA_FOLDER=

# the PROJECT_FOLDER is in the simoms project foldes
PROJECT_FOLDER="$PARENT_FOLDER"
CODE_FOLDER="$CURRENT_FOLDER"
PROJECT_INPUT_DATA_FOLDER="$PROJECT_FOLDER"/project_input_data

## run id variable - run id variable -  has the timestamp
dt=`date +%Y%m%d`
tm=`date +%H%M%S`
RUN_ID=RUN_"$dt"_"$tm"

## set a variable with the name of the directory the output (and interim data) will be placed, and then create the folders
# TO DO set the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# INTERIM_DATA_FOLDER="$PROJECT_FOLDER"/"$RUN_ID"/interim_data
# TO DO create the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# mkdir -p "$INTERIM_DATA_FOLDER" 
OUTPUT_DATA_FOLDER="$PROJECT_FOLDER"/output_pooled_rabt_blastx_ensemlbe_gtf_mouse"_"$RUN_ID
SAVED_CODE="$OUTPUT_DATA_FOLDER"/saved_code
## (use -p option below) mkdir "$PROJECT_FOLDER"/"$RUN_ID"/output_data
mkdir -p "$OUTPUT_DATA_FOLDER" # the -p option will create any leading directories that do not already exist.
# no interim data folder for this project # mkdir -p "$INTERIM_DATA_FOLDER"

# create a run log file in the output folder
# This file can be used to capture any log information you need in your script.
RUN_LOG_FILE="$OUTPUT_DATA_FOLDER"/"$RUN_ID"_LOG.txt
echo begin script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 
echo "The qsub job name: $PBS_JOBNAME" at `date` >> "$RUN_LOG_FILE"
echo "The qsub job id: $PBS_JOBID" at `date` >> "$RUN_LOG_FILE"
echo "The project folder is: $PROJECT_FOLDER" >> "$RUN_LOG_FILE"
echo "The code folder is: $CODE_FOLDER" >> "$RUN_LOG_FILE"

#################################################################################################
### BEGIN THE REAL WORK NOW THAT THE INITIAL HOUSE KEEPING IS DONE ##############################
#################################################################################################


#
#  HEY, this is a bash function, please also see the code below the function
#
function DOBLAST {

	# inputs to this function
	# COMB_GTF, path to gtf from cuffcompare with RPKM added by homer
	# SUB_OUTPUT_FLDR, path to output folder
	# FILE_PREFIX, file prefix for output files
	# GREP_VAR, either 'class_code "u"' or 'class_code "i"'  (either to find intergenic or intragenic exons)

	mkdir $SUB_OUTPUT_FLDR
	
	module load bedtools/2.25.0
	module load R/3.2.3
	
	# filter the combined gtf file from cuffcompare, using grep, to get all the rows with either 'class_code "u"' or 'class_code "i"' 
	EXONS_GTF="$SUB_OUTPUT_FLDR"/"$FILE_PREFIX"_transcripts.gtf
	grep "$GREP_VAR" $COMB_GTF > $EXONS_GTF
	
	

	#	use bedops to convert the intergenic .gtf file to an intergenic .bed file
	module load bedops/2.4.2
	EXONS_BED="$SUB_OUTPUT_FLDR"/"$FILE_PREFIX"_transcripts.bed
	INTERGENIC_BED="$SUB_OUTPUT_FLDR"/intergenic_transcripts.bed
	# the command below does:
	# uses the gtf2bed command to convert the gtf file to a bed file
	# uses awk to replace the contents of column 10 with just a "."
	# uses sort to sort the results
	# uses uniq to get only the unique rows
	# uses sortBed to sort it as a properly sorted bed file
	gtf2bed < $EXONS_GTF | awk -F $'\t'  'BEGIN {OFS=FS} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,".",$11}' | sort | uniq | sortBed > $EXONS_BED

	# We have found that homer can calculate a slightly different rpkm value for the same
	# exon (that is, chrom, start, end and strand  are the same) for different transcripts
	# we consider those as "duplicates" and call the R code below to eliminate the duplicates
	# and chose the first of the duplicate rows as the one to keep.  We do this without 
	# consideration of the rpkm value because we have seen the values are only slightly different.
	EXONS_BED_NO_DUPS="$SUB_OUTPUT_FLDR"/intergenic_transcripts_no_dups.bed
	Rscript keep_one_exon_with_multi_rpkm_values.R input_file="$EXONS_BED" output_file="$EXONS_BED_NO_DUPS"
	rm "$EXONS_BED"
	mv "$EXONS_BED_NO_DUPS" "$EXONS_BED"

	#######  Start of Filter by repeats ####################################################
	repeats_short_chrom_names="/home/groups/simons/Joe/mm_blast_project/input_data/repeats/get_mm9_repeat_masker-20170103_short_chrom_names.bed"

	# Do a bedtools intersect and keep the intersections for additional analysis
	bedtools intersect -a $EXONS_BED -b $repeats_short_chrom_names -wa  > "$SUB_OUTPUT_FLDR"/all_cuffcompare_intergenic_repeat_masker_intersects.bed
	echo "Cuffdiff Intergenic exons that overlap with a repeat element:" >> "$RUN_LOG_FILE"
	sort -u "$SUB_OUTPUT_FLDR"/all_cuffcompare_intergenic_repeat_masker_intersects.bed | wc -l >> "$RUN_LOG_FILE"

	# now eliminate those cuffcompare intergenic loci with overlap with repeat masker
	EXONS_BED_REPEAT_FILTER="$SUB_OUTPUT_FLDR"/intergenic_transcripts_filtered_for_repeats.bed
	bedtools intersect -v -a $EXONS_BED -b $repeats_short_chrom_names  > $EXONS_BED_REPEAT_FILTER
	####### Done with filtering for repeats  #####


	####### Start filtering for "Yale Pseudo60" Pseudo Genes #####
	pseudo_genes_short_chrom_names="/home/groups/simons/Joe/mm_blast_project/input_data/pseudo_genes/Yale_pseudo60-20170103_short_chrom_names.bed"

	# Do a bedtools intersect and keep the intersections for additional analysis
	bedtools intersect -a $EXONS_BED_REPEAT_FILTER -b $pseudo_genes_short_chrom_names -wa  > "$SUB_OUTPUT_FLDR"/pseudo_genes_intersect.bed
	echo "Cuffdiff Intergenic exons that overlap with a pseudo gene after filtering by repeat elements:" >> "$RUN_LOG_FILE"
	sort -u "$SUB_OUTPUT_FLDR"/pseudo_genes_intersect.bed | wc -l >> "$RUN_LOG_FILE"

	# now eliminate those cuffcompare intergenic loci with overlap with a pseudo gene
	EXONS_BED_REPEAT_AND_PSEUDO_GENE_FILTER="$SUB_OUTPUT_FLDR"/intergenic_transcripts_filtered_for_repeats_and_pseudo_genes.bed
	bedtools intersect -v -a $EXONS_BED_REPEAT_FILTER -b $pseudo_genes_short_chrom_names  > $EXONS_BED_REPEAT_AND_PSEUDO_GENE_FILTER
	####### Done with filtering for pseudo genes  #####


	####### Now filter by the RPKM values supplied by HOMER
	EXONS_BED_REPEAT_PSEUDO_GENE_AND_RPKM_FILTER="$SUB_OUTPUT_FLDR"/intergenic_transcripts_filtered_for_repeats_pseudo_genes_and_rpkm.bed
	awk -F"\t" '($11 > 1.0) {print}' $EXONS_BED_REPEAT_AND_PSEUDO_GENE_FILTER > $EXONS_BED_REPEAT_PSEUDO_GENE_AND_RPKM_FILTER

	# in column 4, add all the data the will be preserved in the fasta file
	TMP="$SUB_OUTPUT_FLDR"/tmp
	awk -F $'\t'  'BEGIN {OFS=FS} {$4 = $1":"$2"-"$3":"$4":"$11":"$6};1' $EXONS_BED_REPEAT_PSEUDO_GENE_AND_RPKM_FILTER >$TMP && mv $TMP $EXONS_BED_REPEAT_PSEUDO_GENE_AND_RPKM_FILTER

	# use bedtool getfasta to create a search fasta file of the intergenic regions seleted
	GENOME_FA=/home/groups/simons/Joe/ensembl37_NCBIM37_genome/genome.fa
	SEARCH_FA="$SUB_OUTPUT_FLDR"/blast_serch.fa
	bedtools getfasta -name -fi $GENOME_FA -bed $EXONS_BED_REPEAT_PSEUDO_GENE_AND_RPKM_FILTER -fo $SEARCH_FA

	# for this test we are getting the protein files using the below commands ...
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.2.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.3.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.4.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.5.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.6.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.7.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.8.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.9.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.10.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.11.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.12.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.13.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.14.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.15.protein.faa.gz
	# wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.16.protein.faa.gz
	# gunzip *.faa
	# cat human.1.protein.faa > human.protein.cat.file.faa
	# cat human.2.protein.faa >> human.protein.cat.file.faa
	# cat human.3.protein.faa >> human.protein.cat.file.faa
	# cat human.4.protein.faa >> human.protein.cat.file.faa
	# cat human.5.protein.faa >> human.protein.cat.file.faa
	# cat human.6.protein.faa >> human.protein.cat.file.faa
	# cat human.7.protein.faa >> human.protein.cat.file.faa
	# cat human.8.protein.faa >> human.protein.cat.file.faa
	# cat human.9.protein.faa >> human.protein.cat.file.faa
	# cat human.10.protein.faa >> human.protein.cat.file.faa
	# cat human.11.protein.faa >> human.protein.cat.file.faa
	# cat human.12.protein.faa >> human.protein.cat.file.faa
	# cat human.13.protein.faa >> human.protein.cat.file.faa
	# cat human.14.protein.faa >> human.protein.cat.file.faa
	# cat human.15.protein.faa >> human.protein.cat.file.faa
	# cat human.16.protein.faa >> human.protein.cat.file.faa
	# 

	# make a blast DB of the protein sequences in the protein .faa file
	module load blast+/2.3.0
	FASTA_IN=/home/groups/simons/Joe/mm_blast_project/input_data/ftp_human_protein_faa_files_from_ncbi/human.protein.cat.file.faa
	DB_OUT=/home/groups/simons/Joe/mm_blast_project/input_data/ftp_human_protein_faa_files_from_ncbi/human.protein.cat.file.db
	# ran this once already - no need to run again # makeblastdb -in $FASTA_IN -input_type fasta -dbtype prot -out $DB_OUT -taxid 9606 -parse_seqids

	# now use blast to get some results
	BLASTX_OUT="$SUB_OUTPUT_FLDR"/blast_results_fmt1.txt
	# no need to run this # blastx -db $DB_OUT -query $SEARCH_FA -out $BLASTX_OUT -outfmt 1

	# now use blast to get some results
	# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
	BLASTX_OUT="$SUB_OUTPUT_FLDR"/blast_results_fmt6_ftp_db.txt
	blastx -db $DB_OUT -query $SEARCH_FA -out $BLASTX_OUT -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid evalue bitscore pident length mismatch gapopen qstart qend sstart send salltitles"

	# split qseqid (which is the coordinates all as chrom:start-end) into 3 separate fields
	# first split on the first occurance of ':' (write it to tmp, and them move tmp back to the original file name)
	TMP="$SUB_OUTPUT_FLDR"/tmp
	awk '{ sub(/:/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
	# second split on the first occurance of "-" (write it to tmp, and them move tmp back to the original file name)
	awk '{ sub(/-/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
	# third splint on next occurance of :
	awk '{ sub(/:/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
	# fourth splint on next occurance of :
	awk '{ sub(/:/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
	# fifth splint on next occurance of :
	awk '{ sub(/:/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT

	# now add header

	echo "# script name = $my_script_name" > $TMP
	echo "# project folder = $PROJECT_FOLDER" >> $TMP
	echo "# code folder = $CODE_FOLDER" >> $TMP
	echo "# run id = $RUN_ID" >> $TMP
	echo "# output folder = $SUB_OUTPUT_FLDR" >> $TMP
	VALUE=`wc -l $EXONS_BED`
	echo "# Intergenic Transfrags (exons) per cuffcompare $VALUE" >> TMP
	VALUE=`wc -l $EXONS_BED_REPEAT_FILTER`
	echo "# Transfrags with NO overlap with repeat coordinates $VALUE" >> TMP
	VALUE=`wc -l $EXONS_BED_REPEAT_AND_PSEUDO_GENE_FILTER`
	echo "# And with NO overlap with Yale Pseudo60 coordinates $VALUE" >> TMP
	VALUE=`wc -l $EXONS_BED_REPEAT_PSEUDO_GENE_AND_RPKM_FILTER`
	echo "# Transfrags who belong to a transcript with an average across samples of  > 1 rpkm $VALUE" >> TMP
	echo "#" >> $TMP
	echo -e chrom"\t"start"\t"end"\t"name"\t"rpkm"\t"strand"\t"sseqid"\t"evalue"\t"bitscore"\t"pident"\t"length"\t"mismatch"\t"gapopen"\t"qstart"\t"qend"\t"sstart"\t"send"\t"salltitles >> $TMP

	TMP2="$SUB_OUTPUT_FLDR"/tmp2

	cat $TMP $BLASTX_OUT > $TMP2 && mv $TMP2 $BLASTX_OUT
	
}	
	
# A30C
# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_mm/input_data/RPKM/rabt_pooled_output/A30C_w_rpkm_counts.gtf
################# A30C_intERgenic ###########################
# set grep var to find intERgenic exons
GREP_VAR='class_code "u"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="A30C_intERgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUTPUT_FLDR="$OUTPUT_DATA_FOLDER"/A30C_intERgenic
# call the blast function
DOBLAST

# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_mm/input_data/RPKM/rabt_pooled_output/A30C_w_rpkm_counts.gtf
################# A30C_intRAgenic ###########################
# set grep var to find intRAgenic exons
GREP_VAR='class_code "i"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="A30C_intRAgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUTPUT_FLDR="$OUTPUT_DATA_FOLDER"/A30C_intRAgenic
# call the blast function
DOBLAST

# FC30C
# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_mm/input_data/RPKM/rabt_pooled_output/FC30C_w_rpkm_counts.gtf
################# FC30C_intERgenic ###########################
# set grep var to find intERgenic exons
GREP_VAR='class_code "u"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="FC30C_intERgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUTPUT_FLDR="$OUTPUT_DATA_FOLDER"/FC30C_intERgenic
# call the blast function
DOBLAST

# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_mm/input_data/RPKM/rabt_pooled_output/FC30C_w_rpkm_counts.gtf
################# FC30C_intRAgenic ###########################
# set grep var to find intRAgenic exons
GREP_VAR='class_code "i"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="FC30C_intRAgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUTPUT_FLDR="$OUTPUT_DATA_FOLDER"/FC30C_intRAgenic
# call the blast function
DOBLAST

# H30C
# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_mm/input_data/RPKM/rabt_pooled_output/H30C_w_rpkm_counts.gtf
################# H30C_intERgenic ###########################
# set grep var to find intERgenic exons
GREP_VAR='class_code "u"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="H30C_intERgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUTPUT_FLDR="$OUTPUT_DATA_FOLDER"/H30C_intERgenic
# call the blast function
DOBLAST

# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_mm/input_data/RPKM/rabt_pooled_output/H30C_w_rpkm_counts.gtf
################# H30C_intRAgenic ###########################
# set grep var to find intRAgenic exons
GREP_VAR='class_code "i"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="H30C_intRAgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUTPUT_FLDR="$OUTPUT_DATA_FOLDER"/H30C_intRAgenic
# call the blast function
DOBLAST


#####################################################################################
### END OF THE REAL WORK - DO the FINAL HOUSE KEEPING  ##############################
#####################################################################################

# copy the contents of the current folder (with this script and other code) to the saved code folder
cp -R "$CODE_FOLDER" "$SAVED_CODE"
echo end script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 


