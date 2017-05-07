#!/bin/bash
#PBS -A simons
#PBS -l walltime=168:00:00
#PBS -l nodes=1:ppn=12
#PBS -N Mouse_Cufflinks
#PBS -o Mouse_Cufflinks.out
#PBS -e Mouse_Cufflinks.err
#PBS -m abe
#PBS -M jmtroy2@igb.illinois.edu

# this script is meant to be run in the same folder as the .tophat output folders
# use the below qsub command to run
# qsub -S /bin/bash mouse_cufflinks_ncbim37_a_30_c_RABT.sh

cd $PBS_O_WORKDIR

echo `pwd`

module load samtools/0.1.19
module load tophat2/2.0.8
module load bowtie2/2.1.0 
module load cufflinks/2.1.1

gtf=/home/groups/simons/Joe/ensembl37_NCBIM37_genome/genes.gtf
genome=/home/groups/simons/Joe/ensembl37_NCBIM37_genome/genome

echo "Start of Cufflinks"

for file in ./A_30_C*.tophat;
do
        if [ -d "$file" ]; then
	                BASENM=$(basename $file)
			BASENM=${BASENM%%.tophat}
			cufflinks -p 12 -o $BASENM.cufflinks_rabt --GTF-guide $gtf  $file/accepted_hits.bam
	continue
	fi
done

echo "End of Cufflinks"
