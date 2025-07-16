#!/bin/bash
#SBATCH -J Pic_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"
for histName in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
	 picard MarkDuplicates I=${projPath}/05Alignment/bam/${histName}_bowtie2_filter_sorted.bam \
	 O=${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.bam \
	 M=${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.matrix \
	 ASO=coordinate REMOVE_DUPLICATES=true >${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.log 2>&1
  }
done