#!/bin/bash
#SBATCH -J SamV_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"
cores=8

for histName in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
  	 samtools view -h -b -@ 8 -f 3 -F 12 -F 256 -q 20 ${projPath}/05Alignment/sam/${histName}_bowtie2.sam > ${projPath}/05Alignment/bam/${histName}_bowtie2_filter.bam
  }
done
