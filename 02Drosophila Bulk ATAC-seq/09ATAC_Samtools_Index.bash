#!/bin/bash
#SBATCH -J SamIn_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=60G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"
for histName in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
    samtools view -h -@ 12 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.bam | grep -v 'mitochondrion_genome' | samtools view -b -@ 12 -o ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam 

for histName in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
	samtools index -@ 12 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam
  }
done