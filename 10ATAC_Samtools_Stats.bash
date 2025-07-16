#!/bin/bash
#SBATCH -J Stats_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"
for histName in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
     samtools stats -@ 8 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam > ${projPath}/05Alignment/Final_QC/Mapping_Stats/${histName}_bowtie2_filter_sorted_rmdup_rmmito.mapping_stats
  }
done