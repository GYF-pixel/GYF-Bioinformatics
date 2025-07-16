#!/bin/bash
#SBATCH -J bCov_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"

for histName in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_raw.bw
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_BPM_20bin.bw --binSize 20 --normalizeUsing BPM
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_RPKM_20bin.bw --binSize 20 --normalizeUsing RPKM
  }
done
