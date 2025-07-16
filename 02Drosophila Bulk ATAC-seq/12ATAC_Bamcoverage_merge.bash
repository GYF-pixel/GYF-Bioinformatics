#!/bin/bash
#SBATCH -J merg_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"

samtools merge 7_bowtie2_filter_sorted_rmdup_rmmito.bam 7A_bowtie2_filter_sorted_rmdup_rmmito.bam 7B_bowtie2_filter_sorted_rmdup_rmmito.bam 7C_bowtie2_filter_sorted_rmdup_rmmito.bam &
samtools merge Ras403_bowtie2_filter_sorted_rmdup_rmmito.bam Ras403A_bowtie2_filter_sorted_rmdup_rmmito.bam Ras403B_bowtie2_filter_sorted_rmdup_rmmito.bam Ras403C_bowtie2_filter_sorted_rmdup_rmmito.bam &
samtools merge F5d_bowtie2_filter_sorted_rmdup_rmmito.bam F5dA_bowtie2_filter_sorted_rmdup_rmmito.bam F5dB_bowtie2_filter_sorted_rmdup_rmmito.bam F5dC_bowtie2_filter_sorted_rmdup_rmmito.bam &
samtools merge F9d_bowtie2_filter_sorted_rmdup_rmmito.bam F9dA_bowtie2_filter_sorted_rmdup_rmmito.bam F9dB_bowtie2_filter_sorted_rmdup_rmmito.bam F9dC_bowtie2_filter_sorted_rmdup_rmmito.bam &
samtools merge F13d_bowtie2_filter_sorted_rmdup_rmmito.bam F13dA_bowtie2_filter_sorted_rmdup_rmmito.bam F13dB_bowtie2_filter_sorted_rmdup_rmmito.bam F13dC_bowtie2_filter_sorted_rmdup_rmmito.bam &

 for histName in Ras403 7 F13d F9d F5d;do
   {
     samtools index $projPath/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_raw.bw
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_BPM_20bin.bw --binSize 20 --normalizeUsing BPM
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_RPKM_20bin.bw --binSize 20 --normalizeUsing RPKM
   }
done
