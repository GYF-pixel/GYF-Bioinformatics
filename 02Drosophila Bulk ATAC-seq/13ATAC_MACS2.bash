#!/bin/bash
#SBATCH -J Macs_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"

for histName in 7A 7B 7C 7 Ras403A Ras403B Ras403C Ras403 F5dA F5dB F5dC F5d F9dA F9dB F9dC F9d F13dA F13dB F13dC F13d;do
  {
macs2 callpeak -t ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam \
      -g dm -f BAMPE --nomodel --shift -100 --extsize 200 -n macs2_${histName}_peak_q0.05 --outdir $projPath/07MACS2 -q 0.05 --keep-dup all 2>${projPath}/07MACS2/macs2Peak_${histName}_summary.txt
    }
done
