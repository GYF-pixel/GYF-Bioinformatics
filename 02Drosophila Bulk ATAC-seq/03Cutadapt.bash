#!/bin/bash
#SBATCH -J ctad_ATAC
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=30720
#SBATCH -c 8
for i in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
cutadapt -j 8 --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 25 \
-a CTGTCTCTTATACACATC \
-A CTGTCTCTTATACACATC \
-o /storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq/03Cutadapt/${i}_R1_cutadapt.fq.gz \
-p /storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq/03Cutadapt/${i}_R2_cutadapt.fq.gz \
/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq/01Rawdata/${i}_R1.fq.gz \
/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq/01Rawdata/${i}_R2.fq.gz > /storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq/03Cutadapt/${i}_cutadapt_infor.log 2>&1
  }
done

