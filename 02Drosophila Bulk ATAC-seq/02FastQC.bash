#!/bin/bash
#SBATCH -J QC_ATAC
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=20480
#SBATCH -c 8
fastqc -o /storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq/02FastQC /storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq/01Rawdata/*.fq.gz