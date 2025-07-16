#!/bin/bash
#SBATCH -J bwt_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=160G
#SBATCH -c 20
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"
ref="/storage/maxianjueLab/guoyifan/species_reference/DroMelanogaster/Drosophila_melanogaster.BDGP6.32"
cores=20

for histName in 7A 7B 7C Ras403A Ras403B Ras403C F5dA F5dB F5dC F9dA F9dB F9dC F13dA F13dB F13dC;do
  {
		bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 0 -X 2000 --no-unal -p ${cores} -x ${ref} \
		-1 ${projPath}/03Cutadapt/${histName}_R1_cutadapt.fq.gz \
		-2 ${projPath}/03Cutadapt/${histName}_R2_cutadapt.fq.gz \
		-S ${projPath}/05Alignment/sam/${histName}_bowtie2.sam \
		> ${projPath}/05Alignment/sam/bowtie2_summary/${histName}_bowtie2_infor.log 2>&1 
  }
done

