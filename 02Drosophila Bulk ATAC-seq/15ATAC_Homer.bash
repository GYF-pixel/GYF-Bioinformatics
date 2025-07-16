#!/bin/bash
#SBATCH -J Homer
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"
ref="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"

for histName in 7A 7B 7C 7 Ras403A Ras403B Ras403C Ras403 F5dA F5dB F5dC F5d F9dA F9dB F9dC F9d F13dA F13dB F13dC F13d;do
{
findMotifsGenome.pl ${projPath}/09Homer/macs2_${histName}_peak_q0.05_peaks_homer.tmp ${ref} ${projPath}/09Homer/${histName} -size 200 -len 8,10,12,15,18 -p 12
  }
done

