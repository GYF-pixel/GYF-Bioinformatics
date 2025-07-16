#!/bin/bash
#SBATCH -J depTSS_GYF
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=200G
#SBATCH -c 20
projPath="/storage/maxianjueLab/guoyifan/Ras_Identification/ATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"

for histName in 7A 7B 7C 7 Ras403A Ras403B Ras403C Ras403 F5dA F5dB F5dC F5d F9dA F9dB F9dC F9d F13dA F13dB F13dC F13d;do
{
     computeMatrix reference-point -S ${projPath}/06BigWig/${histName}_normalized_RPKM_20bin.bw \
							   -p 20 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_data.gz
     plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#925E9F' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_heatmap.pdf
    }
done

