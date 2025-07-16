
#ChIPseeker

#load packages
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(clusterProfiler)
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot

#work dir
getwd()
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\04Fry_specific\\02Fry ATAC-seq")

#load ref GTF file downloaded from Ensemble
txdb <- makeTxDbFromGFF("D:\\OneDrive - 西湖大学\\04Codes\\00RefGenome\\Drosophila_melanogaster\\Drosophila_melanogaster.BDGP6.32.109.gtf",
                        format="gtf")    #GTF and GFF3 can be used
						
#If interested, you can use the following methods to explore what content txdb contains
keytypes(txdb)    
keys(txdb)

#1. Read in a single summary file (summits)
peaks <- readPeakFile("macs2_F5d_peak_q0.05_summits.bed")

#2. Structural annotation
peakAnno <- annotatePeak(peaks,
                         TxDb=txdb,
                         tssRegion=c(-2000, 2000))

#3. After the annotation is completed, visualize and choose from multiple images
#3.1
plotAnnoBar(peakAnno)
#3.2
plotDistToTSS(peakAnno)
#3.3
vennpie(peakAnno)
#3.4
plotAnnoPie(peakAnno)
#3.5
#install.packages("ggupset")
library(ggupset)
upsetplot(peakAnno)
#3.6
#install.packages("ggimage")
library(ggimage)
upsetplot(peakAnno, vennpie=TRUE)
 
#4. Output files
#Finally, convert our annotation results into a data box for easy viewing
df <- as.data.frame(peakAnno)
#Extract the annotated genes (column 14) for subsequent functional analysis
gene <- df[,14]

write.csv(gene, "17 ChIPseeker_annotated_genes.csv")

