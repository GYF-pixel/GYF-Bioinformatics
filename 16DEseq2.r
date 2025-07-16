
#ChIPseeker

#load packages
library(methods)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(clusterProfiler)
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomicAlignments)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot
library(rtracklayer)
library(ChIPQC)
library(DiffBind)
library(Rsubread)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(tracktables)

#work dir
getwd()
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\04Fry_specific\\02Fry ATAC-seq\\07MACS2")

#load peaks called by MACS2
peaks <- dir("./", pattern = "*.narrowPeak", full.names = TRUE)
peaks

#Using the apply function to obtain the peaks matrix
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
# myPeaks is a list object

allPeaksSet_nR <- reduce(unlist(GRangesList(myPeaks)))
overlap <- list()
for (i in 1:length(myPeaks)) {
  overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}

# create a matrix of presence/absence of these peaks over each sample.
overlapMatrix <- do.call(cbind, overlap)
colnames(overlapMatrix) <- basename(peaks)
mcols(allPeaksSet_nR) <- overlapMatrix
allPeaksSet_nR[1:2, ]

#Ensure that peaks are expressed in at least two samples
occurrences <- rowSums(as.data.frame(elementMetadata(nrToCount)))
nrToCount <- nrToCount[occurrences >= 2, ]
nrToCount

#Use bam file to form myCounts
bamsToCount <- dir("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\04Fry_specific\\02Fry ATAC-seq\\05Bam", full.names = TRUE, pattern = "*.\\.bam$")
myCounts <- summarizeOverlaps(nrToCount, bamsToCount, singleEnd = FALSE)
colnames(myCounts) <- c("7A","7B","7C",
										"Ras403A","Ras403B","Ras403C",
										"F5dA","F5dB","F5dC","F9dA",
										"F9dB","F9dC","F13dA",
										"F13dB","F13dC")
										
#check myCounts
temp <- assay(myCounts)
dim(temp)
head(temp)

#DESeq2

#group information
Group <- factor(c(rep("WT",3), rep("Ras",3), rep("RasFry_5d",3), rep("RasFry_9d",3), rep("RasFry_13d",3)))
Group

metaData <- data.frame(Group, row.names = colnames(myCounts))
metaData

#Build a DESeq2 object
atacDDS <- DESeqDataSetFromMatrix(assay(myCounts), metaData, ~Group, rowRanges = rowRanges(myCounts))
atacDDS <- DESeq(atacDDS)

#Differential accessibility analysis between each two groups
Diffpeak <- results(atacDDS, c("Group", "Ras", "WT"), format = "GRanges")									#control in the later
Diffpeak <- results(atacDDS, c("Group", "RasFry_5d", "Ras"), format = "GRanges")						#control in the later
Diffpeak <- results(atacDDS, c("Group", "RasFry_9d", "RasFry_5d"), format = "GRanges")			#control in the later
Diffpeak <- results(atacDDS, c("Group", "RasFry_13d", "RasFry_9d"), format = "GRanges")			#control in the later

Diffpeak <- Diffpeak[order(Diffpeak$pvalue)]
Diffpeak

#Extract only the open area near the promoter
#upstream and downstream 2kb were taken into account
toOverLap <- promoters(TxDb.Dmelanogaster.UCSC.dm3.ensGene, 2000, 2000)   
Diffpeak <- Diffpeak[(!is.na(Diffpeak$padj) &
                                          Diffpeak$padj < 0.05) & Diffpeak %over% toOverLap, ]

#Generate html report
myReport <- makebedtable(Diffpeak, "Diffpeak.html", getwd())
browseURL(myReport)


