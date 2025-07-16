
#GroupABC_RcisTarget

#load packages
library(Mfuzz)
library(limma)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(devtools)
library(fgsea)
library(SummarizedExperiment)
library(DelayedArray)
library(biomaRt)
library(curl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
library(xlsx)
library(RcisTarget)

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\09WGCNA_GroupABC\\05RcisTarget")
list.files()

#1. Gene list
gene_list <- read.xlsx2(file = "12 RcisTarget_Input.xlsx",sheetIndex = 1,header = T)

geneList1 = gene_list

#2. Gene annotation
mart <- useEnsembl(biomart = "ensembl", dataset = 'dmelanogaster_gene_ensembl',mirror = "uswest")         	#This is useEnsembl, mirror can be equal to "uswest" or "asia"
my_ensembl_gene_id<-gene_list[,i]
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
head(study_symbols$external_gene_name)
a = study_symbols$external_gene_name
write.csv(a, "13 RcisTarget_Input.csv")

geneLists <- list(study_symbols$external_gene_name)

#3. RcisTarget Enrichment
geneList1 = read.xlsx2(file = "12 RcisTarget_Input.xlsx",sheetIndex = 1,header = T)
#geneList1 = read.xlsx2(file = "12 RcisTarget_All_Input.xlsx",sheetIndex = 1,header = T)


for (i in 1:length(colnames(geneList1)))
{
#Annotation
#mart <- useEnsembl(biomart = "ensembl", dataset = 'dmelanogaster_gene_ensembl',mirror = "uswest")         
#my_ensembl_gene_id<-geneList1[,i]
#study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
#head(study_symbols)

#Genelist
geneLists <- list(geneListName=geneList1[,i]) 

#Loading
data(motifAnnotations_dmel)

motifRankings <- importRankings("./dm6_v10_clust.genes_vs_motifs.rankings.feather")				#downloaded from offcial web

#Enrichment
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                               motifAnnot=motifAnnotations)
							   
file_name <- paste('14 RcisTarget_TF_Enrichment',colnames(geneList1)[i],".csv")

write.csv(motifEnrichmentTable_wGenes,file_name,quote=F)
cat("File", file_name, "saved.\n")    ###check save
}

