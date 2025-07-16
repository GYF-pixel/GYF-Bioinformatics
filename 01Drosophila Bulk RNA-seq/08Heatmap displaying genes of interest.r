
#load package
library(ggplot2)
library(tidyverse)
library(devtools)
library(biomaRt)
library(curl)
library(ggpubr)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(ggunchull)
library(scales)

#setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigeneisis\\04Fry_specific\\01Fry RNA-seq\\01Analysis")                    #设置工作目录

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\04Analysis")
list.files()         

#1. Count file arrangement

#input files after batheffect removed
data.filter = read.csv("newdata_filter_remove_pre_batheffect_removed.csv", header = T,row.names = 1)

#condition table
colnames(data.filter)
data.filter <- data.filter %>% dplyr::select(
#wildtype
"WT_79EA","WT_79EB","WT_79EC",
#Ras group
"Ras1A","Ras2A","Ras3A",
"X4.03A","X4.03B","X4.03C",
"Ras82BA","Ras82BB","Ras82BC",
#5d_tumors
"X4.41A1A","X4.41A2A","X4.41A3A",
"X11.231A1A","X11.231A2A","X11.231A3A",
"X6.128A1A","X6.128A2A","X6.128A3A",
"X6.144A1A","X6.144A2A","X6.144A3A",
"X6.147A1A","X6.147A2A","X6.147A3A",
"X7.51A1A","X7.51A2A","X7.51A3A",
"X6.119A1A","X6.119A2A","X6.119A3A",
"X4.47A1A","X4.47A2A","X4.47A3A",
"X4.90A1A","X4.90A2A","X4.90A3A",
"X3.199A1A","X3.199A2A","X3.199A3A",
#9d_tumors
"X4.41B1A","X4.41B2A","X4.41B3A",
"X11.231B1A","X11.231B2A","X11.231B3A",
"X6.128B1A","X6.128B2A","X6.128B3A",
"X6.144B1A","X6.144B2A","X6.144B3A",
"X6.147B1A","X6.147B2A","X6.147B3A",  
"X7.51B1A","X7.51B2A","X7.51B3A",
"X6.119B1A","X6.119B2A","X6.119B3A",
"X4.47B1A","X4.47B2A","X4.47B3A",
"X4.90B1A","X4.90B2A","X4.90B3A",
"X3.199B1A","X3.199B2A","X3.199B3A",  
#13d_tumors
"X4.41C1A","X4.41C2A","X4.41C3A",
"X11.231C1A","X11.231C2A","X11.231C3A",
"X6.128C1A","X6.128C2A","X6.128C3A",
"X6.144C1A","X6.144C2A","X6.144C3A",
"X6.147C1A","X6.147C2A","X6.147C3A",
"X7.51C1A","X7.51C2A","X7.51C3A",
"X6.119C1A","X6.119C2A","X6.119C3A",
"X4.47C1A","X4.47C2A","X4.47C3A",
"X4.90C1A","X4.90C2A","X4.90C3A",
"X3.199C1A","X3.199C2A","X3.199C3A"
)

#group information
combat_Expr = data.filter
colnames(combat_Expr)

group_info = c(rep("WT",3), rep("40ARas",3), rep("Ras79E",3), rep("Ras82B",3), 
rep("Scrib_5d",3), rep("lgl_5d",3), rep("Vps36_5d",3), rep("Syx7_5d",3),rep("Rabex5_5d",3), rep("TSG101_5d",3), rep("fry_5d",3), rep("fmt_5d",3), rep("emei_5d",3), rep("msn_5d",3),
rep("Scrib_9d",3), rep("lgl_9d",3), rep("Vps36_9d",3), rep("Syx7_9d",3),rep("Rabex5_9d",3), rep("TSG101_9d",3), rep("fry_9d",3), rep("fmt_9d",3), rep("emei_9d",3), rep("msn_9d",3),
rep("Scrib_13d",3), rep("lgl_13d",3), rep("Vps36_13d",3), rep("Syx7_13d",3),rep("Rabex5_13d",3), rep("TSG101_13d",3), rep("fry_13d",3), rep("fmt_13d",3), rep("emei_13d",3), rep("msn_13d",3)
)


#2. Heatmap input data preparation

#Data processing gene count matrix
dat=combat_Expr 
dat <- log2(edgeR::cpm(dat)+1)		#log2(CPM+1)
dat[1:4,1:4]
dim(dat)

#The expression values of the same genotype in the same time series were averaged.
avereps_df  <- t(limma::avereps( t(dat) , ID = group_info))

#gene list to display in Heatmap
genes <- read.table(file="genes.txt",header = T)		#gene symbol form

#Gene annotation
mart <- useEnsembl(biomart = "ensembl", dataset = 'dmelanogaster_gene_ensembl',mirror = "asia")         
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'external_gene_name', values = genes, mart = mart)
head(study_symbols)

#Generate heatmap input file
Expression <- avereps_df[study_symbols[,1],]
rownames(study_symbols) <- study_symbols$ensembl_gene_id
rownames(Expression) <- study_symbols[rownames(Expression),"external_gene_name"]
write.csv(Expression,file ="Expression.csv")			


#3. Draw Heatmap

#Set color palette
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

#Output heatmap
input_heatmap_regroup = Expression

heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="",
                 fontsize = 12,
                 scale="row",
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 20,
                 cellheight = 12,cellwidth = 12,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap

ggsave("Part3_Heatmap_genes.pdf", plot = heatmap, width = 25, height = 45) 



