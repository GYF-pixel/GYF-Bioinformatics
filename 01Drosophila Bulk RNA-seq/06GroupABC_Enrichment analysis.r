# KEGG and GO enrichment

#load packages
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(stringr)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(Cairo)
library(enrichplot)
library(xlsx)

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\09WGCNA_GroupABC")
list.files()

#load gene list
gene_list <- read.xlsx2(file = "geneList.xlsx",sheetIndex = 1,header = T)
diff_nameSig = gene_list

#1. GO Enrichment 
head(diff_nameSig)
#keyType = "ENSEMBL"
ALL <- enrichGO(gene = diff_nameSig$ensembl_gene_id, 
                OrgDb = org.Dm.eg.db, 
                keyType = "ENSEMBL",
                ont = 'ALL',
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",  
                qvalueCutoff  = 0.1, readable=T) 
				
write.csv(as.data.frame(ALL@result), file="Part2_GOALL.csv",sep="\t", quote=FALSE)
GOFile = read.csv("Part2_GOALL.csv", header = T, sep="\t")

Goplot <- dotplot(ALL, split = "ONTOLOGY", font.size = 8, showCategory = 10) + 
				facet_grid(ONTOLOGY ~ ., scale = "free") + 
				scale_y_discrete(labels = function(x) str_wrap(x, width =50)) + 
                scale_size(range=c(2, 6))   #set dot size

ggsave("Part2_Goplot.pdf", plot = Goplot, width = 8, height = 8) 
head(ALL,1);dim(ALL)				

#Biological process
sum(ALL$ONTOLOGY=="BP") 
#Cellular component
sum(ALL$ONTOLOGY=="CC") 
#Molecular function
sum(ALL$ONTOLOGY=="MF") 

#2. KEGG analysis
library(stringr)
library(DOSE)
data(geneList, package="DOSE")

columns(org.Dm.eg.db)

gene.df <- bitr(diff_nameSig$ensembl_gene_id, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Dm.eg.db) 
			  
gene.kegg <- bitr_kegg(gene.df$ENTREZID,fromType="ncbi-geneid",
                        toType="kegg",organism='dme')
head(gene.kegg)

options(clusterProfiler.download.method = "auto")

ekegg <- enrichKEGG(gene = gene.kegg$kegg, 
                    organism = "dme", 
                    keyType = "kegg",
                    pvalueCutoff  = 0.05,
                    pAdjustMethod = "BH",  
                    qvalueCutoff  = 0.1)  




