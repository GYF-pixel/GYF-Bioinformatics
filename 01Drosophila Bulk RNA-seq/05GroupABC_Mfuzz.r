
#GroupABC_Mfuzz

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

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\09WGCNA_GroupABC")
list.files()

#1 Gene list
gene_list <- read.xlsx2(file = "08 Mfuzz_GroupA.xlsx",sheetIndex = 1,header = T)
gene_list <- read.xlsx2(file = "08 Mfuzz_GroupB.xlsx",sheetIndex = 1,header = T)
gene_list <- read.xlsx2(file = "08 Mfuzz_GroupC.xlsx",sheetIndex = 1,header = T)

#2 Gene expression count
#data = read.csv("newdata_filter.csv", header = T)
data.filter = read.csv("06 WCGNA_GroupA_Count.csv", header = T,row.names = 1)
data.filter = read.csv("06 WCGNA_GroupB_Count.csv", header = T,row.names = 1)
data.filter = read.csv("06 WCGNA_GroupC_Count.csv", header = T,row.names = 1)

group_info = c(rep("WT",3), rep("Ras",9), rep("5d_Tumors",18), rep("9d_Tumors",18), rep("13d_Tumors",18))			#group A
group_info = c(rep("WT",3), rep("Ras",9), rep("5d_Tumors",6), rep("9d_Tumors",6), rep("13d_Tumors",6))				#group B
group_info = c(rep("WT",3), rep("Ras",9), rep("5d_Tumors",6), rep("9d_Tumors",6), rep("13d_Tumors",6))				#group C

#3 Input
dat=data.filter[gene_list[,1],] ##extract genes in gene list

dat <- log2(edgeR::cpm(dat)+1)
dat[1:4,1:4]
dim(dat)

avereps_df  <- t(limma::avereps( t(dat) , ID = group_info))  ##Take the average of the expression values for the same time series

#4.Mfuzz

#4.1 Filtering---

#Steps such as removing genes with low expression levels or small changes between different time points
#Mfuzz clustering requires an object of ExpressionSet type, so it is necessary to first construct such an object using expression vectors.
eset <- new("ExpressionSet",exprs = avereps_df)

#Processing NA values
eset <- filter.NA(eset, thres = 0.25)
eset <- fill.NA(eset, mode = 'mean')

#Remove genes with small differences between samples based on standard deviation
#The number of genes removed varies among different datasets
eset <- filter.std(eset,min.std=0)
##10818 genes excluded. 
eset


#4.2 Standardisation----

#When clustering, a numerical value is needed to represent the distance between different genes, and the Euclidean distance is used in Mfuzz,
#Due to the fact that the definition of ordinary Euclidean distance does not take into account the differences in scales between different dimensions, standardization is required first
eset <- standardise(eset)

#4.3 Setting of parameters for FCM clustering----
#The clustering algorithm in Mfuzz requires two parameters to be provided,
#The first parameter is the desired number of clusters to be obtained, which we directly specify
#The second parameter is called the fuzzier value, represented by the lowercase letter m, and can be evaluated as the optimal value through a function
c <- 10
m <- mestimate(eset) 						#Evaluate the optimal value of m
set.seed(123)
cl <- mfuzz(eset, c = c, m = m) 		#The result of clustering without setting membership


#4.4 visualise----

#Draw a graph, set the time.labels parameter timeline, and correspond it to the columns in the original gene expression dataset
library(RColorBrewer)
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
pdf('mfuzz_clusters_plot.pdf',height = 5.5,width = 14)

#Mfuzz.flot is the simplest chart
mfuzz.plot(eset,cl,mfrow=c(2,5),
           new.window= FALSE,
           time.labels= colnames(eset) ,
           colo = color.2)
dev.off()
#Mfuzz.dlot2 adjustable graphs

#Do not filter any genes min.mem=0
pdf('mfuzz_clusters_plot_min_mem0.pdf',height = 5.5,width = 14)
mfuzz.plot2(eset,cl,mfrow=c(2,5), min.mem=0,
#           colo="fancy",
           time.labels= colnames(eset),
		   xlab="Time",ylab="Expression changes",
		   x11=FALSE,  	#	If TRUE, a new window will be open for plotting
           ax.col="black",bg = "white",col.axis="black",col.lab="black",
           col.main="black",col.sub="black",col="black",
           Xwidth=9,Xheight=9,
		   ylim.set=c(-2,2),
		   centre=TRUE, # TRUE--lines for cluster centres
		   centre.col="black",centre.lwd=2,
		   single=FALSE  # Integer if a specific cluster is to be plotted, otherwise it should be set to FALSE.
           )
dev.off()


#Filter genes with membership<0.2 min.mem=0.2
pdf('mfuzz_clusters_plot_min_mem02.pdf',height = 5.5,width = 14)
mfuzz.plot2(eset,cl,mfrow=c(2,5), min.mem=0.2,
#           colo="fancy",
           time.labels= colnames(eset),
		   xlab="Time",ylab="Expression changes",
		   x11=FALSE,  	#	If TRUE, a new window will be open for plotting
           ax.col="black",bg = "white",col.axis="black",col.lab="black",
           col.main="black",col.sub="black",col="black",
           Xwidth=9,Xheight=9,
		   ylim.set=c(-2,2),
		   centre=TRUE, # TRUE--lines for cluster centres
		   centre.col="black",centre.lwd=2,
		   single=FALSE  # Integer if a specific cluster is to be plotted, otherwise it should be set to FALSE.
           )
dev.off()

#Filter genes with membership<0.5 min.mem=0.5
pdf('mfuzz_clusters_plot_min_mem05.pdf',height = 5.5,width = 14)
mfuzz.plot2(eset,cl,mfrow=c(2,5), min.mem=0.5,
#           colo="fancy",
           time.labels= colnames(eset),
		   xlab="Time",ylab="Expression changes",
		   x11=FALSE,  	#	If TRUE, a new window will be open for plotting
           ax.col="black",bg = "white",col.axis="black",col.lab="black",
           col.main="black",col.sub="black",col="black",
           Xwidth=9,Xheight=9,
		   ylim.set=c(-2,2),
		   centre=TRUE, # TRUE--lines for cluster centres
		   centre.col="black",centre.lwd=2,
		   single=FALSE  # Integer if a specific cluster is to be plotted, otherwise it should be set to FALSE.
           )
dev.off()


#4.5 glimpse results----

#The complete clustering results are saved in the cl object, and the common operations for this object are as follows

cl$size 							#View the number of genes in each cluster
table(cl$cluster)
cluster_gene_number = table(cl$cluster)
write.csv(cluster_gene_number, "cluster_gene_number.csv")

cl$cluster[cl$cluster == 1] 		#Extract genes from a certain cluster

## cluster cores
# membership values can also indicate the similarity of vectors to each other.
eset
# extracts genes forming the alpha cores of soft clusters
cl.thres <- acore(eset,cl,min.acore=0.2)   #min.acore--minimum membership values of gene belonging to the cluster core.

head(cl.thres[[1]]) #View genes in cluster1 after filtering
length(cl.thres)

cl.thres[[1]]

#4.6 save results----

###Gene annotation
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))

q <-  cl.thres   #Save cl.thres as object q
for (i in 1:length(q))
{
  files <- q[[i]]
  my_ensembl_gene_id<-row.names(files)
  study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
  head(study_symbols)
  ensembl_gene_id<-rownames(files)
  files <- cbind(ensembl_gene_id,files)
  colnames(files)[1]<-c("ensembl_gene_id")
  #The output results of all genes' edgers - diff_2
  diff_name <- merge(files,study_symbols,by="ensembl_gene_id")
  q[[i]] <- diff_name
}

###for loop to save the object q

for (i in 1:length(q))
{
  newdata <- q[[i]]
  file_name <- paste('mfuzz_filter0.2_',i,".csv")
  write.csv(newdata,file_name, row.names = FALSE)
  cat("File", file_name, "saved.\n")    ###check save
}

