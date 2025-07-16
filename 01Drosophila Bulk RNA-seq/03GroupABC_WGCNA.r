
#GroupABC_WGCNA

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
library(SCENIC)
library(scales)
library(WGCNA)
library(reshape2)
library(stringr)

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\09WGCNA_GroupABC")
list.files()
options(stringsAsFactors = FALSE)

#Open multithreading
enableWGCNAThreads()
#Allowing parallel execution with up to 12 working processes.

#1. Gene expression data and phenotype data
#load expression data
data = read.csv("newdata_filter_remove_pre_batheffect_removed.csv", header = T,row.names = 1)
dat <- log2(edgeR::cpm(data)+1)
dat[1:4,1:4]
dim(dat)

#load phenotype data
trait <- read.csv(file="WGCNA Traits Input.csv", header=T, row.names=1)
rownames(trait)		#The first 17 lines are phenotype data
colnames(trait)		#The column names are the same as those of gene expression data


#1.1 Group A
data.filter = as.data.frame(dat)
colnames(data.filter)
data.filter <- data.filter %>% dplyr::select(
"WT_79EA","WT_79EB","WT_79EC",  

"X4.03A","X4.03B","X4.03C",  
"Ras82BA","Ras82BB","Ras82BC",
"Ras1A","Ras2A","Ras3A",

"X4.41A1A","X4.41A2A","X4.41A3A",					#scrib
"X11.231A1A","X11.231A2A","X11.231A3A",  		#lgl
"X6.128A1A","X6.128A2A","X6.128A3A",				#Vps36
"X6.119A1A","X6.119A2A","X6.119A3A",				#fry
"X4.47A1A","X4.47A2A","X4.47A3A",					#fmt
"X4.90A1A","X4.90A2A","X4.90A3A",					#emei

"X4.41B1A","X4.41B2A","X4.41B3A",				
"X11.231B1A","X11.231B2A","X11.231B3A",		
"X6.128B1A","X6.128B2A","X6.128B3A",
"X6.119B1A","X6.119B2A","X6.119B3A",
"X4.47B1A","X4.47B2A","X4.47B3A",
"X4.90B1A","X4.90B2A","X4.90B3A",

"X4.41C1A","X4.41C2A","X4.41C3A",
"X11.231C1A","X11.231C2A","X11.231C3A",
"X6.128C1A","X6.128C2A","X6.128C3A",
"X6.119C1A","X6.119C2A","X6.119C3A", 
"X4.47C1A","X4.47C2A","X4.47C3A",
"X4.90C1A","X4.90C2A","X4.90C3A"
)

colnames(data.filter)					#The column names correspond to the tumor samples of the corresponding group
write.csv(data.filter, "01 WCGNA_GroupA_expression.csv",row.names = T)

trait.filter <- trait[c(1:17), colnames(data.filter)]			#The first 17 lines are phenotype data
write.csv(trait.filter, "01 WCGNA_GroupA_trait.csv",row.names = T)


#1.2 Group B
data.filter = as.data.frame(dat)
colnames(data.filter)
data.filter <- data.filter %>% dplyr::select(
"WT_79EA","WT_79EB","WT_79EC",  

"X4.03A","X4.03B","X4.03C",  
"Ras82BA","Ras82BB","Ras82BC",
"Ras1A","Ras2A","Ras3A",

"X3.199A1A","X3.199A2A","X3.199A3A",				#msn
"X6.144A1A","X6.144A2A","X6.144A3A",				#Syx7

"X3.199B1A","X3.199B2A","X3.199B3A",  
"X6.144B1A","X6.144B2A","X6.144B3A",

"X3.199C1A","X3.199C2A","X3.199C3A",
"X6.144C1A","X6.144C2A","X6.144C3A"

)

colnames(data.filter)
write.csv(data.filter, "01 WCGNA_GroupB_expression.csv",row.names = T)

trait.filter <- trait[c(1:17), colnames(data.filter)]			#The first 17 lines are phenotype data
write.csv(trait.filter, "01 WCGNA_GroupB_trait.csv",row.names = T)


#1.3 Group C
data.filter = as.data.frame(dat)
colnames(data.filter)
data.filter <- data.filter %>% dplyr::select(
"WT_79EA","WT_79EB","WT_79EC",  

"X4.03A","X4.03B","X4.03C",  
"Ras82BA","Ras82BB","Ras82BC",
"Ras1A","Ras2A","Ras3A",

"X6.147A1A","X6.147A2A","X6.147A3A",				#Rabex-5
"X7.51A1A","X7.51A2A","X7.51A3A",					#TSG101

"X6.147B1A","X6.147B2A","X6.147B3A",  
"X7.51B1A","X7.51B2A","X7.51B3A",

"X6.147C1A","X6.147C2A","X6.147C3A",
"X7.51C1A","X7.51C2A","X7.51C3A"
)

colnames(data.filter)
write.csv(data.filter, "01 WCGNA_GroupC_expression.csv",row.names = T)

trait.filter <- trait[c(1:17), colnames(data.filter)]			##The first 17 lines are phenotype data
write.csv(trait.filter, "01 WCGNA_GroupC_trait.csv",row.names = T)



#2. WGCNA

#Read the feature dataset
#dataExpr = read.csv(file="01 WCGNA_GroupA_expression.csv", header=T, row.names=1)
#dataExpr = read.csv(file="01 WCGNA_GroupB_expression.csv", header=T, row.names=1)
dataExpr = read.csv(file="01 WCGNA_GroupC_expression.csv", header=T, row.names=1)

#traits = read.csv(file="01 WCGNA_GroupA_trait.csv", header=T, row.names=1)
#traits = read.csv(file="01 WCGNA_GroupB_trait.csv", header=T, row.names=1)
traits = read.csv(file="01 WCGNA_GroupC_trait.csv", header=T, row.names=1)
#avereps_df  <- t(limma::avereps( t(dat) , ID = group_info))				##Take the average of the expression values for the same time series (Optionally)


#2.1 Condition Setting

#Official recommendation for "signed" or "signed hybrid"
#To be consistent with the original document, no modifications were made
type = "unsigned"


#2.2 Correlation calculation

#Official recommendation: Biweight mid correction&Bicor
# corType: pearson or bicor
#To be consistent with the original document, no modifications were made
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
#When calculating the correlation of binary variables, such as sample trait information,
#When gene expression is heavily dependent on disease status, the following parameters need to be set
maxPOutliers = ifelse(corType=="pearson",1,0.05)

#When associating binary variables of sample traits, set
robustY = ifelse(corType=="pearson",T,F)


#2.2 Data filtering

#Select genes with a median absolute deviation of at least 75% and a MAD greater than 0.01
#After filtering, the computational load will be reduced and some information will be lost
#You can also skip filtering and make MAD greater than 0
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                 max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
				 
#Convert to a matrix with samples in rows and genes in columns
dataExpr <- as.data.frame(t(dataExprVar))

#Detecting missing values
gsg = goodSamplesGenes(dataExpr, verbose = 3)
##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

## [1]  134 2697

head(dataExpr)[,1:8]


#2.3 Soft threshold screening

#Check if there are any outlier samples
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

pdf("01 sampleTree.pdf",width = 20,height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
text(6.5,4)
dev.off()

#The filtering principle of soft threshold is to make the constructed network more in line with the characteristics of scale-free networks.
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
## pickSoftThreshold: will use block size 2697.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 2697 of 2697
##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.1370  0.825          0.412 587.000  5.95e+02  922.0
## 2      2   0.0416 -0.332          0.630 206.000  2.02e+02  443.0
## 3      3   0.2280 -0.747          0.920  91.500  8.43e+01  247.0
## 4      4   0.3910 -1.120          0.908  47.400  4.02e+01  154.0
## 5      5   0.7320 -1.230          0.958  27.400  2.14e+01  102.0
## 6      6   0.8810 -1.490          0.916  17.200  1.22e+01   83.7
## 7      7   0.8940 -1.640          0.869  11.600  7.29e+00   75.4
## 8      8   0.8620 -1.660          0.827   8.250  4.56e+00   69.2
## 9      9   0.8200 -1.600          0.810   6.160  2.97e+00   64.2
## 10    10   0.8390 -1.560          0.855   4.780  2.01e+00   60.1
## 11    12   0.8020 -1.410          0.866   3.160  9.61e-01   53.2
## 12    14   0.8470 -1.340          0.909   2.280  4.84e-01   47.7
## 13    16   0.8850 -1.250          0.932   1.750  2.64e-01   43.1
## 14    18   0.8830 -1.210          0.922   1.400  1.46e-01   39.1
## 15    20   0.9110 -1.180          0.926   1.150  8.35e-02   35.6
## 16    22   0.9160 -1.140          0.927   0.968  5.02e-02   32.6
## 17    24   0.9520 -1.120          0.961   0.828  2.89e-02   29.9
## 18    26   0.9520 -1.120          0.944   0.716  1.77e-02   27.5
## 19    28   0.9380 -1.120          0.922   0.626  1.08e-02   25.4
## 20    30   0.9620 -1.110          0.951   0.551  6.49e-03   23.5			

#Save as image using PDF () and dev. off()
pdf("02 Scale independence and Mean connectivity plot.pdf",width = 8,height = 5)   

par(mfrow = c(1,2))   #Image segmentation into 1 row and 2 columns
cex1 = 0.9 				 #Font size setting

#The horizontal axis represents Soft threshold (power), and the vertical axis represents the evaluation parameters of scale-free networks. The higher the value, the more the network conforms to non scale characteristics
head(sft[["fitIndices"]])
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
	 
#Screening criteria. R-square=0.85
abline(h=0.85,col="red")

#Soft threshold and average connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
	 
abline(h=100,col="red")

dev.off()
	 
#The power value recommended by the system is 3
power = sft$powerEstimate
power

power = 3	#group A
power = 12	#group B
power = 14	#group C


#When the power of an undirected network is less than 15 or the power of a directed network is less than 30, there is no power value that can make
#The scale-free network graph structure R ^ 2 reaches 0.8, and the average connectivity is relatively high, such as above 100, which may be due to
#Some samples differ significantly from others. This may be influenced by batch effects, sample heterogeneity, or experimental conditions
#Due to the significant impact of expression. You can view grouping information and the presence of abnormal samples by drawing sample clusters.
#If this is indeed caused by meaningful biological changes, the following empirical power values can also be used.

#if (is.na(power)){
#  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
#          ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
#          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
#          ifelse(type == "unsigned", 6, 12))       
#          )
#          )
#}


#2.4 Network construction

#2.4.1 One Step Network Construction and Module Detection

#Power: Soft threshold calculated in the previous step
#MaxBlockSize: The maximum number of genes in a module that a computer can process (default is 5000);
#A 4G memory computer can handle 8000-10000 units, a 16GB memory computer can handle 20000 units, and a 32GB memory computer can handle 30000 units
#If computing resources permit, it is best to place it in a block.
# corType: pearson or bicor
#NumericLabels: returns a number instead of a color as the module name, which can then be converted to a color
#SaveTOMs: The most time-consuming calculation, stored for future use
#MergeCutHeight: The threshold for merging modules, the larger the threshold, the fewer modules

exprMat = ("newdata_filter_remove_pre_batheffect_removed_WGCNA")

net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.5,						#group A 0.15; group B 0.4; group C 0.5
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)
					   
##  Calculating module eigengenes block-wise from all genes
##    Flagging genes and samples with too many missing values...
##     ..step 1
##  ..Working on block 1 .
##     TOM calculation: adjacency..
##     ..will use 47 parallel threads.
##      Fraction of slow calculations: 0.000000
##     ..connectivity..
##     ..matrix multiplication (system BLAS)..
##     ..normalization..
##     ..done.
##    ..saving TOM for block 1 into file WGCNA/LiverFemaleClean.txt.tom-block.1.RData
##  ....clustering..
##  ....detecting modules..
##  ....calculating module eigengenes..
##  ....checking kME in modules..
##      ..removing 3 genes from module 1 because their KME is too low.
##      ..removing 5 genes from module 12 because their KME is too low.
##      ..removing 1 genes from module 14 because their KME is too low.
##  ..merging modules that are too close..
##      mergeCloseModules: Merging modules whose distance is less than 0.25
##        Calculating new MEs...

#Sort in descending order based on the number of genes in the module, numbered as' 1- Maximum number of modules'.
#* * 0 (grey) * * represents genes that have not been classified into any modules.  
table(net$colors)

##When merges CutHeight=0.25, the output is a module
##   0    1    2    3    4    5    6    7    8    9   10 
## 647 3549 2212 1160  817  358  214  143   83   45   36 

#2.4.2 Hierarchical clustering tree displays each module

#The gray ones represent genes that have not been classified into modules
# Convert labels to colors for plotting

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

#Plot the dendrogram and the module colors underneath
#If you are not satisfied with the results, you can also recall Blockwise Trees to save computation time

pdf("03 Moduledendrograms plot.pdf",width = 8,height = 4)   
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

#2.4.3 Draw correlation diagrams between modules

# module eigengene, A line graph can be drawn as a display of gene expression trends for each module
MEs = net$MEs

#No need to recalculate, just change the following names
#The official tutorial is recalculated, so there is no need to go through so much trouble at the beginning
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

#Correlation maps between modules obtained by clustering based on gene expression levels
#Set the bottom, left, top, and right margins for marDendro/marHeatmap

pdf("03 Eigengene adjacency heatmap.pdf",width = 6,height = 5)   

plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(0,3,2,7),
                      marHeatmap = c(3,4,2,7), plotDendrograms = T, 
                      xLabelsAngle = 90)

dev.off()		  

#If there is phenotype data, it can also be combined with ME data and plotted together

#MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
#plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
#                      marDendro = c(3,3,2,4),
#                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
#                      xLabelsAngle = 90)


#2.4.4 Phenotypic association analysis
		  
trait <- "01 WCGNA_GroupC_trait.csv"

#Read phenotype data
if(trait != "") {
  traitData <- t(read.csv(file=trait, header=T, row.names=1))  #转置一下row为samplenamte，col为traits
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}

#Module and phenotype data association
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

#Signif means to retain a few decimal places
pdf("04 Module and Trait heatmap.pdf",width = 16,height = 5)   

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.8, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               cex.lab.x = 1, cex.lab.y =1)
			   
dev.off()		  


#2.4.5 Export all genes and hub genes from a specific module

#Check the dimensions and row names to ensure that the rows are samples and the columns are genes

dim(dataExpr)
head(rownames(dataExpr)) 	#Check sample names
head(colnames(dataExpr)) 		#Check gene names

#Create objects corresponding to genes and modules
geneModuleAssignment <- data.frame(Gene = colnames(dataExpr), Module = moduleColors)

#Calculate the module membership relationships of genes in each module
geneModuleMembership = cor(dataExpr, MEs_col, use = "p")

##Extract genes of interest
#Select a module of interest, such as the blue module
interest_module = "blue"

modGenes = (geneModuleAssignment$Module == interest_module)  

table(modGenes)

#Retrieve Hub genes from the module
module = paste0("ME",interest_module)

ModuleCorrespondingGenes = names(dataExpr)[which(modGenes)]

hubGenes = names(dataExpr)[which(modGenes)][order(geneModuleMembership[modGenes, module], decreasing = TRUE)[1:10]]
print(hubGenes)


#2.4.6 Export all genes and hub genes from relative specific module

#Annotate genes
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))

a <- as.data.frame(table(geneModuleAssignment$Module),row.names = 1)
write.csv(a, "04 Module_gene_number.csv",row.names = T)

for (i in row.names(a))
{
interest_module = i
modGenes = (geneModuleAssignment$Module == interest_module) 
table(modGenes)

#########ModuleCorrespondingGenes
ModuleCorrespondingGenes = names(dataExpr)[which(modGenes)]

my_ensembl_gene_id<- ModuleCorrespondingGenes
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)

file_name <- paste('05 ME_',i,"ModuleCorrespondingGenes.csv")
write.csv(study_symbols,file_name, row.names = FALSE)
cat("File", file_name, "saved.\n")    ###reply save successfully

#########hubGenes
module = paste0("ME",interest_module)

hubGenes = names(dataExpr)[which(modGenes)][order(geneModuleMembership[modGenes, module], decreasing = TRUE)[1:10]]
head(hubGenes)

my_ensembl_gene_id<- ModuleCorrespondingGenes
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)

file_name <- paste('05 ME_',i,"ModulehubGenes.csv")
write.csv(study_symbols,file_name, row.names = FALSE)
cat("File", file_name, "saved.\n")    ###reply save successfully

}

