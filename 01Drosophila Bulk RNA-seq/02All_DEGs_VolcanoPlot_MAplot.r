
#load package
library(biomaRt)
library(curl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(scales)

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\04Analysis")
list.files()         

#1. DEGs identification by edgeR

#foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625
padj=0.05

#input files after batheffect removed
data.filter = read.csv("newdata_filter_remove_pre_batheffect_removed.csv", header = T,row.names = 1)
combat_Expr = data.filter

#group information
colnames(combat_Expr)
group_info = c(rep("WT",3), rep("Ras",9), rep("5d_Tumors",30), rep("9d_Tumors",30), rep("13d_Tumors",30))


#DGEList constructor
dge.list.obj <- DGEList(counts = combat_Expr, group = group_info)
dge.list.obj
# Normalization method: "TMM","TMMwsp","RLE","upperquartile","none"
dge.list.obj <- calcNormFactors(dge.list.obj,method = "RLE") # DESeq2, cuffdiff
dge.list.obj$samples
#plotMDS(dge.list.obj)
# make design matrix
design.mat <- model.matrix(~group_info)
# estimate dispersion
dge.list.obj <- estimateDisp(dge.list.obj,design.mat)
dge.list.obj$common.dispersion
dge.list.obj$tagwise.dispersion
	### Equal to the following steps
	#### 1st common dispersion
	#dge.list.obj <- estimateCommonDisp(dge.list.obj)
	#### 2nd tagwise dispersion
	#dge.list.obj <- estimateTagwiseDisp(dge.list.obj)
	#### plot dispersion
	#plotBCV(dge.list.obj, cex = 0.8)
	#### plot var and mean
	#plotMeanVar(dge.list.obj, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

# test with likelihood ratio test
##1.  GLM general linear model
	#fit <- glmFit(dge.list.obj, design.mat)
	#lrt <- glmLRT(fit, coef=2)
	#DEGs.res.lrt <- as.data.frame(topTags(lrt,n=nrow(count_df.filter),sort.by = "logFC"))
##2.  exactTest
#dge.list.res <- exactTest(dge.list.obj, pair = c("WT","Ras")) 
#dge.list.res <- exactTest(dge.list.obj, pair = c("Ras","5d_Tumors")) 
#dge.list.res <- exactTest(dge.list.obj, pair = c("5d_Tumors","9d_Tumors")) 
#dge.list.res <- exactTest(dge.list.obj, pair = c("9d_Tumors","13d_Tumors")) 
dge.list.res <- exactTest(dge.list.obj, pair = c("Ras","13d_Tumors")) #13d_Tumors vs Ras

#topTags Ranking of differential data according to criteria "logFC"
topTags(dge.list.res)
ordered_tags <- topTags(dge.list.res, n=1000000, sort.by = "logFC")
head(ordered_tags)
DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(combat_Expr),sort.by = "logFC"))



#2. Output related results 

#filter valid results
allDiff=DEGs.res[is.na(DEGs.res$FDR)==FALSE,]
diff=allDiff

#gene annotation for gene_symbol
mart <- useEnsembl(biomart = "ensembl", dataset = 'dmelanogaster_gene_ensembl',mirror = "uswest")         	#useEnsembl, mirror can be set as "uswest", "asia", etc.
my_ensembl_gene_id<-row.names(diff)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)

#Integrate gene annotation information with differential analysis results
ensembl_gene_id<-rownames(diff)
diff <- cbind(ensembl_gene_id,diff)
colnames(diff)[1]<-c("ensembl_gene_id")

#edgeR output results for all genes--diff_name
diff_name <- merge(diff,study_symbols,by="ensembl_gene_id")
write.table(diff_name,file="Part1_edgerOut.xls",sep="\t",quote=F)

#edgeR output results for all DEGs--diff_nameSig
diff_nameSig = diff_name[(diff_name$FDR < padj & (diff_name$logFC>foldChange | diff_name$logFC<(-foldChange))),]
write.table(diff_nameSig, file="Part1_diff_nameSig.xls",sep="\t",quote=F)
nrow(diff_nameSig)

#edgeR output results for all Upregulated DEGs--diff_nameUp
diff_nameUp = diff_name[(diff_name$FDR < padj & (diff_name$logFC>foldChange)),]
write.table(diff_nameUp, file="Part1_diff_nameup.xls",sep="\t",quote=F)
nrow(diff_nameUp)

#edgeR output results for all Downregulated DEGs--diff_nameDown
diff_nameDown = diff_name[(diff_name$FDR < padj & (diff_name$logFC<(-foldChange))),]
write.table(diff_nameDown, file="Part1_diff_namedown.xls",sep="\t",quote=F)
nrow(diff_nameDown)

#edgeR output results for normalized counts table of all genes--normalizeExp
normalizeExp=cpm(dge.list.obj)
write.table(normalizeExp,file="Part1_normalizeExp.xls",sep="\t",quote=F,col.names=T)  

#edgeR output results for normalized counts table of DEGs--diffExp
diffExp=normalizeExp[diff_nameSig$ensembl_gene_id,]
write.table(diffExp,file="Part1_diffmRNAExp.xls",sep="\t",quote=F,col.names=T)        
nrow(diffExp)



#3. Volcano Plot

#convert FDR to -1*log10 to enhance the differentiation between expressed genes.
diff_name2 <- diff_name
rownames(diff_name2)=diff_name2[,1]
diff_name2$log10diff_namepadj <- -log10(diff_name2$FDR)

#annotate up or down-regulation information 
diff_name2$group <- "nonsignificance"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC > 0.5849625)]="Up"			#-log10(FDR=0.05)=1.30103 log2(FC=1.5)=0.5849625
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC < -0.5849625)]="Down"		#-log10(FDR=0.05)=1.30103 log2(FC=1.5)=0.5849625
table(diff_name2$group)

#lable top 20 up or down-regulated genes based on FDR
diff_name2$label=""
diff_name2 <- diff_name2[order(diff_name2$FDR),]
diff_name2_upgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="Up")],20)
diff_name2_downgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="Down")],20)
diff_name2_top20genes_FDR <- c(as.character(diff_name2_upgenes_FDR),as.character(diff_name2_downgenes_FDR))
diff_name2$label[match(diff_name2_top20genes_FDR,diff_name2$external_gene_name)] <- diff_name2_top20genes_FDR

#lable top 10 DEGs based on logFC
diff_name2$logFC_abs=abs(diff_name2$logFC)
diff_name2 <- diff_name2[order(diff_name2$logFC_abs),]
diff_name2_upgenes_logFC <- tail(diff_name2$external_gene_name[which(diff_name2$group=="Up")],10)
diff_name2_downgenes_logFC <- head(diff_name2$external_gene_name[which(diff_name2$group=="Down")],10)
diff_name2_top10genes_logFC <- c(as.character(diff_name2_upgenes_logFC),as.character(diff_name2_downgenes_logFC))
diff_name2$label[match(diff_name2_top10genes_logFC, diff_name2$external_gene_name)] <- diff_name2_top10genes_logFC
write.table(diff_name2,file="Part1_diff_name2.xls",sep="\t",quote=F,col.names=T)   

#draw the volcano plot
volcano <- ggscatter(data = diff_name2,x = "logFC",y = "log10diff_namepadj",
	color = "group", 
	palette = c("#2f5688","#BBBBBB","#CC0000"),
	size = 1,
	font.label = c(8, "plain"),
	label = diff_name2$label,
	repel = T,
	xlab="Log2FoldChange",
	ylab="-Log10(FDR)",)+theme_base()+
	geom_hline(yintercept = 1.3,linetype="dashed")+ 
	geom_vline(xintercept = c(-0.5849625,0.5849625),linetype="dashed")
ggsave("Part2_Volcanoplot.pdf", plot = volcano, width = 11, height = 8) 



#4. MAplot

ggplot(diff_name2,aes(x=logCPM,y=logFC)) +
    geom_point(aes(color=group, alpha = group),size = 1.2,
               show.legend = T) +
    geom_hline(yintercept =  c(-0.5849625, 0.5849625),lty=2,lwd = 1) +
    theme_bw(base_size = 12) + 
    ggtitle("Neg vs Pos") + 
    scale_color_manual(values = c("#2f5688","#BBBBBB","#CC0000","darkorange")) +
    scale_size_manual(values = c(1.5,2,1.5,1.5)) +
    scale_alpha_manual(values = c(0.6,1,0.6,0.6)) +
    scale_x_log10() + 
    theme(panel.grid=element_blank(),
	panel.border = element_rect(color = "black", size = 2, fill = NA),  #border line
	axis.text.x = element_text(face = "plain", size = 12, colour = "black"),   #text  face = "italic" ("plain", "italic", "bold", "bold.italic")
	axis.text.y =element_text(face = "plain", size = 12, colour = "black")) +   
    scale_y_continuous(limits = c(-15,15),breaks = c(-10,-5,0,5,10)) +
    xlab("Mean Normalized Counts") + 
    ylab("Log2FoldChange") +
    geom_text_repel(aes(label=NA), color="black",fontface="italic",  #if label gene id, add "label=label"
                    size=4, segment.size=0.5,hjust=0,
                    nudge_x=1, nudge_y = 1)
ggsave(filename="Part2_MAplot.pdf", plot=last_plot(), scale = 1, width = 13, height = 8, units = c("cm"), dpi = 600,device = "pdf") 



