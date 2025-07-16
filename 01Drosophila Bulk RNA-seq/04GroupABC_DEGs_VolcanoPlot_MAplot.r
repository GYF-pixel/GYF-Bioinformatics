
#GroupABC_DEGs_VolcanoPlot_MAplot

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

#0. Arrange data
data = read.csv("newdata_filter_remove_pre_batheffect_removed.csv", header = T,row.names = 1)

#0.1 Group A
data.filter = as.data.frame(data)
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

colnames(data.filter)					#列名是对应组的肿瘤样品
write.csv(data.filter, "06 WCGNA_GroupA_Count.csv",row.names = T)


#0.2 Group B
data.filter = as.data.frame(data)
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
write.csv(data.filter, "06 WCGNA_GroupB_Count.csv",row.names = T)


#0.3 Group C
data.filter = as.data.frame(data)
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
write.csv(data.filter, "06 WCGNA_GroupC_Count.csv",row.names = T)

#1. DEGs identification by edgeR

#load packages
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

setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\09WGCNA_GroupABC")
list.files()


##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625
padj=0.05

#data = read.csv("newdata_filter.csv", header = T)
data.filter = read.csv("06 WCGNA_GroupA_Count.csv", header = T,row.names = 1)
data.filter = read.csv("06 WCGNA_GroupB_Count.csv", header = T,row.names = 1)
data.filter = read.csv("06 WCGNA_GroupC_Count.csv", header = T,row.names = 1)

#group information
combat_Expr = data.filter
colnames(combat_Expr)

group_info = c(rep("WT",3), rep("Ras",9), rep("5d_Tumors",18), rep("9d_Tumors",18), rep("13d_Tumors",18))			#group A
group_info = c(rep("WT",3), rep("Ras",9), rep("5d_Tumors",6), rep("9d_Tumors",6), rep("13d_Tumors",6))				#group B
group_info = c(rep("WT",3), rep("Ras",9), rep("5d_Tumors",6), rep("9d_Tumors",6), rep("13d_Tumors",6))				#group C


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



