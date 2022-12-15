
###1.edgeR分析差异基因DEGs
#source("http://bioconductor.org/biocLite.R")   #source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625  
padj=0.05

##去除重复和低表达的样品
setwd("F:\\RNA_Sequencing_Data\\rnaseq\\my_rnaseq\\5.deg")                    #设置工作目录
library("edgeR")
##手动删除all_gene_count.tsv没有数据对应的行
rt=read.table("all_gene_count.tsv",sep="\t",header=T,check.names=F)                  #改成自己的文件名
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

group=c("Ras","Scrib_Ras","Scrib_Ras","Scrib_Ras","Scrib_Wts","Scrib_Wts","Scrib_Wts")
group=c("WT","O53","O53","Wts_O53","O53_Wts","O53_Wts")
group=c("E75+DHR3","E75+DHR3","E75+DHR3","E75-IR","E75-IR","E75-IR","E75","E75","E75","WT","WT","WT")
group=c("Ras-lgl-GPX","Ras-lgl-GPX","Ras-lgl-GPX","Ras-lgl","Ras-lgl","Ras-lgl","Ras","Ras","Ras")
group=c("R-FB","R-FB","R-FB","RS-FB","RS-FB","RS-FB")
group=c("M6","M6","M6","Ras","Ras","Ras")
group=c("Ptp61f-lar","Ptp61f-lar","Tumor","Tumor","Ras","Ras","W1118","W1118")
group=c("Control","Control","Control","Vhl_polybromo_KD","Vhl_polybromo_KD","Vhl_polybromo_KD")
group=c("4-03","4-03","4-03","119-6","119-6","119-6","119-13","119-13","119-13","141-13","141-13","141-13","189-13","189-13","189-13")
group=c(rep("normal",4),rep("Sporadic_VS", 10),rep("NF2_associated_VS", 17))
group=c("79E","79E","79E","fry-mut","fry-mut","fry-mut")
group=c("119-6","119-6","119-6","119-10","119-10","119-10","119-13","119-13","119-13","189-13","189-13","189-13","141-13","141-13","141-13","189-Larv","189-Larv","189-Larv","4-03-Larv","4-03-Larv","4-03-Larv","4-03","4-03","4-03")
group=c("Q_S-Spz6_IR","Q_S-Spz6_IR","Q_S-Spz6_IR","Q_S","Q_S","Q_S")
group=c("QS11-FB","QS11-FB","QS11-FB","QS8-FB","QS8-FB","QS8-FB","R-FB","R-FB","R-FB","W-FB","W-FB","W-FB")
#group=c(rep("normal",1),rep("tumor",471))                         #按照癌症和正常样品数目修改
design <- model.matrix(~group)
y <- DGEList(counts=data, group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("W-FB","QS8-FB"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值（normalizeExp.txt）
normalizeExp1 <- normalizeExp[-1,]
#在Excel中删除为0的个体
normalizeExp_Non_Zero_Genes <- read.table("normalizeExp_Non_Zero_Genes.txt",sep="\t",header=T,check.names=F)

diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmRNAExp.txt）


###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(diff)
diff <- cbind(ensembl_gene_id,diff)
colnames(diff)[1]<-c("ensembl_gene_id")
diff_name <- merge(diff,study_symbols,by="ensembl_gene_id")

###对取交集的基因"interactions.txt"进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff)
Inter <- read.table("Overlapped vs R-FB.txt",sep="\t",header=F,check.names=F)
head(Inter)
rownames(Inter)<-Inter[,1]
my_ensembl_gene_id<-row.names(Inter)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(Inter)
Inter <- cbind(ensembl_gene_id,Inter)
colnames(Inter)[1]<-c("ensembl_gene_id")
inter_diff_name <- merge(Inter,study_symbols,by="ensembl_gene_id")
head(inter_diff_name )
write.table(inter_diff_name,file="inter_diff_name.txt",sep="\t",quote=F,col.names=F)   



#画火山图
library("ggplot2")
library("ggpubr")
library("ggthemes")
library("pheatmap")
#在画图之前需要将padj转换成-1*log10，这样可以拉开表达基因之间的差距   log10(0.05)=-1.30103 log2(1.5)=0.5849625
diff_name2 <- diff_name
rownames(diff_name2)=diff_name2[,1]
diff_name2$log10diff_namepadj <- -log10(diff_name2$FDR)
diff_name2$group <- "nonsignificance"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC > 0.5849625)]="up-regulated"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC < -0.5849625)]="down-regulated"
table(diff_name2$group)
#Label-FDR最高的20个基因(ensembl_gene_id和external_gene_name可以互换)
diff_name2$label=""
diff_name2 <- diff_name2[order(diff_name2$FDR),]
diff_name2_upgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="up-regulated")],20)
diff_name2_downgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="down-regulated")],20)
diff_name2_top20genes_FDR <- c(as.character(diff_name2_upgenes_FDR),as.character(diff_name2_downgenes_FDR))
diff_name2$label[match(diff_name2_top20genes_FDR,diff_name2$external_gene_name)] <- diff_name2_top20genes_FDR
#Label-logFC最高的10个基因
diff_name2$logFC_abs=abs(diff_name2$logFC)
diff_name2 <- diff_name2[order(diff_name2$logFC_abs),]
diff_name2_upgenes_logFC <- tail(diff_name2$external_gene_name[which(diff_name2$group=="up-regulated")],10)
diff_name2_downgenes_logFC <- head(diff_name2$external_gene_name[which(diff_name2$group=="down-regulated")],10)
diff_name2_top10genes_logFC <- c(as.character(diff_name2_upgenes_logFC),as.character(diff_name2_downgenes_logFC))
diff_name2$label[match(diff_name2_top10genes_logFC, diff_name2$external_gene_name)] <- diff_name2_top10genes_logFC
write.table(diff_name2,file="diff_name2.xls",sep="\t",quote=F,col.names=T)   
#画分界线（1.30是-log10(0.05)的值） log10(0.05)=-1.30103 log2(1.5)=0.5849625
ggscatter(data = diff_name2,x = "logFC",y = "log10diff_namepadj",
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

#画热图
library(pheatmap)
library(RColorBrewer)
#可以用的多种颜色
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

#所有基因的heatmap
TPS <- newData

heatmap=pheatmap(TPS,color = colors,
                 main="TPS",
                 fontsize = 15,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "TPS.pdf"
)

#diffSig基因的heatmap
##Gene_Ensemble_id
diffSig_rownames <- rownames(diffSig)
TPS_diffSig <- as.matrix(TPS[diffSig_rownames,])

##Gene_id
diffSig_rownames <- rownames(diffSig)
TPS_diffSig <- as.matrix(TPS[diffSig_rownames,])
TPS_diffSig2 <- diff_name2[match(rownames(TPS_diffSig),diff_name2$ensembl_gene_id),]
TPS_diffSig2 <- TPS_diffSig2[order(TPS_diffSig2$ensembl_gene_id),]
rownames(TPS_diffSig) <- TPS_diffSig2$external_gene_name  

##画图
heatmap=pheatmap(TPS_diffSig,color = colors,
                 main="TPS_diffSig",
                 fontsize = 15,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "TPS_diffSig.pdf"
)

#diffSig基因的top基因的heatmap
diff_name2 <- diff_name2[order(diff_name2$label),]
diffSig_top_genes_gene_id <- tail(diff_name2,60)
TPS_top_genes <- as.matrix(TPS[rownames(diffSig_top_genes_gene_id),])
#diffSig_top_genes_gene_id$label[c(1,2)]<- c("FBgn0032144","FBgn0036024")
rownames(TPS_top_genes) <- diffSig_top_genes_gene_id$external_gene_name
TPS_top_genes <-  TPS_top_genes[,c("QS11-1A","QS11-2A","QS11-3A","W-FB-1A","W-FB-2A","W-FB-3A")]

heatmap=pheatmap(TPS_top_genes,color = colors,
                 main="TPS_top_genes",
                 fontsize = 15,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "TPS_top_genes.pdf"
)

#需输出为txt根据热图结果手动删除不好的基因
Q <- row.names(TPS_top_genes)
write.table(Q,file="Q.txt",sep="\t",quote=F,col.names=F)
Q1 <- read.table("Q.txt",sep="\t",header=F,check.names=F)
TPS_top_genes_new <- as.matrix(TPS_top_genes[c(Q1[,2]),])

heatmap=pheatmap(TPS_top_genes_new,color = colors,
                 main="TPS_top_genes",
                 fontsize = 15,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "TPS_top_genes_new.pdf"
)




###差异表达基因功能富集分析
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(Cairo)
library(enrichplot)
#get the ENTREZID for the next analysis 以来于clusterProfiler包

sig.gene <- read.table("diffmRNAExp.txt",sep="\t",header=F,check.names=F)
head(sig.gene)
gene<-sig.gene[,1]
head(gene)
gene.df <- bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Dm.eg.db)
head(gene.df)
write.table(gene.df,file="gene_df.txt",sep="\t",quote=F,col.names=F)         #输出差异基因的三种id，这里会比sig.gene数目减少，是因为db里有些没有ensemble_id

#Go classification & #Go enrichment
##分别进行
#GO:CC(细胞组分） #keyType表示差异的基因的索引方式
ego_cc <- enrichGO(gene = gene.df$ENSEMBL,
                   OrgDb      = "org.Dm.eg.db",
                   keyType    = 'ENSEMBL',
                   ont        = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

barplot(ego_cc,showCategory = 20,title="The GO_CC enrichment analysis of up and down DEGs ")+ 
  scale_size(range=c(2, 12))+
  scale_y_discrete(labels=function(ego_cc) str_wrap(ego_cc,width = 25))
 
dotplot(ego_cc,showCategory = 20, title="The GO_CC enrichment analysis of DEGs", orderBy = "x")+
  scale_x_discrete(labels=function(ego_cc) str_wrap(ego_cc,width = 25))
 
#GO:BP(生物过程）
ego_bp <- enrichGO(gene = gene.df$ENSEMBL,
                 OrgDb      = "org.Dm.eg.db",
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

barplot(ego_bp,showCategory = 20,title="The GO_BP enrichment analysis of up and down DEGs")+ 
  scale_size(range=c(2, 12))+
  scale_y_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))
  
dotplot(ego_bp,showCategory = 20, title="The GO_BP enrichment analysis of DEGs")+
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))


#GO:MF(分子功能）  
ego_mf <- enrichGO(gene = gene.df$ENSEMBL,
                 OrgDb      = "org.Dm.eg.db",
                 keyType    = 'ENSEMBL',
                 ont        = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

barplot(ego_mf,showCategory = 20,title="The GO_MF enrichment analysis of up and down DEGs ")+ 
  scale_size(range=c(2, 10))+
  scale_y_discrete(labels=function(ego_mf) str_wrap(ego_mf,width = 25))
  
dotplot(ego_mf,showCategory = 20, title="The GO_MF enrichment analysis of DEGs")+
  scale_size(range=c(2, 10))+
  scale_x_discrete(labels=function(ego_mf) str_wrap(ego_mf,width = 20))                              


##新GO分析 
https://blog.csdn.net/weixin_33174392/article/details/112173971?utm_medium=distribute.pc_relevant_download.none-task-blog-2~default~blogcommendfrombaidu~default-3.nonecase&depth_1-utm_source=distribute.pc_relevant_download.none-task-blog-2~default~blogcommendfrombaidu~default-3.nonecas
##如下
ALL <- enrichGO(gene = gene.df$ENSEMBL, "org.Dm.eg.db", keyType = "ENSEMBL",ont = 'ALL',pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)  #一步到位
write.table(as.data.frame(ALL@result), file="GOALL.txt",sep="\t", quote=FALSE)
GOFile = read.csv("GOALL.txt", header = T, sep="\t")

dotplot(ALL, split = "ONTOLOGY", font.size = 8, showCategory = 10) + 
				facet_grid(ONTOLOGY ~ ., scale = "free") + 
				scale_y_discrete(labels = function(x) str_wrap(x, width =50)) + 
                scale_size(range=c(2, 6))   #设置点的大小
head(ALL,1);dim(ALL)				
sum(ALL$ONTOLOGY=="BP") #Biological process基因产物参与的生物路径或机制
sum(ALL$ONTOLOGY=="CC") #Cellular component基因产物在细胞内外的位置
sum(ALL$ONTOLOGY=="MF") #Molecular function基因产物分子层次的功能

#cnetplot
#对于基因和富集的GO terms之间的对应关系进行展示
#图中灰色的点代表基因，黄色的点代表富集到的GO terms
#如果一个基因位于一个GO Terms下，则将该基因与GO连线
#黄色节点的大小对应富集到的基因个数，默认画top5富集到的GO terms
library(ggnewscale)
cnetplot(ALL,showCategory=5)

FCgenelist <- diff_name2$logFC   #numeric vector
names(FCgenelist) <- as.character(diff_name2$external_gene_name) #named vector
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
head(FCgenelist)
cnetplot(ALL,showCategory=10,foldChange=FCgenelist,circular=TRUE,colorEdge=FALSE) 

##MF的cnetplot
ego_mf_new <- enrichGO(gene = gene.df$SYMBOL,
                 OrgDb      = "org.Dm.eg.db",
                 keyType    = 'ENSEMBL',
                 ont        = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

p1 <- cnetplot(ego_mf_new, showCategory=5,node_label="category", fontsize = 0.3)
p2 <- cnetplot(ego_mf_new, showCategory=5,node_label="gene", fontsize = 0.3) 
p3 <- cnetplot(ego_mf_new, showCategory=5,node_label="all", fontsize = 0.3) 
p4 <- cnetplot(ego_mf_new, showCategory=5,node_label="none", fontsize = 0.3) 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

#upsetplot
library(UpSetR)
library(DOSE)
upsetplot(ego_mf)

#heatplot
heatplot(ALL,foldChange=FCgenelist)

#有向无环图 GO DAG graph
goplot(ego_mf,showCategory=5)


#GO terms关系网络图 Enrichment Map
#对于富集到的GO terms之间的基因重叠关系进行展示
#每个节点是一个富集到的GO term，默认画top30个富集到的GO terms
#节点大小对应该GO terms下富集到的基因个数，节点的颜色对应p.adjust的值，红色小蓝色大
#如果两个GO terms的差异基因存在重叠，说明这两个节点存在overlap关系，用线条连接起来

emapplot(ALL,showCategory=10) 

##KEGG enrichment
https://www.jianshu.com/p/d484003dced5

##ID_Transformation
#install.packages("stringr")
library(stringr)
library(DOSE)
data(geneList, package="DOSE")
gene.kegg <- bitr_kegg(gene.df$ENTREZID,fromType="ncbi-geneid",
                        toType="kegg",organism='dme')
head(gene.kegg)

##在线富集分析
#keyType：one of "kegg",’ncbi-geneid’,’ncib-proteinid’ and ’uniprot’
#use_internal_data=FALSE：logical, use KEGG.db or latest online KEGG data
ekegg <- enrichKEGG(gene.df$ENTREZID, organism='dme',keyType="ncbi-geneid",
                    pvalueCutoff=0.5,pAdjustMethod='BH',
                    minGSSize=10,maxGSSize=500,use_internal_data=F)
					
write.table(as.data.frame(ekegg@result), file="KEGGALL.txt",sep="\t", quote=FALSE)
ekegg[1:30]
ekeggx <- setReadable(ekegg,'org.Dm.eg.db','ENTREZID')


barplot(ekegg,showCategory = 25, title="The KEGG enrichment analysis of up and down DEGs")+
  scale_size(range=c(2, 12))+
  scale_y_discrete(labels=function(ekegg) str_wrap(ekegg,width = 25))
dotplot(ekegg,showCategory = 35, title="The KEGG enrichment analysis of DEGs")+
  scale_size(range=c(2, 9))+
  scale_y_discrete(labels=function(ekegg) str_wrap(ekegg,width = 25))

##本地KEGG分析
##下载KEGG db
remotes::install_github("YuLab-SMU/createKEGGdb")
library(createKEGGdb)
create_kegg_db("dme")
install.packages("KEGG.db_1.0.tar.gz",repos=NULL,type="source")
library(KEGG.db)





##WGCNA
#BiocManager::install("WGCNA")
#BiocManager::install("reshape2")
library(reshape2)
library(stringr)
library(WGCNA)

library(WGCNA)
options(stringsAsFactors = FALSE)
# 指允许R语言程序最大线程运行
allowWGCNAThreads(12)
setwd("F:\\RNA_Sequencing_Data\\rnaseq\\my_rnaseq\\5.deg")                    #设置工作目录

femData = read.table("all_gene_count.tsv",sep="\t",header=T,check.names=F);
rownames(femData) <- femData[,1]
femData <- femData[,-1]
femData <- femData[,c(1:9,22:24)]
datExpr <- femData
##判断有无缺失值
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
##删除缺失值
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = hclust(dist(datExpr), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#dev.off()

traitData = read.csv("traits.txt",header =T, sep="\t")
rownames(traitData) <- traitData[,1]
traitData <- traitData[,-c(1,12)]
traitData <- traitData[colnames(datExpr),]
datTraits = traitData
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap",width = 12, height = 9)



##GSEA

DEGdata <- subset(diff_name2, diff_name2$logFC>foldChange | diff_name2$logFC<(-foldChange) )
DEgenelist <- as.character(DEGdata$external_gene_name)
head(DEgenelist)
#GSEA的基因排序列表
FCgenelist <- diff_name2$logFC   #numeric vector
names(FCgenelist) <- as.character(diff_name2$external_gene_name) #named vector
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
head(FCgenelist)

#MSigDb analysis
#Molecular Signatures Database包含了8种预定义的基因集合。
library(msigdbr)
msigdbr_show_species() #支持的物种
Dm_msigdbr <- msigdbr(species="Drosophila melanogaster")
head(Dm_msigdbr, 2) %>% as.data.frame
DmGO <- msigdbr(species="Drosophila melanogaster",category="C5") %>% 
                  dplyr::select(gs_name, entrez_gene, gene_symbol)
head(DmGO)
## 通用的富集分析
#TERM2GENE=gmt,TERM2NAME=NA 都是两列的数据框
#GENE是基因名,TERM表示GO term编号,NAME表示Description
em <- enricher(DEgenelist,TERM2GENE=DmGO[,c(1,3)])
head(em,1)
em1 <- GSEA(FCgenelist,TERM2GENE=DmGO[,c(1,3)])
head(em1,1)


gene.tx <- bitr(names(FCgenelist),fromType="SYMBOL",toType=c("ENTREZID"),
                OrgDb = org.Dm.eg.db)
colnames(gene.tx)[1] <- "external_gene_name"
gene.tx <- merge(gene.tx,diff_name2,by="external_gene_name")
FCgenelist <- diff_name2$logFC #numeric vector
names(FCgenelist) <- as.character(gene.tx$ENTREZID) #named vector
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order

egseKEGG <- gseKEGG(FCgenelist,organism='dme',keyType="ncbi-geneid",
                    nPerm=1000, minGSSize=10, maxGSSize=500,
                    pvalueCutoff=0.05, pAdjustMethod = "BH")
head(egseKEGG,1);dim(egseKEGG)





##Specific_Pathway_Heatmap
setwd("G:\\OneDrive - 西湖大学\\Ma Lab\\Analysis_for_lab\\Peng Liu\\LP_Analysis\\Tumor_vs_Ras")                    #设置工作目录
#KEGG list transformation
Origin_list <- read.table(file = "longevity.txt",sep="\t",header=F,check.names=F)   #格式为Ilp5的gene list
head(Origin_list)
###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id_list_count_new <- Origin_list[,1]
study_symbols_list_count_new <- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'external_gene_name', values = my_ensembl_gene_id_list_count_new, mart = mart)
head(study_symbols_list_count_new)
write.table(study_symbols_list_count_new,file="list_diffSig.txt",sep="\t",quote=F)
normalizeExp <- read.table(file="normalizeExp.txt",sep="\t",header=T,check.names=F) 
head(normalizeExp)
diffSig <- read.table(file="diffSig.xls",sep="\t",header=T,check.names=F)
head(diffSig)
rownames(normalizeExp) <- normalizeExp[,1]
normalizeExp <- normalizeExp[,-1]
list_count <- normalizeExp [study_symbols_list_count_new[,1],]
write.table(list_count,file="list_count.xls",sep="\t",quote=F) ##输出list_count后，删除NA的行，修改第一列为ID的列名
list_diffSig <- diffSig [study_symbols_list_count_new[,1],]
###重新读取
write.table(list_diffSig,file="list_diffSig.xls",sep="\t",quote=F) ##输出list_diffSig后，删除NA的行，修改第一列为ID的列名
list_count_new <- read.table(file="list_count.xls",sep="\t",header=T,check.names=F)
list_diffSig_new <- read.table(file="list_diffSig.xls",sep="\t",header=T,check.names=F)


#Flybase list transformation
library(xlsx)
Specific_list <- read.table(file = "Overlapped vs R-FB.txt",sep="\t",header=F,check.names=F)   #类似flybase的list， header=F
head(Specific_list)
normalizeExp <- read.table(file="normalizeExp.txt",sep="\t",header=T,check.names=F) 
head(normalizeExp)
diffSig <- read.table(file="diffSig.xls",sep="\t",header=T,check.names=F)
head(diffSig)
rownames(normalizeExp) <- normalizeExp[,1]
normalizeExp <- normalizeExp[,-1]
list_count <- normalizeExp [Specific_list[,1],]
write.table(list_count,file="list_count.xls",sep="\t",quote=F) ##输出list_count后，删除NA的行，修改第一列为ID的列名
list_diffSig <- diffSig [Specific_list[,1],]
write.table(list_diffSig,file="list_diffSig.xls",sep="\t",quote=F) ##输出list_diffSig后，删除NA的行，修改第一列为ID的列名
list_base <- diff_name2 [Specific_list[,1],]
write.table(list_base,file="list_base.xls",sep="\t",quote=F) ##输出list_base后，删除NA的行，修改第一列为ID的列名
###重新读取
list_count_new <- read.table(file="list_count.xls",sep="\t",header=T,check.names=F)
list_diffSig_new <- read.table(file="list_diffSig.xls",sep="\t",header=T,check.names=F)

###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id_list_count_new <- list_count_new[,1]
study_symbols_list_count_new <- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id_list_count_new, mart = mart)
head(study_symbols_list_count_new)
ensembl_gene_id_list_count_new <- list_count_new[,1]
diff <- cbind(ensembl_gene_id,diff)
colnames(list_count_new)[1] <- c("ensembl_gene_id")
list_count_new <- merge(list_count_new,study_symbols_list_count_new,by="ensembl_gene_id")

#输出list中差异显著的基因
list_count_new2 <- list_count_new
rownames(list_count_new2) <- list_count_new2$ensembl_gene_id
list_diffSig_new$external_gene_name <- list_count_new2[list_diffSig_new[,1],]$external_gene_name
write.table(list_diffSig_new, file="list_diffSig_new.xls",sep="\t",quote=F)

#画热图
library(pheatmap)
library(RColorBrewer)

##Gene_Ensemble_id整个pathway 的所有基因
pheatmap_list_count_new <- list_count_new
rownames(pheatmap_list_count_new) <- pheatmap_list_count_new$external_gene_name
colnames(pheatmap_list_count_new)  ##根据列名修改下面的索引结果
pheatmap_list_count_new <- pheatmap_list_count_new[,c("W-FB-1A","W-FB-2A","W-FB-3A","QS8-1A","QS8-2A","QS8-4A","QS11-1A","QS11-2A","QS11-3A")]

##Gene_Ensemble_id整个pathway 的差异基因
pheatmap_list_count_new2 <- pheatmap_list_count_new[list_diffSig_new$external_gene_name,]
colnames(pheatmap_list_count_new2)  ##根据列名修改下面的索引结果
pheatmap_list_count_new <- pheatmap_list_count_new[,c("W-FB-1A","W-FB-2A","W-FB-3A","QS8-1A","QS8-2A","QS8-4A","QS11-1A","QS11-2A","QS11-3A")]

##画图   #记得改热图的名字 选择pheatmap_list_count_new 或者pheatmap_list_count_new2 来决定是整个pathway的基因还是差异基因
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(pheatmap_list_count_new2, color = colors,
                 main="Overlapped vs R-FB",
                 fontsize = 10,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "Overlapped vs R-FB.pdf"
)
##比对热图删除不一致的基因
write.table(pheatmap_list_count_new, file="pheatmap_list_count_new.xls",sep="\t",quote=F) 
##输出pheatmap_list_count_new后，修改第一列为ID的列名
pheatmap_list_count_new2 <- read.table(file="pheatmap_list_count_new.xls",sep="\t",header=T,check.names=F)
rownames(pheatmap_list_count_new2) <- pheatmap_list_count_new2[,1]
pheatmap_list_count_new2 <- pheatmap_list_count_new2[,-1]
#重新画图 #记得改热图的名字
heatmap=pheatmap(pheatmap_list_count_new2, color = colors,
                 main="Hippo target gene",
                 fontsize = 10,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "pheatmap_list_count_new2.pdf"
)

#差异基因的相关性分析
library(psych)
library(corrplot)
#list_diffSig_new <- read.table(file = "list_diffSig_new.xls",sep="\t",header=T,check.names=F)      读取list_diffSig_new.xls修改correlation图
#pheatmap_list_count_new
list_diffSig_new_Correlation <- pheatmap_list_count_new[list_diffSig_new$external_gene_name,]
list_diffSig_new_Correlation <- t(list_diffSig_new_Correlation)
##method参数支持的相关性可选类型除了pearson外，还有spearman和kendall，并用来查看输出结果
list_diffSig_new_CorrelationResults <- corr.test(list_diffSig_new_Correlation, use = "complete", method = "spearman", adjust = "none")
list_diffSig_new_CorrelationResults$r
list_diffSig_new_CorrelationResults$ci  
res_cor <- cor(list_diffSig_new_Correlation)
corrplot(corr=res_cor)
corrplot(corr =res_cor,order = "AOE",type="upper",tl.pos = "d")
corrplot(corr =res_cor,add=TRUE, type="lower", method="number",order="AOE",diag=FALSE,tl.pos="n", cl.pos="n", title = "")





##同样根据最后热图的结果删掉list_diffSig_new.xls中均一性不好的基因，另存为CorelationAnalysis.xls用以展示和后续的相关性分析
#差异基因的相关性分析logFC
Correlation <- read.table(file="CorrelationAnalysis.xls",sep="\t",header=T,check.names=F)
rownames(Correlation) <- Correlation [,1]
Correlation  <- Correlation [,-1]
Corelation

