###多数据集整合分析DeSeq2

#这个一定要看！！！！！！！！！！#
https://blog.csdn.net/yayiling/article/details/113664434
https://www.plob.org/article/14417.html
##https://blingarida.github.io/2019/10/10/Differential%20expression%20analysis/

##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625



setwd("F:\\RNA_Sequencing_Data\\rnaseq\\my_rnaseq\\5.deg")                    #设置工作目录
##手动删除all_gene_count.tsv没有数据对应的行

raw_count <- read.table("all_gene_count.tsv",sep="\t",header=T,check.names=F) #改成自己的文件名
head(raw_count) 
rownames(raw_count)=raw_count[,1]
count_data=raw_count[,2:ncol(raw_count)]

raw_count_1 <- read.table("all_gene_count1.tsv",sep="\t",header=T,check.names=F) #改成自己的文件名
head(raw_count_1) 
rownames(raw_count_1)=raw_count_1[,1]
count_data_1=raw_count_1[,2:ncol(raw_count_1)]

mrna_names <- intersect(rownames(count_data),rownames(count_data_1))
expr <- cbind(count_data[mrna_names,], count_data_1[mrna_names,])
write.table(expr,file="expr.txt",sep="\t",quote=F)

#####检测Batch Effect 批次效应
# t() 转置函数
# dist() 距离函数：按照指定规则求行向量间的距离，因此要转置
dist_mat <- dist(t(expr))
clustering <- hclust(dist_mat) # hclust 的输入结构与 dist 相同！
# 按照batch批次信息聚类
plot(clustering)

##32,12,11分别是提取出的每个样本的数量，样本1命名为batch1，样本2命名为batch2，样本3命名为batch3
###examples：
batch <- paste0("batch",rep(c(1,2,3),c(32,12,11)))
batch
###cases for use：
batch <- paste0("batch",rep(c(1,2),c(3,3)))
batch

##25,7代表样本1中有25个实验组，7个对照组，以此类推   (combat的数据如果为两组则必须为同时存在control and case的情况，三组时，则没必要)
###examples：
tissue <- rep(c("case","control","case","control","case","control"),c(25,7,6,6,6,5))
tissue
###cases for use：
tissue <- rep(c("case","control","case","control"),c(0,3,3,0))
tissue
table(batch,tissue)
mod <- model.matrix(~tissue)
#用Combat()进行处理
library(sva)
expr_batch <- ComBat(dat=expr,batch=batch,mod=mod)
write.table(expr_batch,file="expr_batch.txt",sep="\t",quote=F)
##############
若出现以下报错，则可以删除mod=mod，使用expr_batch <- ComBat(dat=expr,batch=batch)
Found3batches
Adjusting for2covariate(s) or covariate level(s)
Error in ComBat(dat = expr, batch = batch, mod = mod) : 
  At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat

##使用SVA包中的ComBat来移除批次效应,注释信息为官方用法
# ComBat(dat, batch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
# dat: 基因组测量矩阵（探针维度 X 样本数），探针维度例如marker数、基因数.....，例如表达谱矩阵
# batch: 批次协变量，只能传入一个批次协变量！
# mod: 这是一个模式矩阵，里面包含了我们感兴趣的变量！
# par.prior: 基于参数/非参数，默认为基于参数

#####检测消除Batch Effect 批次效应之后的distance
# t() 转置函数
# dist() 距离函数：按照指定规则求行向量间的距离，因此要转置
dist_mat2 <- dist(t(expr_batch))
clustering2 <- hclust(dist_mat2) # hclust 的输入结构与 dist 相同！
# 按照batch批次信息聚类
plot(clustering2)

##heatmap比对前后差异
library(pheatmap)
library(RColorBrewer)
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

heatmap=pheatmap(expr,color = colors,
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
                 filename = "expr.pdf"
)

heatmap=pheatmap(expr_batch,color = colors,
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
                 filename = "expr_batch.pdf"
)

##
foldChange=1   
padj=0.05

#提取需要的列
expr_batch2 = expr[,1:7]
group=c("M6","M6","M6","Ras","Ras","Ras","Ras")
group=c("Scrib_Wts","Scrib_Wts","Scrib_Wts","E75_Scrib_Wts","E75_Scrib_Wts","E75_Scrib_Wts")
#batch信息和group信息
batch <- as.numeric(c("1","1","1","2","2","2"))
condition <- factor(c(rep("control",3),rep("treat",3)), levels = c("control","treat"))
col_data <- data.frame(row.names = colnames(expr_batch2), condition)
col_data$batch = batch

###DESeq2差异分析
library('DESeq2')
dds <- DESeqDataSetFromMatrix(countData = expr_batch2,colData = col_data, design = ~ batch + condition)
dds <- DESeqDataSetFromMatrix(countData = expr_batch2,colData = col_data, design = ~ condition)
nrow(dds)
dds_filter <- dds[ rowSums(counts(dds))>1, ]
dds_out <- DESeq(dds_filter)
res <- results(dds_out)
summary(res)
table(res$padj<0.05)
res_deseq <- res[order(res$padj),]
#一般选取Foldchange值和经过FDR矫正过后的p值，取padj值(p值经过多重校验校正后的值)小于0.05，log2FoldChange大于1的基因作为差异基因集
diff_gene_deseq2 <- subset(res_deseq, padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
res_diff_data <- merge(as.data.frame(res),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(res_diff_data,file = "alpaca_data_new.csv",row.names = F)


# 可以把前面生成的results.xls文件提交到www.ehbio.com/ImageGP绘制火山图
# 或者使用我们的s-plot https://github.com/Tong-Chen/s-plot

sampleA = "treat"
sampleB = "control"
contrastV <- c("conditions", sampleB, sampleA)
res <- results(dds,  contrast=contrastV)
res
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleA]

if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
} else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA
head(baseMeanA)

baseB <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleB]
if (is.vector(baseB)){
        baseMeanB <- as.data.frame(baseB)
} else {
        baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB
head(baseMeanB)

res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
head(res)

# 增加ID信息
res <- cbind(ID=rownames(res), as.data.frame(res))
res$baseMean <- rowMeans(cbind(baseA, baseB))
 
# 校正后p-value为NA的复制为1
res$padj[is.na(res$padj)] <- 1
 
# 按pvalue排序, 把差异大的基因放前面
res <- res[order(res$pvalue),]
head(res)

###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(mat)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(mat)
mat <- cbind(ensembl_gene_id, mat)
colnames(mat)[1]<-c("ensembl_gene_id")
diff_name <- merge(mat,study_symbols,by="ensembl_gene_id")


logCounts <- log2(res$baseMean+1)
logFC <- res$log2FoldChange
FDR <- res$padj
#png(filename=paste(file_base, "Volcano.png", sep="."))
plot(logFC, -1*log10(FDR), col=ifelse(FDR<=0.01, "red", "black"),
 xlab="logFC", ylab="-1*log1o(FDR)", main="Volcano plot", pch=".")
#dev.off()

##PCA出图
rld <- rlog(dds, blind = FALSE)
plotPCA(rld,intgroup=c("condition"))


library("genefilter")
library("pheatmap")
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),100)
mat  <- assay(rld)[ topVarGene, ]

###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(mat)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(mat)
mat <- cbind(ensembl_gene_id, mat)
colnames(mat)[1]<-c("ensembl_gene_id")
diff_name <- merge(mat,study_symbols,by="ensembl_gene_id")

##输出差异基因的热图
newData <- diff_name
row.names(newData) <- newData$external_gene_name
newData2 <- newData[,2:(nrow(col_data)+1)]
dimnames=list(rownames(newData2),colnames(newData2))
data=matrix(as.numeric(as.matrix(newData2)),nrow=nrow(newData2),dimnames=dimnames)

colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(data,color = colors,
                 main="TPS_TOP_diffSig",
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
                 filename = "TPS_TOP_diffSig.pdf"
)


pheatmap(mat, annotation_col=col_data$condition)
pheatmap(mat,cluster_row=T,scale="row", annotation_col=col_data) res0.5 <- results(dds, contrast = c("condition","Basal","LP"),alpha=0.05)
#另一种绘图方式
mat  <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=col_data)




