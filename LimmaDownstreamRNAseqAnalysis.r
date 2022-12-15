
library(GEOquery)
setwd("F:\\RNA_Sequencing_Data\\NF2\\PNAS")


####################################################################################################
##MicroarrayData导入数据和测序平台
geo2 <- getGEO(filename = "GSE108524_series_matrix.txt.gz", destdir = "GSE108524", getGPL = FALSE)
head(geo2@assayData[["exprs"]])
gpl1 <- getGEO(GEO = "GPL17586")
gpl1 <- read.table("GPL17586-45144.txt",sep="\t",header=T,check.names=F) )
#pData(): 返回与实验相关的表型数据或元数据，即phenotypic data的缩写，以data.frame的格式呈现
colnames(pData(geo2))
#看一下元数据的“title”列信息，告诉我们样本信息、处理信息、以及重复信息，factor格式数据：
pData(geo2)$title
#exprs(): Retrieve expression data from eSets.
head(exprs(geo2))
# expression matrix 表达矩阵
expr_cel <- exprs(geo2)
expr_cel <- as.data.frame(expr_cel)
colnames(expr_cel) ##查看样品名
#对得到的表达矩阵操作判断需不需要进行log2的处理
ex <- expr_cel
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

####注释probe
#https://www.jianshu.com/p/d89a25d43549

anno=gpl1@dataTable@table
colnames(anno)
anno[1:4,1:4]
probe2symbol=anno[,c(1,grep("gene_assignment",colnames(anno),ignore.case = T))]
head(probe2symbol)
#去掉未匹配到的情况
probe2symbol=probe2symbol[!(probe2symbol[,2] %in% c("",NA,"---")),]

##可以看到当一个探针匹配到多个基因时，是用///分隔；用//分隔同一基因的不同角度注释，其中Symbol ID一般位于第二条信息。需要用一些R语言技巧提取出来。
#如果出现一个探针组对应多个基因，使用“///”分割，可以选择保留第一个，删除其他的冗余信息
i <- sapply(probe2symbol, is.factor) # Change factor to character
probe2symbol[i] <- lapply(probe2symbol[i], as.character)
probe2symbol[,2] <- data.frame(sapply(probe2symbol[,2], function(x) unlist(strsplit(x,'///'))[1]),stringsAsFactors = F)[,1]
probe2symbol[,2] <- data.frame(sapply(probe2symbol[,2], function(x) unlist(strsplit(x,'//'))[2]),stringsAsFactors = F)[,1]
#去掉注释为"---"的结果
probe2symbol[,2] <- trimws(probe2symbol[,2])
probe2symbol=probe2symbol[!(probe2symbol[,2] %in% c("",NA,"---")),]
rownames(probe2symbol) <- probe2symbol[,1]

#添加symbol列
expr_cel$ID <- rownames(expr_cel)
sameprobe=intersect(expr_cel$ID,probe2symbol$ID )
expr_cel <- expr_cel[sameprobe,1:(ncol(expr_cel)-1)]
dim(expr_cel)


#如果多个探针对应一个基因 需要把重复的探针去掉，计算每一个探针在各个样本中平均值，留下最大的那个
expr_cel$gene_assignment <- probe2symbol[rownames(expr_cel),2]
rowMeans <- apply(expr_cel, 1, function(x) mean(as.numeric(x), na.rm = T))  #apply(m, 2, mean) 是将mean函数应用与矩阵 m 的列上，参数 2 指列（1 指行）
expr_cel <- expr_cel[order(rowMeans, decreasing =T),]
expr_cel <- expr_cel[!duplicated(expr_cel[, dim(expr_cel)[2]]),]  #expr_cel的最后一列为gene name
expr_cel <- na.omit(expr_cel)
dim(expr_cel)

rownames(expr_cel) <- expr_cel[,dim(expr_cel)[2]]
expr_cel <- expr_cel[,-(dim(expr_cel)[2])]
head(expr_cel)

### 保存数据，备下一步分析用
save(geo2, gpl1, expr_cel, file = "./geo2.Rdata")


###相关性分析
library(ggplot2)
#默认pearson
#p <-cor(expr_cel,method="pearson")
#spearman相关系数
p <-cor(expr_cel,method="spearman")
#kendall相关系数
#p <-cor(expr_cel,method="kendall")

#加cluster
head(p)
cor_heatmap <- heatmap(p,scale="none")
ggsave("Corelation.pdf", width = 6, height = 6)
#不想加cluster
#h<-heatmap(p,Rowv = NA, Colv = NA, scale="none")


#导入分组信息
group=c(rep("normal",4),rep("Sporadic_VS", 10),rep("NF2_associated_VS", 17))
names(group) <- colnames(expr_cel)
group <- as.data.frame(group)
head(group)

#PCA分析
expr_df <-as.data.frame(t(expr_cel))#转置
#R语言常用PCA分析函数有prcomp与princomp, 二者主要是计算方法的差别，建议采用prcomp(SVD方法)
pca.results <- prcomp(expr_df, center = TRUE, scale. = FALSE)
#定义足够多的颜色，用于展示分组
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

#画图
#install.packages("ggimage")
#install.packages("devtools")
#install.packages("qrcode")
#library(devtools) #未发布至CRAN 需要devtools安装
#devtools::install_github('fawda123/ggord')
#devtools::install_github("GuangchuangYu/yyplot")
#如果下载yyplot不顺利，可以直接调用文件夹里的geom_ord_ellipse.R
#source('./geom_ord_ellipse.R') 
#library(yyplot)

library(ggimage)
library(ggplot2)
library(plyr)
library(ggord)
library(yyplot)
#用ggord画基本PCA图
ggord(pca.results, grp_in = group$group, repel=TRUE,
      ellipse = FALSE, #不显示置信区间背景色
      size = 2, #样本的点大小
      alpha=0.5, #设置点为半透明，出现叠加的效果
      #如果用自定义的颜色，就运行下面这行
      cols = mycol[1:length(unique(group$group))],
      arrow = NULL,txt = NULL) + #不画箭头和箭头上的文字
  theme(panel.grid =element_blank()) + #去除网格线
  
  #用yyplot添加置信区间圆圈
  geom_ord_ellipse(ellipse_pro = .95, #设置置信区间
                   size=1.5, #线的粗细
                   lty=1 ) #实线

#保存到pdf文件
ggsave("PCA.pdf", width = 6, height = 6)


#DGE for microarray by limma
library('ggplot2')
library('limma')
####limma只能两两比较，所以这里重新分组并设置group_information
View(group)
exprSet <- expr_cel[,c(1:14)]
group_information <- as.data.frame(group[colnames(exprSet),])
colnames(group_information) <- c("group")
design <- model.matrix(~0+factor(group_information$group))#把group设置成一个model matrix#
colnames(design)=levels(factor(group_information$group))
rownames(design)=colnames(exprSet)
design

fit <- lmFit(exprSet,design)##线性拟合
##根据分组信息修改下放的"normal"    "Sporadic_VS"   "NF2_associated_VS"等
cont.matrix<-makeContrasts(Sporadic_VS_VS-normal,levels = design)
fit2=contrasts.fit(fit,cont.matrix)##用对比模型进行差值计算
fit2 <- eBayes(fit2)  ##贝叶斯检验
##eBayes() with trend=TRUE
tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH",sort.by="B",resort.by="M") 
nrDEG = na.omit(tempOutput)

write.csv(nrDEG, "limmaOut.csv")

#筛选有显著差异的基因  adj.P.Val意思是矫正后P值
foldChange=0.5849625 #fold change=1意思是差异是两倍log2(1.5)=0.5849625
pvalue =0.01 #0.05

diff <- nrDEG
diffSig = diff[(diff$P.Val < pvalue & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
#write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
write.csv(diffSig, "diffSig.csv")

#把上调和下调分别输入up和down两个文件
diffUp = diff[(diff$P.Val < pvalue & (diff$logFC>foldChange)),]#foldchange>0是上调，foldchange<0是下调#
#write.table(diffUp, file="up.xls",sep="\t",quote=F)#把上调和下调分别输入up和down两个文件#
write.csv(diffUp, "diffUp.csv")
diffDown = diff[(diff$P.Val< pvalue & (diff$logFC<(-foldChange))),]
#write.table(diffDown, file="down.xls",sep="\t",quote=F)
write.csv(diffDown, "diffDown.csv")


##############################################################################################################
#Limma处理RNAseq Data
library('ggplot2')
library('limma')

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
expr <- count_data
write.table(expr,file="expr.txt",sep="\t",quote=F)
#对得到的表达矩阵操作判断需不需要进行log2的处理 
#接着按照文档的说明以及limma包的习惯，我们需要对count进行标准化以及转化为log2的值，这里标准化的方法为TMM，使用edgeR里面的calcNormFactors函数即可
library(edgeR)
dge <- DGEList(counts = expr)
dim(dge)
#Filtering  Keep genes with total counts more than 50
 A <- rowSums(dge$counts)
 isexpr <- A > 50
dge <- dge[isexpr, ]
dim(dge)
dge <- calcNormFactors(dge)  #Apply TMM normalization
dge 
plotMDS(dge)


###当不同批次的数据中同时存在normal和control的时候，用这个方法
batch <- factor(c(1,1,1,1,1,2,2,2,2))
group=factor(c(rep("Scrib_Wts",3),rep("Ras_Scrib",2),rep("Scrib_Wts",2),rep("Ras_Scrib",2)), levels=c("Scrib_Wts","Ras_Scrib"))
design <- model.matrix(~batch+group)
v <- voom(dge,design,plot=TRUE) #Apply voom to convert the read counts to log2-cpm with associated weights:
fit <- lmFit(v, design)
fit.de <- eBayes(fit, robust=TRUE)
topTable(fit.de, coef="Ras_Scrib")

###当不同批次的数据中并没有同时存在normal和control的时候，用这个方法
#整体的batch信息,以下两种方法都可以
#batch <- paste0("batch",rep(c(1,2,3),c(3,3,3)))
batch <- factor(c(1,1,1,1,1,1,2,2,2))
group=factor(c(rep("Scrib_Wts",3),rep("Ras_Scrib",3),rep("E75_Scrib_Wts",3)), levels=c("Scrib_Wts","Ras_Scrib","E75_Scrib_Wts"))
batch <- factor(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2))
group=factor(c(rep("82B_eye",3),rep("82B_fat",3),rep("Scrib_eye",3),rep("Scrib_fat",3),rep("TL6_Scrib_eye",3)), levels=c("82B_eye","82B_fat","Scrib_eye","Scrib_fat","TL6_Scrib_eye"))

design <- model.matrix(~0+group)
v <- voom(dge,design,plot=TRUE) #Apply voom to convert the read counts to log2-cpm with associated weights:

#####检测Batch Effect 批次效应
ex_b_limma <-  v[["E"]]
# t() 转置函数
# dist() 距离函数：按照指定规则求行向量间的距离，因此要转置
#dist_mat <- dist(t(expr))
#clustering <- hclust(dist_mat) # hclust 的输入结构与 dist 相同！
# 按照batch批次信息聚类
#plot(clustering)

## 使用 limma 的 removeBatchEffect 函数,design 并非必要,design must be a numeric matrix
##v[["E"]]是voom对象v中已经取了log2之后的表达矩阵
#ex_b_limma <- removeBatchEffect(v[["E"]], batch = batch)
#clustering2 <- hclust(dist(t(ex_b_limma)))
#plot(clustering2)
#write.table(ex_b_limma,file="ex_b_limma.txt",sep="\t",quote=F)
#ex_b_limma <- removeBatchEffect(expr,batch = batch, design=design)


#DGE by limma
library('ggplot2')
library('limma')

####limma两两比较，所以这里重新分组并设置group_information
View(group)
exprSet <- ex_b_limma[,c(1:3,7:9)]
group=factor(c(rep("82B_eye",3),rep("82B_fat",3),rep("Scrib_eye",3),rep("Scrib_fat",3),rep("TL6_Scrib_eye",3)), levels=c("82B_eye","82B_fat","Scrib_eye","Scrib_fat","TL6_Scrib_eye"))
group_information <- as.data.frame(group)
rownames(group_information)  <- colnames(exprSet)
group_information
design2 <- model.matrix(~0+factor(group_information$group))#把group设置成一个model matrix#
colnames(design2)=levels(factor(group_information$group))
rownames(design2)=colnames(exprSet)
design2
fit <- lmFit(exprSet,design2)##线性拟合
##根据分组信息修改下放的"normal"    "Sporadic_VS"   "NF2_associated_VS"等
##cont.matrix<-makeContrasts(Sporadic_VS-normal,levels = design2)
##E75_Scrib_Wts 比 Scrib_Wts为E75_Scrib_Wts-Scrib_Wts

cont.matrix<-makeContrasts(E75_Scrib_Wts-Scrib_Wts,levels = design2)
fit2=contrasts.fit(fit,cont.matrix)##用对比模型进行差值计算
fit2 <- eBayes(fit2)  ##贝叶斯检验
##eBayes() with trend=TRUE
tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH",sort.by="B",resort.by="M") 
nrDEG = na.omit(tempOutput)

write.csv(nrDEG, "limmaOut.csv")


######################################################
#多组同时比较：   （一定要看！！！！）
# > contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
# > fit2 <- contrasts.fit(fit, contrast.matrix)
# > fit2 <- eBayes(fit2)
# A list of top genes differential expressed in group2 versus group1 can be obtained from
# > topTable(fit2, coef=1, adjust="BH")
######################################################
exprSet <- ex_b_limma[,c(1:9)]
group=factor(c(rep("Scrib_Wts",3),rep("Ras_Scrib",3),rep("E75_Scrib_Wts",3)), levels=c("Scrib_Wts","Ras_Scrib","E75_Scrib_Wts"))

exprSet <- ex_b_limma[,c(1:15)]
group=factor(c(rep("E82B_eye",3),rep("E82B_fat",3),rep("Scrib_eye",3),rep("Scrib_fat",3),rep("TL6_Scrib_eye",3)), levels=c("E82B_eye","E82B_fat","Scrib_eye","Scrib_fat","TL6_Scrib_eye"))
group_information <- as.data.frame(group)
rownames(group_information)  <- colnames(exprSet)
group_information  ###查看分组是否正确

design2 <- model.matrix(~0+factor(group_information$group))#把group设置成一个model matrix#
colnames(design2)=levels(factor(group_information$group))
rownames(design2)=colnames(exprSet)
design2   ###查看experiment design是否正确

fit <- lmFit(exprSet,design2)##线性拟合
cont.matrix<-makeContrasts(Ras_Scrib-Scrib_Wts, E75_Scrib_Wts-Scrib_Wts, E75_Scrib_Wts-Ras_Scrib,levels = design2)
cont.matrix<-makeContrasts(Scrib_fat-E82B_fat, Scrib_eye-E82B_eye, TL6_Scrib_eye-Scrib_eye, TL6_Scrib_eye-E82B_eye,levels = design2)
fit2=contrasts.fit(fit,cont.matrix)##用对比模型进行差值计算
fit2 <- eBayes(fit2)  ##贝叶斯检验

##eBayes() with trend=TRUE  
##使用coef=1来制定查看第几组的结果如 coef=1 为E75_Scrib_Wts-Scrib_Wts，coef=2 为E75_Scrib_Wts-Ras_Scrib，并通过设置cols_selected来标记

tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH",sort.by="B",resort.by="M") 
nrDEG = na.omit(tempOutput)
write.csv(nrDEG, "limmaOut.csv")

tempOutput = topTable(fit2,coef=2,n=Inf,adjust="BH",sort.by="B",resort.by="M") 
nrDEG = na.omit(tempOutput)
write.csv(nrDEG, "limmaOut.csv")

tempOutput = topTable(fit2,coef=3,n=Inf,adjust="BH",sort.by="B",resort.by="M") 
nrDEG = na.omit(tempOutput)
write.csv(nrDEG, "limmaOut.csv")

tempOutput = topTable(fit2,coef=4,n=Inf,adjust="BH",sort.by="B",resort.by="M") 
nrDEG = na.omit(tempOutput)
write.csv(nrDEG, "limmaOut.csv")

#筛选有显著差异的基因  adj.P.Val意思是矫正后P值
foldChange=0.5849625 #fold change=1意思是差异是两倍log2(1.5)=0.5849625
pvalue =0.01 #0.05
FDR =0.05

###对基因进行注释-获取gene_symbol     
library("biomaRt")
library("curl")
diff <- nrDEG
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(diff)
diff <- cbind(ensembl_gene_id,diff)
colnames(diff)[1]<-c("ensembl_gene_id")
diff_name <- merge(diff,study_symbols,by="ensembl_gene_id")
#删除FBti的所有信息
diff_name <- diff_name[order(diff_name$ensembl_gene_id),]
diff_name[c(12618:12622),]   #显示FBti在哪里
diff_name <- diff_name[c(1:9941),]   


diff <- diff_name
rownames(diff) <- diff$external_gene_name
write.csv(diff, "diff_name.csv")
diffSig = diff[(diff$adj.P.Val < FDR & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
write.csv(diffSig, "diffSig.csv")
#把上调和下调分别输入up和down两个文件
diffUp = diff[(diff$adj.P.Val < FDR & (diff$logFC>foldChange)),]#foldchange>0是上调，foldchange<0是下调#
#write.table(diffUp, file="up.xls",sep="\t",quote=F)#把上调和下调分别输入up和down两个文件#
write.csv(diffUp, "diffUp.csv")
diffDown = diff[(diff$adj.P.Val < FDR & (diff$logFC<(-foldChange))),]
#write.table(diffDown, file="down.xls",sep="\t",quote=F)
write.csv(diffDown, "diffDown.csv")

#画火山图
library("ggplot2")
library("ggpubr")
library("ggthemes")
library("pheatmap")
#在画图之前需要将padj转换成-1*log10，这样可以拉开表达基因之间的差距   log10(0.05)=-1.30103 log2(1.5)=0.5849625
diff_name2 <- diff
diff_name2$log10diff_namepadj <- -log10(diff_name2$adj.P.Val)
diff_name2$external_gene_name <- rownames(diff_name2)
diff_name2$group <- "nonsignificance"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC > 0.5849625)]="up-regulated"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC < -0.5849625)]="down-regulated"
table(diff_name2$group)
#Label-FDR最高的20个基因(ensembl_gene_id和external_gene_name可以互换)
diff_name2$label=""
diff_name2 <- diff_name2[order(diff_name2$adj.P.Val),]
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
#画分界线（1.30是-log10(0.05)的值） log10(0.05)=-1.30103 log2(1.5)=0.5849625 为横轴纵轴线xlim=c(-40,40),	ylim=c(-0,8)
ggscatter(data = diff_name2,x = "logFC",y = "log10diff_namepadj",
	color = "group", 
	palette = c("#2f5688","#BBBBBB","#CC0000"),
	size = 1,
	font.label = c(8, "plain"),
	label = diff_name2$label,
	repel = T,
	xlab="Log2FoldChange",
	ylab="-Log10(FDR)",
)+theme_base()+
	geom_hline(yintercept = 1.3,linetype="dashed")+ 
	geom_vline(xintercept = c(-0.5849625,0.5849625),linetype="dashed")


##plotMA 用DESeq2画 原文链接：https://blog.csdn.net/leo12354/article/details/107689656
library('DESeq2')
library('apeglm')
mycounts <- expr[,c(1:3,7:9)]
condition <- factor(c(rep("control",3),rep("treat",3)), levels = c("control","treat"))
condition

colData <- data.frame(row.names=colnames(mycounts), condition)
colData

dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds

res = results(dds, contrast=c("condition", "control", "treat"))
res = res[order(res$pvalue),]
summary(res)

plotMA(res,ylim=c(-2,2))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=6, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

res_order<-res[order(row.names(res)),]
res = res_order
res.shrink <- lfcShrink(dds,contrast = c("condition","treat","control"), res=res)
#但结果提示Error in lfcShrink(dds, contrast = c("condition", "treat", "control"),  :   type='apeglm' shrinkage only for use with 'coef'，目前不确定原因，故采取了下边两种方式
#1. 
resApe <-lfcShrink(dds, coef=2,type="apeglm")
#2. 
resLFC <- lfcShrink(dds, coef=2)
plotMA(resApe, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
#DESeq2提供了一个plotCounts()函数来查看某一个感兴趣的gene在组间的差别。counts会根据groups分组。
d1=plotCounts(dds, gene="FBgn0000568", intgroup="condition", returnData=TRUE) 
ggplot(d1,aes(condition, count)) + geom_boxplot(aes(fill=condition)) + scale_y_log10() + ggtitle("Eip75B")


###差异表达基因功能富集分析
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Dm.eg.db)
library(ggplot2)
library(stringr)
library(AnnotationDbi)
library(Cairo)
library(enrichplot)

#Drosophila的Go
sig.gene <- diffSig$external_gene_name
gene.df <- bitr(sig.gene, fromType = "SYMBOL", 
              toType = c("ENSEMBL","ENTREZID"),
              OrgDb = org.Dm.eg.db)
head(gene.df)
write.table(gene.df,file="gene_df.txt",sep="\t",quote=F,col.names=F)         #输出差异基因的三种id，这里会比sig.gene数目减少，是因为db里有些没有ensemble_id

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





#人的GO
sig.gene <- diff_name2$external_gene_name
gene.df <- bitr(sig.gene, fromType = "SYMBOL", 
              toType = c("ENSEMBL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(gene.df)
write.table(gene.df,file="gene_df.txt",sep="\t",quote=F,col.names=F)         #输出差异基因的三种id，这里会比sig.gene数目减少，是因为db里有些没有ensemble_id

ALL <- enrichGO(gene = gene.df$ENSEMBL, OrgDb      = "org.Hs.eg.db", keyType = "ENSEMBL",ont = 'ALL',pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)  #一步到位
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

#GO:CC(细胞组分） #keyType表示差异的基因的索引方式
ego_cc <- enrichGO(gene = gene.df$ENSEMBL,
                   OrgDb      = "org.Hs.eg.db",
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
                 OrgDb      = "org.Hs.eg.db",
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

barplot(ego_bp,showCategory = 20, fontsize=8, title="The GO_BP enrichment analysis of up and down DEGs")+ 
  scale_size(range=c(5, 12))+
  scale_y_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 55))
  
dotplot(ego_bp,showCategory = 20, title="The GO_BP enrichment analysis of DEGs")+
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))


#GO:MF(分子功能）  
ego_mf <- enrichGO(gene = gene.df$ENSEMBL,
                 OrgDb      = "org.Hs.eg.db",
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

##KEGG
##在线富集分析
#keyType：one of "kegg",’ncbi-geneid’,’ncib-proteinid’ and ’uniprot’
#use_internal_data=FALSE：logical, use KEGG.db or latest online KEGG data
#organism='dme'   'hsa' 
ekegg <- enrichKEGG(gene.df$ENTREZID, organism='dme',keyType="ncbi-geneid",
                    pvalueCutoff=0.5,pAdjustMethod='BH',qvalueCutoff=0.2,
                    minGSSize=10,maxGSSize=500,use_internal_data=F)
ekegg[1:30]
write.table(as.data.frame(ekegg@result), file="KEGGALL.txt",sep="\t", quote=FALSE)
ekeggx <- setReadable(ekegg,'org.Dm.eg.db','ENTREZID')

barplot(ekegg,showCategory = 30, title="The KEGG enrichment analysis of up and down DEGs")+
  scale_size(range=c(2, 12))+
  scale_y_discrete(labels=function(ekegg) str_wrap(ekegg,width = 30))
  
dotplot(ekegg,showCategory = 25, title="The KEGG enrichment analysis of DEGs")+
  scale_size(range=c(2, 9))+
  scale_y_discrete(labels=function(ekegg) str_wrap(ekegg,width = 35))
  
  
  

#diffSig基因的heatmap
TPS <- ex_b_limma
colnames(TPS)
##Gene_Ensemble_id  to  Gene_id
#cols_selected <- c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")
#cols_selected <- c("Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")
#cols_selected <- c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B")
cols_selected <- c("82B-eye-1A","82B-eye-2A","82B-eye-3A","TL6-scrib-1A","TL6-scrib-2A","TL6-scrib-3A")
diffSig_rownames <- diffSig$ensembl_gene_id
TPS_diffSig <- as.data.frame(TPS[diffSig_rownames,])
rownames(TPS_diffSig) <- diffSig$external_gene_name 
TPS_diffSig <- TPS_diffSig[,cols_selected]

##画图TPS_diffSig
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(TPS_diffSig,color = colors,
                 main="TPS_diffSig",
                 fontsize = 15,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "TPS_diffSig.pdf"
)
##画图TPS_diffUp
diffUp_rownames <- diffUp$ensembl_gene_id
TPS_diffUp <- as.data.frame(TPS[diffUp_rownames,])
rownames(TPS_diffUp) <- diffUp$external_gene_name 
TPS_diffUp <- TPS_diffUp[,cols_selected]

colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(TPS_diffUp,color = colors,
                 main="TPS_diffUp",
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
                 filename = "TPS_diffUp.pdf"
)
##画图TPS_diffDown
diffDown_rownames <- diffDown$ensembl_gene_id
TPS_diffDown <- as.data.frame(TPS[diffDown_rownames,])
rownames(TPS_diffDown) <- diffDown$external_gene_name 
TPS_diffDown <- TPS_diffDown[,cols_selected]

colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(TPS_diffDown,color = colors,
                 main="TPS_diffDown",
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
                 filename = "TPS_diffDown.pdf"
)
  
#Flybase list transformation list 中差异显著的基因
library(xlsx)
Specific_list <- read.table(file = "FlyBase_IDs_TNFα-Eiger Signaling Pathway.txt",sep="\t",header=F,check.names=F)   #类似flybase的list， header=F
head(Specific_list)
diffSig <- read.table(file="diffSig.xls",sep="\t",quote = "",header=T,check.names=F)
head(diffSig)
diffSig_name <- intersect(rownames(ex_b_limma),diffSig$ensembl_gene_id)
normalizeExp <- ex_b_limma[diffSig_name,]
normalizeExp <- normalizeExp[,colnames(exprSet)]
normalizeExp <- as.data.frame(normalizeExp)
list_count <- normalizeExp[Specific_list[,1],]
list_count <- na.omit(list_count)
write.table(list_count,file="list_count.xls",sep="\t",quote=F) ##输出list_count后，删除NA的行，修改第一列为ID的列名
diff_in_list <- intersect(diffSig$ensembl_gene_id,Specific_list[,1]) 
rownames(diffSig) <- diffSig$ensembl_gene_id
list_diffSig <- diffSig[diff_in_list,]
write.table(list_diffSig,file="list_diffSig.xls",sep="\t",quote=F) ##输出list_diffSig后，删除NA的行，修改第一列为ID的列名

#画热图
library(pheatmap)
library(RColorBrewer)

##Gene_Ensemble_id
pheatmap_list_count_new <- list_count
pheatmap_list_count_new2 <- list_diffSig
rownames(pheatmap_list_count_new) <- pheatmap_list_count_new2[rownames(pheatmap_list_count_new),]$external_gene_name
colnames(pheatmap_list_count_new)  ##根据列名修改下面的索引结果
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
pheatmap_list_count_new <- pheatmap_list_count_new[,c("82B-Fat-1A","82B-Fat-2A","82B-Fat-3A","scrib-Fat-1A","scrib-Fat-2A","scrib-Fat-3A")]
##画图   #记得改热图的名字
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(pheatmap_list_count_new, color = colors,
                 main="TNFα-Eiger Signaling",
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
                 filename = "TNFα-Eiger Signaling.pdf"
)
  



###Symbol List
library(xlsx)
normalizeExp <- ex_b_limma[,colnames(exprSet)]
normalizeExp <- as.data.frame(normalizeExp)

Specific_list <- read.table(file = "signaling.txt",sep="\t",header=F,check.names=F)   #类似flybase的list， header=F
head(Specific_list)
#diffSig = diff[(diff$adj.P.Val < FDR & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
diffSig <- read.table(file="diffSig.xls",sep="\t",quote = "",header=T,check.names=F)
head(diffSig)
diffSig_nameID <- diffSig[Specific_list[,1],]$external_gene_name
diffSig_nameFly <- diffSig[Specific_list[,1],]$ensembl_gene_id
list_count <- normalizeExp[diffSig_nameFly,]
list_count <- na.omit(list_count)
write.table(list_count,file="list_count.xls",sep="\t",quote=F) ##输出list_count后，删除NA的行，修改第一列为ID的列名
list_diffSig <- diffSig[diffSig_nameID,]
write.table(list_diffSig,file="list_diffSig.xls",sep="\t",quote=F) ##输出list_diffSig后，删除NA的行，修改第一列为ID的列名

#画热图
library(pheatmap)
library(RColorBrewer)

##Gene_Ensemble_id
diffSig2 <- diffSig
rownames(diffSig2) <- diffSig$ensembl_gene_id
pheatmap_list_count_new <- list_count
rownames(pheatmap_list_count_new) <- diffSig2[rownames(list_count),]$external_gene_name
colnames(pheatmap_list_count_new)  ##根据列名修改下面的索引结果
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B")]
pheatmap_list_count_new <- pheatmap_list_count_new[,c("82B-eye-1A","82B-eye-2A","82B-eye-3A","TL6-scrib-1A","TL6-scrib-2A","TL6-scrib-3A")]
##画图   #记得改热图的名字
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(pheatmap_list_count_new, color = colors,
                 main="KEGG signaling",
                 fontsize = 10,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = F,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "KEGG signaling.pdf"
)
  


#GradientChangedGenes--Heatmap 
###FlybaseID
library(xlsx)
normalizeExp <- ex_b_limma[,colnames(exprSet)]
normalizeExp <- as.data.frame(normalizeExp)

Specific_list <- read.table(file = "GradientChangedGenes.txt",sep="\t",header=F,check.names=F)   #类似flybase的list， header=F
head(Specific_list)
#diffSig = diff[(diff$adj.P.Val < FDR & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-Specific_list[,1]
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-Specific_list[,1]
colnames(Specific_list)[1]<-c("ensembl_gene_id")
Specific_list_names <- study_symbols


diffSig_nameID <- Specific_list_names$external_gene_name
diffSig_nameFly <- Specific_list_names$ensembl_gene_id
list_count <- normalizeExp[diffSig_nameFly,]
list_count <- na.omit(list_count)
rownames(Specific_list_names) <- Specific_list_names[,1]
rownames(list_count) <- Specific_list_names[rownames(list_count),]$external_gene_name
write.table(list_count,file="list_count.xls",sep="\t",quote=F) ##输出list_count后，删除NA的行，修改第一列为ID的列名

#画热图
library(pheatmap)
library(RColorBrewer)

##Gene_Ensemble_id
pheatmap_list_count_new <- list_count
colnames(pheatmap_list_count_new)  ##根据列名修改下面的索引结果
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
pheatmap_list_count_new <- pheatmap_list_count_new[,c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
pheatmap_list_count_new <- pheatmap_list_count_new[,c("82B-Fat-1A","82B-Fat-2A","82B-Fat-3A","scrib-Fat-1A","scrib-Fat-2A","scrib-Fat-3A")]
##画图   #记得改热图的名字
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(pheatmap_list_count_new, color = colors,
                 main="Ion Channel Activity",
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
                 filename = "Ion Channel Activity.pdf"
)
  
  

###GeneSymbolID
library(xlsx)
normalizeExp <- ex_b_limma[,colnames(exprSet)]
normalizeExp <- as.data.frame(normalizeExp)

Specific_list <- read.table(file = "CoChangedGenes.txt",sep="\t",header=F,check.names=F)   #类似flybase的list， header=F
head(Specific_list)

library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_external_gene_name <- Specific_list[,1]
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'external_gene_name', values = my_external_gene_name, mart = mart)
head(study_symbols)

colnames(Specific_list)[1]<-c("external_gene_name")
Specific_list_names <- study_symbols
rownames(Specific_list_names) <- Specific_list_names[,2]
Specific_list_names <- Specific_list_names[my_external_gene_name,]

diffSig_nameID <- Specific_list_names$external_gene_name
diffSig_nameFly <- Specific_list_names$ensembl_gene_id
list_count <- normalizeExp[diffSig_nameFly,]
list_count <- na.omit(list_count)
rownames(Specific_list_names) <- Specific_list_names[,1]
rownames(list_count) <- Specific_list_names[rownames(list_count),]$external_gene_name
write.table(list_count,file="list_count.xls",sep="\t",quote=F) ##输出list_count后，删除NA的行，修改第一列为ID的列名

#画热图
library(pheatmap)
library(RColorBrewer)

##Gene_Ensemble_id
pheatmap_list_count_new <- list_count
colnames(pheatmap_list_count_new)  ##根据列名修改下面的索引结果
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
#pheatmap_list_count_new <- pheatmap_list_count_new[,c("Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]
pheatmap_list_count_new <- pheatmap_list_count_new[,c("Scrib_Wts-1B","Scrib_Wts-2B","Scrib_Wts-3B","Ras_Scrib-1B","Ras_Scrib-2B","Ras_Scrib-3B","E75-scrib-wts-1A","E75-scrib-wts-2A","E75-scrib-wts-3A")]

##画图   #记得改热图的名字
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heatmap=pheatmap(pheatmap_list_count_new, color = colors,
                 main="Signaling Genes",
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
                 filename = "Signaling Genes.pdf"
)
  
  
  
  
  
  
  
  
  
  
  
#PATHVIEW
#BiocManager::install("pathview")
#BiocManager::install("gage")
#BiocManager::install("gageData")
#install.packages("dplyr")
library(gage)
library(pathview)
library(gageData)
library(dplyr)
#加载数据
data(kegg.sets.mm)
data(sigmet.idx.mm)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs =  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)

#提取logFC值并与id相连
foldchanges = diff_name2$logFC
names(foldchanges)= gene.df$ENTREZID
head(foldchanges)
#开始pathway分析，获取结果
keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)
lapply(keggres, head)
# Look at both up (greater), down (less), and statatistics.

#得到分析的结果输出前30个的pathview结果
 keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
       tibble::as_tibble() %>% 
       filter(row_number()<=120) %>% 
       .$id %>% 
       as.character()
	   
keggrespathways

keggresids = substr(keggrespathways, start=1, stop=120)
keggresids
#通过pathview包中的pathway()函数画图。下面写一个函数，这样好循环画出上面产生的前10个通路图。
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))











