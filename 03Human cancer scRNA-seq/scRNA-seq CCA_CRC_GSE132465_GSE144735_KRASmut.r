

##Step 1. Reading Files
##Resolving the difference between left-sided and right-sided colorectal cancer by single-cell sequencing
##https://insight.jci.org/articles/view/152616
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188711
##GSE188711
##
library(Seurat)
library(dplyr)
library(ggpubr)
library(harmony)
getwd()
setwd("F:\\05Human Relavance")
setwd("J:\\05Human Relavance")
list.files()
list.files("./05CRC/GSE188711")

GSM5688706 = Read10X("./05CRC/GSE188711/GSM5688706")
GSM5688707 = Read10X("./05CRC/GSE188711/GSM5688707")
GSM5688708 = Read10X("./05CRC/GSE188711/GSM5688708")
GSM5688709 = Read10X("./05CRC/GSE188711/GSM5688709")
GSM5688710 = Read10X("./05CRC/GSE188711/GSM5688710")
GSM5688711 = Read10X("./05CRC/GSE188711/GSM5688711")


GSM5688706 = CreateSeuratObject(counts = GSM5688706, project = "GSM5688706", min.features = 200)
GSM5688707 = CreateSeuratObject(counts = GSM5688707, project = "GSM5688707", min.features = 200)
GSM5688708 = CreateSeuratObject(counts = GSM5688708, project = "GSM5688708", min.features = 200)
GSM5688709 = CreateSeuratObject(counts = GSM5688709, project = "GSM5688709", min.features = 200)
GSM5688710 = CreateSeuratObject(counts = GSM5688710, project = "GSM5688710", min.features = 200)
GSM5688711 = CreateSeuratObject(counts = GSM5688711, project = "GSM5688711", min.features = 200)


GSE188711 = merge(GSM5688706, y=c(GSM5688707, GSM5688708, GSM5688709, GSM5688710, GSM5688711))
data = GSE188711
data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-") 
data[["percent.rb"]] = PercentageFeatureSet(data, pattern = "^RP[SL]")

violin <- VlnPlot(data,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  pt.size = 0, 
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 20, height = 6)

data2 = subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1500 & nCount_RNA < 40000 & percent.mt < 15) 
saveRDS(data2, file="CRC_Raw_GSE188711.Rds")                     #6 samples

violin <- VlnPlot(data2,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  pt.size = 0, 
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 20, height = 6)


##Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer
##https://www.nature.com/articles/s41588-020-0636-z
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735
##GSE144735
##
library(Seurat)
library(dplyr)
library(ggpubr)
library(harmony)
library(data.table)
getwd()
setwd("F:\\05Human Relavance")
list.files()
list.files("./05CRC/GSE144735")

# Clear all objects in the environment
rm(list = ls())
# Set options to prevent the read string from being converted to a factor type
options(stringsAsFactors = F)

# Define data directory path
dir = './05CRC/GSE144735/'

# Read count files in. gz format using the fread function in the data.table package
ct <- fread('./05CRC/GSE144735/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz', data.table = F)
print(ct[1:4, 1:4])

#Set the row name as the value of the first column, which are usually gene identifiers
rownames(ct) <- ct[, 1]

#Delete the first column as it has already been set as a row name
ct <- ct[, -1]

#Create Seurat object and use the read count data
data <- CreateSeuratObject(counts = ct)

data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-") 
data[["percent.rb"]] = PercentageFeatureSet(data, pattern = "^RP[SL]")

data2 = subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1500 & nCount_RNA < 40000 & percent.mt < 15) 
saveRDS(data2, file="CRC_Raw_GSE144735.Rds")                     #6 samples

##Isolate KRAS-mut Patients
Cell.sub <- subset(data2@meta.data, orig.ident %in% c("KUL21-T","KUL30-T","KUL31-T"))
table(Cell.sub$orig.ident)

scRNAsub <- subset(data2, cells=row.names(Cell.sub))
table(scRNAsub@meta.data$orig.ident)
saveRDS(scRNAsub, file="CRC_Raw_GSE144735_KRASmut.Rds")                     #3 samples

##GSE132465
##
library(Seurat)
library(dplyr)
library(ggpubr)
library(harmony)
library(data.table)
getwd()
setwd("F:\\05Human Relavance")
list.files()
list.files("./05CRC/GSE132465")

#Clear all objects in the environment
rm(list = ls())
#Set options to prevent the read string from being converted to a factor type
options(stringsAsFactors = F)

#Define data directory path
dir = './05CRC/GSE132465/'

#Read count files in. gz format using the fread function in the data.table package
ct <- fread('./05CRC/GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz', data.table = F)

#Print and read the first 4 rows and 4 columns of data
print(ct[1:4, 1:4])

#Set the row name as the value of the first column, which are usually gene identifiers
rownames(ct) <- ct[, 1]

#Delete the first column as it has already been set as a row name
ct <- ct[, -1]

#Create Seurat object and use the read count data
data <- CreateSeuratObject(counts = ct)

data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-") 
data[["percent.rb"]] = PercentageFeatureSet(data, pattern = "^RP[SL]")

data2 = subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1500 & nCount_RNA < 40000 & percent.mt < 15) 
saveRDS(data2, file="CRC_Raw_GSE132465.Rds")                     #6 samples

##Isolate KRAS-mut Patients
Cell.sub <- subset(data2@meta.data, orig.ident %in% c("SMC01-T","SMC04-T","SMC05-T","SMC06-T","SMC07-T","SMC14-T","SMC19-T","SMC20-T","SMC23-T","SMC24-T"))
table(Cell.sub$orig.ident)

scRNAsub <- subset(data2, cells=row.names(Cell.sub))
table(scRNAsub@meta.data$orig.ident)

saveRDS(scRNAsub, file="CRC_Raw_GSE132465_KRASmut.Rds")                     #10 samples


##Step 2. Data Integration

#GSE188711 <- readRDS("CRC_Raw_GSE188711.Rds")
GSE144735 <- readRDS("CRC_Raw_GSE144735.Rds")
GSE132465 <- readRDS("CRC_Raw_GSE132465.Rds")

#KRASmut Only
GSE144735 <- readRDS("CRC_Raw_GSE144735_KRASmut.Rds")
GSE132465 <- readRDS("CRC_Raw_GSE132465_KRASmut.Rds")

#DefaultAssay(GSE188711) <- "RNA"
DefaultAssay(GSE144735) <- "RNA"
DefaultAssay(GSE132465) <- "RNA"

GSE.anchors <- FindIntegrationAnchors(object.list = list(GSE132465, GSE144735), anchor.features = 2000, dims = 1:50)
GSE.combined <- IntegrateData(anchorset = GSE.anchors, dims = 1:50)
pbmc <-  GSE.combined

#pbmc <- NormalizeData(pbmc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2500)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:40
pbmc <- RunUMAP(pbmc, dims=pc.num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
pbmc = FindClusters(pbmc,resolution = 1.0)

saveRDS(pbmc, file="CCA_CRC_GSE132465_GSE144735.Rds") 
saveRDS(pbmc, file="CCA_CRC_GSE132465_GSE144735_KRASmut.Rds") 



##Step 3. tSNEplot
library(devtools)
library(DropletUtils)
library(Seurat)
library(scran)
library(scDblFinder)
library(scater)
library(DoubletFinder)
library(scDblFinder)
library(ggplot2)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
#library(scCustomize)
library(viridis)
library(RColorBrewer)
library(gridExtra)
set.seed(1234)
getwd()
setwd("J:\\05Human Relavance")


pbmc <- readRDS("CCA_CRC_GSE132465_GSE144735_KRASmut.Rds")

##Check the default assay
DefaultAssay(pbmc)

##Check how many clusters there are
levels(pbmc)

##Set the number of color combinations for the corresponding cluster
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(30)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(30)

show_col(col5)

##Normal UMAP dimensionality reduction data results
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

pc.num=1:40
pbmc <- RunTSNE(pbmc, dims=pc.num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
pbmc = FindClusters(pbmc,resolution = 0.5)

##Check how many clusters there are
levels(pbmc)
##Set the number of color combinations for the corresponding cluster
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(20)
show_col(col5)

UMPplot_label <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_TSNEplot_label.pdf", plot = plot_grid(UMPplot_label), width = 9, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(legend.position = "none", panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_TSNEplot_label_nolegend.pdf", plot = plot_grid(UMPplot_unlabel), width = 8.5, height = 8)

UMPplot_label_nolegend <- DimPlot(pbmc, reduction = "tsne", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_TSNEplot_unlabel.pdf", plot = plot_grid(UMPplot_label_nolegend), width = 9, height = 8)
UMPplot_label_nolegend <- DimPlot(pbmc, reduction = "tsne", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(legend.position = "none", panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_TSNEplot_unlabel_nolegend.pdf", plot = plot_grid(UMPplot_label_nolegend), width = 8.5, height = 8)

##Step 4. Cell Type Annotation
#Marker genes for each cell type/or genes that need to be displayed
DefaultAssay(pbmc) <- "RNA"

pbmc <- NormalizeData(pbmc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)

T_cells_genes <- c("CD3D", "CD3E", "CD8A", "CD4", "CD2")
B_cells_genes <- c("CD19", "CD79A", "MS4A1", "VPREB3")
Plasma_cells_genes <- c("SDC1","MZB1","IGHG1")
NK_cells_genes <- c("FGFBP2", "FCG3RA", "CX3CR1","NCAM1","GNLY")
Myeloid_cells_genes <- c("PTPRC", "CD14", "AIF1", "TYROBP", "CD163","CD11B")
Epithelial_cells_genes <- c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")   # epi or tumor 
Glial_cells_genes <- c("S100B","CDH2")
Mast_cells_genes <- c("IL1RL1", "KIT", "MS4A2", "CPA3", "CST3")

#Mesenchymal_stromal_cells
	Fibroblasts_genes <- c("ACTA2", "COL1A1", "PDGFRB", "DCN","FAP")
	Endothelial_cells_genes <- c("PECAM1", "VWF","CD31", "CD34")

#Myeloid_cells
Myeloid_cells_genes <- c("PTPRC", "CD14", "AIF1", "TYROBP", "CD163","CD11b")
	Monocytes_genes <- c("CD14", "CD16","CD11b","LYVE1")
	Macrophages_genes <- c("CD68", "CD163", "CD14","CD11b")
	Dendritic_cells_genes <- c("CD1C","CD141")
	Neutrophils_genes <- c("CD11b","CD15","CD66b","LYZ","MPO","ELANE","CSF3R")
	Eosinophils_genes <- c("CD11b","CD163","EG2","MBP","EOS")
	Basophils_genes <- c("CD123","CD203c","FCER1A","HRH1","CMA1")
	Mast_cells_genes <- c("IL1RL1", "KIT", "MS4A2", "CPA3", "CST3")

#Epithelial cells
Epithelial_cells_genes <- c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")   # epi or tumor 
	Intestinal_stem_cells_genes <- c()
	Paneth_cells_genes <- c()
	Enterocytes_genes <- c()
	Goblet_cells_genes <- c()
	Enteroendocrine_cells_genes <- c()
 

#The genes are stored as a list, which is the input of DotPlot to achieve the facet effect
features <- list("T_cell" = T_cells_genes,
                 "B_cell" = B_cells_genes,
				 "Plasma_cell" = Plasma_cells_genes,
				 "NK_cell" = NK_cells_genes,
				 "Myeloid_cell" = Myeloid_cells_genes,
				 "Epithelial_cell" = Epithelial_cells_genes,
				 "Glial_cell" = Glial_cells_genes,
				 "Fibroblast" = Fibroblasts_genes,
				 "Endothelial_cell" = Endothelial_cells_genes,
				 "Mast_cells" = Mast_cells_genes)
				 
#Rapid annotation of cell types
Adotplot <- DotPlot(object = pbmc, features=features)
Adotplot 

##set celltype
levels(pbmc)
levels(pbmc@active.ident)
head(pbmc@meta.data)

pbmc@meta.data$celltype='Unknown cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(0, 1, 3, 11, 15, 16),'celltype'] = 'T/NK cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(8),'celltype'] = 'B cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(6),'celltype'] = 'Plasma cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(2,4),'celltype'] = 'Myeloid cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(7, 9, 10, 13, 17, 18, 19),'celltype'] = 'Epithelial cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(5, 14),'celltype'] = 'Fibroblasts'
pbmc@meta.data[pbmc$seurat_clusters %in% c(12),'celltype'] = 'Endothelial cells'

head(pbmc@meta.data)
table(pbmc@meta.data[["celltype"]])


saveRDS(pbmc, file="CCA_CRC_GSE132465_GSE144735_KRASmut_Annotated.Rds") 


##Step 5. Cell Annotation Draw Dotplot
pbmc <- readRDS("CCA_CRC_GSE132465_GSE144735_KRASmut_Annotated.Rds")

T_cells_genes <- c("CD3D", "CD3E", "CD2", "CD8A", "CD4")
B_cells_genes <- c("CD19",  "MS4A1", "VPREB3","CD79A")
Plasma_cells_genes <- c("MZB1","SDC1","IGHG1")
NK_cells_genes <- c("GNLY","NKG2D","KIRs","PRF1","GZMB","IFNG")
Myeloid_cells_genes <- c("CD14", "AIF1", "TYROBP", "CD163","CD11B")
Epithelial_cells_genes <- c("EPCAM", "KRT19", "CD24", "PROM1")   # epi or tumor 
Glial_cells_genes <- c("S100B","CDH2","SOX10","RET","PHOX2B","ASCL1","CHAT")
Mast_cells_genes <- c("IL1RL1", "KIT", "MS4A2", "CPA3","TPSAB1","IL4","IL13")

#Mesenchymal_stromal_cells
	Fibroblasts_genes <- c("ACTA2", "COL1A1", "PDGFRB", "DCN","FAP")
	Endothelial_cells_genes <- c("PECAM1", "VWF","CD31", "CD34")

#The genes are stored as a list, which is the input of DotPlot to achieve the facet effect
features <- list("T_cell" = T_cells_genes,
                 "NK_cell" = NK_cells_genes,
                 "B_cell" = B_cells_genes,
				 "Plasma_cell" = Plasma_cells_genes,
				 "Myeloid_cell" = Myeloid_cells_genes,
				 "Epithelial_cell" = Epithelial_cells_genes,
				 "Fibroblast" = Fibroblasts_genes,
				 "Endothelial_cell" = Endothelial_cells_genes
				 )

##Set cell type
levels(pbmc)
levels(pbmc@active.ident)

pbmc$celltype2 <- pbmc@active.ident

new.cluster.ids <- c(
							"0"="T/NK cells",
							"1"="T/NK cells",
							"2"="Myeloid cells",
							"3"="T/NK cells",
							"4"="Myeloid cells",
							"5"="Fibroblasts",
							"6"="Plasma cells",
							"7"="Epithelial cells",
							"8"="B cells",
							"9"="Epithelial cells",
							"10"="Epithelial cells",
							"11"="T/NK cells",
							"12"="Endothelial cells",
							"13"="Epithelial cells",
							"14"="Fibroblasts",
							"15"="T/NK cells",							
							"16"="T/NK cells",
							"17"="Epithelial cells",
							"18"="Epithelial cells",
							"19"="Epithelial cells"
							)

pbmc <- RenameIdents(pbmc, new.cluster.ids)    

levels(pbmc) 

#Set the display order of cell types
levels(pbmc)  <- c("T/NK cells","B cells","Plasma cells","Myeloid cells","Epithelial cells","Fibroblasts","Endothelial cells")


#Graphic modification
Adotplot <- DotPlot(object = pbmc, features=features)&
  theme_bw()& #Set theme
  geom_point(shape=21, aes(size=pct.exp),stroke=1)& 		#Add borders to points
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = 'black', size = 12, face = "bold"),
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        strip.background = element_blank(), 	#Remove the background of the sectional image
        strip.text = element_blank(),					#Remove the text from the sectional diagram
        plot.margin=unit(c(1, 1, 1, 1),'cm'),		#Narrowing down the drawing scope
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),		#Change border size
        panel.spacing = unit(0.12, "cm"),			#Modify the spacing of the sectional drawing
        # legend.frame = element_rect(colour = "black"),											#Modify the legend border
        # legend.ticks = element_line(colour = "black", linewidth  = 0),
        legend.key.width = unit(0.3, "cm"),																#Modify the width of the legend
        legend.key.height = unit(0.5, "cm"),																#Modify the height of the legend
        legend.title = element_text(color = 'black', face = "bold", size=9))& 				#Modify the font size of the legend title
#  scale_color_gradientn(colours = colorRampPalette(c("navy","white","firebrick3"))(100))& 	#Modify continuous legend color
scale_color_gradientn(colours = colorRampPalette(c("white","firebrick3"))(100))& 
  labs(tag = "Cell type annotation")& 			#Add tags
  theme(plot.tag.position = c(0.3, 1.05), 		#Set the label position and adjust it yourself
        plot.tag = element_text(size = 12,face = "bold"))& 						#Set label text size
  guides(size=guide_legend(title="Proportion of\nexpressing cells"), 	#Set the title of the legend size
         colour=guide_colorbar(title="Average\nexpression"))					#Set the title of the legend color bar

Adotplot

#45 marker genes 14 cm
ggsave("P2_Marker_gene_Adotplot.pdf",plot=Adotplot,width = 12, height = 3.5)


##Check how many clusters there are
levels(pbmc)

##Set the number of color combinations for the corresponding cluster
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(30)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(7)
show_col(col5)


UMPplot_label <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "celltype", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P3_TSNEplot_Annotated_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "celltype", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(legend.position = "none", panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P3_TSNEplot_Annotated_label_nolegend.pdf", plot = plot_grid(UMPplot_unlabel), width = 8.5, height = 8)


#marker gene
FeaturePlot(pbmc,'NRAS',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NRAS")
ggsave("P5_pbmcTSNE1_NRAS.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'KRAS',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P5_pbmcTSNE1_KRAS.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'HRAS',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("HRAS")
ggsave("P5_pbmcTSNE1_HRAS.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'KRT19',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRT19")
ggsave("P5_pbmcTSNE1_KRT19.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'EPCAM',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("EPCAM ")
ggsave("P5_pbmcTSNE1_EPCAM .pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'VIM', reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("VIM")
ggsave("P5_pbmcTSNE1_VIM.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'CD44',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CD44")
ggsave("P5_pbmcTSNE1_CD44.pdf", width = 5.5, height = 5)

FeaturePlot(pbmc,'NOTCH1',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH1")
ggsave("P5_pbmcTSNE1_NOTCH1.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'NOTCH2',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH2")
ggsave("P5_pbmcTSNE1_NOTCH2.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'NOTCH3',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH3")
ggsave("P5_pbmcTSNE1_NOTCH3.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'NOTCH4',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH4")
ggsave("P5_pbmcTSNE1_NOTCH4.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'REL',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("REL")
ggsave("P5_pbmcTSNE1_REL.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'RELA',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("RELA")
ggsave("P5_pbmcTSNE1_RELA.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'NFKB1',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NFKB1")
ggsave("P5_pbmcTSNE1_NFKB1.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'YAP1',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("YAP1")
ggsave("P5_pbmcTSNE1_YAP1.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'TAZ',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("TAZ")
ggsave("P5_pbmcTSNE1_TAZ.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'TEAD1',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("TEAD1")
ggsave("P5_pbmcTSNE1_TEAD1.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'JUN',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("JUN")
ggsave("P5_pbmcTSNE1_JUN.pdf", width = 5.5, height = 5)
FeaturePlot(pbmc,'FOS',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("FOS")
ggsave("P5_pbmcTSNE1_FOS.pdf", width = 5.5, height = 5)



#Step 5. Extract specific cell subtypes
pbmc <- readRDS("CCA_CRC_GSE132465_GSE144735_KRASmut_Annotated.Rds")

Cell.sub <- subset(pbmc@meta.data, celltype %in% c("Epithelial cells"))
table(Cell.sub$orig.ident)

scRNAsub <- subset(pbmc, cells=row.names(Cell.sub))
table(scRNAsub@meta.data$celltype)

##Re clustering analysis
DefaultAssay(scRNAsub) 
DefaultAssay(scRNAsub) <- "integrated"

scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2500)
scRNAsub <- ScaleData(scRNAsub, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub), npcs = 60, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:60
#scRNAsub <- RunUMAP(scRNAsub, dims=pc.num)
scRNAsub <- RunTSNE(scRNAsub, dims=pc.num)
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num)
scRNAsub = FindClusters(scRNAsub,resolution = 0.7)

##查看有多少个cluster
levels(scRNAsub)

##Set the number of color combinations for the corresponding cluster
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(15)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(13)

show_col(col5)


##Normal UMAP dimensionality reduction data results
UMPplot_scRNAsub_label <- DimPlot(scRNAsub, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_UMPplot_scRNAsub_label.pdf", plot = plot_grid(UMPplot_scRNAsub_label), width = 6.5, height = 6)
UMPplot_scRNAsub_unlabel <- DimPlot(scRNAsub, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_UMPplot_scRNAsub_unlabel.pdf", plot = plot_grid(UMPplot_scRNAsub_unlabel), width = 6.5, height = 6)

TSNEplot_scRNAsub_label <- DimPlot(scRNAsub, reduction = "tsne", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_TSNEplot_scRNAsub_label.pdf", plot = plot_grid(TSNEplot_scRNAsub_label), width = 8, height = 7)
TSNEplot_scRNAsub_unlabel <- DimPlot(scRNAsub, reduction = "tsne", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_TSNEplot_scRNAsub_unlabel.pdf", plot = plot_grid(TSNEplot_scRNAsub_unlabel), width = 8, height = 7)


DefaultAssay(scRNAsub) <- "RNA"

FeaturePlot <- FeaturePlot(scRNAsub,'KRAS', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P7_scRNAsubUMAP1_KRAS.pdf", plot = plot_grid(FeaturePlot), width = 4.5, height = 4)
FeaturePlot <- FeaturePlot(scRNAsub,'NRAS', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NRAS")
ggsave("P7_scRNAsubUMAP1_NRAS.pdf", plot = plot_grid(FeaturePlot), width = 4.5, height = 4)
FeaturePlot <- FeaturePlot(scRNAsub,'HRAS', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("HRAS")
ggsave("P7_scRNAsubUMAP1_HRAS.pdf", plot = plot_grid(FeaturePlot), width = 4.5, height = 4)


FeaturePlot(scRNAsub,'NOTCH1', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH1")
ggsave("P7_scRNAsubUMAP1_NOTCH1.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'NOTCH2', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH2")
ggsave("P7_scRNAsubUMAP1_NOTCH2.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'NOTCH3', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH3")
ggsave("P7_scRNAsubUMAP1_NOTCH3.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'NOTCH4', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NOTCH4")
ggsave("P7_scRNAsubUMAP1_NOTCH4.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'REL', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("REL")
ggsave("P7_scRNAsubUMAP1_REL.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'RELA', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("RELA")
ggsave("P7_scRNAsubUMAP1_RELA.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'NFKB1', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NFKB1")
ggsave("P7_scRNAsubUMAP1_NFKB1.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'YAP1', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("YAP1")
ggsave("P7_scRNAsubUMAP1_YAP1.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'TAZ', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("TAZ")
ggsave("P7_scRNAsubUMAP1_TAZ.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'TEAD1', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("TEAD1")
ggsave("P7_scRNAsubUMAP1_TEAD1.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'JUN', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("JUN")
ggsave("P7_scRNAsubUMAP1_JUN.pdf", width = 4.5, height = 4)
FeaturePlot(scRNAsub,'FOS', reduction = "tsne", order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("FOS")
ggsave("P7_scRNAsubUMAP1_FOS.pdf", width = 4.5, height = 4)

saveRDS(scRNAsub, file="CCA_CRC_GSE132465_GSE144735_KRASmut_Annotated_Epithelial.Rds") 


#Step 7. Addmodulescore
library(Seurat)
library(GSVA)

scRNAsub <- readRDS(file="CCA_CRC_GSE132465_GSE144735_KRASmut_Annotated_Epithelial.Rds")
pbmc <- scRNAsub

##Select RasFry related genes
RasFry_genes <- c(
#JNK
"HRAS", "KRAS", "NRAS",  "RAF1", "BRAF", "ARAF", "MAP3K1", "MAP3K2", "MAP3K3", "MAP3K5","MAP2K4", "MAP2K7","MAPK8", "MAPK9", "MAPK10","JUN", "FOS", "ATF2", "ELK1", "DUSP1", "DUSP5", "DUSP6",
#
"NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4","JAG1", "JAG2", "DLL1", "DLL3", "DLL4","MIB1", "MIB2", "ADAM10", "ADAM17","RBPJ", "HES1", "HES5", "HEY1", "HEY2", "HEYL","MYC", "CCND1", "HES6", "HES7",
#TLR
"TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TLR10", "MYD88", "TICAM1", "TICAM2", "TIRAP",  "IRAK1", "IRAK4", "TBK1", "CHUK", "IKBKB", "IKBKG", "NFKB1", "RELA", "RELB", "NFKB2", "REL", "IRF3", "IRF7",  "IL1B", "TNF", "IFNB1",
#HippoInactivation	   
"YAP1", "TAZ", "TEAD1", "TEAD2", "TEAD3", "TEAD4","CCN2", "CCN1", "ANKRD1", "AREG", "EREG", "FZD7", "JUN", "MYC", "BIRC5","CCNE1","CCND1","LAMA5","IGF2","LIN28B","AXIN2","PPP2R2A"
)

gene_sets2 <- as.data.frame(RasFry_genes)
pbmc <- AddModuleScore(pbmc, features = gene_sets2,name = "RasFry_score")

FeaturePlot(pbmc,'RasFry_score1',reduction = "tsne",order = TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("RasFry_score1")
ggsave("P8_pbmcTSNE_RasFry_score.pdf", width = 4.5, height = 4)


#Step 8. Export Meta file to draw correlation analysis

Meta <- as.data.frame(pbmc@meta.data)
write.csv(Meta, file = "P9_pbmc_metadata.csv",row.names = T)

a <- Meta[,c("RAS_score1", "TLR_score1", "NOTCH_score1", "JNK_score1", "MAPK_score1", "HIPPO_inactivation_score1","Glycolysis_score1")]


library(ggplot2)
library(ggpubr)
library(ggExtra)

gene1="RAS_score1"             			
gene2="Glycolysis_score1"              

x=Meta$RAS_score1
y=Meta$TLR_score1

df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value

p1=ggplot(df1, aes(x, y)) + 
      xlab(gene1)+ylab(gene2)+
      geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))

p2

saveRDS(pbmc, file="CCA_CRC_GSE132465_GSE144735_KRASmut_Annotated_Epithelial_Scored.Rds")


#Step 9. Export Ras highly activated clusters to draw correlation analysis
library(devtools)
library(DropletUtils)
library(Seurat)
library(scran)
library(scDblFinder)
library(scater)
library(DoubletFinder)
library(scDblFinder)
library(ggplot2)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(reshape)
library(ggpubr)
set.seed(1234)
getwd()
setwd("I:\\05Human Relavance")

pbmc <- readRDS("CCA_CRC_GSE132465_GSE144735_KRASmut_Annotated_Epithelial_Scored.Rds")

levels(pbmc)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(13)


VlnPlot(
    pbmc,raster = FALSE, pt.size = 0.1, cols = col5,
    features = c("RAS_score1"), 
    group.by = "seurat_clusters"  
) + theme_classic()+ 
	scale_y_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2)) + 
	theme(plot.title = element_text(size=20,hjust = 0.5),  
			axis.line = element_line(color = "black",linewidth = 1),
			text = element_text (size = 18),
			axis.text = element_text (color = "black", size = 18))
ggsave("P10_pbmcVlnPlot_RAS_signaling_score.pdf", width = 8, height = 4)  

VlnPlot(
    pbmc,raster = FALSE, pt.size = 0, cols = col5,
    features = c("RAS_score1"), 
    group.by = "seurat_clusters"  
) + theme_classic()+ 
	scale_y_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2)) + 
	theme(plot.title = element_text(size=20,hjust = 0.5),  
			axis.line = element_line(color = "black",linewidth = 1),
			text = element_text (size = 18),
			axis.text = element_text (color = "black", size = 18))
ggsave("P10_pbmcVlnPlot_RAS_signaling_score_no_dots.pdf", width = 8, height = 4)  


VlnPlot(
    pbmc,raster = FALSE, pt.size = 0.1, ncol = 3, cols = col5,
    features = c("RAS_score1", "MAPK_score1", "JNK_score1", "TLR_score1", "NOTCH_score1","HIPPO_inactivation_score1"), 
    group.by = "seurat_clusters",   
)
ggsave("P10_pbmcVlnPlot_signaling_score.pdf", width = 15, height = 6)


#Warning: 'ncol' is ignored with 'stack' is TRUE
#Warning: 'y.max' is ignored when 'stack' is TRUE



#Step 9.1 Ras Signaling correlations with Other Signaling

clusters <- pbmc@meta.data$seurat_clusters
#Select the corresponding genes
signaling_scores <- pbmc@meta.data [["RAS_score1"]]
#Calculate the average value of a specific signaling pathway

cluster_means <- tapply(signaling_scores, clusters, mean)
cluster_means1 <- as.data.frame(cluster_means)
cluster_means1$cluster <- rownames(cluster_means1)

cluster_means2 <- cluster_means1[order(-cluster_means1$cluster_means),]

print(cluster_means)
print(cluster_means2)

#Extract specific cell subtypes
Cell.sub <- subset(pbmc@meta.data, seurat_clusters %in% c("17", "18", "4"))

table(Cell.sub$orig.ident)

scRNAsub <- subset(pbmc, cells=row.names(Cell.sub))
table(scRNAsub@meta.data$seurat_clusters)
meta <- scRNAsub@meta.data

meta <- meta[,c("seurat_clusters","RAS_score1", "RasFry_score1", "MAPK_score1", "JNK_score1", "TLR_score1", "NOTCH_score1","HIPPO_inactivation_score1")]

library(reshape)
library(ggpubr)


a <- c("RasFry_score1", "MAPK_score1", "JNK_score1", "TLR_score1", "NOTCH_score1","HIPPO_inactivation_score1")
names(a) <- c(1:6)

for (i in 1:6) {

meta_signaling_wide <- meta[,c("seurat_clusters","RAS_score1", a[i])]
file_name <- paste('P10_Cor_RAS_score1', a[i],".csv")
write.csv(meta_signaling_wide,file_name, row.names = FALSE)
cat("File", file_name, "saved.\n")    

input <- meta_signaling_wide

file_name2 <- paste('P11_scRNAsub_Cor_RAS_score_', a[i],".pdf")

print(colnames(input)[3])
colnames(input)[3] = "NewSignaling"

Corplot <- ggplot(input, mapping=aes(x=NewSignaling, y=RAS_score1))+  
  geom_point(size=1)+
  geom_smooth(method = 'lm',
              formula = 'y ~ x',
              se=T,
              lwd=1,
              color = "#9f0000", 
              fill = "lightgrey")+
  stat_cor(method='spearman',
           label.x = -0.1, 
           label.y = 0.7, 
           label.sep = "\n",
           size=5,p.accuracy = 0.001)+
 # labs(title='Correlationships of Different Signaling Pathways')+
		   theme_classic()+ 
		   scale_y_continuous(limits = c(-0.3, 0.9), breaks = seq(-0.3, 0.9, 0.3)) + 
		   #scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2) + 
  theme(plot.title = element_text(size=20,hjust = 0.5),  
			axis.line = element_line(color = "black",linewidth = 1),
			text = element_text (size = 18),
			axis.text = element_text (color = "black", size = 18))

ggsave(filename= file_name2, plot= Corplot, width = 4.5, height = 4.5)

cat("File", file_name2, "saved.\n")   
}


#####Step 9.2 RasFry scoring correlations with Other Cancer Progression Biomarkers

meta <- pbmc@meta.data
meta <- meta[,c("seurat_clusters","RasFry_score1", "Cell_Proliferation_score1", "EMT_genes1", "Cell_Migration_genes1", "Immune_Evasion_genes1","Glycolysis_score1")]

library(reshape)
library(ggpubr)

a <- c("Cell_Proliferation_score1", "EMT_genes1", "Cell_Migration_genes1", "Immune_Evasion_genes1","Glycolysis_score1")
names(a) <- c(1:5)

for (i in 1:5) {

###########################Step 9.2.1 all cells

meta_signaling_wide <- meta[,c("seurat_clusters","RasFry_score1", a[i])]
file_name <- paste('P12_Cor_RasFry_score1', a[i],".csv")
write.csv(meta_signaling_wide,file_name, row.names = FALSE)
cat("File", file_name, "saved.\n")   

input <- meta_signaling_wide

file_name2 <- paste('P13_pbmc_Cor_RasFry_score_', a[i],".pdf")

print(colnames(input)[3])
colnames(input)[3] = "NewSignaling"

Corplot <- ggplot(input, mapping=aes(x=NewSignaling, y=RasFry_score1))+  
  geom_point(size=1)+
  geom_smooth(method = 'lm',
              formula = 'y ~ x',
              se=T,
              lwd=1,
              color = "#9f0000", 
              fill = "lightgrey")+
  stat_cor(method='spearman',
           label.x = 0.0, 
           label.y = 0.4, 
           label.sep = "\n",
           size=5,p.accuracy = 0.001)+
 # labs(title='Correlationships of Different Signaling Pathways')+
		   theme_classic()+ 
		   scale_y_continuous(limits = c(-0.3, 0.5), breaks = seq(-0.3, 0.5, 0.2)) +
		   #scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2) +
  theme(plot.title = element_text(size=20,hjust = 0.5),  
			axis.line = element_line(color = "black",linewidth = 1),
			text = element_text (size = 18),
			axis.text = element_text (color = "black", size = 18))

ggsave(filename= file_name2, plot= Corplot, width = 4.5, height = 4.5)

cat("File", file_name2, "saved.\n")  


###########################Step 9.2.2 Ras activated cells

meta_scRNAsub <- subset(meta, seurat_clusters %in% c("5", "6", "3"))

meta_signaling_wide <- meta_scRNAsub[,c("seurat_clusters","RasFry_score1", a[i])]
file_name <- paste('P14_scRNAsub_Cor_RasFry_score1', a[i],".csv")
write.csv(meta_signaling_wide,file_name, row.names = FALSE)
cat("File", file_name, "saved.\n")   

input <- meta_signaling_wide

file_name2 <- paste('P15_scRNAsub_Cor_RasFry_score_', a[i],".pdf")

print(colnames(input)[3])
colnames(input)[3] = "NewSignaling"

Corplot <- ggplot(input, mapping=aes(x=NewSignaling, y=RasFry_score1))+ 
  geom_point(size=1)+
  geom_smooth(method = 'lm',
              formula = 'y ~ x',
              se=T,
              lwd=1,
              color = "#9f0000", 
              fill = "lightgrey")+
  stat_cor(method='spearman',
           label.x = 0.0, 
           label.y = 0.4, 
           label.sep = "\n",
           size=5,p.accuracy = 0.001)+
 # labs(title='Correlationships of Different Signaling Pathways')+
		   theme_classic()+ 
		   scale_y_continuous(limits = c(-0.3, 0.5), breaks = seq(-0.3, 0.5, 0.2)) + 
		   #scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2) + 
  theme(plot.title = element_text(size=20,hjust = 0.5), 
			axis.line = element_line(color = "black",linewidth = 1),
			text = element_text (size = 18),
			axis.text = element_text (color = "black", size = 18))

ggsave(filename= file_name2, plot= Corplot, width = 4.5, height = 4.5)

cat("File", file_name2, "saved.\n")   
}



