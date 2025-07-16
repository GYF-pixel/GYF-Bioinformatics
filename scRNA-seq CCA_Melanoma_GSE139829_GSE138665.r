

##Step 1. Reading Files

#Single-cell RNA sequencing reveals intratumoral heterogeneity in primary uveal melanomas and identifies HES6 as a driver of the metastatic disease
#https://www.nature.com/articles/s41418-020-00730-7
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138665
##GSE138665
library(Seurat)
library(dplyr)
library(ggpubr)
library(harmony)
getwd()
setwd("F:\\05Human Relavance")
list.files()
list.files("./01Melanoma/GSE138665_RAW")

LH16_3814 = Read10X_h5(".\\01Melanoma\\GSE138665_RAW\\GSM4107899_LH16.3814_raw_gene_bc_matrices_h5.h5")
LH17_364 = Read10X_h5(".\\01Melanoma\\GSE138665_RAW\\GSM4107900_LH17.364_raw_gene_bc_matrices_h5.h5")
LH17_530 = Read10X_h5(".\\01Melanoma\\GSE138665_RAW\\GSM4107901_LH17.530_raw_gene_bc_matrices_h5.h5")
LH17_3222 = Read10X_h5(".\\01Melanoma\\GSE138665_RAW\\GSM4107902_LH17.3222_raw_gene_bc_matrices_h5.h5")
LH17_3554 = Read10X_h5(".\\01Melanoma\\GSE138665_RAW\\GSM4107903_LH17.3554_raw_gene_bc_matrices_h5.h5")
LH18_277 = Read10X_h5(".\\01Melanoma\\GSE138665_RAW\\GSM4107904_LH18.277_raw_gene_bc_matrices_h5.h5")

LH16_3814 = CreateSeuratObject(counts = LH16_3814, project = "LH16_3814", min.features = 100)
LH17_364 = CreateSeuratObject(counts = LH17_364, project = "LH17_364", min.features = 100)
LH17_530 = CreateSeuratObject(counts = LH17_530, project = "LH17_530", min.features = 100)
LH17_3222 = CreateSeuratObject(counts = LH17_3222, project = "LH17_3222", min.features = 100)
LH17_3554 = CreateSeuratObject(counts = LH17_3554, project = "LH17_3554", min.features = 100)
LH18_277 = CreateSeuratObject(counts = LH18_277, project = "LH18_277", min.features =100)


GSE138665  = merge(LH16_3814,y=c(LH17_364,LH17_530,LH17_3222,LH17_3554,LH18_277))
data = GSE138665
data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-") 
data[["percent.rb"]] = PercentageFeatureSet(data, pattern = "^RP[SL]")
data2 = subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10) 
saveRDS(data2, file="Melanoma_Raw_GSE138665.Rds")                     #6 samples



##Single-cell analysis reveals new evolutionary complexity in uveal melanoma
##https://www.nature.com/articles/s41467-019-14256-1
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139829
##GSE139829
list.files("./01Melanoma/GSE139829_RAW/")

GSM4147091_BSSR0022 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147091_BSSR0022")
GSM4147092_UMM041L = Read10X("./01Melanoma/GSE139829_RAW/GSM4147092_UMM041L")
GSM4147093_UMM059 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147093_UMM059")
GSM4147094_UMM061 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147094_UMM061")
GSM4147095_UMM062 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147095_UMM062")
GSM4147096_UMM063 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147096_UMM063")
GSM4147097_UMM064 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147097_UMM064")
GSM4147098_UMM065 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147098_UMM065")
GSM4147099_UMM066 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147099_UMM066")
GSM4147100_UMM067L = Read10X("./01Melanoma/GSE139829_RAW/GSM4147100_UMM067L")
GSM4147101_UMM069 = Read10X("./01Melanoma/GSE139829_RAW/GSM4147101_UMM069")


GSM4147091_BSSR0022 = CreateSeuratObject(counts = GSM4147091_BSSR0022, project = "GSM4147091_BSSR0022", min.features = 200)
GSM4147092_UMM041L = CreateSeuratObject(counts = GSM4147092_UMM041L, project = "GSM4147092_UMM041L", min.features = 200)
GSM4147093_UMM059 = CreateSeuratObject(counts = GSM4147093_UMM059, project = "GSM4147093_UMM059", min.features = 200)
GSM4147094_UMM061 = CreateSeuratObject(counts = GSM4147094_UMM061, project = "GSM4147094_UMM061", min.features = 200)
GSM4147095_UMM062 = CreateSeuratObject(counts = GSM4147095_UMM062, project = "GSM4147095_UMM062", min.features = 200)
GSM4147096_UMM063 = CreateSeuratObject(counts = GSM4147096_UMM063, project = "GSM4147096_UMM063", min.features = 200)
GSM4147097_UMM064 = CreateSeuratObject(counts = GSM4147097_UMM064, project = "GSM4147097_UMM064", min.features = 200)
GSM4147098_UMM065 = CreateSeuratObject(counts = GSM4147098_UMM065, project = "GSM4147098_UMM065", min.features = 200)
GSM4147099_UMM066 = CreateSeuratObject(counts = GSM4147099_UMM066, project = "GSM4147099_UMM066", min.features = 200)
GSM4147100_UMM067L = CreateSeuratObject(counts = GSM4147100_UMM067L, project = "GSM4147100_UMM067L", min.features = 200)
GSM4147101_UMM069 = CreateSeuratObject(counts = GSM4147101_UMM069, project = "GSM4147101_UMM069", min.features = 200)

GSE139829 = merge(GSM4147091_BSSR0022,y=c(GSM4147092_UMM041L, GSM4147093_UMM059, GSM4147094_UMM061, 
																							GSM4147095_UMM062, GSM4147096_UMM063, GSM4147097_UMM064,
																							GSM4147098_UMM065, GSM4147099_UMM066, GSM4147100_UMM067L, GSM4147101_UMM069))
data = GSE139829
data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-") 
data[["percent.rb"]] = PercentageFeatureSet(data, pattern = "^RP[SL]")
data2 = subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10) 
saveRDS(data2, file="Melanoma_Raw_GSE139829.Rds")                     #11 samples


#Melanoma-2018-Cell-Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma
#https://www.cell.com/cell/fulltext/S0092-8674(18)31394-1?cid=tw%26p#sec-4
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575
##GSE120575
library(Seurat)
library(dplyr)
library(ggpubr)
library(harmony)
library(data.table)
library(gridExtra)
getwd()
setwd("F:\\05Human Relavance")
list.files()
list.files("./01Melanoma/GSE120575")


#Clear all objects in the environment
rm(list = ls())
#Set options to prevent the read string from being converted to a factor type
options(stringsAsFactors = F)

#Read data
file_path <- "./01Melanoma/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"

	 geo_acc <- "GSE120575"
	 datadir <- "./01Melanoma/GSE120575"
   ## Load expression matrix and metadata
    exp.mat <- read.delim(sprintf("./01Melanoma/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
        datadir, geo_acc), header = F, sep = "\t")
    genes <- exp.mat[c(-1, -2), 1]
    cells <- as.vector(t(exp.mat[1, 2:16292]))
    samples <- as.factor(t(exp.mat[2, 2:16292]))

    exp.mat <- exp.mat[c(-1, -2), 2:16292]
    colnames(exp.mat) <- cells
    rownames(exp.mat) <- genes

    meta <- read.delim(sprintf("./01Melanoma/GSE120575/GSE120575_patient_ID_single_cells.txt.gz",
        datadir, geo_acc), header = T, sep = "\t", skip = 19, nrows = 16291)
    meta <- meta[, 1:7]

    treat <- factor(ifelse(grepl("Post", samples), "Post", "Pre"))
    response <- factor(meta$characteristics..response)
    therapy <- factor(meta$characteristics..therapy)

    ## Create Seurat object and add meta data
    query.object <- CreateSeuratObject(counts = exp.mat, project = "SadeFeldman",
        min.cells = 10)
    rm(exp.mat)
    query.object@meta.data$Sample <- samples
    query.object@meta.data$Time <- treat
    query.object@meta.data$Response <- response
    query.object@meta.data$Therapy <- therapy

	saveRDS(query.object, file="Melanoma_Raw_GSE120575.Rds")                     #11 samples





##Step 2. Data Integration
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
#library(scCustomize) # 需要Seurat版本4.3.0
library(viridis)
library(RColorBrewer)
library(gridExtra)
set.seed(1234)
getwd()
#setwd("J:\\05Human Relavance")
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\05Human Relavance")

GSE139829 <- readRDS("Melanoma_Raw_GSE139829.Rds")
GSE138665 <- readRDS("Melanoma_Raw_GSE138665.Rds")

DefaultAssay(GSE139829) <- "RNA"
DefaultAssay(GSE138665) <- "RNA"

GSE.anchors <- FindIntegrationAnchors(object.list = list(GSE139829, GSE138665), anchor.features = 2000, dims = 1:50)
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

saveRDS(pbmc, file="CCA_Melanoma_GSE139829_GSE138665.Rds") 


##Step 3. tSNEplot
pbmc <- readRDS("CCA_Melanoma_GSE139829_GSE138665.Rds")

##Check the default assay
DefaultAssay(pbmc)

##Check how many clusters there are
levels(pbmc)

##Set the number of color combinations for the corresponding cluster
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(30)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(54)

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
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(40)
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

T_cells_genes <- c("TRAC","CD3D", "CD69","CD3E", "CD8A", "CD8B", "CD4", "CD2")
NK_cells_genes <- c("FGFBP2", "FCG3RA", "CX3CR1","NCAM1","GNLY","NKG7","KLRD1")
B_cells_genes <- c("CD19", "CD20", "CD79A", "MS4A1", "VPREB3", "IGHM", "IGLL1")
Plasma_cells_genes <- c("SDC1","MZB1","IGHG1","JCHAIN")
Myeloid_cells_genes <- c("C1QA", "CPA3","CD68","CD14", "AIF1", "TYROBP", "CD163","CD11B")
Melanoma_cells_genes <- c("MITF", "PMEL", "MLANA", "TYR", "CDK4", "BCL2", "RAB27A") 
Photoreceptor_cell_genes <- c("RCVRN","CLU","PTGDS","RBP1","FRZB","TTR")

#Mesenchymal_stromal_cells
	Fibroblasts_genes <- c("ACTA2", "COL1A1", "COL1A2", "PDGFRB", "PDGFRA", "THY1", "DCN","FAP")
	Endothelial_cells_genes <- c("FN1","RETN","PECAM1", "VWF","CD31", "CD34","CLDN5", "COL4A1")

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
				 "NK_cell" = NK_cells_genes,
				 "Plasma_cell" = Plasma_cells_genes,
				 "Myeloid_cell" = Myeloid_cells_genes,
				 "Melanoma_cell" = Melanoma_cells_genes,
				 "Photoreceptor_cell" = Photoreceptor_cell_genes,
				 "Fibroblast" = Fibroblasts_genes,
				 "Endothelial_cell" = Endothelial_cells_genes
)
				 
#Rapid annotation of cell types
Adotplot <- DotPlot(object = pbmc, features=features)
Adotplot 

##set celltype
levels(pbmc)
levels(pbmc@active.ident)
head(pbmc@meta.data)

pbmc@meta.data$celltype='Malignant tumor cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(3, 8, 13),'celltype'] = 'T/NK cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(18, 34),'celltype'] = 'B/Plasma cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(4, 14),'celltype'] = 'Myeloid cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(36),'celltype'] = 'Photoreceptor cells'

#pbmc@meta.data[pbmc$seurat_clusters %in% c(20),'celltype'] = 'Glial cells'
pbmc@meta.data[pbmc$seurat_clusters %in% c(26),'celltype'] = 'Fibroblasts'
pbmc@meta.data[pbmc$seurat_clusters %in% c(32),'celltype'] = 'Endothelial cells'

head(pbmc@meta.data)
table(pbmc@meta.data[["celltype"]])


saveRDS(pbmc, file="CCA_Melanoma_GSE139829_GSE138665_Annotated.Rds") 


##Step 5. Cell Annotation Draw Dotplot
pbmc <- readRDS("CCA_Melanoma_GSE139829_GSE138665_Annotated.Rds")

T_cells_genes <- c("CD3D", "CD69","CD3E", "CD8A", "CD8B", "CD2", "CD4")
NK_cells_genes <- c("NCAM1","NKG7","KLRD1","GNLY")
B_cells_genes <- c("CD19", "CD20", "MS4A1", "VPREB3", "IGHM", "CD79A")
Plasma_cells_genes <- c("SDC1","MZB1","IGHG1","JCHAIN")
Myeloid_cells_genes <- c("C1QA", "CD68","CD14", "AIF1", "TYROBP", "CD163","CD11B")
Melanoma_cells_genes <- c("MITF", "PMEL", "MLANA", "CDK4", "TYR", "BCL2") 
Photoreceptor_cell_genes <- c("RCVRN","RBP3", "RGS9", "RHO","PDE6B","NRL")

#Mesenchymal_stromal_cells
	Fibroblasts_genes <- c("ACTA2", "COL1A1", "COL1A2", "PDGFRB","THY1", "DCN","FN1")
	Endothelial_cells_genes <- c("COL4A1","PECAM1", "VWF","CD31", "CD34","CLDN5")


#The genes are stored as a list, which is the input of DotPlot to achieve the facet effect
features <- list("T_cell" = T_cells_genes,
				 "NK_cell" = NK_cells_genes,
                 "B_cell" = B_cells_genes,
				 "Plasma_cell" = Plasma_cells_genes,
				 "Myeloid_cell" = Myeloid_cells_genes,
				 "Melanoma_cell" = Melanoma_cells_genes,
				 "Photoreceptor_cell" = Photoreceptor_cell_genes,
				 "Fibroblast" = Fibroblasts_genes,
				 "Endothelial_cell" = Endothelial_cells_genes
)

##set celltype
levels(pbmc)
levels(pbmc@active.ident)

pbmc$celltype2 <- pbmc@active.ident

new.cluster.ids <- c(
							"0"="Malignant tumor cells",
							"1"="Malignant tumor cells",
							"2"="Malignant tumor cells",
							"3"="T/NK cells",
							"4"="Myeloid cells",
							"5"="Malignant tumor cells",
							"6"="Malignant tumor cells",
							"7"="Malignant tumor cells",
							"8"="T/NK cells",
							"9"="Malignant tumor cells",
							"10"="Malignant tumor cells",
							"11"="Malignant tumor cells",
							"12"="Malignant tumor cells",
							"13"="T/NK cells",
							"14"="Myeloid cells",
							"15"="Malignant tumor cells",							
							"16"="Malignant tumor cells",
							"17"="Malignant tumor cells",
							"18"="B/Plasma cells",
							"19"="Malignant tumor cells",
							"20"="Malignant tumor cells",
							"21"="Malignant tumor cells",
							"22"="Malignant tumor cells",
							"23"="Malignant tumor cells",
							"24"="Malignant tumor cells",
							"25"="Malignant tumor cells",
							"26"="Fibroblasts",
							"27"="Malignant tumor cells",
							"28"="Malignant tumor cells",
							"29"="Malignant tumor cells",
							"30"="Malignant tumor cells",
							"31"="Malignant tumor cells",
							"32"="Endothelial cells",
							"33"="Malignant tumor cells",
							"34"="B/Plasma cells",
							"35"="Malignant tumor cells",							
							"36"="Photoreceptor cells",
							"37"="Malignant tumor cells",
							"38"="Malignant tumor cells",
							"39"="Malignant tumor cells"
							)

pbmc <- RenameIdents(pbmc, new.cluster.ids)    

levels(pbmc) 

#Set the display order of cell types
levels(pbmc)  <- c("T/NK cells","B/Plasma cells","Myeloid cells","Malignant tumor cells","Photoreceptor cells","Fibroblasts", "Endothelial cells")


#Graphic modification
Adotplot <- DotPlot(object = pbmc, features=features)&
  theme_bw()& 
  geom_point(shape=21, aes(size=pct.exp),stroke=1)& 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = 'black', size = 12, face = "bold"),
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.margin=unit(c(1, 1, 1, 1),'cm'),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),
        panel.spacing = unit(0.12, "cm"),
        # legend.frame = element_rect(colour = "black"),
        # legend.ticks = element_line(colour = "black", linewidth  = 0),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title = element_text(color = 'black', face = "bold", size=9))& 
#  scale_color_gradientn(colours = colorRampPalette(c("navy","white","firebrick3"))(100))& 
scale_color_gradientn(colours = colorRampPalette(c("white","firebrick3"))(100))& 
  labs(tag = "Cell type annotation")& 
  theme(plot.tag.position = c(0.3, 1.05), 
        plot.tag = element_text(size = 12,face = "bold"))& 
  guides(size=guide_legend(title="Proportion of\nexpressing cells"), 
         colour=guide_colorbar(title="Average\nexpression"))

Adotplot

#45 marker genes 14 cm
ggsave("P2_Marker_gene_Adotplot.pdf",plot=Adotplot,width = 13, height = 3.5)


##Check how many clusters there are
levels(pbmc)

##Set the number of color combinations for the corresponding cluster
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(30)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(12)

show_col(col5)


UMPplot_label <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "celltype", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P3_TSNEplot_Annotated_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "celltype", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(legend.position = "none", panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P3_TSNEplot_Annotated_label_nolegend.pdf", plot = plot_grid(UMPplot_unlabel), width = 8.5, height = 8)


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




#Step 6. Extract specific cell subtypes
pbmc <- readRDS("CCA_Melanoma_GSE139829_GSE138665_Annotated.Rds")
#Extract specific cell subtypes
Cell.sub <- subset(pbmc@meta.data, celltype %in% c("Malignant tumor cells"))
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
scRNAsub = FindClusters(scRNAsub,resolution = 0.3)

##Check how many clusters there are
levels(scRNAsub)

##Set the number of color combinations for the corresponding cluster
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(15)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(37)

show_col(col5)

##Normal UMAP dimensionality reduction data results
UMPplot_scRNAsub_label <- DimPlot(scRNAsub, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_UMPplot_scRNAsub_label.pdf", plot = plot_grid(UMPplot_scRNAsub_label), width = 6.5, height = 6)
UMPplot_scRNAsub_unlabel <- DimPlot(scRNAsub, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_UMPplot_scRNAsub_unlabel.pdf", plot = plot_grid(UMPplot_scRNAsub_unlabel), width = 6.5, height = 6)

TSNEplot_scRNAsub_label <- DimPlot(scRNAsub, reduction = "tsne", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_TSNEplot_scRNAsub_label.pdf", plot = plot_grid(TSNEplot_scRNAsub_label), width = 8.5, height = 7)
TSNEplot_scRNAsub_unlabel <- DimPlot(scRNAsub, reduction = "tsne", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P4_TSNEplot_scRNAsub_unlabel.pdf", plot = plot_grid(TSNEplot_scRNAsub_unlabel), width = 8.5, height = 7)


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

saveRDS(scRNAsub, file="CCA_Melanoma_GSE139829_GSE138665_Annotated_Malignant.Rds") 


##############################Step 7. Addmodulescore
library(Seurat)
library(GSVA)

scRNAsub <- readRDS(file="CCA_Melanoma_GSE139829_GSE138665_Annotated_Malignant.Rds")
pbmc <- scRNAsub

##Select RasFry related genes
RasFry_genes <- c(
#JNK
"HRAS", "KRAS", "NRAS",  "RAF1", "BRAF", "ARAF", "MAP3K1", "MAP3K2", "MAP3K3", "MAP3K5","MAP2K4", "MAP2K7","MAPK8", "MAPK9", "MAPK10","JUN", "FOS", "ATF2", "ELK1", "DUSP1", "DUSP5", "DUSP6",
#NOTCH
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


##Step 8. Export Meta file to draw correlation analysis
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

saveRDS(pbmc, file="CCA_Melanoma_GSE139829_GSE138665_Annotated_Malignant_Scored.Rds")


##Step 9. Export Ras highly activated clusters to draw correlation analysis
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
set.seed(1234)
getwd()

setwd("G:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\05Human Relavance")

pbmc <- readRDS("CCA_Melanoma_GSE139829_GSE138665_Annotated_Malignant_Scored.Rds")

levels(pbmc)
col5 <- colorRampPalette(brewer.pal(12,"Set3"))(37)

VlnPlot(
    pbmc,raster = FALSE, pt.size = 0.1, cols = col5,
    features = c("RAS_score1"), 
    group.by = "seurat_clusters"  
) + theme_classic()+ 
	scale_y_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2)) + 
	theme(legend.position = "none",plot.title = element_text(size=20,hjust = 0.5),  
			axis.line = element_line(color = "black",linewidth = 1),
			text = element_text (size = 18),
			axis.text = element_text (color = "black", size = 18))
ggsave("P10_pbmcVlnPlot_RAS_signaling_score.pdf", width = 16, height = 4)  

VlnPlot(
    pbmc,raster = FALSE, pt.size = 00, cols = col5,
    features = c("RAS_score1"), 
    group.by = "seurat_clusters"  
) + theme_classic()+ 
	scale_y_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2)) + 
	theme(legend.position = "none",plot.title = element_text(size=20,hjust = 0.5),  
			axis.line = element_line(color = "black",linewidth = 1),
			text = element_text (size = 18),
			axis.text = element_text (color = "black", size = 18))
ggsave("P10_pbmcVlnPlot_RAS_signaling_score_no_dots.pdf", width = 16, height = 4)  


VlnPlot(
    pbmc,raster = FALSE, pt.size = 0.1, y.max = 1.5, ncol = 3, cols = col5,
    features = c("RAS_score1", "MAPK_score1", "JNK_score1", "TLR_score1", "NOTCH_score1","HIPPO_inactivation_score1"), 
    group.by = "seurat_clusters",   
    log = TRUE
)
ggsave("P10_pbmcVlnPlot_signaling_score.pdf", width = 30, height = 6)

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


