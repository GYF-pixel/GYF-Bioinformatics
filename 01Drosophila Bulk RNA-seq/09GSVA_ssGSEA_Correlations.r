
#GSVA_ssGSEA_Malignancy Correlations

#load packages
library(ggplot2) 
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase) 
library(dplyr) 
library(data.table) 
library(biomaRt)	

#work dir
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\05Analysis per Tumor\\Part 2 GSVA&ssGSEA")
list.files()


#1. Signaling Pathway List
#geneset_name <- c("BMP","EGFR","FGFR","Hedgehog","Hippo","Imd","InR","JAK_STAT","Notch","TNFalpha_Eiger","Toll","Wnt")
geneset_name <- read.table(file = "./00GeneList/FlyBase_IDs_Signaling Pahtway Name.txt",header = T)
rownames(geneset_name) <- geneset_name[,1]
genelist <- split(geneset_name,rownames(geneset_name))  


##1.1 Signaling Pathway All
for (i in 1: nrow(geneset_name)) {
									file_input = 	paste('./00GeneList/FlyBase_IDs', geneset_name[i,], "Signaling Pathway All.txt")
									print(file_input)
									Input = read.table(file = file_input, header = F)
									genelist[[i]] <- Input[,1]

}

##1.2 Signaling Pathway Positive Regulators
for (i in 1: nrow(geneset_name)) {
									file_input = 	paste('./00GeneList/FlyBase_IDs', geneset_name[i,], "Signaling Pathway Positive Regulators.txt")
									print(file_input)
									Input = read.table(file = file_input, header = F)
									genelist[[i]] <- Input[,1]

}

##1.3 Signaling Pathway Negative Regulators
for (i in 1: nrow(geneset_name)) {
									file_input = 	paste('./00GeneList/FlyBase_IDs', geneset_name[i,], "Signaling Pathway Negative Regulators.txt")
									print(file_input)
									Input = read.table(file = file_input, header = F)
									genelist[[i]] <- Input[,1]

}


#2. Reading Expression Profile reads the original gene Count file
#data = read.csv(file = "newdata_filter.csv", header = T,row.names = 1)
data = read.csv(file = "newdata_filter_remove_pre_batheffect_removed.csv", header = T,row.names = 1)

dat <- log2(edgeR::cpm(data)+1)
dat[1:4,1:4]
dim(dat)


#3. GSVA
#Standardized data uses Gaussian distribution, i.e. log2 transformed data # tpm after log uses Gaussian distribution
#Use Poisson distribution for count values
#If it is uncertain whether gene set S will be enriched at the top or bottom, using method two may not be appropriate. 
#In the GSVA package, select methods through the parameter mix.diff. If mix.diff=T, choose method two, and if mix.diff=F, choose method one.

expr = dat
genesets4gsva = genelist
expr_geneset <- gsva(expr = as.matrix(expr), 			#It cannot be data.frame
                     gset.idx.list = genesets4gsva,
                     method="gsva",
                     kcdf="Gaussian", 								#Use Gaussian distribution for TPM after logging
                     parallel.sz=10 									#Multi threading
                     )
write.csv(expr_geneset, "GSVA_All Samples All.csv", row.names = TRUE)
write.csv(expr_geneset, "GSVA_All Samples Positive Regulators.csv", row.names = TRUE)
write.csv(expr_geneset, "GSVA_All Samples Negative Regulators.csv", row.names = TRUE)

#4. ssGSEA

expr_geneset_ssGSEA <- gsva(expr = as.matrix(expr), 	#It cannot be data.frame
                     gset.idx.list = genesets4gsva,
                     method="ssgsea",
                     parallel.sz=10 											#Multi threading
                     )
write.csv(expr_geneset_ssGSEA, "ssGSEA_All Samples All.csv", row.names = TRUE)
write.csv(expr_geneset_ssGSEA, "ssGSEA_All Samples Positive Regulators.csv", row.names = TRUE)
write.csv(expr_geneset_ssGSEA, "ssGSEA_All Samples Negative Regulators.csv", row.names = TRUE)



#5.  For the for loop of Step 3&Step 4 for per Tumor
getwd()
countfile_name <- c("3-199","4-41","4-47","4-90","6-119","6-128","6-144","6-147","7-51","11-231")
length(countfile_name)
genesets4gsva = genelist

for (i in 1: length(countfile_name)) {
						countfile = paste('./00GeneCountFile/newdata_filter_batcheffect_removed', countfile_name[i], ".csv")
						data = read.csv(file = countfile, header = T,row.names = 1)
						dat <- log2(edgeR::cpm(data)+1)
						dat[1:4,1:4]
						dim(dat)
						expr = dat

#GSVA scoring
expr_geneset <- gsva(expr = as.matrix(expr), 			#It cannot be data.frame
                     gset.idx.list = genesets4gsva,
                     method="gsva",
                     kcdf="Gaussian",								#Use Gaussian distribution for TPM after logging
                     parallel.sz=10 									#Multi threading
                     )

  savingfile_name <- paste('GSVA_',countfile_name[i],"Positive Regulators.csv")
  write.csv(expr_geneset, savingfile_name, row.names = TRUE)
  cat("File", savingfile_name, "saved.\n")    				###check save

#ssGSEA	scoring
expr_geneset_ssGSEA <- gsva(expr = as.matrix(expr), 	#It cannot be data.frame
                     gset.idx.list = genesets4gsva,
                     method="ssgsea",
                     parallel.sz=10 									#Multi threading
                     )

  savingfile_name_ssGSEA <- paste('ssGSEA_',countfile_name[i],"Positive Regulators.csv")
  write.csv(expr_geneset_ssGSEA, savingfile_name_ssGSEA, row.names = TRUE)
  cat("File", savingfile_name_ssGSEA, "saved.\n")    ###check save
}





#6. GSVA heatmap for all samples

#load packages
library(ggplot2) 
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase) 
library(dplyr) 
library(data.table) 
library(biomaRt)	
library(curl)	
library(pheatmap)
library(RColorBrewer)

#work dir 
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\05Analysis per Tumor\\Part 2 GSVA&ssGSEA")
list.files()

#color palette
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

#Retrieve GSVA score
normalizeExp=read.csv(file="GSVA_All Samples All.csv",header=T,row.names = 1)
normalizeExp=read.csv(file="GSVA_All Samples Positive Regulators.csv",header=T,row.names = 1)
normalizeExp=read.csv(file="GSVA_All Samples Negative Regulators.csv",header=T,row.names = 1)
head(normalizeExp)

colnames(normalizeExp)


input_heatmap_regroup = normalizeExp[,c(
"WT_79EA","WT_79EB","WT_79EC",

"Ras1A","Ras2A","Ras3A",
"X4.03A","X4.03B","X4.03C",
"Ras82BA","Ras82BB","Ras82BC",

"X4.41A1A","X4.41A2A","X4.41A3A",
"X11.231A1A","X11.231A2A","X11.231A3A",
"X6.128A1A","X6.128A2A","X6.128A3A",
"X6.144A1A","X6.144A2A","X6.144A3A",
"X6.147A1A","X6.147A2A","X6.147A3A",
"X7.51A1A","X7.51A2A","X7.51A3A",
"X6.119A1A","X6.119A2A","X6.119A3A",
"X4.47A1A","X4.47A2A","X4.47A3A",
"X4.90A1A","X4.90A2A","X4.90A3A",
"X3.199A1A","X3.199A2A","X3.199A3A",

"X4.41B1A","X4.41B2A","X4.41B3A",
"X11.231B1A","X11.231B2A","X11.231B3A",
"X6.128B1A","X6.128B2A","X6.128B3A",
"X6.144B1A","X6.144B2A","X6.144B3A",
"X6.147B1A","X6.147B2A","X6.147B3A",  
"X7.51B1A","X7.51B2A","X7.51B3A",
"X6.119B1A","X6.119B2A","X6.119B3A",
"X4.47B1A","X4.47B2A","X4.47B3A",
"X4.90B1A","X4.90B2A","X4.90B3A",
"X3.199B1A","X3.199B2A","X3.199B3A",  


"X4.41C1A","X4.41C2A","X4.41C3A",
"X11.231C1A","X11.231C2A","X11.231C3A",
"X6.128C1A","X6.128C2A","X6.128C3A",
"X6.144C1A","X6.144C2A","X6.144C3A",
"X6.147C1A","X6.147C2A","X6.147C3A",
"X7.51C1A","X7.51C2A","X7.51C3A",
"X6.119C1A","X6.119C2A","X6.119C3A",
"X4.47C1A","X4.47C2A","X4.47C3A",
"X4.90C1A","X4.90C2A","X4.90C3A",
"X3.199C1A","X3.199C2A","X3.199C3A"
)]

group_info = c(rep("WT",3), rep("40ARas",3), rep("Ras79E",3), rep("Ras82B",3), 
						rep("5d_4-41",3), rep("5d_11-231",3), rep("5d_6-128",3),rep("5d_6-144",3), rep("5d_6-147",3), rep("5d_7-51",3), rep("5d_6-119",3), rep("5d_4-47",3), rep("5d_4-90",3), rep("5d_3-199",3),
						rep("9d_4-41",3), rep("9d_11-231",3), rep("9d_6-128",3),rep("9d_6-144",3), rep("9d_6-147",3), rep("9d_7-51",3), rep("9d_6-119",3), rep("9d_4-47",3), rep("9d_4-90",3), rep("9d_3-199",3),
						rep("13d_4-41",3), rep("13d_11-231",3), rep("13d_6-128",3),rep("13d_6-144",3), rep("13d_6-147",3), rep("13d_7-51",3), rep("13d_6-119",3), rep("13d_4-47",3), rep("13d_4-90",3), rep("13d_3-199",3)
						)


#Take the average of the same repeated samples
input_heatmap_regroup_avereps_df  <- t(limma::avereps( t(input_heatmap_regroup) , ID = group_info))

#input_heatmap_regroup
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="",
                 fontsize = 15,
                 scale="none",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 6,cellwidth = 8,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap
ggsave("Part3_Heatmap_genes_negative.pdf", plot = heatmap, width = 18, height = 8) 



#input_heatmap_regroup_avereps_df
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
heatmap=pheatmap(input_heatmap_regroup_avereps_df,color = colors,
                 main="",
                 fontsize = 12,
                 scale="row",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 12,cellwidth = 18,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap

ggsave("Part3_Heatmap_genes_regroup_avereps_df_negative.pdf", plot = heatmap, width = 15, height = 6) 



#7. GSVA heatmap for each samples

#load packages
library(ggplot2) 
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase) 
library(dplyr) 
library(data.table) 
library(biomaRt)	
library(curl)	
library(pheatmap)
library(RColorBrewer)

#work dir 
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\05Analysis per Tumor\\Part 2 GSVA&ssGSEA")
list.files()

#color palette
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

getwd()
countfile_name <- c("3-199","4-41","4-47","4-90","6-119","6-128","6-144","6-147","7-51","11-231")

for (i in 1: length(countfile_name)) {
						#调取GSVA score
						GSVA_file = paste('./GSVA_', countfile_name[i], "Positive Regulators.csv")
						normalizeExp = read.csv(file = GSVA_file, header = T,row.names = 1)
						print(colnames(normalizeExp))
						input_heatmap_regroup = normalizeExp
###绘图input_heatmap_regroup
heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="",
                 fontsize = 12,
                 scale="none",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 15,cellwidth = 18,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap
savingfile_name_GSVA_heatmap <- paste('GSVA_',countfile_name[i],"Positive Regulators_Heatmap.pdf")
ggsave(filename = savingfile_name_GSVA_heatmap, plot = heatmap, width = 15, height = 6) 

						
						
###绘图input_heatmap_regroup_avereps_df
						group_info = c(rep("WT",3), rep("Ras",3), rep("5d_T",3), rep("9d_T",3), rep("13d_T",3))
						input_heatmap_regroup_avereps_df  <- t(limma::avereps( t(input_heatmap_regroup) , ID = group_info))##对相同重复样品取平均
heatmap2=pheatmap(input_heatmap_regroup_avereps_df,color = colors,
                 main="",
                 fontsize = 12,
                 scale="none",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 15,cellwidth = 18,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T, angle_col = c("90")
)
heatmap2
savingfile_name_GSVA_heatmap_avereps <- paste('GSVA_',countfile_name[i],"Positive Regulators_Heatmap_avereps.pdf")
ggsave(filename = savingfile_name_GSVA_heatmap_avereps, plot = heatmap2, width = 15, height = 6) 

}


#8. ssGSEA heatmap for all samples

#load packages
library(ggplot2) 
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase) 
library(dplyr) 
library(data.table) 
library(biomaRt)	
library(curl)	
library(pheatmap)
library(RColorBrewer)

#work dir 
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\05Analysis per Tumor\\Part 2 GSVA&ssGSEA")
list.files()

#color palette
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

#调取ssGSEA score
normalizeExp=read.csv(file="ssGSEA_All Samples All.csv",header=T,row.names = 1)
normalizeExp=read.csv(file="ssGSEA_All Samples Positive Regulators.csv",header=T,row.names = 1)
normalizeExp=read.csv(file="ssGSEA_All Samples Negative Regulators.csv",header=T,row.names = 1)
head(normalizeExp)

colnames(normalizeExp)

#Sample order
input_heatmap_regroup = normalizeExp[,c(
"WT_79EA","WT_79EB","WT_79EC",

"Ras1A","Ras2A","Ras3A",
"X4.03A","X4.03B","X4.03C",
"Ras82BA","Ras82BB","Ras82BC",

"X4.41A1A","X4.41A2A","X4.41A3A",
"X11.231A1A","X11.231A2A","X11.231A3A",
"X6.128A1A","X6.128A2A","X6.128A3A",
"X6.144A1A","X6.144A2A","X6.144A3A",
"X6.147A1A","X6.147A2A","X6.147A3A",
"X7.51A1A","X7.51A2A","X7.51A3A",
"X6.119A1A","X6.119A2A","X6.119A3A",
"X4.47A1A","X4.47A2A","X4.47A3A",
"X4.90A1A","X4.90A2A","X4.90A3A",
"X3.199A1A","X3.199A2A","X3.199A3A",

"X4.41B1A","X4.41B2A","X4.41B3A",
"X11.231B1A","X11.231B2A","X11.231B3A",
"X6.128B1A","X6.128B2A","X6.128B3A",
"X6.144B1A","X6.144B2A","X6.144B3A",
"X6.147B1A","X6.147B2A","X6.147B3A",  
"X7.51B1A","X7.51B2A","X7.51B3A",
"X6.119B1A","X6.119B2A","X6.119B3A",
"X4.47B1A","X4.47B2A","X4.47B3A",
"X4.90B1A","X4.90B2A","X4.90B3A",
"X3.199B1A","X3.199B2A","X3.199B3A",  


"X4.41C1A","X4.41C2A","X4.41C3A",
"X11.231C1A","X11.231C2A","X11.231C3A",
"X6.128C1A","X6.128C2A","X6.128C3A",
"X6.144C1A","X6.144C2A","X6.144C3A",
"X6.147C1A","X6.147C2A","X6.147C3A",
"X7.51C1A","X7.51C2A","X7.51C3A",
"X6.119C1A","X6.119C2A","X6.119C3A",
"X4.47C1A","X4.47C2A","X4.47C3A",
"X4.90C1A","X4.90C2A","X4.90C3A",
"X3.199C1A","X3.199C2A","X3.199C3A"
)]

group_info = c(rep("WT",3), rep("40ARas",3), rep("Ras79E",3), rep("Ras82B",3), 
						rep("5d_4-41",3), rep("5d_11-231",3), rep("5d_6-128",3),rep("5d_6-144",3), rep("5d_6-147",3), rep("5d_7-51",3), rep("5d_6-119",3), rep("5d_4-47",3), rep("5d_4-90",3), rep("5d_3-199",3),
						rep("9d_4-41",3), rep("9d_11-231",3), rep("9d_6-128",3),rep("9d_6-144",3), rep("9d_6-147",3), rep("9d_7-51",3), rep("9d_6-119",3), rep("9d_4-47",3), rep("9d_4-90",3), rep("9d_3-199",3),
						rep("13d_4-41",3), rep("13d_11-231",3), rep("13d_6-128",3),rep("13d_6-144",3), rep("13d_6-147",3), rep("13d_7-51",3), rep("13d_6-119",3), rep("13d_4-47",3), rep("13d_4-90",3), rep("13d_3-199",3)
						)

input_heatmap_regroup_avereps_df  <- t(limma::avereps( t(input_heatmap_regroup) , ID = group_info))


#input_heatmap_regroup
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="",
                 fontsize = 15,
                 scale="none",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 6,cellwidth = 8,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap
ggsave("Part3_Heatmap_genes_ssGSEA.pdf", plot = heatmap, width = 18, height = 8) 

#input_heatmap_regroup_avereps_df
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
heatmap=pheatmap(input_heatmap_regroup_avereps_df,color = colors,
                 main="",
                 fontsize = 12,
                 scale="none",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 12,cellwidth = 18,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap

ggsave("Part3_Heatmap_genes_regroup_avereps_df_ssGSEA.pdf", plot = heatmap, width = 15, height = 6) 




#9. ssGSEA heatmap for each sample

#load packages
library(ggplot2) 
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase) 
library(dplyr) 
library(data.table) 
library(biomaRt)	
library(curl)	
library(pheatmap)
library(RColorBrewer)

#work dir 
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\05Analysis per Tumor\\Part 2 GSVA&ssGSEA")
list.files()

#color palette
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

getwd()
countfile_name <- c("3-199","4-41","4-47","4-90","6-119","6-128","6-144","6-147","7-51","11-231")


for (i in 1: length(countfile_name)) {
						#调取ssGSEA score
						ssGSEA_file = paste('./ssGSEA_', countfile_name[i], "Positive Regulators.csv")
						normalizeExp = read.csv(file = ssGSEA_file, header = T,row.names = 1)
						print(colnames(normalizeExp))
						input_heatmap_regroup = normalizeExp
###绘图input_heatmap_regroup
heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="",
                 fontsize = 12,
                 scale="row",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 15,cellwidth = 18,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap
savingfile_name_ssGSEA_heatmap <- paste('ssGSEA_',countfile_name[i],"Positive Regulators_Heatmap_scalebyrow.pdf")
ggsave(filename = savingfile_name_ssGSEA_heatmap, plot = heatmap, width = 15, height = 6) 

						
						
###绘图input_heatmap_regroup_avereps_df
						group_info = c(rep("WT",3), rep("Ras",3), rep("5d_T",3), rep("9d_T",3), rep("13d_T",3))
						input_heatmap_regroup_avereps_df  <- t(limma::avereps( t(input_heatmap_regroup) , ID = group_info))##对相同重复样品取平均
heatmap2=pheatmap(input_heatmap_regroup_avereps_df,color = colors,
                 main="",
                 fontsize = 12,
                 scale="row",				#No Scale  'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 15,cellwidth = 18,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T, angle_col = c("90")
)
heatmap2
savingfile_name_ssGSEA_heatmap_avereps <- paste('ssGSEA_',countfile_name[i],"Positive Regulators_Heatmap_avereps_scalebyrow.pdf")
ggsave(filename = savingfile_name_ssGSEA_heatmap_avereps, plot = heatmap2, width = 15, height = 6) 

}



#10. The correlations of signaling pathway activitiy 

#load packages
library(ggplot2) 
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase) 
library(dplyr) 
library(data.table) 
library(biomaRt)	
library(curl)	
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(limma)
library(ggplot2)
library(ggprism)
library(corrplot)
library(psych)
library(reshape2)

#work dir
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\01Multitumor RNA-seq\\05Analysis per Tumor\\Part 2 GSVA&ssGSEA\\03GSVA&ssGSEA Correlated with Phenotypic data")
list.files()


####10.1 GSVA cor WCGCNA Traits Data
getwd()
countfile_name <- c("3-199","4-41","4-47","4-90","6-119","6-128","6-144","6-147","7-51","11-231")

for (i in 1: length(countfile_name)) {
						#Load phenotype data
						WGCNA_Traits_input = read.csv(file = "01 WGCNA Traits Input.csv", header = T,row.names = 1)	
						
						WGCNA_file = paste('./WGCNA_', countfile_name[i], "All.csv")
						WGCNA_input = read.csv(file = WGCNA_file, header = T,row.names = 1)	
						print(colnames(WGCNA_input))
						
						WGCNA_input2 =  WGCNA_Traits_input[c("Transparency Ratio Summary","Transparent Area Summary",
																							  "GFP Region to Whole Larva Body Ratio","GFP Region Size",
																							  "Disc Tissue Volume","GFP Region Volume",
																							  "VNC Invasion Stage1","VNC Invasion Stage2","VNC Invasion Stage3","VNC Invasion Stage4",
																							  "VNC Invasion Stage1(%)","VNC Invasion Stage2(%)","VNC Invasion Stage3(%)","VNC Invasion Stage4(%)"),c(colnames(WGCNA_input))] #只取前17行
						##Save the corresponding phenotype file
						savingfile_name <- paste('WGCNA_',countfile_name[i],"All.csv")
						write.csv(WGCNA_input2, savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ### check save
						
						#Load GSVA score
						GSVA_file = paste('./GSVA_', countfile_name[i], "All.csv")
						GSVA_input = read.csv(file = GSVA_file, header = T,row.names = 1)
						print(colnames(GSVA_input))

						#Correlation analysis
						#If the adjust parameter is not specified, the function will default to performing Holm correction
						#Many studies do not add this parameter and mistakenly believe that it has not been corrected, which eventually leads to problems
						#corr_matrix <- corr.test(x, y, method = 'spearman')
						#If the p-value is not corrected, please specify adjust='none '
						#corr_matrix <- corr.test(x, y, method = 'spearman', adjust = 'none')
						#If you want to correct the p-value, please specify the specific parameters of adjust, such as Benjamini correction
						corr_matrix <- corr.test(t(GSVA_input), t(WGCNA_input2), method = 'spearman', adjust = 'BH')
						
						##save R value 
						savingfile_name <- paste('GSVA_cor_WGCNA_',countfile_name[i],"R value.csv")
						write.csv(corr_matrix[["r"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						##save P value 
						savingfile_name <- paste('GSVA_cor_WGCNA_',countfile_name[i],"P value.csv")
						write.csv(corr_matrix[["p"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						##save P.adjust value
						savingfile_name <- paste('GSVA_cor_WGCNA_',countfile_name[i],"P.adjust value.csv")
						write.csv(corr_matrix[["p.adj"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")   ###check save
						
						##Organize the data format into bubble heat maps for easy drawing
						#Use melt() to convert data boxes from wide format to long format 
						
						##Convert P value to long_fata
						long_df <- melt(corr_matrix[["p"]])
						colnames(long_df) <- c("SIgnaling Pahtway","Phenotypic Data","P value")
						##save P value 
						savingfile_name <- paste('GSVA_cor_WGCNA_',countfile_name[i],"P value_long.csv")
						write.csv(long_df, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")   ###check save
						
						##transfer R value to long_format						
						long_df2 <- melt(corr_matrix[["r"]])
						colnames(long_df2) <- c("SIgnaling Pahtway","Phenotypic Data","R value")
						
						##save R value in long_format	
						savingfile_name <- paste('GSVA_cor_WGCNA_',countfile_name[i],"R value_long.csv")
						write.csv(long_df2, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##Merge two files
						long_df3 <- cbind(long_df,long_df2)  		#The rows of matrices a and b in cbind (a, b) must be the same and sorted in the same order
						long_df3 <- long_df3[,c("SIgnaling Pahtway","Phenotypic Data","P value","R value")]
						long_df3$log10Pvalue = -log10(long_df3[,"P value"])
						savingfile_name <- paste('GSVA_cor_WGCNA_',countfile_name[i],"R value&P value_long.csv")
						write.csv(long_df3, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
}



####11.2 ssGSEA cor WCGCNA Traits Data
getwd()
countfile_name <- c("3-199","4-41","4-47","4-90","6-119","6-128","6-144","6-147","7-51","11-231")

for (i in 1: length(countfile_name)) {
						#Load phenotype data
						WGCNA_Traits_input = read.csv(file = "01 WGCNA Traits Input.csv", header = T,row.names = 1)	
						
						WGCNA_file = paste('./WGCNA_', countfile_name[i], "All.csv")
						WGCNA_input = read.csv(file = WGCNA_file, header = T,row.names = 1)	
						print(colnames(WGCNA_input))
						
						WGCNA_input2 =  WGCNA_Traits_input[c("Transparency Ratio Summary","Transparent Area Summary",
																							  "GFP Region to Whole Larva Body Ratio","GFP Region Size",
																							  "Disc Tissue Volume","GFP Region Volume",
																							  "VNC Invasion Stage1","VNC Invasion Stage2","VNC Invasion Stage3","VNC Invasion Stage4",
																							  "VNC Invasion Stage1(%)","VNC Invasion Stage2(%)","VNC Invasion Stage3(%)","VNC Invasion Stage4(%)"),c(colnames(WGCNA_input))] #只取前17行
						##Save the corresponding phenotype file
						savingfile_name <- paste('WGCNA_',countfile_name[i],"All.csv")
						write.csv(WGCNA_input2, savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						#Load GSVA score
						ssGSEA_file = paste('./ssGSEA_', countfile_name[i], "All.csv")
						ssGSEA_input = read.csv(file = ssGSEA_file, header = T,row.names = 1)
						print(colnames(ssGSEA_input))

						#Correlation analysis
						#If the adjust parameter is not specified, the function will default to performing Holm correction
						#Many studies do not add this parameter and mistakenly believe that it has not been corrected, which eventually leads to problems
						#corr_matrix <- corr.test(x, y, method = 'spearman')
						#If the p-value is not corrected, please specify adjust='none '
						#corr_matrix <- corr.test(x, y, method = 'spearman', adjust = 'none')
						#If you want to correct the p-value, please specify the specific parameters of adjust, such as Benjamini correction
						corr_matrix <- corr.test(t(ssGSEA_input), t(WGCNA_input2), method = 'spearman', adjust = 'BH')
						
						##save R value 
						savingfile_name <- paste('ssGSEA_cor_WGCNA_',countfile_name[i],"R value.csv")
						write.csv(corr_matrix[["r"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						##save P value 
						savingfile_name <- paste('ssGSEA_cor_WGCNA_',countfile_name[i],"P value.csv")
						write.csv(corr_matrix[["p"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n") 	###check save
						##save P.adjust value 
						savingfile_name <- paste('ssGSEA_cor_WGCNA_',countfile_name[i],"P.adjust value.csv")
						write.csv(corr_matrix[["p.adj"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##Organize the data format into bubble heat maps for easy drawing
						#Use melt() to convert data boxes from wide format to long format 
						
						##Convert P value to long_fata
						long_df <- melt(corr_matrix[["p"]])
						colnames(long_df) <- c("SIgnaling Pahtway","Phenotypic Data","P value")
						
						##save P value 
						savingfile_name <- paste('ssGSEA_cor_WGCNA_',countfile_name[i],"P value_long.csv")
						write.csv(long_df, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##transfer R value to long_format						
						long_df2 <- melt(corr_matrix[["r"]])
						colnames(long_df2) <- c("SIgnaling Pahtway","Phenotypic Data","R value")
						
						##save R value in long_format	
						savingfile_name <- paste('ssGSEA_cor_WGCNA_',countfile_name[i],"R value_long.csv")
						write.csv(long_df2, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##Merge two files
						long_df3 <- cbind(long_df,long_df2)  		#The rows of matrices a and b in cbind (a, b) must be the same and sorted in the same order
						long_df3 <- long_df3[,c("SIgnaling Pahtway","Phenotypic Data","P value","R value")]
						long_df3$log10Pvalue = -log10(long_df3[,"P value"])
						savingfile_name <- paste('ssGSEA_cor_WGCNA_',countfile_name[i],"R value&P value_long.csv")
						write.csv(long_df3, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")   ###check save

}



####11.3 ssGSEA_Positive cor WCGCNA Traits Data

getwd()
countfile_name <- c("3-199","4-41","4-47","4-90","6-119","6-128","6-144","6-147","7-51","11-231")

for (i in 1: length(countfile_name)) {
						#Load phenotype data
						WGCNA_Traits_input = read.csv(file = "01 WGCNA Traits Input.csv", header = T,row.names = 1)	
						
						WGCNA_file = paste('./WGCNA_', countfile_name[i], "All.csv")
						WGCNA_input = read.csv(file = WGCNA_file, header = T,row.names = 1)	
						print(colnames(WGCNA_input))
						
						WGCNA_input2 =  WGCNA_Traits_input[c("Transparency Ratio Summary","Transparent Area Summary",
																							  "GFP Region to Whole Larva Body Ratio","GFP Region Size",
																							  "Disc Tissue Volume","GFP Region Volume",
																							  "VNC Invasion Stage1","VNC Invasion Stage2","VNC Invasion Stage3","VNC Invasion Stage4",
																							  "VNC Invasion Stage1(%)","VNC Invasion Stage2(%)","VNC Invasion Stage3(%)","VNC Invasion Stage4(%)"),c(colnames(WGCNA_input))] #只取前17行
						##Save the corresponding phenotype file
						savingfile_name <- paste('WGCNA_',countfile_name[i],"All.csv")
						write.csv(WGCNA_input2, savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						#load GSVA score
						ssGSEA_file = paste('./ssGSEA_', countfile_name[i], "Positive Regulators.csv")
						ssGSEA_input = read.csv(file = ssGSEA_file, header = T,row.names = 1)
						print(colnames(ssGSEA_input))

						#Correlation analysis
						#If the adjust parameter is not specified, the function will default to performing Holm correction
						#Many studies do not add this parameter and mistakenly believe that it has not been corrected, which eventually leads to problems
						#corr_matrix <- corr.test(x, y, method = 'spearman')
						#If the p-value is not corrected, please specify adjust='none '
						#corr_matrix <- corr.test(x, y, method = 'spearman', adjust = 'none')
						#If you want to correct the p-value, please specify the specific parameters of adjust, such as Benjamini correction
						corr_matrix <- corr.test(t(ssGSEA_input), t(WGCNA_input2), method = 'spearman', adjust = 'BH')
						
						##save R value
						savingfile_name <- paste('ssGSEA_Pos_cor_WGCNA_',countfile_name[i],"R value.csv")
						write.csv(corr_matrix[["r"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##save P value
						savingfile_name <- paste('ssGSEA_Pos_cor_WGCNA_',countfile_name[i],"P value.csv")
						write.csv(corr_matrix[["p"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##save P.adjust value
						savingfile_name <- paste('ssGSEA_Pos_cor_WGCNA_',countfile_name[i],"P.adjust value.csv")
						write.csv(corr_matrix[["p.adj"]], savingfile_name, row.names = TRUE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##Organize the data format into bubble heat maps for easy drawing
						#Use melt() to convert data boxes from wide format to long format 
						
						##Convert P value to long_fata
						long_df <- melt(corr_matrix[["p"]])
						colnames(long_df) <- c("SIgnaling Pahtway","Phenotypic Data","P value")
						
						##save P value 
						savingfile_name <- paste('ssGSEA_Pos_cor_WGCNA_',countfile_name[i],"P value_long.csv")
						write.csv(long_df, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")   ###check save
						
						##transfer R value in long_format							
						long_df2 <- melt(corr_matrix[["r"]])
						colnames(long_df2) <- c("SIgnaling Pahtway","Phenotypic Data","R value")
						
						##save R value in long_format	
						savingfile_name <- paste('ssGSEA_Pos_cor_WGCNA_',countfile_name[i],"R value_long.csv")
						write.csv(long_df2, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")    ###check save
						
						##Merge two files
						long_df3 <- cbind(long_df,long_df2)  		#The rows of matrices a and b in cbind (a, b) must be the same and sorted in the same order
						long_df3 <- long_df3[,c("SIgnaling Pahtway","Phenotypic Data","P value","R value")]
						long_df3$log10Pvalue = -log10(long_df3[,"P value"])
						savingfile_name <- paste('ssGSEA_Pos_cor_WGCNA_',countfile_name[i],"R_P_value_long.csv")    		#CSV files have limit length of file name
						write.csv(long_df3, savingfile_name, row.names = FALSE)
						cat("File", savingfile_name, "saved.\n")    ###check save

}


