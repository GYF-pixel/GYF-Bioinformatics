
#Robust linear regression model

#1. Transformation between long and wide data

#load packages
library(reshape2)
library(xlsx)

#work dir
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\Linear model")
list.files()

# Read all numeric data files of phenotypic malignancies

#01Tumor-mediated Cachexia
data <- read.xlsx2(file = "01Tumor-mediated Cachexia_Transparency.xlsx", sheetIndex = 1)
data <- read.xlsx2(file = "01Tumor-mediated Cachexia_Transparent Ratio.xlsx", sheetIndex = 1)
#02Larva Tumor Burden
data <- read.xlsx2(file = "02Larva Tumor Burden_Final_GFP Region.xlsx", sheetIndex = 1)
data <- read.xlsx2(file = "02Larva Tumor Burden_Final_Ratio.xlsx", sheetIndex = 1)
#03Tumor volume
data <- read.xlsx2(file = "03Tumor volume_Final_Disc.xlsx", sheetIndex = 1)
data <- read.xlsx2(file = "03Tumor volume_Final_GFP.xlsx", sheetIndex = 1)

#All columns are measure.vars
melted_data <- melt(data, measure.vars = colnames(data), variable.name = "Variable", value.name = "Value")

#Output long data
write.xlsx(melted_data, file="01Tumor-mediated Cachexia_Transparency_long.xlsx", col.names = T)
write.xlsx(melted_data, file="01Tumor-mediated Cachexia_Transparent Ratio_long.xlsx", col.names = T)

write.xlsx(melted_data, file="02Larva Tumor Burden_Final_GFP Region_long.xlsx", col.names = T)
write.xlsx(melted_data, file="02Larva Tumor Burden_Final_Ratio_long.xlsx", col.names = T)

write.xlsx(melted_data, file="03Tumor volume_Final_Disc_long.xlsx", col.names = T)
write.xlsx(melted_data, file="03Tumor volume_Final_GFP_long.xlsx", col.names = T)



#2. Loading long data for plotting

install.packages("car") #Install the car package to calculate the variance inflation factor (VIF).
install.packages("glmnet") #Install the glmnet package for LASSO regression and ridge regression.
install.packages("ggpmisc") #Insert formula

#load packages
library(car) 
library(reshape2)
library(xlsx)
library(ggplot2)
library(glmnet) 
library(ggpmisc) #Insert regression formula
library(MASS) #use rlm for robust linear regression
library(viridis)
library(RColorBrewer)
library(ggsci)

#work dir
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\Linear model")
list.files()

#file list
name.list <- c("01Tumor-mediated Cachexia_Transparency_long.xlsx",					#scale_y_continuous(limits = c(0, 28), breaks = seq(0, 28, 7))
						"01Tumor-mediated Cachexia_Transparent Ratio_long.xlsx",		#scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 10))
						"02Larva Tumor Burden_Final_GFP Region_long.xlsx",					#scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 4))
						"02Larva Tumor Burden_Final_Ratio_long.xlsx",							#scale_y_continuous(limits = c(0, 57), breaks = seq(0, 57, 10))
						"03Tumor volume_Final_Disc_long.xlsx",										#scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 5))
						"03Tumor volume_Final_GFP_long.xlsx"    									#scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 5))
						)
						
length(name.list)
#Set parameter i, subsequently use name.list[i] to determine the order of data to be read
i = 1			#1 to  6

#load file
data <- read.xlsx2(file = name.list[i], sheetIndex = 1)
data$Time = as.numeric(data$Time)
data$Value = as.numeric(data$Value)
head(data)



#3. Formal Drawing

#Set theme parameters
mytheme <- theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(hjust = 0.5, size = 14, color = 'black'), 
          axis.text.y = element_text(hjust = 0.5, size = 14, color = 'black'),
          axis.title.y = element_text(size = 14), 
          axis.title.x = element_text(size = 14), 
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.position = "right",
          legend.background = element_blank())
		  
#Set color palette
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(20)
set.seed(1234)

#method parameter defaults to the loess method (fitting)
#method = "gam" requires loading the mgcv package and specifying formula = y ~ s(x). For large datasets, use formula = y ~ s(x, bs = "cs")
#method = "lm" performs linear fitting. method = "rlm" is similar to lm but employs a robust fitting algorithm to reduce the influence of outliers; it requires the MASS package.

#3.1 Average: Regression equation based on all Ras-related tumors
ggplt <- ggplot(data = data, 
                mapping = aes(x=Time, y=Value)) +
    #geom_point(size = 2, alpha = 0.5, position = "jitter") + #alpha represents transparency
    geom_jitter(size = 2, alpha = 0.5) + 			#The jitter geom is a convenient shortcut for geom_point(position = "jitter"), the continuous points were separated and some random noise was added.
    geom_smooth(method = 'rlm',					#using rlm
                se=T,												#Add confidence interval, default is T
                lwd=1.5,											#Line width
                #color = "#9f0000",							#Fit curve color
                fill = "lightgrey", 								#Confidence interval color
                fullrange=TRUE) + 							#Show total length
    scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 5))  + 				#Set Y-axis size
    scale_x_continuous(limits = c(3.5, 15.5), breaks = seq(1.5, 17.5, 2)) + 		#Set X-axis size
    mytheme +      
    scale_fill_npg() +									
    scale_color_npg() +							
	#scale_color_manual(values = c('#FF0000','#FF7F00','#FFFF00','#00FF00')) + 	#Adjust the plot by manually adding specific colors
    theme(panel.border = element_rect(color = "black",linewidth = 1)) +						#Add border
	#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 	#Remove gridlines
    stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
                 formula = y ~ x,  parse = TRUE,
                 size = 5, 																									#Formula font size
                 label.x = 0.05, 																							#The position of the formula in the figure, the ratio between 0-1
                 label.y = 0.95)
ggplt 

#save files
file_name2 = paste("P1_rlm_Average",name.list[i],".pdf")
ggsave(filename= file_name2, plot= ggplt, width = 8, height = 4.5)


#3.2 Single: Regression equation based on each Ras-related tumors
ggplt <- ggplot(data = data, 
                mapping = aes(x=Time, y=Value, 
                              color=Genotype, 				#颜色标记
                              fill=Genotype, 			#添加置信区间
                )) + 
	#geom_point(size = 2, alpha = 0.5, position = "jitter") + #alpha represents transparency
    geom_jitter(size = 2, alpha = 0.5) + 			#The jitter geom is a convenient shortcut for geom_point(position = "jitter"), the continuous points were separated and some random noise was added.
    geom_smooth(method = 'rlm',					#using rlm
                se=T,												#Add confidence interval, default is T
                lwd=1.5,											#Line width
                #color = "#9f0000",							#Fit curve color
                fill = "lightgrey", 								#Confidence interval color
                fullrange=TRUE,  							#Show total length
				aes(color=Genotype)) + 
    scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 5))  + 				#Set Y-axis size
    scale_x_continuous(limits = c(3.5, 15.5), breaks = seq(1.5, 17.5, 2)) + 		#Set X-axis size
    mytheme +      
    scale_fill_npg() +									
    scale_color_npg() +							
	#scale_color_manual(values = c('#FF0000','#FF7F00','#FFFF00','#00FF00')) + 	#Adjust the plot by manually adding specific colors
    theme(panel.border = element_rect(color = "black",linewidth = 1)) +						#Add border
	#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 	#Remove gridlines
    stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
                 formula = y ~ x,  parse = TRUE,
                 size = 5, 																									#Formula font size
                 label.x = 0.05,  																							#The position of the formula in the figure, the ratio between 0-1
                 label.y = seq(from = 0.15, to = 0.95, length.out = 10))
ggplt 

#save files
file_name2 = paste("P1_rlm_Single",name.list[i],".pdf")
ggsave(filename= file_name2, plot= ggplt, width = 8, height = 4.5)



#4. PCA analysis

#load package
library(RColorBrewer)
library(ggsci)
library(reshape2)
library(xlsx)
library(pheatmap)
library(FactoMineR)
library(ggplot2)			#ggplot2 draws a two-dimensional scatter plot
library(ggrepel)			#Identify the sample name and use the extension package ggrepel of ggplot2 to complete it

#work dir
setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\Linear model")
list.files()

#load data files
data = read.table("04PCA_Input.txt", header = T) #using Intercept, slope, and AdjR2 of robust linear regression model of each Ras tumor model as input
head(data)			#check data fomat

#                        Ras_emei Ras_fmt Ras_fry Ras_lgl Ras_msn Ras_Rabex Ras_scrib Ras_Syx7 Ras_TSG101 Ras_Vps36
#Cachexia_area_Intercept     2.940   1.940   5.340   5.220  -0.886     7.960      6.46    0.766      3.450     7.350
#Cachexia_area_slope         0.972   0.824   0.498   0.674   1.340     0.373      0.71    1.050      0.603     0.443
#Cachexia_area_AdjR2         0.350   0.450   0.220   0.360   0.600     0.140      0.27    0.540      0.290     0.100
#Cachexia_ratio_Intercept   13.100  16.000  24.600  18.100   3.900    25.100     18.60    8.130     11.900    19.900
#Cachexia_ratio_slope        2.460   1.440   0.917   1.670   3.260     0.965      1.72    2.510      1.940     1.660
#Cachexia_ratio_AdjR2        0.500   0.250   0.130   0.340   0.590     0.180      0.24    0.570      0.330     0.160

row.names(data)
data.filter = data

#4.1 PCA analysis
gene <- t(data.filter)
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
gene.pca 					#View specific data analysis results
plot(gene.pca)  			#PCA diagram


#4.2 contributions of principle components
gene.pca$var$contrib

#save file
a <- gene.pca$var$contrib
file_name = "05PCA contributions of the individuals.csv"
write.csv(a, file_name, row.names = TRUE)

input_heatmap_regroup = a
colnames(input_heatmap_regroup) = c("PC1", "PC2")

#color palette
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

#draw heatmap
heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="",
                 fontsize = 12,
                 scale="none",											#scale argument shoud take values: 'none', 'row' or 'column'
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = F, cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 20,
                 cellheight = 12,cellwidth = 12,
                 cutree_row=3,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap

ggsave("P3_Heatmap_PC_contributions.pdf", plot = heatmap, width = 5, height = 6) 


#4.3 Draw PCA projection in two dimensions
#Extract the coordinates of the sample in the first two axes of PCA
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
pca_sample			#show results

#Extract the contribution of the first two axes of PCA
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )
pca_eig1
pca_eig2

#Add grouping information
colnames(data.filter)
group_info = colnames(data.filter)

group_info2 <- as.data.frame(group_info)
rownames(group_info2) <- colnames(data.filter)
group_info2$samples <- rownames(group_info2)

pca_sample <- cbind(pca_sample, group_info2)
pca_sample  #The graphical data contains sample coordinates and grouping information

#color palette
col5 <- pal_npg(palette = c("nrc"))(10)

#Draw PC1-PC2 two dimensional projection
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
					geom_point(aes(color = group_info), size = 5) +	#Draw a two-dimensional scatter plot based on the sample coordinates
					scale_color_manual(values = col5) +  					#Custom color
					theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',linewidth = 2), legend.key = element_rect(fill = 'transparent'), #Remove background and gridlines, set border thickness
							text = element_text (size = 18),axis.text = element_text (size = 18)) +  #Set font size
					labs(x =  paste('PC1:', pca_eig1, '%'), y = paste('PC2:', pca_eig2, '%'), color = '')  #Add PCA axis contribution to the coordinate axis title
p

#label sample id
p2 <- p + 
geom_text_repel(aes(label = samples), size = 3, show.legend = FALSE, box.padding = unit(0.5, 'lines'))
p2

ggsave("P3_PCA_orignial_label.pdf", plot = p2, width = 7, height = 5)


#4.4 K-means clustering

#load packages
library(cluster)
library(factoextra)

#K-means clustering with 3 centers
result <- kmeans(pca_sample[,1:2], 3)
result

#PCA results + K-means clustering
P <- fviz_cluster(object = result, 
							data = pca_sample[,1:2],
							palette =c( "#E64B35FF","#43CD80", "#2E9FDF" ),
							main = "K-means with k = 3",													#Main title
							repel = TRUE,
							ggtheme = theme_bw() + 														
							theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',linewidth = 2), 
										legend.position = "right",												#The legend is placed at the bottom
										legend.key = element_rect(fill = 'transparent'), 		#Remove background and gridlines, set border thickness
										text = element_text (size = 18),axis.text = element_text (size = 18))
							)

P

#save files
ggsave("P3_PCA_Kmeans_label.pdf", plot = P, width = 6.5, height = 5.5)

