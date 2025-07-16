#Mime
https://github.com/l-magnificence/Mime

# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

if (!requireNamespace("fastAdaboost", quietly = TRUE))
  devtools::install_github("souravc83/fastAdaboost")

if (!requireNamespace("Mime", quietly = TRUE))
  devtools::install_github("l-magnificence/Mime")
  


##Step 1. Data downloading
https://github.com/l-magnificence/Mime

#Data Downloading https://xenabrowser.net/datapages/

#COADREAD
https://xenabrowser.net/datapages/?cohort=TCGA%20Colon%20and%20Rectal%20Cancer%20(COADREAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

https://xenabrowser.net/datapages/?cohort=TCGA%20Colon%20Cancer%20(COAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
https://xenabrowser.net/datapages/?cohort=TCGA%20Rectal%20Cancer%20(READ)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

#LUADLUSC
https://xenabrowser.net/datapages/?cohort=TCGA%20Lung%20Adenocarcinoma%20(LUAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
https://xenabrowser.net/datapages/?cohort=TCGA%20Lung%20Squamous%20Cell%20Carcinoma%20(LUSC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

#PAAD
https://xenabrowser.net/datapages/?cohort=TCGA%20Pancreatic%20Cancer%20(PAAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

#SKCM
https://xenabrowser.net/datapages/?cohort=TCGA%20Melanoma%20(SKCM)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443


##Step 2. Prepare for input data

library(data.table)  
library(dplyr)       
library(tidyverse)   
library(limma)
library(caret)
library(Mime1)

#setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\08Machine Learning")
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\08Machine Learning")
list.files()

#Step 2.1 Read RNAseq datasets

#表达矩阵读取
#readCount <- fread("./Oringinal data/TCGA-LUAD.star_tpm.tsv.gz", header = T, sep = '\t', data.table = F)
#readCount <- fread("./Oringinal data/TCGA-LUSC.star_tpm.tsv.gz", header = T, sep = '\t', data.table = F)
readCount <- fread("./Oringinal data/TCGA-COAD.star_tpm.tsv.gz", header = T, sep = '\t', data.table = F)
#readCount <- fread("./Oringinal data/TCGA-READ.star_tpm.tsv.gz", header = T, sep = '\t', data.table = F)
#readCount <- fread("./Oringinal data/TCGA-PAAD.star_tpm.tsv.gz", header = T, sep = '\t', data.table = F)
#readCount <- fread("./Oringinal data/TCGA-SKCM.star_tpm.tsv.gz", header = T, sep = '\t', data.table = F)

#ID Annotation
ID <- fread("./Oringinal data/gencode.v36.annotation.gtf.gene.probemap", header = T, sep = '\t', data.table = F)
#Gene Annotation
TCGA_gset <- readCount %>%
  inner_join(ID, by = c("Ensembl_ID" = "id")) %>%
  dplyr::select(gene, starts_with("TCGA") )
#Remove repeat genes & Average gene expression
TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$gene))
#Save files
#write.csv(TCGA_gset,'TCGA_LUAD_TPM.csv')
#write.csv(TCGA_gset,'TCGA_LUSC_TPM.csv')
#write.csv(TCGA_gset,'TCGA_COAD_TPM.csv')
#write.csv(TCGA_gset,'TCGA_READ_TPM.csv')
#write.csv(TCGA_gset,'TCGA_PAAD_TPM.csv')
#write.csv(TCGA_gset,'TCGA_SKCM_TPM.csv')

#Step 2.2 Read Survival datasets
#readSurvival <- read.table(file="./Oringinal data/TCGA-LUAD.survival.tsv", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
#readSurvival <- read.table(file="./Oringinal data/TCGA-LUSC.survival.tsv", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
readSurvival <- read.table(file="./Oringinal data/TCGA-COAD.survival.tsv", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
#readSurvival <- read.table(file="./Oringinal data/TCGA-READ.survival.tsv", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
#readSurvival <- read.table(file="./Oringinal data/TCGA-PAAD.survival.tsv", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
#readSurvival <- read.table(file="./Oringinal data/TCGA-SKCM.survival.tsv", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
readSurvival$ID = rownames(readSurvival)

#Step 2.3 Merge two information
readCount <- t(TCGA_gset)
final_samples <- intersect(rownames(readCount), rownames(readSurvival))

readCount_Survival = cbind(readSurvival[final_samples, c("ID","OS.time","OS")], readCount[final_samples, ])

#rename ID
readCount_Survival[,"ID"] <- substring(readCount_Survival[,"ID"],1,12) %>% gsub("-",".",.)


#Step 2.4 Seperate Training datasets and Test datasets
inTrain <- createDataPartition(y = readCount_Survival[,3], p = 0.7, list = F)		#Using survival time to sepereate data

train <- readCount_Survival[inTrain, ]
train = as.data.frame(train)
rownames(train) = make.names(rownames(train))

test <- readCount_Survival[-inTrain, ]
test = as.data.frame(test)
rownames(test) = make.names(rownames(test))

#Step 2.5 Save the datasets
input_data <- list(object1 = train, object2 = test)

#save(input_data, file = "./Prepared data/TCGA-LUAD_MLinput.RData")
#save(input_data, file = "./Prepared data/TCGA-LUSC_MLinput.RData")
save(input_data, file = "./Prepared data/TCGA-COAD_MLinput.RData")
#save(input_data, file = "./Prepared data/TCGA-READ_MLinput.RData")
#save(input_data, file = "./Prepared data/TCGA-PAAD_MLinput.RData")
#save(input_data, file = "./Prepared data/TCGA-SKCM_MLinput.RData")


#Step 3 Machine Learning
#BiocManager::install("snowfall")
#install.packages("quadprog")

#load packages
library(data.table)  
library(dplyr)       
library(tidyverse)   
library(limma)
library(caret)
library(Mime1)

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\08Machine Learning")
list.files()

#Step 3.1 Read genelist
#genelist <- read.table(file = "./Prepared data/Candidate_gene_list.txt", header = T)					#NFkB + Hippo + Notch + MAPK
genelist <- read.table(file = "./Prepared data/Candidate_gene_list2.txt", header = T)				#18 TSs + JAK/STAT + NFkB + Hippo + Notch + MAPK
genelist2 = genelist$ID

#Step 3.2 Load input for machine learning

#load("./Prepared data/TCGA-LUAD_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-LUSC_MLinput.Rdata")    #input_data
load("./Prepared data/TCGA-COAD_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-READ_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-PAAD_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-SKCM_MLinput.Rdata")    #input_data
input_data[["object1"]][1:5,1:5]


#Step 3.3 ML.Dev.Prog.Sig()
#ML. Dev. Prog. Sig() provides three modes: all, single, and double. Means using all ten algorithms and combinations. Single means using only one of the ten algorithms. Double means combining two algorithms.
#In most cases, we usually use the all mode to analyze data. If set to (default), univariate Cox regression will be performed first between the genes provided in the training dataset to screen for prognostic variables, which will then be used to construct the model.

res <- ML.Dev.Prog.Sig(train_data = input_data$object1,
                     list_train_vali_Data = input_data,
                     unicox.filter.for.candi = T,
                     unicox_p_cutoff = 0.05,
                     candidate_genes = genelist[,"ID"],
                     mode = 'all',nodesize =5,seed = 5201314 )

#save(res, file = "./TCGA-LUAD_MLres.RData")
#save(res, file = "./TCGA-LUSC_MLres.RData")
save(res, file = "./Data0.7_genelist2_Res/TCGA-COAD_MLres.RData")
#save(res, file = "./TCGA-READ_MLres.RData")
#save(res, file = "./Data0.7_genelist2_Res/TCGA-PAAD_MLres.RData")
#save(res, file = "./TCGA-SKCM_MLres.RData")


#Step 4 Loading Results

library(data.table)  
library(dplyr)       
library(tidyverse)   
library(limma)
library(caret)
library(Mime1)

#setwd("J:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\08Machine Learning")
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\08Machine Learning")
list.files()

#Step 4.1  Loading  input data 
#load("./Prepared data/TCGA-LUAD_MLinput.Rdata")    #input_data
load("./Prepared data/TCGA-COAD_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-PAAD_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-SKCM_MLinput.Rdata")    #input_data
list_train_vali_Data = input_data

#Step 4.2  Loading  input results 
#load("./Data0.7_genelist2_Res/TCGA-LUAD_MLres.Rdata")    #input_data
load("./Data0.7_genelist2_Res/TCGA-COAD_MLres.Rdata")    #input_data
#load("./Data0.7_genelist2_Res/TCGA-PAAD_MLres.Rdata")    #input_data
#load("./Data0.7_genelist2_Res/TCGA-SKCM_MLres.Rdata")    #input_data


#Step 4.3  Outputs
Unicox_genes <- res[["Sig.genes"]]
write.csv(Unicox_genes,'P0_Unicox_genes.csv')

cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35)
ggsave("P1_Overall_117_Cindex.pdf", width = 7, height = 12)

#Draw the C-index of specific models in different datasets
cindex_dis_select(res,
                  model="StepCox[forward] + plsRcox",
                  order= names(list_train_vali_Data))

#Draw patient survival curves based on risk scores calculated using specific models on different datasets
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + plsRcox",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)


#Calculate the AUC score for each model
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["object1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                            auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["object1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                            auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["object1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                            auc_cal_method="KM")


#Step 5.2  1/3/5 year AUC
#Draw 1, 3, and 5-year AUC for specific models in different datasets using all.auc.1y as an example
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
ggsave("P2_Overall_117_Cindex_AUC1.pdf", width = 7, height = 12)

#Draw 1, 3, and 5-year AUC for specific models in different datasets using all.auc.3y as an example
auc_dis_all(all.auc.3y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=3)
ggsave("P2_Overall_117_Cindex_AUC3.pdf", width = 7, height = 12)

#Draw 1, 3, and 5-year AUC for specific models in different datasets using all.auc.5y as an example
auc_dis_all(all.auc.5y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=5)
ggsave("P2_Overall_117_Cindex_AUC5.pdf", width = 7, height = 12)


#Step 5.3  1/3/5 year Machine Learning Model Survival

#Draw ROC for specific models in different datasets using all.auc.1y as an example
roc_vis(all.auc.1y,
        model_name = "StepCox[forward] + RSF",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1) + theme(panel.border = element_rect(color = "black",linewidth = 2), 
									text = element_text (size = 16, color = 'black'), 
									axis.text = element_text (size = 16, color = 'black',hjust = 0.5)
									)
ggsave("P3_Overall_117_Cindex_AUC1_ROC_top.pdf", width = 6.5, height = 6.5)

#Draw ROC for specific models in different datasets using all.auc.3y as an example
roc_vis(all.auc.3y,
        model_name = "StepCox[forward] + RSF",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=3) + theme(panel.border = element_rect(color = "black",linewidth = 2), 
									text = element_text (size = 16, color = 'black'), 
									axis.text = element_text (size = 16, color = 'black',hjust = 0.5)
									)
ggsave("P3_Overall_117_Cindex_AUC3_ROC_top.pdf", width = 6.5, height = 6.5)

#Draw ROC for specific models in different datasets using all.auc.5y as an example
roc_vis(all.auc.5y,
        model_name = "StepCox[both] + RSF",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=5) + theme(panel.border = element_rect(color = "black",linewidth = 2), 
									text = element_text (size = 16, color = 'black'), 
									axis.text = element_text (size = 16, color = 'black',hjust = 0.5)
									)
ggsave("P3_Overall_117_Cindex_AUC5_ROC_top.pdf", width = 6.5, height = 6.5)


#Step 5.4 All Survival 
#SKCM: Top 5 the same
selected_algorithms <- c("StepCox[forward] + RSF", "StepCox[both] + RSF","StepCox[backward] + RSF","RSF","Lasso + RSF")
#PAAD:  Top 3 the same
selected_algorithms <- c("StepCox[forward] + RSF", "RSF","Lasso + RSF")
#COAD: Top 5 the same
selected_algorithms <- c("StepCox[forward] + RSF", "StepCox[both] + RSF","StepCox[backward] + RSF")
#LUAD: Top 5 the same
selected_algorithms <- c("StepCox[forward] + RSF", "StepCox[both] + RSF","StepCox[backward] + RSF","RSF","Lasso + RSF")


for (j in c(1:length(selected_algorithms))) {
selected_model <- selected_algorithms[j]

survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = selected_model, dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)

#file name
file_name = paste('P4_Survival_117_Cindex_All_top', selected_model,".pdf")  	#"P4_Survival_117_Cindex_All_top.pdf"
ggsave(filename= file_name, width = 10, height = 5.5)
cat("File", file_name, "saved.\n")    ###Check save
}


#Step 5. Selection of core genes

#load packages
library(data.table)  
library(dplyr)       
library(tidyverse)   
library(limma)
library(caret)
library(Mime1)

#work dir
setwd("F:\\02The Drosophila Multitumor Model Reveals Gene Regulatory Networks Driving Tumorigenesis\\12New Part\\08Machine Learning")
list.files()

#Step 5.1  Loading  input data 
#load("./Prepared data/TCGA-LUAD_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-COAD_MLinput.Rdata")    #input_data
#load("./Prepared data/TCGA-PAAD_MLinput.Rdata")    #input_data
load("./Prepared data/TCGA-SKCM_MLinput.Rdata")    #input_data
list_train_vali_Data = input_data

#Step 5.2 Read genelist
#genelist <- read.table(file = "./Prepared data/Candidate_gene_list.txt", header = T)					#NFkB + Hippo + Notch + MAPK
genelist <- read.table(file = "./Prepared data/Candidate_gene_list2.txt", header = T)				#18 TSs + JAK/STAT + NFkB + Hippo + Notch + MAPK
genelist2 = genelist$ID


res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$object1,
                                            candidate_genes = genelist2,
                                            mode = "all",nodesize =5,seed = 5201314 )


#Visualization
core_feature_select(res.feature.all)

#Draw the contribution of genes screened through different methods
core_feature_rank(res.feature.all, top=20)


