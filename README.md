# Cross-Species Insights from ART-D to Uncover Evolutionarily Conserved Oncogenic Mechanisms

## **Summary:**  
Cancer arises from oncogenic clones, yet the dynamic mechanisms driving their stepwise evolution toward malignancy remain incompletely understood. Here, we establish the Atlas of *Ras*-driven Tumors in *Drosophila* (ART-D), a systematic, cross-species platform that dissects the molecular and phenotypic trajectories of tumorigenesis across ten genetically defined RasV12-driven models. By integrating longitudinal phenotypic profiling, we define three conserved stages of tumor development—initiation, promotion, and progression—distinguished by distinct shifts in tumor burden and tumor-induced cachexia. Transcriptomic analysis reveals stage-specific signaling rewiring: early tumorigenesis is characterized by co-activation of JAK/STAT, NF-κB/Toll, and MAPK pathways, whereas malignant progression is driven by Notch hyperactivation and Hippo pathway inactivation. Through integrative multi-omics and machine learning, we uncover an evolutionarily conserved pathogenic network coordinating JNK, NF-κB/Toll, Notch, and Hippo signaling, which we functionally validate across species. ART-D serves as a transformative resource bridging Drosophila genetics and human cancer biology, offering a robust framework for decoding conserved oncogenic principles and identifying of stage-specific vulnerabilities in RAS-driven cancers.

## **Keywords:**  
*RAS*; Tumorigenesis; Tumor malignancy; Notch; NK-κB/Toll; Hippo; *fry*; *Drosophila*; machine learning 

## **Multi-omics data analysis encompasses:**  
In this study, the longitudial phenotypic data profiling identified the stage-specific features of the Atlas of *Ras*-driven Tumors in *Drosophila* (ART-D), revealing 1) Initial stage (0–5.5 days AEL): Marked by tumor cell overgrowth. 2) Promotion stage (5.5–9.5 days AEL): Characterized by uncontrolled proliferation, onset of invasion, and early cachexia. 3) Progression stage (9.5–13.5 days AEL): Distinguished by severe cachexia and aggressive tissue invasion. The phenotype-dependent tumor subtyping identified three subgroups of *Ras*-driven Tumors (DRTs): Group A comprises *RasV12* tumors with deficiencies in *scrib*, *l(2)gl*, *Vps36*, *fry*, *emei*, or *fmt*; Group B consists of tumors lacking *Syx7* or *msn*; and Group C includes tumors deficient in *TSG101* or *Rabex-5*. 

To delineat the underlying mechanisms and transformative value of ART-D, we perfomed the following data analysis.

### 1. *Drosophila* Bulk RNA-seq
A transcriptomic atlas comprising *Drosophila* *WT* eye-antennal discs, *RasV12* benign tumors, and ten genetically defined *Drosophila* RasV12-driven tumors (DRTs) harboring deficiencies in *scrib*, *l(2)gl*, *Vps36*, *Syx7*, *Rabex-5*, *TSG101*, *fmt*, *emei*, *fry*, or *msn* at three different stages (initial/promotion/progression stage). The integrative analysis of phenotypic data and transcriptomic dynamic changes delineated the **Transcriptome Dynamics of Drosophila Ras-driven Tumor Subtypes**, revealing **a two-stage signaling paradigm governing tumor progression**.

### 2. *Drosophila* Bulk ATAC-seq
The integrative analysis of multi-omcs data sets of bulk RNA-seq and bulk ATAC-seq comprising *Drosophila* WT eye-antennal discs, Ras benign tumors, and RasV12-driven fry-deficient tumors at three different stages (initial/promotion/progression stage) implicated **JNK, NF-κB/Toll, Hippo, and Notch signaling as cooperative drivers** of tumorigenesis in *Ras*-activated, *fry*-deficient contexts.

### 3. Human cancer WGS/WES, bulk RNA-seq, Survival, and scRNA-seq datasets
The integrative pan-cancer analysis of WGS/WES, bulk RNA-seq, survival time of cBioProtal/GTEx/TCGA samples validated the clinical relevance of *Drosophila* tumor suppressors in human cancers with *RAS* activation;   
The analysis of scRNA-seq data from LUAD, PADC, COAD and SKCM patients consolidated the idea of JNK, NF-κB/Toll, Notch, and Hippo pathways co-operation to be **an evolutionarily conserved mechanism to drive tumorigenesis**.

### 4. Human cancer machine learning modeling
To explore the diagnostic value of tumor suppressors in DRTs for human cancers, we performed **117 combinations of machine learning algorithms to develop prognostic prediction models** tailored to different human KRAS/NRAS cancers. 

## **For more information, please read our paper:**  
Cross-Species Insights from ART-D to Uncover Evolutionarily Conserved Oncogenic Mechanisms    
https://doi.org/10.1101/2025.10.20.683390

The sequencing data can be downloaded freely, and codes developed for this study are provided.
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1209517/    
        
The gene count matrix after batch effect removal has been provided by    
newdata_filter_remove_pre_batheffect_removed.csv   
