'''
# Citation: Under submission
# Email: qususu@cibr.ac.cn
# Seurat v4.0.3 in R v4.1.0
# This file xxxx
'''


https://github.com/vshanka23/The-Drosophila-Brain-on-Cocaine-at-Single-Cell-Resolution/blob/master/Rcode_for_analyses.R
https://github.com/czbiohub/scell_lung_adenocarcinoma/blob/master/scripts/01_Import_data_and_metadata.Rmd
https://github.com/XinJason/Cross-kingdom-synthetic-microbiota



这个文件原来是：d100_scRNAseq_flydrug_merge12.R。为了发表到github现在在删减内容；


library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(dplyr)
library(rlang)
library(TopKLists)
library(readxl)
library(openxlsx)
library(xlsx)
library(Vennerable)
library(clustree) 
library(future)
plan("multiprocess",workers=9)
library(DoubletFinder)
options(future.globals.maxSize= 53687091200000000)

##FASTQ generation, demultiplexing and alignment were analyzed by the Basecall software and Cell Ranger v6.0.1, next all data processed via Seurat v4.0.3 in R v4.1.0.

##############################################################################
##Read data in & Create Seurat objects

#Read data in
A_R1_Data<-Read10X(data.dir="/home/zhangli_lab/qususu/DATA/projects/prj_scRNAseq_flydrug/CIBR-SUSU-10X-dmel-20210719/A/outs/filtered_feature_bc_matrix")
C_R1_Data<-Read10X(data.dir="/home/zhangli_lab/qususu/DATA/projects/prj_scRNAseq_flydrug/CIBR-SUSU-10X-dmel-20210719/C/outs/filtered_feature_bc_matrix")
M_R1_Data<-Read10X(data.dir="/home/zhangli_lab/qususu/DATA/projects/prj_scRNAseq_flydrug/CIBR-SUSU-10X-dmel-20210719/M/outs/filtered_feature_bc_matrix")
A_R2_Data<-Read10X(data.dir="/home/zhangli_lab/qususu/DATA/projects/prj_scRNAseq_flydrug/CIBR-SUSU-10X-dmel-20210719/A-0811/outs/filtered_feature_bc_matrix")
C_R2_Data<-Read10X(data.dir="/home/zhangli_lab/qususu/DATA/projects/prj_scRNAseq_flydrug/CIBR-SUSU-10X-dmel-20210719/C-0811/outs/filtered_feature_bc_matrix")
M_R2_Data<-Read10X(data.dir="/home/zhangli_lab/qususu/DATA/projects/prj_scRNAseq_flydrug/CIBR-SUSU-10X-dmel-20210719/M-0811/outs/filtered_feature_bc_matrix")

#Create Seurat objects
A_R1 <- CreateSeuratObject(counts=A_R1_Data,project="Male_ATX",min.cells=5)
C_R1 <- CreateSeuratObject(counts=C_R1_Data,project="Male_Control",min.cells=5)
M_R1 <- CreateSeuratObject(counts=M_R1_Data,project="Male_MPH",min.cells=5)
A_R2 <- CreateSeuratObject(counts=A_R2_Data,project="Male_ATX",min.cells=5)
C_R2 <- CreateSeuratObject(counts=C_R2_Data,project="Male_Control",min.cells=5)
M_R2 <- CreateSeuratObject(counts=M_R2_Data,project="Male_MPH",min.cells=5)

#Set identities for condition
A_R1$stim <- "ATX"
C_R1$stim <- "Control"
M_R1$stim <- "MPH"
A_R2$stim <- "ATX"
C_R2$stim <- "Control"
M_R2$stim <- "MPH"

#Set identities for gender and condition
A_R1$gender_stim <- "male_ATX"
C_R1$gender_stim <- "male_Control"
M_R1$gender_stim <- "male_MPH"
A_R2$gender_stim <- "male_ATX"
C_R2$gender_stim <- "male_Control"
M_R2$gender_stim <- "male_MPH"

#Set identities for sample
A_R1$sample_id <- "male_ATX_R1"
C_R1$sample_id <- "male_Control_R1"
M_R1$sample_id <- "male_MPH_R1"
A_R2$sample_id <- "male_ATX_R2"
C_R2$sample_id <- "male_Control_R2"
M_R2$sample_id <- "male_MPH_R2"

##############################################################################
##Pre-processing: quality control

#Calculate mitochondrial genes ratio
A_R1[["percent.mt"]] <- PercentageFeatureSet(A_R1, pattern = "^mt:")
C_R1[["percent.mt"]] <- PercentageFeatureSet(C_R1, pattern = "^mt:")
M_R1[["percent.mt"]] <- PercentageFeatureSet(M_R1, pattern = "^mt:")
A_R2[["percent.mt"]] <- PercentageFeatureSet(A_R2, pattern = "^mt:")
C_R2[["percent.mt"]] <- PercentageFeatureSet(C_R2, pattern = "^mt:")
M_R2[["percent.mt"]] <- PercentageFeatureSet(M_R2, pattern = "^mt:")

#Plots to see the distributions of the mitochondrial genes of each samples
sampleName=A_R1
plot1 <- FeatureScatter(sampleName, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sampleName, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1) 
#ggsave("A_R1_pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5)


#Remove spurious features (low and high)
A_R1.f <- subset(A_R1,subset=nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt <15) #15554
C_R1.f <- subset(C_R1,subset=nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt <15) #14708
M_R1.f <- subset(M_R1,subset=nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt <15) #17983
A_R2.f <- subset(A_R2,subset=nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt <15) #11495
C_R2.f <- subset(C_R2,subset=nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt <15) #14016
M_R2.f <- subset(M_R2,subset=nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt <15) #16011

#Normalize
AR1.f <- SCTransform(A_R1.f,verbose = FALSE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)
CR1.f <- SCTransform(C_R1.f,verbose = FALSE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)
MR1.f <- SCTransform(M_R1.f,verbose = FALSE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)
AR2.f <- SCTransform(A_R2.f,verbose = FALSE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)
CR2.f <- SCTransform(C_R2.f,verbose = FALSE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)
MR2.f <- SCTransform(M_R2.f,verbose = FALSE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)


# Remove doublets by DoubletFinder
AR1.f <- RunPCA(object = AR1.f, verbose = FALSE)
pcaplot <- ElbowPlot(AR1.f)
ggsave("AR1f_ElbowPlot_d1003.pdf", plot = pcaplot, width = 12, height = 5)
AR1.f <- RunUMAP(AR1.f,reduction = "pca", dims = 1:13)

CR1.f <- RunPCA(object = CR1.f, verbose = FALSE)
pcaplot <- ElbowPlot(CR1.f)
ggsave("CR1f_ElbowPlot_d1003.pdf", plot = pcaplot, width = 12, height = 5)
CR1.f <- RunUMAP(CR1.f,reduction = "pca", dims = 1:13)

MR1.f <- RunPCA(object = MR1.f, verbose = FALSE)
pcaplot <- ElbowPlot(MR1.f)
ggsave("MR1f_ElbowPlot_d1003.pdf", plot = pcaplot, width = 12, height = 5)
MR1.f <- RunUMAP(MR1.f,reduction = "pca", dims = 1:13)

AR2.f <- RunPCA(object = AR2.f, verbose = FALSE)
pcaplot <- ElbowPlot(AR2.f)
ggsave("AR2f_ElbowPlot_d1003.pdf", plot = pcaplot, width = 12, height = 5)
AR2.f <- RunUMAP(AR2.f,reduction = "pca", dims = 1:13)

CR2.f <- RunPCA(object = CR2.f, verbose = FALSE)
pcaplot <- ElbowPlot(CR2.f)
ggsave("CR2f_ElbowPlot_d1003.pdf", plot = pcaplot, width = 12, height = 5)
CR2.f <- RunUMAP(CR2.f,reduction = "pca", dims = 1:13)

MR2.f <- RunPCA(object = MR2.f, verbose = FALSE)
pcaplot <- ElbowPlot(MR2.f)
ggsave("MR2f_ElbowPlot_d1003.pdf", plot = pcaplot, width = 12, height = 5)
MR2.f <- RunUMAP(MR2.f,reduction = "pca", dims = 1:13)


#这是一个测试最佳参数的过程，运行速度慢
sweep.res.list.AR1.f <- paramSweep_v3(AR1.f, PCs = 1:13, sct = T)
sweep.stats.AR1.f <- summarizeSweep(sweep.res.list.AR1.f, GT = FALSE)  
bcmvn.AR1.f <- find.pK(sweep.stats.AR1.f) #可以看到最佳参数的点
pK_bcmvn.AR1.f <- bcmvn.AR1.f$pK[which.max(bcmvn.AR1.f$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
#0.27 win

sweep.res.list.CR1.f <- paramSweep_v3(CR1.f, PCs = 1:13, sct = T)
sweep.stats.CR1.f <- summarizeSweep(sweep.res.list.CR1.f, GT = FALSE)  
bcmvn.CR1.f <- find.pK(sweep.stats.CR1.f) #可以看到最佳参数的点
#显示为NULL，很奇怪
pK_bcmvn.CR1.f <- bcmvn.CR1.f$pK[which.max(bcmvn.CR1.f$BCmetric)] %>% as.character() %>% as.numeric()
#0.01 win

sweep.res.list.MR1.f <- paramSweep_v3(MR1.f, PCs = 1:13, sct = T)
sweep.stats.MR1.f <- summarizeSweep(sweep.res.list.MR1.f, GT = FALSE)  
bcmvn.MR1.f <- find.pK(sweep.stats.MR1.f) #可以看到最佳参数的点
#显示为NULL，很奇怪
pK_bcmvn.MR1.f <- bcmvn.MR1.f$pK[which.max(bcmvn.MR1.f$BCmetric)] %>% as.character() %>% as.numeric()
#0.08 win

sweep.res.list.AR2.f <- paramSweep_v3(AR2.f, PCs = 1:13, sct = T)
sweep.stats.AR2.f <- summarizeSweep(sweep.res.list.AR2.f, GT = FALSE)  
bcmvn.AR2.f <- find.pK(sweep.stats.AR2.f) #可以看到最佳参数的点
pK_bcmvn.AR2.f <- bcmvn.AR2.f$pK[which.max(bcmvn.AR2.f$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
#0.02 win

sweep.res.list.CR2.f <- paramSweep_v3(CR2.f, PCs = 1:13, sct = T)
sweep.stats.CR2.f <- summarizeSweep(sweep.res.list.CR2.f, GT = FALSE)  
bcmvn.CR2.f <- find.pK(sweep.stats.CR2.f) #可以看到最佳参数的点
#显示为NULL，很奇怪
pK_bcmvn.CR2.f <- bcmvn.CR2.f$pK[which.max(bcmvn.CR2.f$BCmetric)] %>% as.character() %>% as.numeric()
#0.1 win

sweep.res.list.MR2.f <- paramSweep_v3(MR2.f, PCs = 1:13, sct = T)
sweep.stats.MR2.f <- summarizeSweep(sweep.res.list.MR2.f, GT = FALSE)  
bcmvn.MR2.f <- find.pK(sweep.stats.MR2.f) #可以看到最佳参数的点
#显示为NULL，很奇怪
pK_bcmvn.MR2.f <- bcmvn.MR2.f$pK[which.max(bcmvn.MR2.f$BCmetric)] %>% as.character() %>% as.numeric()
#0.29 win





# 10x: https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
#每个样本的load细胞数量不太一样：第一批按照1w进行捕获；第二批按照9000捕获。但是实际load细胞数量都比较高。
DoubletRate1 = 0.080    # 直接查表，10000细胞对应的doublets rate是～8.0%
DoubletRate2 = 0.072


# Replication 1
annotations.AR1.f <- AR1.f@meta.data$ClusteringResults
homotypic.prop.AR1.f <- modelHomotypic(annotations.AR1.f)
nExp_poi.AR1.f <- round(DoubletRate1*nrow(AR1.f@meta.data)) 
nExp_poi.adj.AR1.f <- round(nExp_poi.AR1.f*(1-homotypic.prop.AR1.f))

AR1.f <- doubletFinder_v3(AR1.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.AR1.f, nExp = nExp_poi.adj.AR1.f, reuse.pANN = F, sct = T)
#DimPlot(AR1.f,reduction="umap",group.by="DF.classifications_0.25_0.27_1244")
# name of the DF prediction can change, so extract the correct column name.
DF.name.AR1.f = colnames(AR1.f@meta.data)[grepl("DF.classification", colnames(AR1.f@meta.data))]
doubletsPlot <- cowplot::plot_grid(ncol = 2, DimPlot(AR1.f, group.by = "orig.ident") + NoAxes(), DimPlot(AR1.f, group.by = DF.name.AR1.f) + NoAxes())
ggsave("AR1f_doubletsPlot_d1003.pdf", plot = doubletsPlot, width = 12, height = 5)

AR1.f <- doubletFinder_v3(AR1.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.AR1.f, nExp = nExp_poi.adj.AR1.f, reuse.pANN = "pANN_0.25_0.27_1244", sct = T)
AR1.f = AR1.f[, AR1.f@meta.data[, DF.name.AR1.f] == "Singlet"]
#AR1.f <- subset(AR1.f, cells=rownames(AR1.f@meta.data)[which(AR1.f@meta.data$DF.classifications_0.25_0.27_1244 == "Singlet")])
table(AR1.f@meta.data$orig.ident)
#14310



annotations.CR1.f <- CR1.f@meta.data$ClusteringResults
homotypic.prop.CR1.f <- modelHomotypic(annotations.CR1.f)
nExp_poi.CR1.f <- round(DoubletRate1*nrow(CR1.f@meta.data)) 
nExp_poi.adj.CR1.f <- round(nExp_poi.CR1.f*(1-homotypic.prop.CR1.f))

CR1.f <- doubletFinder_v3(CR1.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.CR1.f, nExp = nExp_poi.adj.CR1.f, reuse.pANN = F, sct = T)
#DimPlot(CR1.f,reduction="umap",group.by="DF.classifications_0.25_0.01_1177")
# name of the DF prediction can change, so extract the correct column name.
DF.name.CR1.f = colnames(CR1.f@meta.data)[grepl("DF.classification", colnames(CR1.f@meta.data))]
doubletsPlot <- cowplot::plot_grid(ncol = 2, DimPlot(CR1.f, group.by = "orig.ident") + NoAxes(), DimPlot(CR1.f, group.by = DF.name.CR1.f) + NoAxes())
ggsave("CR1f_doubletsPlot_d1003.pdf", plot = doubletsPlot, width = 12, height = 5)

CR1.f <- doubletFinder_v3(CR1.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.CR1.f, nExp = nExp_poi.adj.CR1.f, reuse.pANN = "pANN_0.25_0.01_1177", sct = T)
CR1.f = CR1.f[, CR1.f@meta.data[, DF.name.CR1.f] == "Singlet"]
#CR1.f <- subset(CR1.f, cells=rownames(CR1.f@meta.data)[which(CR1.f@meta.data$DF.classifications_0.25_0.01_1177 == "Singlet")])
table(CR1.f@meta.data$orig.ident)
#13531


annotations.MR1.f <- MR1.f@meta.data$ClusteringResults
homotypic.prop.MR1.f <- modelHomotypic(annotations.MR1.f)
nExp_poi.MR1.f <- round(DoubletRate1*nrow(MR1.f@meta.data)) 
nExp_poi.adj.MR1.f <- round(nExp_poi.MR1.f*(1-homotypic.prop.MR1.f))

MR1.f <- doubletFinder_v3(MR1.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.MR1.f, nExp = nExp_poi.adj.MR1.f, reuse.pANN = F, sct = T)
#DimPlot(MR1.f,reduction="umap",group.by="DF.classifications_0.25_0.08_1439")
# name of the DF prediction can change, so extract the correct column name.
DF.name.MR1.f = colnames(MR1.f@meta.data)[grepl("DF.classification", colnames(MR1.f@meta.data))]
doubletsPlot <- cowplot::plot_grid(ncol = 2, DimPlot(MR1.f, group.by = "orig.ident") + NoAxes(), DimPlot(MR1.f, group.by = DF.name.MR1.f) + NoAxes())
ggsave("MR1f_doubletsPlot_d1003.pdf", plot = doubletsPlot, width = 12, height = 5)

MR1.f <- doubletFinder_v3(MR1.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.MR1.f, nExp = nExp_poi.adj.MR1.f, reuse.pANN = "pANN_0.25_0.08_1439", sct = T)
MR1.f = MR1.f[, MR1.f@meta.data[, DF.name.MR1.f] == "Singlet"]
#MR1.f <- subset(MR1.f, cells=rownames(MR1.f@meta.data)[which(MR1.f@meta.data$DF.classifications_0.25_0.08_1439 == "Singlet")])
table(MR1.f@meta.data$orig.ident)
#16544

# Replication 2
annotations.AR2.f <- AR2.f@meta.data$ClusteringResults
homotypic.prop.AR2.f <- modelHomotypic(annotations.AR2.f)
nExp_poi.AR2.f <- round(DoubletRate2*nrow(AR2.f@meta.data)) 
nExp_poi.adj.AR2.f <- round(nExp_poi.AR2.f*(1-homotypic.prop.AR2.f))

AR2.f <- doubletFinder_v3(AR2.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.AR2.f, nExp = nExp_poi.adj.AR2.f, reuse.pANN = F, sct = T)
#DimPlot(AR2.f,reduction="umap",group.by="DF.classifications_0.25_0.02_828")
# name of the DF prediction can change, so extract the correct column name.
DF.name.AR2.f = colnames(AR2.f@meta.data)[grepl("DF.classification", colnames(AR2.f@meta.data))]
doubletsPlot <- cowplot::plot_grid(ncol = 2, DimPlot(AR2.f, group.by = "orig.ident") + NoAxes(), DimPlot(AR2.f, group.by = DF.name.AR2.f) + NoAxes())
ggsave("AR2f_doubletsPlot_d1003.pdf", plot = doubletsPlot, width = 12, height = 5)

AR2.f <- doubletFinder_v3(AR2.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.AR2.f, nExp = nExp_poi.adj.AR2.f, reuse.pANN = "pANN_0.25_0.02_828", sct = T)
AR2.f = AR2.f[, AR2.f@meta.data[, DF.name.AR2.f] == "Singlet"]
#AR2.f <- subset(AR2.f, cells=rownames(AR2.f@meta.data)[which(AR2.f@meta.data$DF.classifications_0.25_0.02_828 == "Singlet")])
table(AR2.f@meta.data$orig.ident)
#10667


annotations.CR2.f <- CR2.f@meta.data$ClusteringResults
homotypic.prop.CR2.f <- modelHomotypic(annotations.CR2.f)
nExp_poi.CR2.f <- round(DoubletRate2*nrow(CR2.f@meta.data)) 
nExp_poi.adj.CR2.f <- round(nExp_poi.CR2.f*(1-homotypic.prop.CR2.f))


CR2.f <- doubletFinder_v3(CR2.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.CR2.f, nExp = nExp_poi.adj.CR2.f, reuse.pANN = F, sct = T)
#DimPlot(CR2.f,reduction="umap",group.by="DF.classifications_0.25_0.1_1009")
# name of the DF prediction can change, so extract the correct column name.
DF.name.CR2.f = colnames(CR2.f@meta.data)[grepl("DF.classification", colnames(CR2.f@meta.data))]
doubletsPlot <- cowplot::plot_grid(ncol = 2, DimPlot(CR2.f, group.by = "orig.ident") + NoAxes(), DimPlot(CR2.f, group.by = DF.name.CR2.f) + NoAxes())
ggsave("CR2f_doubletsPlot_d1003.pdf", plot = doubletsPlot, width = 12, height = 5)

CR2.f <- doubletFinder_v3(CR2.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.CR2.f, nExp = nExp_poi.adj.CR2.f, reuse.pANN = "pANN_0.25_0.1_1009", sct = T)
CR2.f = CR2.f[, CR2.f@meta.data[, DF.name.CR2.f] == "Singlet"]
#CR2.f <- subset(CR2.f, cells=rownames(CR2.f@meta.data)[which(CR2.f@meta.data$DF.classifications_0.25_0.1_1009 == "Singlet")])
table(CR2.f@meta.data$orig.ident)
#13007


annotations.MR2.f <- MR2.f@meta.data$ClusteringResults
homotypic.prop.MR2.f <- modelHomotypic(annotations.MR2.f)
nExp_poi.MR2.f <- round(DoubletRate2*nrow(MR2.f@meta.data)) 
nExp_poi.adj.MR2.f <- round(nExp_poi.MR2.f*(1-homotypic.prop.MR2.f))


MR2.f <- doubletFinder_v3(MR2.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.MR2.f, nExp = nExp_poi.adj.MR2.f, reuse.pANN = F, sct = T)
#DimPlot(MR2.f,reduction="umap",group.by="DF.classifications_0.25_0.29_1153")
# name of the DF prediction can change, so extract the correct column name.
DF.name.MR2.f = colnames(MR2.f@meta.data)[grepl("DF.classification", colnames(MR2.f@meta.data))]
doubletsPlot <- cowplot::plot_grid(ncol = 2, DimPlot(MR2.f, group.by = "orig.ident") + NoAxes(), DimPlot(MR2.f, group.by = DF.name.MR2.f) + NoAxes())
ggsave("MR2f_doubletsPlot_d1003.pdf", plot = doubletsPlot, width = 12, height = 5)

MR2.f <- doubletFinder_v3(MR2.f, PCs = 1:13, pN = 0.25, pK = pK_bcmvn.MR2.f, nExp = nExp_poi.adj.MR2.f, reuse.pANN = "pANN_0.25_0.29_1153", sct = T)
MR2.f = MR2.f[, MR2.f@meta.data[, DF.name.MR2.f] == "Singlet"]
#MR2.f <- subset(MR2.f, cells=rownames(MR2.f@meta.data)[which(MR2.f@meta.data$DF.classifications_0.25_0.29_1153 == "Singlet")])
table(MR2.f@meta.data$orig.ident)
#14858


##############################################################################
##Integrate all samples after pre-processing

#Integration & cluster
Brain.features <- SelectIntegrationFeatures(object.list=c(AR1.f,CR1.f,MR1.f,AR2.f,CR2.f,MR2.f),nfeatures=3000)
Brain.features
Brain.list <- PrepSCTIntegration(object.list = c(AR1.f,CR1.f,MR1.f,AR2.f,CR2.f,MR2.f),anchor.features = Brain.features,verbose = FALSE)
##Identify anchors #并行很快
Brain.anchors <- FindIntegrationAnchors(object.list = Brain.list, normalization.method = "SCT", anchor.features = Brain.features, verbose = FALSE)   
Brain.integrated <- IntegrateData(anchorset = Brain.anchors, normalization.method = "SCT", verbose = FALSE)

Brain.integrated <- RunPCA(object = Brain.integrated, verbose = FALSE)
pcaplot <- ElbowPlot(Brain.integrated)
ggsave("1_filter_merge12_ElbowPlot2_d1003_re.pdf", plot = pcaplot, width = 10, height = 10)

Brain.integrated <- RunUMAP(Brain.integrated, reduction = "pca", dims = 1:15)
Brain.integrated <- FindNeighbors(Brain.integrated,reduction = "pca", dims = 1:15)
Brain.integrated.clusters <- FindClusters(object = Brain.integrated,resolution = c(seq(.0,1.2,.1)))
#clusTree <- clustree(Brain.integrated.clusters@meta.data, prefix = "integrated_snn_res.")
#ggsave("1_filter_merge12_clusTree_pc15_d1003_re.pdf", plot = clusTree, width = 14, height = 14)

Brain.integrated <- FindClusters(Brain.integrated, resolution = 0.1) #28 clusters

Brain.integrated$Celltype <- Idents(Brain.integrated)
plots <- DimPlot(Brain.integrated,reduction="umap",group.by=c("stim","Celltype","sample_id"),combine=FALSE)
plots <- lapply(X = plots, FUN = function(x) + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2), byrow= TRUE, override.aes = list(size=3)))
combined <- CombinePlots(plots)

#ggsave("1_filter_merge12_pc15_dimPlot_d1003_re.pdf", plot = combined, width = 18, height = 14)


#Find top 10 markers genes of clusters 
DefaultAssay(Brain.integrated) <- "RNA"
Brain.integrated <- NormalizeData(Brain.integrated, verbose = FALSE)

Brain.integrated_markers_all_4 <- FindAllMarkers(Brain.integrated, min.pct=0.25, logfc.threshold = 0.25, only.pos = TRUE, assay="SCT",slot="data",test.use="MAST")

Brain.integrated_cluster_markers10_4 <- Brain.integrated_markers_all_4 %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) %>% print(n=3*37) 
write.csv(Brain.integrated_cluster_markers10_4, "1_filter_merge12_cluster_markers_top10_SCTdata_MAST.csv")


#Identify cell types, rename clusters according the top marker genes & known genes of the tissue--------
new.cluster.ids.Combined <- c("0-ClusterA", "1-OL1", "2-EB1", "3-PNs1", "4-ClusterB", "5-ClusterC", "6-OL2", "7-Glia", "8-PNs2", "9-MBKC_a", "10-PNs3", "11-PNs4", "12-EB2", "13-ClusterD", "14-ClusterE", "15-OL3", "16-PNs5", "17-OL4", "18-MBKC_b", "19-OL5", "20-Monoamines", "21-OL6", "22-OL7", "23-OL8", "24-ClusterF", "25-PNs6", "26-ClusterG", "27-ClusterH")
names(new.cluster.ids.Combined)<-levels(Brain.integrated)
Brain.integrated<-RenameIdents(Brain.integrated,new.cluster.ids.Combined)

dimPlot.rename <- DimPlot(Brain.integrated,label = TRUE,pt.size = 1,reduction="umap",label.size=3, label.box = T)
ggsave("1_filter_merge12_renameClusters_1st.pdf", plot = dimPlot.rename, width =14, height = 10)
##############################################################################