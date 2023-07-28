'''

# Citation: Under Revision
# Email: qususu@cibr.ac.cn
# Date: 2023-05-11
# Brief introduction:
    Multiple-level analyses were conducted after getting the primary clusters as shown in the Methods. Main steps were processed by using R v4.1.0, 
    such as DEG analysis, re-cluster of some clusters and so on. Other steps including pathway analysis, cell-cell communication analysis and potential 
    gene analysis were analyzed easily by using other softwares which are not introduced here. Users can repeat follow our manuscripts.
    
# Main steps were processed by using R v4.1.0 as follows:
    1.DEG analysis of drug treatment group compared with control based on 28 clusters seperately;
    2.Re-clustering of specific primary clusters:
        1)re-cluster of 20-Monoamines
        2)re-cluster of 7-Glia
        
'''

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
library(gtools)
library(RColorBrewer)

####################################################################################################
####################################### Perform DEG analysis #######################################

#Create new identities for DE analysis downstream
Brain.integrated$celltype.stim <- paste(Idents(Brain.integrated),Brain.integrated$stim,sep="_")
Brain.integrated$cluster_ids <- Idents(Brain.integrated)
Idents(Brain.integrated) <- "celltype.stim"

#ATX with control
Brain_ATX_merge12_DE_male_stim_C0 <- FindMarkers(Brain.integrated,ident.1 = "0-ClusterA_ATX",ident.2 = "0-ClusterA_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C1 <- FindMarkers(Brain.integrated,ident.1 = "1-OL1_ATX",ident.2 = "1-OL1_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C2 <- FindMarkers(Brain.integrated,ident.1 = "2-EB1_ATX",ident.2 = "2-EB1_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C3 <- FindMarkers(Brain.integrated,ident.1 = "3-PNs1_ATX",ident.2 = "3-PNs1_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C4 <- FindMarkers(Brain.integrated,ident.1 = "4-ClusterB_ATX",ident.2 = "4-ClusterB_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C5 <- FindMarkers(Brain.integrated,ident.1 = "5-ClusterC_ATX",ident.2 = "5-ClusterC_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C6 <- FindMarkers(Brain.integrated,ident.1 = "6-OL2_ATX",ident.2 = "6-OL2_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C7 <- FindMarkers(Brain.integrated,ident.1 = "7-Glia_ATX",ident.2 = "7-Glia_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C8 <- FindMarkers(Brain.integrated,ident.1 = "8-PNs2_ATX",ident.2 = "8-PNs2_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C9 <- FindMarkers(Brain.integrated,ident.1 = "9-MBKC_a_ATX",ident.2 = "9-MBKC_a_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C10 <- FindMarkers(Brain.integrated,ident.1 = "10-PNs3_ATX",ident.2 = "10-PNs3_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C11 <- FindMarkers(Brain.integrated,ident.1 = "11-PNs4_ATX",ident.2 = "11-PNs4_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C12 <- FindMarkers(Brain.integrated,ident.1 = "12-EB2_ATX",ident.2 = "12-EB2_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C13 <- FindMarkers(Brain.integrated,ident.1 = "13-ClusterD_ATX",ident.2 = "13-ClusterD_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C14 <- FindMarkers(Brain.integrated,ident.1 = "14-ClusterE_ATX",ident.2 = "14-ClusterE_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C15 <- FindMarkers(Brain.integrated,ident.1 = "15-OL3_ATX",ident.2 = "15-OL3_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C16 <- FindMarkers(Brain.integrated,ident.1 = "16-PNs5_ATX",ident.2 = "16-PNs5_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C17 <- FindMarkers(Brain.integrated,ident.1 = "17-OL4_ATX",ident.2 = "17-OL4_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C18 <- FindMarkers(Brain.integrated,ident.1 = "18-MBKC_b_ATX",ident.2 = "18-MBKC_b_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C19 <- FindMarkers(Brain.integrated,ident.1 = "19-OL5_ATX",ident.2 = "19-OL5_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C20 <- FindMarkers(Brain.integrated,ident.1 = "20-Monoamines_ATX",ident.2 = "20-Monoamines_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C21 <- FindMarkers(Brain.integrated,ident.1 = "21-OL6_ATX",ident.2 = "21-OL6_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C22 <- FindMarkers(Brain.integrated,ident.1 = "22-OL7_ATX",ident.2 = "22-OL7_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C23 <- FindMarkers(Brain.integrated,ident.1 = "23-OL8_ATX",ident.2 = "23-OL8_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C24 <- FindMarkers(Brain.integrated,ident.1 = "24-ClusterF_ATX",ident.2 = "24-ClusterF_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C25 <- FindMarkers(Brain.integrated,ident.1 = "25-PNs6_ATX",ident.2 = "25-PNs6_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C26 <- FindMarkers(Brain.integrated,ident.1 = "26-ClusterG_ATX",ident.2 = "26-ClusterG_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_ATX_merge12_DE_male_stim_C27 <- FindMarkers(Brain.integrated,ident.1 = "27-ClusterH_ATX",ident.2 = "27-ClusterH_Control",test.use = "MAST", assay = "SCT", slot = "data")

Brain_MPH_merge12_DE_male_stim_C0 <- FindMarkers(Brain.integrated,ident.1 = "0-ClusterA_MPH",ident.2 = "0-ClusterA_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C1 <- FindMarkers(Brain.integrated,ident.1 = "1-OL1_MPH",ident.2 = "1-OL1_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C2 <- FindMarkers(Brain.integrated,ident.1 = "2-EB1_MPH",ident.2 = "2-EB1_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C3 <- FindMarkers(Brain.integrated,ident.1 = "3-PNs1_MPH",ident.2 = "3-PNs1_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C4 <- FindMarkers(Brain.integrated,ident.1 = "4-ClusterB_MPH",ident.2 = "4-ClusterB_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C5 <- FindMarkers(Brain.integrated,ident.1 = "5-ClusterC_MPH",ident.2 = "5-ClusterC_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C6 <- FindMarkers(Brain.integrated,ident.1 = "6-OL2_MPH",ident.2 = "6-OL2_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C7 <- FindMarkers(Brain.integrated,ident.1 = "7-Glia_MPH",ident.2 = "7-Glia_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C8 <- FindMarkers(Brain.integrated,ident.1 = "8-PNs2_MPH",ident.2 = "8-PNs2_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C9 <- FindMarkers(Brain.integrated,ident.1 = "9-MBKC_a_MPH",ident.2 = "9-MBKC_a_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C10 <- FindMarkers(Brain.integrated,ident.1 = "10-PNs3_MPH",ident.2 = "10-PNs3_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C11 <- FindMarkers(Brain.integrated,ident.1 = "11-PNs4_MPH",ident.2 = "11-PNs4_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C12 <- FindMarkers(Brain.integrated,ident.1 = "12-EB2_MPH",ident.2 = "12-EB2_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C13 <- FindMarkers(Brain.integrated,ident.1 = "13-ClusterD_MPH",ident.2 = "13-ClusterD_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C14 <- FindMarkers(Brain.integrated,ident.1 = "14-ClusterE_MPH",ident.2 = "14-ClusterE_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C15 <- FindMarkers(Brain.integrated,ident.1 = "15-OL3_MPH",ident.2 = "15-OL3_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C16 <- FindMarkers(Brain.integrated,ident.1 = "16-PNs5_MPH",ident.2 = "16-PNs5_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C17 <- FindMarkers(Brain.integrated,ident.1 = "17-OL4_MPH",ident.2 = "17-OL4_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C18 <- FindMarkers(Brain.integrated,ident.1 = "18-MBKC_b_MPH",ident.2 = "18-MBKC_b_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C19 <- FindMarkers(Brain.integrated,ident.1 = "19-OL5_MPH",ident.2 = "19-OL5_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C20 <- FindMarkers(Brain.integrated,ident.1 = "20-Monoamines_MPH",ident.2 = "20-Monoamines_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C21 <- FindMarkers(Brain.integrated,ident.1 = "21-OL6_MPH",ident.2 = "21-OL6_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C22 <- FindMarkers(Brain.integrated,ident.1 = "22-OL7_MPH",ident.2 = "22-OL7_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C23 <- FindMarkers(Brain.integrated,ident.1 = "23-OL8_MPH",ident.2 = "23-OL8_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C24 <- FindMarkers(Brain.integrated,ident.1 = "24-ClusterF_MPH",ident.2 = "24-ClusterF_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C25 <- FindMarkers(Brain.integrated,ident.1 = "25-PNs6_MPH",ident.2 = "25-PNs6_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C26 <- FindMarkers(Brain.integrated,ident.1 = "26-ClusterG_MPH",ident.2 = "26-ClusterG_Control",test.use = "MAST", assay = "SCT", slot = "data")
Brain_MPH_merge12_DE_male_stim_C27 <- FindMarkers(Brain.integrated,ident.1 = "27-ClusterH_MPH",ident.2 = "27-ClusterH_Control",test.use = "MAST", assay = "SCT", slot = "data")

#Export results to excel
ATX_merge12_DE_list_SCT_data <- list("ATX_merge12_C0_SCT_data" = Brain_ATX_merge12_DE_male_stim_C0,"ATX_merge12_C1"=Brain_ATX_merge12_DE_male_stim_C1,"ATX_merge12_C2"=Brain_ATX_merge12_DE_male_stim_C2,"ATX_merge12_C3"=Brain_ATX_merge12_DE_male_stim_C3,"ATX_merge12_C4"=Brain_ATX_merge12_DE_male_stim_C4,"ATX_merge12_C5"=Brain_ATX_merge12_DE_male_stim_C5,"ATX_merge12_C6"=Brain_ATX_merge12_DE_male_stim_C6,"ATX_merge12_C7"=Brain_ATX_merge12_DE_male_stim_C7,"ATX_merge12_C8"=Brain_ATX_merge12_DE_male_stim_C8,"ATX_merge12_C9"=Brain_ATX_merge12_DE_male_stim_C9,"ATX_merge12_C10"=Brain_ATX_merge12_DE_male_stim_C10,"ATX_merge12_C11"=Brain_ATX_merge12_DE_male_stim_C11,"ATX_merge12_C12"=Brain_ATX_merge12_DE_male_stim_C12,"ATX_merge12_C13"=Brain_ATX_merge12_DE_male_stim_C13,"ATX_merge12_C14"=Brain_ATX_merge12_DE_male_stim_C14,"ATX_merge12_15"=Brain_ATX_merge12_DE_male_stim_C15,"ATX_merge12_C16"=Brain_ATX_merge12_DE_male_stim_C16,"ATX_merge12_C17"=Brain_ATX_merge12_DE_male_stim_C17,"ATX_merge12_C18"=Brain_ATX_merge12_DE_male_stim_C18,"ATX_merge12_C19"=Brain_ATX_merge12_DE_male_stim_C19,"ATX_merge12_C20"=Brain_ATX_merge12_DE_male_stim_C20,"ATX_merge12_C21"=Brain_ATX_merge12_DE_male_stim_C21,"ATX_merge12_C22"=Brain_ATX_merge12_DE_male_stim_C22,"ATX_merge12_C23"=Brain_ATX_merge12_DE_male_stim_C23,"ATX_merge12_C24"=Brain_ATX_merge12_DE_male_stim_C24,"ATX_merge12_C25"=Brain_ATX_merge12_DE_male_stim_C25,"ATX_merge12_C26"=Brain_ATX_merge12_DE_male_stim_C26,"ATX_merge12_C27"=Brain_ATX_merge12_DE_male_stim_C27)
openxlsx::write.xlsx(ATX_merge12_DE_list_SCT_data,file="Brain_ATX_merge12_d1003_re_res0.1_DE_male_SCT_data.xlsx",rowNames=TRUE)

MPH_merge12_DE_list_SCT_data <- list("MPH_merge12_C0_SCT_data" = Brain_MPH_merge12_DE_male_stim_C0,"MPH_merge12_C1"=Brain_MPH_merge12_DE_male_stim_C1,"MPH_merge12_C2"=Brain_MPH_merge12_DE_male_stim_C2,"MPH_merge12_C3"=Brain_MPH_merge12_DE_male_stim_C3,"MPH_merge12_C4"=Brain_MPH_merge12_DE_male_stim_C4,"MPH_merge12_C5"=Brain_MPH_merge12_DE_male_stim_C5,"MPH_merge12_C6"=Brain_MPH_merge12_DE_male_stim_C6,"MPH_merge12_C7"=Brain_MPH_merge12_DE_male_stim_C7,"MPH_merge12_C8"=Brain_MPH_merge12_DE_male_stim_C8,"MPH_merge12_C9"=Brain_MPH_merge12_DE_male_stim_C9,"MPH_merge12_C10"=Brain_MPH_merge12_DE_male_stim_C10,"MPH_merge12_C11"=Brain_MPH_merge12_DE_male_stim_C11,"MPH_merge12_C12"=Brain_MPH_merge12_DE_male_stim_C12,"MPH_merge12_C13"=Brain_MPH_merge12_DE_male_stim_C13,"MPH_merge12_C14"=Brain_MPH_merge12_DE_male_stim_C14,"MPH_merge12_15"=Brain_MPH_merge12_DE_male_stim_C15,"MPH_merge12_C16"=Brain_MPH_merge12_DE_male_stim_C16,"MPH_merge12_C17"=Brain_MPH_merge12_DE_male_stim_C17,"MPH_merge12_C18"=Brain_MPH_merge12_DE_male_stim_C18,"MPH_merge12_C19"=Brain_MPH_merge12_DE_male_stim_C19,"MPH_merge12_C20"=Brain_MPH_merge12_DE_male_stim_C20,"MPH_merge12_C21"=Brain_MPH_merge12_DE_male_stim_C21,"MPH_merge12_C22"=Brain_MPH_merge12_DE_male_stim_C22,"MPH_merge12_C23"=Brain_MPH_merge12_DE_male_stim_C23,"MPH_merge12_C24"=Brain_MPH_merge12_DE_male_stim_C24,"MPH_merge12_C25"=Brain_MPH_merge12_DE_male_stim_C25,"MPH_merge12_C26"=Brain_MPH_merge12_DE_male_stim_C26,"MPH_merge12_C27"=Brain_MPH_merge12_DE_male_stim_C27)
openxlsx::write.xlsx(MPH_merge12_DE_list_SCT_data,file="Brain_MPH_merge12_d1003_re_res0.1_DE_male_SCT_data.xlsx",rowNames=TRUE)

####################################################################################################
############################# Re-clustering of specific primary clusters ###########################

scRNA=Brain.integrated
DefaultAssay(scRNA)<-"integrated"

#re-cluster of 20-Monoamines-----------------------------
Monoamines=subset(scRNA,cluster_ids==c("20-Monoamines"))
Monoamines <- FindVariableFeatures(Monoamines,verbose = FALSE,selection.method = "vst")
Monoamines<-ScaleData(Monoamines)
Monoamines <- RunPCA(object = Monoamines, verbose = FALSE,assay="integrated")
ElbowPlot(Monoamines)
Monoamines<- RunUMAP(Monoamines, reduction = "pca", dims = 1:15)
Monoamines<-RunTSNE(Monoamines,dims = 1:15)
Monoamines<-FindNeighbors(Monoamines,reduction = "pca", dims = 1:15)
Monoamines.clusters <- FindClusters(Monoamines,resolution = c(seq(0,1.2,.1)))
clustree(Monoamines.clusters@meta.data, prefix = "integrated_snn_res.")
Monoamines <- FindClusters(Monoamines, resolution = 0.3) 

Monoamines$Celltype <- Idents(Monoamines)
plots <- DimPlot(Monoamines,reduction="umap",group.by=c("stim","Celltype","sample_id"),combine=FALSE)
plots <- lapply(X = plots, FUN = function(x) + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2), byrow= TRUE, override.aes = list(size=3)))
combined.Monoamines <- CombinePlots(plots)

cols=c('grey', 'lightskyblue1', '#CCCC99', '#20B2AA', 'grey', 'grey', 'lightskyblue1', '#B088FF', '#20B2AA', '#FFD700', '#20B2AA', '#20B2AA', '#CCCC99', 'grey', 'grey', 'lightskyblue1', '#20B2AA', 'lightskyblue1', '#FFD700', 'lightskyblue1', '#FF69B4', 'lightskyblue1', 'lightskyblue1', 'lightskyblue1', 'grey', '#20B2AA', 'grey', 'grey')
dimPlot<-DimPlot(scRNA,label = TRUE,pt.size = 0.5,reduction="umap",label.size=4, label.box = F,label.color="#474747",group.by="cluster_ids",cols=cols)+NoLegend() +labs(x = "UMAP_1", y = "UMAP_2",title = " ")
VmatPlot<-dimPlot+FeaturePlot(scRNA,feature=c("Vmat"),pt.size=0.5)
ggsave("./manuscriptsPics/recheck0712/VmatPlot_allsamples_20220804.pdf", plot = VmatPlot, width =14, height = 7,dpi=600)

#re-cluster of 7-Glia --------------------------------------
Glias=subset(scRNA,cluster_ids==c("7-Glia"))
#in total 2859 cells

Glias <- FindVariableFeatures(Glias,verbose = FALSE,selection.method = "vst")
Glias<-ScaleData(Glias)
Glias <- RunPCA(object = Glias, verbose = FALSE,assay="integrated")
ElbowPlot(Glias)
Glias<- RunUMAP(Glias, reduction = "pca", dims = 1:20)
Glias<-RunTSNE(Glias,dims = 1:20)
Glias<-FindNeighbors(Glias,reduction = "pca", dims = 1:20)

Glias.clusters <- FindClusters(Glias,resolution = c(seq(0,0.2,.02)))
#clustree(Glias.clusters@meta.data, prefix = "integrated_snn_res.")
Glias <- FindClusters(Glias, resolution = 0.04) 

Glias$Celltype <- Idents(Glias)
plots <- DimPlot(Glias,reduction="umap",group.by=c("stim","Celltype","sample_id"),combine=FALSE)
plots <- lapply(X = plots, FUN = function(x) + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2), byrow= TRUE, override.aes = list(size=3)))
combined.Glias <- CombinePlots(plots)
ggsave("./Glias/DimPlot_subGlias.pdf", plot = combined.Glias, width = 14, height = 10)

#Find marker genes 
Glias$celltype.stim <- paste(Idents(Glias),Glias$stim,sep="_")
Glias_markers_all_SCTMAST <- FindAllMarkers(Glias, min.pct=0.25, logfc.threshold = 0.25, only.pos = TRUE, assay="SCT",slot="data",test.use="MAST")
Glias_cluster_markers10_SCTMAST <- Glias_markers_all_SCTMAST %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) %>% print(n=3*37) 
write.csv(Glias_cluster_markers10_SCTMAST, "./Glias/subGlias_cluster_markers_top10_SCTMAST.csv")

new.cluster.ids <- c("Ensheathing Glia", "Astrocyte-like Glia", "Neuron:nrv3/Dop2R", "Surface Glia", "Cortex Glia")
names(new.cluster.ids) <- levels(Glias)
Glias <- RenameIdents(Glias,new.cluster.ids)

DimPlot(Glias,label = TRUE,pt.size = 1,reduction="umap",label.size=5, label.box = T)

#plot to see the top10 marker genes
gliatop10=c("CG42613", "CG9657", "Nepl5", "CAH7", "axo", "Ugt35B1", "CG8369", "trol", "CG16743", "CG9377", 
            "CG9394", "Eaat1", "alrm", "Gat", "CG42342", "CG43693", "wun2", "speck", "bbc", "Gs2", 
            "Cys", "CG3168", "CG6910", "CG14687", "Fas2", "SPARC", "CG9336", "CG9691", "Inos", "regucalcin", 
            "fabp", "CG2852", "CG9686", "CG3270", "CG40470", "wrapper", "Obp44a", "MtnA", "CG15201", "apolpp")

DotPlot(Glias,features=gliatop10,assay="SCT",cols=c("blue","red"),scale=T,dot.scale=8,scale.by="size",idents=c("Ensheathing Glia", "Astrocyte-like Glia", "Surface Glia", "Cortex Glia"))+RotatedAxis()+ labs(y="Identites",x="Marker genes",title = "Classification of glial cell-types")

# DEG analysis of the sub-glial cells
Glias$celltype.stim <- paste(Idents(Glias),Glias$stim,sep="_")
DefaultAssay(Glias) <-"SCT"
Glias$cluster_ids <- Idents(Glias)
Idents(Glias) <- "celltype.stim"

ATX_merge12_DE_subCluster_C0_Ensheathing <- FindMarkers(Glias,ident.1 = "Ensheathing_ATX",ident.2 = "Ensheathing_Control",test.use = "MAST", assay = "SCT", slot = "data")
ATX_merge12_DE_subCluster_C1_Astrocyte-like <- FindMarkers(Glias,ident.1 = "Astrocyte-like_ATX",ident.2 = "Astrocyte-like_Control",test.use = "MAST", assay = "SCT", slot = "data")
ATX_merge12_DE_subCluster_C2_Neuron <- FindMarkers(Glias,ident.1 = "Neuron:nrv3/Dop2R_ATX",ident.2 = "Neuron:nrv3/Dop2R_Control",test.use = "MAST", assay = "SCT", slot = "data")
ATX_merge12_DE_subCluster_C3_Surface <- FindMarkers(Glias,ident.1 = "Surface Glia_ATX",ident.2 = "Surface Glia_Control",test.use = "MAST", assay = "SCT", slot = "data")
ATX_merge12_DE_subCluster_C4_Cortex <- FindMarkers(Glias,ident.1 = "Cortex Glia_ATX",ident.2 = "Cortex Glia_Control",test.use = "MAST", assay = "SCT", slot = "data")

ATX_merge12_DE_list <- list("ATX_merge12_C0_Ensheathing" = ATX_merge12_DE_subCluster_C0_Ensheathing,"ATX_merge12_C1_AstrocyteNeuropi"=ATX_merge12_DE_subCluster_C1_Astrocyte-like,"ATX_merge12_C2_Neuron"=ATX_merge12_DE_subCluster_C2_Neuron,"ATX_merge12_C3_Surface"=ATX_merge12_DE_subCluster_C3_Surface,"ATX_merge12_C4_Cortex"=ATX_merge12_DE_subCluster_C4_Cortex)
openxlsx::write.xlsx(ATX_merge12_DE_list,file="./Glias/ATX_merge12_DE_subglial.xlsx",rowNames=TRUE)


MPH_merge12_DE_subCluster_C0_Ensheathing <- FindMarkers(Glias,ident.1 = "Ensheathing_MPH",ident.2 = "Ensheathing_Control",test.use = "MAST", assay = "SCT", slot = "data")
MPH_merge12_DE_subCluster_C1_Astrocyte-like <- FindMarkers(Glias,ident.1 = "Astrocyte-like_MPH",ident.2 = "Astrocyte-like_Control",test.use = "MAST", assay = "SCT", slot = "data")
MPH_merge12_DE_subCluster_C2_Neuron <- FindMarkers(Glias,ident.1 = "Neuron:nrv3/Dop2R_MPH",ident.2 = "Neuron:nrv3/Dop2R_Control",test.use = "MAST", assay = "SCT", slot = "data")
MPH_merge12_DE_subCluster_C3_Surface <- FindMarkers(Glias,ident.1 = "Surface Glia_MPH",ident.2 = "Surface Glia_Control",test.use = "MAST", assay = "SCT", slot = "data")
MPH_merge12_DE_subCluster_C4_Cortex <- FindMarkers(Glias,ident.1 = "Cortex Glia_MPH",ident.2 = "Cortex Glia_Control",test.use = "MAST", assay = "SCT", slot = "data")

MPH_merge12_DE_list <- list("MPH_merge12_C0_Ensheathing" = MPH_merge12_DE_subCluster_C0_Ensheathing,"MPH_merge12_C1_AstrocyteNeuropi"=MPH_merge12_DE_subCluster_C1_Astrocyte-like,"MPH_merge12_C2_Neuron"=MPH_merge12_DE_subCluster_C2_Neuron,"MPH_merge12_C3_Surface"=MPH_merge12_DE_subCluster_C3_Surface,"MPH_merge12_C4_Cortex"=MPH_merge12_DE_subCluster_C4_Cortex)
openxlsx::write.xlsx(MPH_merge12_DE_list,file="./Glias/MPH_merge12_DE_subglial.xlsx",rowNames=TRUE)

####################################################################################################
################################ cell proportions of subglial cells ################################

onlyGlias=subset(Glias,idents=c("Ensheathing Glia", "Astrocyte-like Glia","Surface Glia", "Cortex Glia"))
DimPlot(onlyGlias,label = TRUE,pt.size = 0.5,reduction="umap",label.size=5, label.box = T,label.color="#474747")+labs(x = "UMAP1", y = "UMAP2",title = " ") +NoLegend()

cell.prop<-as.data.frame(prop.table(table(onlyGlias$Celltype,onlyGlias$sample_id),2))
colnames(cell.prop)<-c("cluster","treat","proportion")
prop.plot2 <-ggplot(cell.prop,aes(treat,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="#474747", size=1, linetype="solid"),
        axis.text.y=element_text(size=13,color="black"),
        axis.text.x=element_text(size=13,color="black"),
        axis.text.x.bottom = element_text( angle =90),
        axis.ticks.length=unit(0.5,'cm'),
        axis.title.x=element_text(size=15,color="black"),
        axis.title.y=element_text(size=15,color="black"))+guides(fill=guide_legend(title=NULL))  

ggsave("./subclusters/Glias/subGlias_cellprop_4subs_reps_20220731.pdf", plot = prop.plot2, width =6, height = 8,dpi=600)
####################################################################################################
####################################################################################################
