#Packages required
library(Seurat)
library(dplyr)
library(HGNChelper)
library(openxlsx)

#Read CellRanger Output into SeuratObject for all 8 datasets
mdx8<-Read10X("MDX8")
mdx5<-Read10X("MDX5")
mdx3<-Read10X("MDX3")
mdx1<-Read10X("MDX1")
mdx8<-CreateSeuratObject(counts=mdx8,min.cells=3,min.features=200)
mdx1<-CreateSeuratObject(counts=mdx1,min.cells=3,min.features=200)
mdx3<-CreateSeuratObject(counts=mdx3,min.cells=3,min.features=200)
mdx5<-CreateSeuratObject(counts=mdx5,min.cells=3,min.features=200)

bl67<-Read10X("BL6_7")
bl66<-Read10X("BL6_6")
bl64<-Read10X("BL6_4")
bl63<-Read10X("BL6_3")
bl67<-CreateSeuratObject(counts=bl67,min.cells=3,min.features=200)
bl66<-CreateSeuratObject(counts=bl66,min.cells=3,min.features=200)
bl64<-CreateSeuratObject(counts=bl64,min.cells=3,min.features=200)
bl63<-CreateSeuratObject(counts=bl63,min.cells=3,min.features=200)


#Add percent mitochondrial gene metadata
mdx8[["percent.mt"]] <- PercentageFeatureSet(mdx8, pattern = "^mt-")
mdx5[["percent.mt"]] <- PercentageFeatureSet(mdx5, pattern = "^mt-")
mdx3[["percent.mt"]] <- PercentageFeatureSet(mdx3, pattern = "^mt-")
mdx1[["percent.mt"]] <- PercentageFeatureSet(mdx1, pattern = "^mt-")

bl67[["percent.mt"]] <- PercentageFeatureSet(bl67, pattern = "^mt-")
bl66[["percent.mt"]] <- PercentageFeatureSet(bl66, pattern = "^mt-")
bl64[["percent.mt"]] <- PercentageFeatureSet(bl64, pattern = "^mt-")
bl63[["percent.mt"]] <- PercentageFeatureSet(bl63, pattern = "^mt-")


#Filter out dying cells
mdx8 <- subset(mdx8, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
mdx5 <- subset(mdx5, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
mdx3 <- subset(mdx3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
mdx1 <- subset(mdx1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

bl67 <- subset(bl67, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
bl66 <- subset(bl66, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
bl64 <- subset(bl64, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
bl63 <- subset(bl63, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

#Normalize the data and find variable features for each dataset
#Add metadata for individual dataset name and overall condition for later use
mdx1<-NormalizeData(mdx1)
mdx1<-FindVariableFeatures(mdx1,selection.method="vst",nfeatures=2000)

cond_mdx1<-rep("MDX1",length(colnames(mdx1)))
MDX<-rep("MDX",length(colnames(mdx1)))
mdx1<-AddMetaData(object=mdx1, metadata = cond_mdx1, col.name="Condition")
mdx1<-AddMetaData(object=mdx1, metadata = MDX, col.name="Overall_Condition")

mdx3<-NormalizeData(mdx3)
mdx3<-FindVariableFeatures(mdx3,selection.method="vst",nfeatures=2000)

cond_mdx3<-rep("MDX3",length(colnames(mdx3)))
MDX<-rep("MDX",length(colnames(mdx3)))
mdx3<-AddMetaData(object=mdx3, metadata = cond_mdx3, col.name="Condition")
mdx3<-AddMetaData(object=mdx3, metadata = MDX, col.name="Overall_Condition")


mdx5<-NormalizeData(mdx5)
mdx5<-FindVariableFeatures(mdx5,selection.method="vst",nfeatures=2000)

cond_mdx5<-rep("MDX5",length(colnames(mdx5)))
MDX<-rep("MDX",length(colnames(mdx5)))
mdx5<-AddMetaData(object=mdx5, metadata = cond_mdx5, col.name="Condition")
mdx5<-AddMetaData(object=mdx5, metadata = MDX, col.name="Overall_Condition")


mdx8<-NormalizeData(mdx8)
mdx8<-FindVariableFeatures(mdx8,selection.method="vst",nfeatures=2000)

cond_mdx8<-rep("MDX8",length(colnames(mdx8)))
MDX<-rep("MDX",length(colnames(mdx8)))
mdx8<-AddMetaData(object=mdx8, metadata = cond_mdx8, col.name="Condition")
mdx8<-AddMetaData(object=mdx8, metadata = MDX, col.name="Overall_Condition")


bl63<-NormalizeData(bl63)
bl63<-FindVariableFeatures(bl63,selection.method="vst",nfeatures=2000)

cond_bl63<-rep("BL63",length(colnames(bl63)))
BL6<-rep("BL6",length(colnames(bl63)))
bl63<-AddMetaData(object=bl63, metadata = cond_bl63, col.name="Condition")
bl63<-AddMetaData(object=bl63, metadata = BL6, col.name="Overall_Condition")

bl63<-NormalizeData(bl63)
bl63<-FindVariableFeatures(bl63,selection.method="vst",nfeatures=2000)

cond_bl64<-rep("BL64",length(colnames(bl64)))
BL6<-rep("BL6",length(colnames(bl64)))
bl64<-AddMetaData(object=bl64, metadata = cond_bl64, col.name="Condition")
bl64<-AddMetaData(object=bl64, metadata = BL6, col.name="Overall_Condition")

bl66<-NormalizeData(bl66)
bl66<-FindVariableFeatures(bl66,selection.method="vst",nfeatures=2000)

cond_bl66<-rep("BL66",length(colnames(bl66)))
BL6<-rep("BL6",length(colnames(bl66)))
bl66<-AddMetaData(object=bl66, metadata = cond_bl66, col.name="Condition")
bl66<-AddMetaData(object=bl66, metadata = BL6, col.name="Overall_Condition")

bl67<-NormalizeData(bl67)
bl67<-FindVariableFeatures(bl67,selection.method="vst",nfeatures=2000)

cond_bl67<-rep("BL67",length(colnames(bl67)))
BL6<-rep("BL6",length(colnames(bl67)))
bl67<-AddMetaData(object=bl67, metadata = cond_bl67, col.name="Condition")
bl67<-AddMetaData(object=bl67, metadata = BL6, col.name="Overall_Condition")

#Selection of integration features for integration
features<-SelectIntegrationFeatures(object.list=c(mdx1,mdx3,mdx5,mdx8,bl63,bl64,bl66,bl67))

#Integrate data
mdx.bl6.anchors<-FindIntegrationAnchors(object.list=c(mdx1,mdx3,mdx5,mdx8,bl63,bl64,bl66,bl67),anchor.features=features)
mdx.bl6.combined<-IntegrateData(anchorset=mdx.bl6.anchors)

DefaultAssay(mdx.bl6.combined) <- "integrated"

#Scale Data, Run PCA and UMAP dimension reduction, and cluster using default Seurat pipeline
mdx.bl6.combined <- ScaleData(mdx.bl6.combined,verbose=FALSE)
mdx.bl6.combined <- RunPCA(mdx.bl6.combined, npcs = 75, verbose=FALSE)


## Code used for JackStraw analysis, no need to run more than once
mdx.bl6.combined <- JackStraw(mdx.bl6.combined, num.replicate = 100,dims=75)
mdx.bl6.combined <- ScoreJackStraw(mdx.bl6.combined, dims = 1:75)

# png(file="Jackstraw_Plot_NewParameters.png", width=2550, height=1100)
# JackStrawPlot(mdx.bl6.combined, dims = 1:75)
# dev.off()

mdx.bl6.combined <- RunUMAP(mdx.bl6.combined, reduction="pca", dims=1:69)
mdx.bl6.combined <- FindNeighbors(mdx.bl6.combined, reduction="pca", dims=1:69)
mdx.bl6.combined <- FindClusters(mdx.bl6.combined, resolution = 0.5)

x<-data.frame(cluster=0:24,cell_type=c("Neurons (Kcnq3,Lmo7,Rasgef1b)","Neurons (Meg3,Rian,Atp1b1)","Neurons (Rplp0,Rps23,Rps11)","Somatostatin GABAergic neurons (Npy,SstC130073E24Rik)","Astrocytes (Ptprz1,Slc1a3,Cst3)","Neurons (Tmsb4x,Mapt,Dpysl3)","Neurons (Chd7,Lzts1,Slc17a6)","Neurons-ECM (Prox1,Sema5a,Adamts18)","Mixture of GABAergic and glutamatergic neurons (Tshz2,Gria2,Ptprd)","Neurons (Nrxn3,Tcf4,Son)","Neurons (Adarb2,Ntng1,Cnr1","Neurons (Mfap4,Top2a,Gm3764","Neurons (Ntm,Mef2c,Acvr2a)","Cajal-Retzius cells (Snhg11,Cacna2d2,Reln)","Bergmann glial cells (Slc1a2,Adgrv1,Gdf10)","Neurons (mt-Co3,mt-Nd2,mt-Co2)","Microglial Cells (Ptgds,Col3a1,Col1a1)","Astrocytes (Wls,Rmst,Slco1c1)","Choroid Plexus Cells (Ttr,Igfbp2,Trpm3)","Oligodendrocyte Precursor Cells (Pdgfra,Cspg4,Serpine2)","Glutamatergic neurons or neuroblasts (Zfhx3,Zic1,Zic4)","Neurons (Spp1,Apoe,Ctsl)","Microglial cells (Hspa4l,Ccdc39,Myb)","Ependymal cells (Vtn,Cald1,Igfbp7)","Mixture of perivascular and choroid-plexus like cells (Gm12840,Anxa3,Mgp)"))

#Assign cell type annotation as metadata to SeuratObject
mdx.bl6.combined@meta.data$Cell_Type = ""
for(j in unique(x$cluster)){
  cl_type = x[x$cluster==j,]; 
  mdx.bl6.combined@meta.data$Cell_Type[mdx.bl6.combined@meta.data$seurat_clusters == j] = as.character(cl_type$cell_type[1])
}


y<-data.frame(cluster=0:24,cell_type=c("Neu (Kcnq3,Lmo7,Rasgef1b)","Neu (Meg3,Rian,Atp1b1)","Neu (Rplp0,Rps23,Rps11)","Sst GABA neu (Npy,SstC130073E24Rik)","Astro (Ptprz1,Slc1a3,Cst3)","Neu (Tmsb4x,Mapt,Dpysl3)","Neu (Chd7,Lzts1,Slc17a6)","Neu-ECM (Prox1,Sema5a,Adamts18)","GABA + glu neu (Tshz2,Gria2,Ptprd)","Neu (Nrxn3,Tcf4,Son)","Neu (Adarb2,Ntng1,Cnr1","Neu (Mfap4,Top2a,Gm3764","Neu (Ntm,Mef2c,Acvr2a)","C-R cells (Snhg11,Cacna2d2,Reln)","B-G cells (Slc1a2,Adgrv1,Gdf10)","Neu (mt-Co3,mt-Nd2,mt-Co2)","Microglia (Ptgds,Col3a1,Col1a1)","Astro (Wls,Rmst,Slco1c1)","C-P Cells (Ttr,Igfbp2,Trpm3)","OPCs (Pdgfra,Cspg4,Serpine2)","Glu neu + neuroblasts (Zfhx3,Zic1,Zic4)","Neu (Spp1,Apoe,Ctsl)","Microglia (Hspa4l,Ccdc39,Myb)","Ependymal cells (Vtn,Cald1,Igfbp7)","Perivascular and C-P like cells (Gm12840,Anxa3,Mgp)"))

#Assign cell type annotation as metadata to SeuratObject
mdx.bl6.combined@meta.data$Cell_Type_Shortened = ""
for(j in unique(x$cluster)){
  cl_type = y[y$cluster==j,]; 
  mdx.bl6.combined@meta.data$Cell_Type_Shortened[mdx.bl6.combined@meta.data$seurat_clusters == j] = as.character(cl_type$cell_type[1])
}


png(file="UMAP_Cluster_GeneMarkers.png", width=1500, height=600)
DimPlot(mdx.bl6.combined, reduction = "umap", split.by = "Overall_Condition", group.by = 'Cell_Type')        
dev.off()

# saveRDS(mdx.bl6.combined, file = "RN_Clustering_61_PCs_CellMarker_NewParams.rds")
