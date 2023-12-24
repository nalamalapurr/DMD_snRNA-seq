# packages needed
library("dplyr")
library("Seurat")
library("magrittr")
library("clusterProfiler")
library("org.Mm.eg.db")
library("enrichplot")
library("ReactomePA")
library("DOSE")
library("cowplot")
library("tidyverse")

#mdx.bl6.combined<-readRDS("RN_Clustering_61_PCs_CellMarker_NewParams.rds")

DefaultAssay(mdx.bl6.combined)<-"RNA"

x<-data.frame(cluster=0:24,cell_type=c("Neurons (Kcnq3,Lmo7,Rasgef1b)","Neurons (Meg3,Rian,Atp1b1)","Neurons (Rplp0,Rps23,Rps11)","Somatostatin GABAergic neurons (Npy,SstC130073E24Rik)","Astrocytes (Ptprz1,Slc1a3,Cst3)","Neurons (Tmsb4x,Mapt,Dpysl3)","Neurons (Chd7,Lzts1,Slc17a6)","Neurons-ECM (Prox1,Sema5a,Adamts18)","Mixture of GABAergic and glutamatergic neurons (Tshz2,Gria2,Ptprd)","Neurons (Nrxn3,Tcf4,Son)","Neurons (Adarb2,Ntng1,Cnr1","Neurons (Mfap4,Top2a,Gm3764","Neurons (Ntm,Mef2c,Acvr2a)","Cajal-Retzius cells (Snhg11,Cacna2d2,Reln)","Bergmann glial cells (Slc1a2,Adgrv1,Gdf10)","Neurons (mt-Co3,mt-Nd2,mt-Co2)","Microglial Cells (Ptgds,Col3a1,Col1a1)","Astrocytes (Wls,Rmst,Slco1c1)","Choroid Plexus Cells (Ttr,Igfbp2,Trpm3)","Oligodendrocyte Precursor Cells (Pdgfra,Cspg4,Serpine2)","Glutamatergic neurons or neuroblasts (Zfhx3,Zic1,Zic4)","Neurons (Spp1,Apoe,Ctsl)","Microglial cells (Hspa4l,Ccdc39,Myb)","Ependymal cells (Vtn,Cald1,Igfbp7)","Mixture of perivascular and choroid-plexus like cells (Gm12840,Anxa3,Mgp)"))

# run differential expression testing for each cluster to identify upregulated genes in C57BL/6
bl6.markers<-FindMarkers(mdx.bl6.combined,ident.1="BL6",group.by="Overall_Condition",subset.ident=levels(mdx.bl6.combined)[1], only.pos=TRUE)
df.bl6<-head(x=bl6.markers,n=1001L)
df.bl6.lim<-df.bl6[df.bl6$avg_log2FC > 0.25 & df.bl6$p_val_adj<0.05,]

df.bl6.genes<-data.frame("genes"=rownames(df.bl6.lim),"cluster"=rep(subset(x, cluster == 0)$cell_type,length(rownames(df.bl6.lim))))

for(i in 1:24){
bl6.markers<-FindMarkers(mdx.bl6.combined,ident.1="BL6",group.by="Overall_Condition",subset.ident=levels(mdx.bl6.combined)[i+1], only.pos=TRUE)
df.bl6<-head(x=bl6.markers,n=1001L)
df.bl6.lim<-df.bl6[df.bl6$avg_log2FC > 0.25 & df.bl6$p_val_adj<0.05,]

df.bl6.genes.add<-data.frame("genes"=rownames(df.bl6.lim),"cluster"=rep(subset(x, cluster == i)$cell_type,length(rownames(df.bl6.lim))))

df.bl6.genes<-rbind(df.bl6.genes,df.bl6.genes.add)
}

df.bl6.genes.split<-split(df.bl6.genes$genes,f=df.bl6.genes$cluster)

geneid.bl6 <- df.bl6.genes.split %>% map(~{
  
  gene.df <- select(org.Mm.eg.db,
                    keys = .x,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

# GO Term enrichment for each cluster
go.bp.bl6 <- geneid.bl6 %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.cc.bl6 <- geneid.bl6 %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.mf.bl6 <- geneid.bl6 %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

# KEGG Pathway enrichment for each cluster
kegg.bl6 <- geneid.bl6 %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'mmu'
  )
  return(eKEGG)
})


barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {
  ## use *height* to satisy barplot generic definition
  ## actually here is an enrichResult object.
  colorBy <- color
  
  if(is.null(object)){
    df <- data.frame(matrix(ncol=9,nrow=0))
    colnames<-c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
    colnames(df)<-colnames
  }
  else{
    df <- fortify(object, showCategory = showCategory, by = x)
    df$Description <- gsub(' - .*','',df$Description)
    df$p.adjust <- -log10(df$p.adjust)
  }
  #df <- df[c(1:3,9:12,15,16),]
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } 
  else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 FDR') + ylab(NULL)
}

lapply(1:length(geneid.bl6), function(x){
  name <- paste0(sub("\\_.*", "", names(geneid.bl6)[[x]]),"_BL6")
  names <- names(geneid.bl6)[[x]]
  g1 = barplotTerm(go.bp.bl6[[x]], showCategory = 20, title = paste0(names, " GO_BP"), color = 'blue', x = 'p.adjust')
  g2 = barplotTerm(go.cc.bl6[[x]], showCategory = 20, title = paste0(names, " GO_CC"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(go.mf.bl6[[x]], showCategory = 20, title = paste0(names, " GO_MF"), color = 'blue', x = 'p.adjust')
  g4 = barplotTerm(kegg.bl6[[x]], showCategory = 20, title = paste0(names, " KEGG"), color = 'blue', x = 'p.adjust')
  png(file=paste0("GO_KEGG_",name,"_RN_barplot.png"), width=1500, height=600)
  print(plot_grid(g1,g2,g3,g4,nrow=2,rel_widths=c(1,1,1,1),align='h'))
  dev.off()
})

# run differential expression testing for each cluster to identify upregulated genes in mdx52
mdx.markers<-FindMarkers(mdx.bl6.combined,ident.1="MDX",group.by="Overall_Condition",subset.ident=levels(mdx.bl6.combined)[1], only.pos=TRUE)
df.mdx<-head(x=mdx.markers,n=1001L)
df.mdx.lim<-df.mdx[df.mdx$avg_log2FC > 0.25 & df.mdx$p_val_adj<0.05,]

df.mdx.genes<-data.frame("genes"=rownames(df.mdx.lim),"cluster"=rep(subset(x, cluster == 0)$cell_type,length(rownames(df.mdx.lim))))


for(i in 1:24){
mdx.markers<-FindMarkers(mdx.bl6.combined,ident.1="MDX",group.by="Overall_Condition",subset.ident=levels(mdx.bl6.combined)[i+1], only.pos=TRUE)
df.mdx<-head(x=mdx.markers,n=1001L)
df.mdx.lim<-df.mdx[df.bl6$avg_log2FC > 0.25 & df.mdx$p_val_adj<0.05,]

df.mdx.genes.add<-data.frame("genes"=rownames(df.mdx.lim),"cluster"=rep(subset(x, cluster == i)$cell_type,length(rownames(df.mdx.lim))))

df.mdx.genes<-rbind(df.mdx.genes,df.mdx.genes.add)
}

df.mdx.genes.split<-split(df.mdx.genes$genes,f=df.mdx.genes$cluster)

geneid.mdx <- df.mdx.genes.split %>% map(~{
  
  gene.df <- select(org.Mm.eg.db,
                    keys = .x,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

# GO term enrichment for each cluster
go.bp.mdx <- geneid.mdx %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.cc.mdx <- geneid.mdx %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.mf.mdx <- geneid.mdx %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

# KEGG pathway enrichment for each cluster
kegg.mdx <- geneid.mdx %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'mmu'
  )
  return(eKEGG)
})

lapply(1:length(geneid.mdx), function(x){
  
  name <- paste0(sub("\\_.*", "", names(geneid.mdx)[[x]]),"_MDX")
  names <- names(geneid.mdx)[[x]] 
  g1 = barplotTerm(go.bp.mdx[[x]], showCategory = 20, title = paste0(names, " GO_BP"), color = 'blue', x = 'p.adjust')
  g2 = barplotTerm(go.cc.mdx[[x]], showCategory = 20, title = paste0(names, " GO_CC"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(go.mf.mdx[[x]], showCategory = 20, title = paste0(names, " GO_MF"), color = 'blue', x = 'p.adjust')
  g4 = barplotTerm(kegg.mdx[[x]], showCategory = 20, title = paste0(names, " KEGG"), color = 'blue', x = 'p.adjust')
  png(file=paste0("GO_KEGG_",name,"_RN_barplot.png"), width=1500, height=600)
  print(plot_grid(g1,g2,g3,g4,nrow=2,rel_widths=c(1,1,1,1),align='h'))
  dev.off()
})

#Print genes associated with each GO term/KEGG pathway
lapply(1:length(geneid.bl6), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.bl6)[[x]]),"_BL6")
	go.bp.bl6.dt<-data.frame("GO Term" = go.bp.bl6[[x]]@result$Description, "Contributing Genes" = go.bp.bl6[[x]]@result$geneID)
	write.csv(go.bp.bl6.dt, paste0(name,"_GO_BP_RN_Genes.csv"))
})

lapply(1:length(geneid.bl6), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.bl6)[[x]]),"_BL6")
	go.cc.bl6.dt<-data.frame("GO Term" = go.cc.bl6[[x]]@result$Description, "Contributing Genes" = go.cc.bl6[[x]]@result$geneID)
	write.csv(go.cc.bl6.dt, paste0(name,"_GO_CC_RN_Genes.csv"))
})

lapply(1:length(geneid.bl6), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.bl6)[[x]]),"_BL6")
	go.mf.bl6.dt<-data.frame("GO Term" = go.mf.bl6[[x]]@result$Description, "Contributing Genes" = go.mf.bl6[[x]]@result$geneID)
	write.csv(go.mf.bl6.dt, paste0(name,"_GO_MF_RN_Genes.csv"))
})

lapply(1:length(geneid.bl6), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.bl6)[[x]]),"_BL6")
	kegg.bl6.dt<-data.frame("GO Term" = kegg.bl6[[x]]@result$Description, "Contributing Genes" = kegg.bl6[[x]]@result$geneID)
	write.csv(kegg.bl6.dt, paste0(name,"_KEGG_RN_Genes.csv"))
})






lapply(1:length(geneid.mdx), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.mdx)[[x]]),"_MDX")
	go.bp.mdx.dt<-data.frame("GO Term" = go.bp.mdx[[x]]@result$Description, "Contributing Genes" = go.bp.mdx[[x]]@result$geneID)
	write.csv(go.bp.mdx.dt, paste0(name,"_GO_BP_RN_Genes.csv"))
})

lapply(1:length(geneid.mdx), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.mdx)[[x]]),"_MDX")
	go.cc.mdx.dt<-data.frame("GO Term" = go.cc.mdx[[x]]@result$Description, "Contributing Genes" = go.cc.mdx[[x]]@result$geneID)
	write.csv(go.cc.mdx.dt, paste0(name,"_GO_CC_RN_Genes.csv"))
})

lapply(1:length(geneid.mdx), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.mdx)[[x]]),"_MDX")
	go.mf.mdx.dt<-data.frame("GO Term" = go.mf.mdx[[x]]@result$Description, "Contributing Genes" = go.mf.mdx[[x]]@result$geneID)
	write.csv(go.mf.mdx.dt, paste0(name,"_GO_MF_RN_Genes.csv"))
})

lapply(1:length(geneid.mdx), function(x){
	name <- paste0(sub("\\_.*", "", names(geneid.mdx)[[x]]),"_MDX")
	kegg.mdx.dt<-data.frame("GO Term" = kegg.mdx[[x]]@result$Description, "Contributing Genes" = kegg.mdx[[x]]@result$geneID)
	write.csv(kegg.mdx.dt, paste0(name,"_KEGG_RN_Genes.csv"))
})

# Run differential expression testing for upregulated genes in the entire C57BL/6 dataset
bl6.overall.markers<-FindMarkers(mdx.bl6.combined,ident.1="BL6",group.by="Overall_Condition", only.pos=TRUE)
df.bl6.overall<-head(x=bl6.overall.markers,n=1001L)
df.bl6.overall.lim<-df.bl6.overall[df.bl6.overall$avg_log2FC > 0.25 & df.bl6.overall$p_val_adj<0.05,]

df.bl6.overall.genes<-data.frame("genes"=rownames(df.bl6.overall.lim))

geneid.bl6.overall <- df.bl6.overall.genes %>% map(~{
  
  gene.df <- select(org.Mm.eg.db,
                    keys = .x,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

# Overall GO term enrichment
go.bp.bl6.overall <- geneid.bl6.overall %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.cc.bl6.overall <- geneid.bl6.overall %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.mf.bl6.overall <- geneid.bl6.overall %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

# Overall KEGG pathway enrichment
kegg.bl6.overall <- geneid.bl6.overall %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'mmu'
  )
  return(eKEGG)
})


names <- "BL6_Overall_Genes" 
g1 = barplotTerm(go.bp.bl6.overall[[1]], showCategory = 20, title = paste0(names, " GO_BP"), color = 'blue', x = 'p.adjust')
g2 = barplotTerm(go.cc.bl6.overall[[1]], showCategory = 20, title = paste0(names, " GO_CC"), color = 'blue', x = 'p.adjust')
g3 = barplotTerm(go.mf.bl6.overall[[1]], showCategory = 20, title = paste0(names, " GO_MF"), color = 'blue', x = 'p.adjust')
g4 = barplotTerm(kegg.bl6.overall[[1]], showCategory = 20, title = paste0(names, " KEGG"), color = 'blue', x = 'p.adjust')
png(file=paste0("GO_KEGG_BL6_Overall_Genes_RN_barplot.png"), width=1500, height=600)
print(plot_grid(g1,g2,g3,g4,nrow=2,rel_widths=c(1,1,1,1),align='h'))
dev.off()


# Run differential expression testing for upregulated genes in the entire mdx52 dataset
mdx.overall.markers<-FindMarkers(mdx.bl6.combined,ident.1="MDX",group.by="Overall_Condition", only.pos=TRUE)
df.mdx.overall<-head(x=mdx.overall.markers,n=1201L)
df.mdx.overall.lim<-df.mdx.overall[df.mdx.overall$avg_log2FC > 0.25 & df.mdx.overall$p_val_adj<0.05,]

df.mdx.overall.genes<-data.frame("genes"=rownames(df.mdx.overall.lim))

geneid.mdx.overall <- df.mdx.overall.genes %>% map(~{
  
  gene.df <- select(org.Mm.eg.db,
                    keys = .x,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

# Overall GO term enrichment
go.bp.mdx.overall <- geneid.mdx.overall %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.cc.mdx.overall <- geneid.mdx.overall %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

go.mf.mdx.overall <- geneid.mdx.overall %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Mm.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
 
})

# Overall KEGG pathway enrichment
kegg.mdx.overall <- geneid.mdx.overall %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'mmu'
  )
  return(eKEGG)
})

# Plot Overall GO Term and KEGG pathway enrichment results
names <- "MDX_Overall_Genes" 
g1 = barplotTerm(go.bp.mdx.overall[[1]], showCategory = 20, title = paste0(names, " GO_BP"), color = 'red', x = 'p.adjust')
g2 = barplotTerm(go.cc.mdx.overall[[1]], showCategory = 20, title = paste0(names, " GO_CC"), color = 'red', x = 'p.adjust')
g3 = barplotTerm(go.mf.mdx.overall[[1]], showCategory = 20, title = paste0(names, " GO_MF"), color = 'red', x = 'p.adjust')
g4 = barplotTerm(kegg.mdx.overall[[1]], showCategory = 20, title = paste0(names, " KEGG"), color = 'red', x = 'p.adjust')
png(file=paste0("GO_KEGG_MDX_Overall_Genes_RN_barplot.png"), width=1500, height=600)
print(plot_grid(g1,g2,g3,g4,nrow=2,rel_widths=c(1,1,1,1),align='h'))
dev.off()


names <- "MDX_Overall_Genes" 
g1 = barplotTerm(go.bp.mdx.overall[[1]], showCategory = 20, title = paste0(names, " GO_BP"), color = 'red', x = 'p.adjust')
g4 = barplotTerm(kegg.mdx.overall[[1]], showCategory = 20, title = paste0(names, " KEGG"), color = 'red', x = 'p.adjust')
png(file=paste0("GO_KEGG_MDX_Overall_Genes_RN_barplot_RED.png"), width=1500, height=600)
print(plot_grid(g1,g4,nrow=1,rel_widths=c(1,1),align='h'))
dev.off()
