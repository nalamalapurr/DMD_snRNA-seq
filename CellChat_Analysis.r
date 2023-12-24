#Packages needed
library("Seurat")
library("CellChat")
library("patchwork")
library("dplyr")
library("tidyverse")
library("ggplot2")

options(future.globals.maxSize = 6e9)

#mdx.bl6.combined<-readRDS("RN_Clustering_61_PCs_CellMarker_NewParams.rds")

DefaultAssay(mdx.bl6.combined)<-"RNA"

Idents(object = mdx.bl6.combined) <- "Cell_Type_Shortened"

mdx.combined<-subset(x=mdx.bl6.combined, subset = Overall_Condition == "MDX")

cellChat.mdx <- createCellChat(object = mdx.combined, group.by = "Cell_Type", assay = "RNA")

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB

cellChat.mdx@DB <- CellChatDB.use

cellChat.mdx <- subsetData(cellChat.mdx)
future::plan("multicore", workers = 4)

cellChat.mdx <- identifyOverExpressedGenes(cellChat.mdx)
cellChat.mdx <- identifyOverExpressedInteractions(cellChat.mdx)
cellChat.mdx <- projectData(cellChat.mdx, PPI.mouse)

cellChat.mdx <- computeCommunProb(cellChat.mdx, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellChat.mdx <- filterCommunication(cellChat.mdx, min.cells = 10)

cellChat.mdx <- computeCommunProbPathway(cellChat.mdx)

cellChat.mdx <- aggregateNet(cellChat.mdx)

groupSize <- as.numeric(table(cellChat.mdx@idents))
par(mfrow = c(1,2), xpd=TRUE)


mat <- cellChat.mdx@net$weight

bl6.combined<-subset(x=mdx.bl6.combined, subset = Overall_Condition == "BL6")

cellChat.bl6 <- createCellChat(object = bl6.combined, group.by = "Cell_Type", assay = "RNA")

cellChat.bl6@DB <- CellChatDB.use

cellChat.bl6 <- subsetData(cellChat.bl6)
future::plan("multicore", workers = 4)

cellChat.bl6 <- identifyOverExpressedGenes(cellChat.bl6)
cellChat.bl6 <- identifyOverExpressedInteractions(cellChat.bl6)
cellChat.bl6 <- projectData(cellChat.bl6, PPI.mouse)

cellChat.bl6 <- computeCommunProb(cellChat.bl6, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellChat.bl6 <- filterCommunication(cellChat.bl6, min.cells = 10)

cellChat.bl6 <- computeCommunProbPathway(cellChat.bl6)

cellChat.bl6 <- aggregateNet(cellChat.bl6)

groupSize <- as.numeric(table(cellChat.bl6@idents))
par(mfrow = c(1,2), xpd=TRUE)

mat <- cellChat.bl6@net$weight

#Merge the two cellChat objects and create a heatmap comparing cell-cell interactions between conditions
object.list <- list(BL6 = cellChat.bl6, MDX = cellChat.mdx)

mergedCellChat<-mergeCellChat(object.list, add.names = names(object.list))

gg1_n <- compareInteractions(mergedCellChat, show.legend = F, group = c(1,2),group.facet="Cell_Type")
gg2_n <- compareInteractions(mergedCellChat, show.legend = F, group = c(1,2), measure = "weight")

png(file="CellChat_Interactions.png", width=1000, height=1000)
gg1_n+gg2_n
dev.off()

png(file="CellChat_Interactions.png", width=1000, height=1000)
gg1_n
dev.off()

png(file="CellChat_Interactions_Weighted.png", width=1000, height=1000)
gg2_n
dev.off()

gg1_m<-netVisual_heatmap(mergedCellChat)
gg2_m<-netVisual_heatmap(mergedCellChat, measure = "weight")

png(file="CellChat_netHeatmap.png", width=1000, height=1000)
gg1_m+gg2_m
dev.off()

pos.dataset = "MDX"
features.name = pos.dataset
mergedCellChat <- identifyOverExpressedGenes(mergedCellChat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

net <- netMappingDEG(mergedCellChat, features.name = features.name)
net.up <- subsetCommunication(mergedCellChat, net = net, datasets = "MDX",ligand.logFC = 0.1, receptor.logFC = 0.1)
net.down <- subsetCommunication(mergedCellChat, net = net, datasets = "BL6",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, mergedCellChat)
gene.down <- extractGeneSubsetFromPair(net.down, mergedCellChat)


pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(mergedCellChat, pairLR.use = pairLR.use.up, sources.use = 1:21, targets.use = 1:21, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(mergedCellChat, pairLR.use = pairLR.use.down, sources.use = 1:21, targets.use = 1:21, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object

png(file="CellChat_netBubble.png", width=5000, height=1000)
gg1 + gg2
dev.off()


png(file="CellChat_chord_mdx_up.png", width=2500, height=2500)
netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = net.up, lab.cex = 3.0, small.gap = 3.5)
dev.off()

png(file="CellChat_chord_mdx_down.png", width=2500, height=2500)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = net.down, lab.cex = 2.0, small.gap = 3.5)
dev.off()

png(file="CellChat_enrich_plot_mdx_up.png", width=1000, height=1000)
computeEnrichmentScore(net.up, species = 'mouse')
dev.off()

png(file="CellChat_enrich_plot_mdx_down.png", width=1000, height=1000)
computeEnrichmentScore(net.down, species = 'mouse')
dev.off()
