#Description: Analysis of Nippostrongylus Enterocytes by Cell Chat 
# Main Question: How are Enterocytes from Nippo infected mice different from comparable Naive enterocytes and other comparable Nippostrongylus enterocytes?
#Figures Figure 5E and Figure 5F Figure 5D
##Load the Libraries----
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(cowplot)
library(scCustomize)
library(forcats)
library(ggbreak)
library(presto)
library(ggforce)
library(glmGamPoi)
library(future)
library(dsb)
library(RColorBrewer)
library(monocle3)
library(gt)
library(patchwork)
#devtools::install_github("jinworks/CellChat")
library(CellChat)
library(NMF)
library(circlize)
library(ComplexHeatmap)
library(SeuratDisk)
options(future.globals.maxSize = 5000 * 1024^2)



#directory set up----

getwd() #/home/hartandrew
setwd("/path/to/analysis/directory")
out <- "/path/to/analysis/directory/Output"
version <- "/path/to/analysis/directory/Version"
seurat <- "/path/to/analysis/directory/Seurat_Files"
images <- "/path/to/analysis/directory/Images"
CSV <- "/path/to/analysis/directory/CSV"


#Color Picking-----
colors8 <- c("#980339","#980339", "#379FFB", "#379FFB", "#379FFB", "#D6A206", "#D6A206","#027102","#027102","#6313B7" )


#Load the Objects ----
library(qs2)
Ileum <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2_clustered_annotated_Ileum.qs2")
Ileum[["RNA"]] <-  as(Ileum[["RNA"]], "Assay5")
unique(Ileum$FinestCellType)
unique(Ileum$FineCellType)
unique(Ileum$secondlevel)
Ileum$CellType <- Ileum$secondlevel


Idents(Ileum) <- "SampleType"
unique(Ileum$SampleType)

Ileum  <- subset(x = Ileum, idents = c("Ileum Lamina Propria","Ileum Intestinal Epithelial Cells"  ))

#Data parsing----
Idents(Ileum) <- "InfectionStatus"

Naive_Ileum <- subset(x = Ileum, idents = "Naive")
Nippostrongylus_Ileum <- subset(x = Ileum, idents = "Nippostrongylus")



#Create input Files and CellChat DB----

Ileum.Nippostrongylus.input<- Nippostrongylus_Ileum[["RNA"]]$data
Ileum.naive.input<- Naive_Ileum[["RNA"]]$data


labels <-Nippostrongylus_Ileum$CellType
Ileum.Nippostrongylus.meta <- data.frame(labels = labels, row.names = names(labels)) 
Ileum.Nippostrongylus.meta$samples <- Nippostrongylus_Ileum$Mouse
Ileum.Nippostrongylus.cellchat <- createCellChat(object = Ileum.Nippostrongylus.input, meta = Ileum.Nippostrongylus.meta, group.by = "labels")

labels <-Naive_Ileum$CellType
Ileum.naive.meta <- data.frame(labels = labels, row.names = names(labels)) 
Ileum.naive.meta$samples <- Naive_Ileum$Mouse
Ileum.naive.cellchat <- createCellChat(object = Ileum.naive.input, meta = Ileum.naive.meta, group.by = "labels")

#Set the database to Use
CellChatDB <-  CellChatDB.mouse
showDatabaseCategory(CellChatDB) #visualize the categories of interactions

Ileum.naive.cellchat@DB <- CellChatDB
Ileum.Nippostrongylus.cellchat@DB <- CellChatDB


#Run Analysis Pt 1----
Ileum.Nippostrongylus.cellchat <- subsetData(Ileum.Nippostrongylus.cellchat)
Ileum.Nippostrongylus.cellchat <- identifyOverExpressedGenes(Ileum.Nippostrongylus.cellchat)
Ileum.Nippostrongylus.cellchat <- identifyOverExpressedInteractions(Ileum.Nippostrongylus.cellchat)
# The number of highly variable ligand-receptor pairs used for signaling inference is 2055

Ileum.naive.cellchat <- subsetData(Ileum.naive.cellchat)
Ileum.naive.cellchat <- identifyOverExpressedGenes(Ileum.naive.cellchat)
Ileum.naive.cellchat <- identifyOverExpressedInteractions(Ileum.naive.cellchat)
# The number of highly variable ligand-receptor pairs used for signaling inference is 1752



#Filter 
future::plan("default")
library(BiocParallel)
SnowParam(workers = 60, type = "SOCK")
Ileum.Nippostrongylus.cellchat <- computeCommunProb(Ileum.Nippostrongylus.cellchat, type = "triMean", nboot = 50)
Ileum.Nippostrongylus.cellchat <- filterCommunication(Ileum.Nippostrongylus.cellchat, min.cells = 4)

Ileum.naive.cellchat <- computeCommunProb(Ileum.naive.cellchat, type = "triMean", nboot = 50)
Ileum.naive.cellchat <- filterCommunication(Ileum.naive.cellchat, min.cells = 4)



#Infer the cell-cell communication at a signaling pathway level
#Calculate the aggregated cell-cell communication network
Ileum.Nippostrongylus.cellchat <- computeCommunProbPathway(Ileum.Nippostrongylus.cellchat)
Ileum.Nippostrongylus.cellchat <- aggregateNet(Ileum.Nippostrongylus.cellchat)

Ileum.naive.cellchat <- computeCommunProbPathway(Ileum.naive.cellchat)
Ileum.naive.cellchat <- aggregateNet(Ileum.naive.cellchat)


#Compute network centrality scores
Ileum.Nippostrongylus.cellchat <- netAnalysis_computeCentrality(Ileum.Nippostrongylus.cellchat, slot.name = "netP")
Ileum.naive.cellchat <- netAnalysis_computeCentrality(Ileum.naive.cellchat, slot.name = "netP")

#Visualize Interactions Circle Plot ----
groupSize <- as.numeric(table(Ileum.Nippostrongylus.cellchat@idents))
netVisual_circle(Ileum.Nippostrongylus.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

groupSize <- as.numeric(table(Ileum.naive.cellchat@idents))
netVisual_circle(Ileum.naive.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")



Ileum.Nippostrongylus.cellchat.df <- subsetCommunication(Ileum.Nippostrongylus.cellchat) 
Ileum.naive.cellchat.df <- subsetCommunication(Ileum.naive.cellchat) 


#Combine to perform Comparisons----
#cellchat object list
object.list <- list(Naive_Ileum = Ileum.naive.cellchat, Nippostrongylus_Ileum = Ileum.Nippostrongylus.cellchat)
save(object.list, file = paste0(out,"cellchat_object.list_unfiltered_NippostrongylusAnalysis.RData"))
#I will create the object list that is not restrained by overlapping cell types
object.list2 <- list(Naive_Ileum = Ileum.naive.cellchat, Nippostrongylus_Ileum = Ileum.Nippostrongylus.cellchat)

#Identify CellTypes to keep - If you include cell types which are not present in different groups, you cannot perform the comparisons

b <-as.character(sort(unique(Ileum.Nippostrongylus.cellchat.df$source)))
c<- as.character(sort(unique(Ileum.naive.cellchat.df$source)))

#Identify the cell types tht all MLN samples share
overlap <- Reduce(intersect, list( b, c))

object.list2[[1]] <- subsetCellChat(object.list2[[1]], idents.use = overlap)
object.list2[[2]] <- subsetCellChat(object.list2[[2]], idents.use = overlap)


cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
cellchat2 <- mergeCellChat(object.list2, add.names = names(object.list2), cell.prefix = TRUE)
cellchat2 <- liftCellChat(cellchat2)
# Users can now export the merged CellChat object and the list of the two separate objects for later use

save(object.list2, file = paste0(out,"cellchat_object.list_filtered_NippostrongylusAnalysis.RData"))




#Begin Visualization ----
pathways.show <- c("IL10", "IFN-II")
netVisual_individual(object.list[[1]], targets.use = c("Stem Cells", "Transit Amplifying Cells", "Enterocytes"), signaling = pathways.show,  layout = "circle")
netVisual_individual(object.list[[2]], targets.use = c("Stem Cells", "Transit Amplifying Cells", "Enterocytes"), signaling = pathways.show,  layout = "circle")

# Total interactions
gg1 <- compareInteractions(cellchat2, show.legend = F, group = c(1, 2),  size.text = 14, x.lab.rot = T, angle.x = 45)+ scale_fill_manual(values = colors8)
gg2 <- compareInteractions(cellchat2, show.legend = F, group = c(1, 2), measure = "weight",  size.text = 14,  x.lab.rot = T, angle.x = 45)+ scale_fill_manual(values = colors8)
gg1 + scale_fill_manual(values = colors8  , breaks = c( "Naive_Ileum", "Nippostrongylus_Ileum" ))
gg1 + gg2

#Heatmap 
gg1 <- netVisual_heatmap(cellchat2, comparison = c("Naive_Ileum", "Nippostrongylus_Ileum"), title.name = "Nippostrongylus vs Naive Ileum")
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat2,comparison = c("Naive_Ileum","Nippostrongylus_Ileum"), title.name = "Nippostrongylus vs Naive Ileum", measure = "weight")
gg1 + gg2




#To better control the node size and edge weights of the inferred networks across different datasets, CellChat computes the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.
#Look at the number of interactions coming from each cell type in each condition -  a good way to visualize comparisons
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfcol= c(3,1))
for (i in 1:3) {
  plot <- netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F,
                           top = 0.1, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfcol= c(1,1))
netVisual_circle(object.list[["Naive_Ileum"]]@net$count, weight.scale = T, label.edge= F,
                 top = 0.1, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)["Naive_Ileum"]))


# Aggregate Scatter ----
# in the below code you can separately graph for each cell type its ingoing and outgoing relative number of interactions and compare across plots
object.list[[1]]@netP$pathways # access pathways
levels(object.list[[1]]@idents) # Access cell types
num.link <- lapply(object.list[1:2], function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

# Combine all values into a single vector
num.link.vec <- unlist(num.link)

# Now compute min and max safely
weight.MinMax <- c(min(num.link.vec), max(num.link.vec))
gg <- list()
for (i in 1:2) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list2[[i]], title = names(object.list2)[i], weight.MinMax = weight.MinMax) + scale_y_continuous(limits = c(0,50)) + scale_x_continuous(limits = c(0,75))
}
patchwork::wrap_plots(plots = gg) 






library(reticulate)
use_python("/home/hartandrew/.conda/envs/sccellfie_env/bin/python", required = TRUE)
ptm = Sys.time()
set.seed(40)
cellchat2 <- computeNetSimilarityPairwise(cellchat2, type = "functional")
cellchat2 <- netEmbedding(cellchat2, type = "functional")
cellchat2 <- netClustering(cellchat2, type = "functional")


netVisual_embeddingPairwise(cellchat2, type = "functional", label.size = 3.5)

unique(Ileum$secondlevel)



gg1 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  targets.use = "Enterocytes",
               cutoff.pvalue = 0.01,
               stacked = T, do.stat = TRUE, comparison = c(2,1))
gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  targets.use = "Enterocytes", stacked = F, do.stat = TRUE, comparison = c(2,1))
gg1 + gg2
svg(paste0(images, "/Late_enterocyte_targets_Ranknet_Nippostrongylus_v_Naive.svg"), 
    width = 10, height = 10)
print(gg1 + gg2)
dev.off()



gg1$data <- subset(gg1$data, pvalues < 0.05)
gg1
svg(paste0(images, "/Late_enterocyte_targets_Ranknet_Nippostrongylus_v_Naive_p_0.05.svg"), 
    width = 10, height = 10)
print(gg1 + gg2)
dev.off()


gg1 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  targets.use = "Enterocytes",
               cutoff.pvalue = 0.01, slot.name = "net",
               stacked = T, do.stat = TRUE, comparison = c(2,1))
gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  targets.use = "Enterocytes", stacked = F, do.stat = TRUE, comparison = c(2,1))
gg1 + gg2
svg(paste0(images, "/Late_enterocyte_targets_Ranknet_Nippostrongylus_v_Naive_net.svg"), 
    width = 10, height = 10)
print(gg1 + gg2)
dev.off()

gg1$data <- subset(gg1$data, pvalues < 0.05)
gg1
svg(paste0(images, "/Late_enterocyte_targets_Ranknet_Nippostrongylus_v_Naive_p_0.05_net.svg"), 
    width = 10, height = 12)
print(gg1 )
dev.off()

netVisual_chord_gene(object.list2[[2]], 
                     sources.use = NULL,
                     signaling = "IL4",
                     targets.use = c("Enterocytes"), lab.cex = 0.5,legend.pos.y = 30)

svg(paste0(images, "/enterocyte_targets_Il4_Nippostrongylus_circle.svg"), 
    width = 10, height = 12)

netVisual_individual(
  object.list2[[2]],
  signaling = "IL4",
  sources.use = NULL,
  targets.use = c("Enterocytes"),
  edge.width.max = 90,
  layout = "circle" 

)
dev.off()

PAirLR <- data.frame(interaction_name = c(""))
netVisual_aggregate(object.list2[[2]], 
                    signaling = "IL4",
                    sources.use = NULL,
                    targets.use = c("Enterocytes"),
                    layout = "circle" ,
                    show.legend = TRUE)

netVisual_aggregate(object.list2[[2]], 
                    signaling = "27HC",
                    sources.use = NULL,
                    targets.use = c("Enterocytes"),
                    layout = "circle" ,
                    show.legend = TRUE)


netVisual_chord_cell(object.list2[[2]], 
                    signaling = "IL4",
                    sources.use = NULL,
                    targets.use = c("Enterocytes"),
                    #layout = "circle" ,
                    show.legend = TRUE)
#Figure 5E
svg(paste0(images, "/enterocyte_targets_Il4_Nippostrongylus_chord_gene.svg"), 
    width = 10, height = 10)
netVisual_chord_gene(
  object.list2[[2]],
  signaling = "IL4",
  sources.use = NULL,
  targets.use = c("Enterocytes"),
  legend.pos.y = 30,
  title.name = "IL4 Signaling: Specific L-R Pairs to Enterocytes"
  # By default, netVisual_chord_gene colors ribbons by the L-R pair
)
dev.off()
netVisual_aggregate(object.list2[[2]], 
                    signaling = "27HC",
                    sources.use = NULL,
                    targets.use = c("Enterocytes"),
                    layout = "circle" ,
                    show.legend = TRUE)

netVisual_chord_gene(
  object.list2[[2]],
  signaling = c("IL4", "27HC"),
  sources.use = NULL,
  targets.use = c("Enterocytes"),
  legend.pos.y = 30,
  title.name = "IL4 Signaling: Specific L-R Pairs to Enterocytes"
  # By default, netVisual_chord_gene colors ribbons by the L-R pair
)

object.list2[[2]]@netP$pathways$IL4

gg1 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  sources.use = "Enterocytes",
               cutoff.pvalue = 0.01,
               stacked = T, do.stat = TRUE, comparison = c(2,1))
gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  sources.use = "Enterocytes", stacked = F, do.stat = TRUE, comparison = c(2,1))
gg1 + gg2
svg(paste0(images, "/Late_enterocyte_sources_Ranknet_Nippostrongylus_v_Naive.svg"), 
    width = 10, height = 10)
print(gg1 + gg2)
dev.off()



gg1$data <- subset(gg1$data, pvalues < 0.05)
gg1
svg(paste0(images, "/Late_enterocyte_sources_Ranknet_Nippostrongylus_v_Naive_p_0.05.svg"), 
    width = 10, height = 10)
print(gg1 + gg2)
dev.off()


gg1 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  sources.use = "Enterocytes",
               cutoff.pvalue = 0.01, slot.name = "net",
               stacked = T, do.stat = TRUE, comparison = c(2,1))
gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  sources.use = "Enterocytes", stacked = F, do.stat = TRUE, comparison = c(2,1))
gg1 + gg2
svg(paste0(images, "/Late_enterocyte_sources_Ranknet_Nippostrongylus_v_Naive_net.svg"), 
    width = 10, height = 10)
print(gg1 + gg2)
dev.off()

gg1$data <- subset(gg1$data, pvalues < 0.05)
gg1
svg(paste0(images, "/Late_enterocyte_sources_Ranknet_Nippostrongylus_v_Naive_p_0.05_net.svg"), 
    width = 10, height = 12)
print(gg1 )
dev.off()

# Figure 5F
object.list2[[2]]@netP$pathways
svg(paste0(images, "/enterocyte_targets_FGF_Nippostrongylus_chord_gene.svg"), 
    width = 10, height = 10)
netVisual_chord_gene(
  object.list2[[2]],
  signaling = c("FGF", "CX3C"),
  sources.use = c("Enterocytes"),
  legend.pos.y = 30,
  title.name = "FGF Signaling: Specific L-R Pairs to Enterocytes"
  # By default, netVisual_chord_gene colors ribbons by the L-R pair
)
dev.off()


















