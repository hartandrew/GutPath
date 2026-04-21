# Description. This script is the work of an iterative process of manually annotating cell types as they appeared in the data. The final version of this script includes mapping previously annotated cells onto the current objects and fine tuning the cell labels. 
---
title: "Analysis_2_Annotation"
author: "AndrewHart"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

##Load the Libraries----
```{r, warning=FALSE, error=FALSE}
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
library(SingleR)
library(celldex)
library(harmony)
#library(Azimuth)
library(glmGamPoi)
library(future)
library(dsb)
library(SeuratDisk)
library(RColorBrewer)
library(monocle3)
library(gt)
library(SeuratWrappers)
library(qs2)
options(future.globals.maxSize = 20000 * 1024^2)
```

#directory set up
```{r}
setwd("/path/to/analysis/directory")
out <- "/path/to/analysis/directory/Output"
version <- "/path/to/analysis/directory/Version"
seurat <- "/path/to/analysis/directory/Seurat_Files"
images <- "/path/to/analysis/directory/Images"
CSV <- "/path/to/analysis/directory/CSV"
```


```{r}
colors8 <- c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
colors9 <- c("#332288", "#88CCEE", "#117733","#44AA99", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
```


#This Script is iterative. Notes about how the annotations were performed can be found in sections of the script that ARE NOT EVALUATED here.

#Load the Objects 
```{r}
MLN <- qs2::qs_read("/path/to/analysis/directory/Seurat_Files/Analysis1_MLN_clustered_regressed.qs2")
Ileum <-  qs2::qs_read("/path/to/analysis/directory/Seurat_Files/Analysis1_Ileum_clustered_regressed.qs2")
#Load objects which have the appropriate annotations and fine clustering previously performed 
MLN_old <-  qs_read("/path/to/analysis/directory/Seurat_Files/MLN_CellxGene_FineFilt2_Feb2025.qs2")
Ileum_old <-  qs2::qs_read("/path/to/analysis/directory/Seurat_Files/Jan2025_Completed_annotation_gsea_Ileum_object.qs2")
```


```{r}
# Convert back to Seurat v5 compatible assays
MLN[["RNA"]] <-  as(MLN[["RNA"]], "Assay5")
MLN[["ADT"]] <-  as(MLN[["ADT"]], "Assay5")
Ileum[["RNA"]] <-  as(Ileum[["RNA"]], "Assay5")
Ileum[["ADT"]] <-  as(Ileum[["ADT"]], "Assay5")
```

# Examine the UMAP and Clusters
```{r}
DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = c( "WNN_clusters_cc"), label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = c("Phase"), label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = c( "Monocle_clusters_cc"), label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = c( "InfectionStatus"), label.size = 2, label = TRUE, raster = F) 
DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = c( "WNN_clusters_cc"), combine = TRUE, label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = c(  "Monocle_clusters_cc"), combine = TRUE, label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = c( "Phase"), combine = TRUE, label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = c( "InfectionStatus"), combine = TRUE, label.size = 2, label = TRUE, raster = F) 
FeaturePlot(Ileum, reduction = "wnn.umap_cc", "Saa1") #I know from past iterations that Saa1 is a defining gene of the yersinia enriched pathogen enterocytes
``` 

#Labeling the Epithelial compartment. 
```{r}
# Haber, A. L.,  O., … Regev, A. (2017). https://doi.org/10.1038/nature24489  and  #Pardy, R. D., Walzer, K. A.,  Hunter, C. A. (2024). https://doi.org/10.1371/journal.ppat.1011820
#Intestinal Epithelial Stem cells 

DefaultAssay(Ileum) <- "RNA"
Ileum <- JoinLayers(Ileum)

stem_features <- list(c( 'Lgr5', 'Ascl2', 'Slc12a2', 'Axin2','Olfm4', 'Gkn3'))
entero_features <- list(c( 'Chga','Chgb','Tac1','Tph1', 'Neurog3'))
tuft_features <- list(c(  'Dclk1',  'Gnat3',  'Plcb2',  'Trpm5',   'Il25',  'Tslp'))
goblet_features <- list(c( 'Clca1', 'Zg16',  'Agr2',   'Muc2',   'Fcgbp'))
entericGlia <- list(c("S100b", "Plp1", "Sox10", "Gfap", "Fabp7"))
paneth_features <- list(c(  'Defa20',  'Defa22',  'Defa21',   'Gm15308',   'Defa17',  'Clps',   'Defa23',   'Gm15299',  'Lyz1',   'Mptx2'))
Ileum <- AddModuleScore(  object = Ileum,  features = paneth_features,  ctrl = 100,  name = 'Paneth_Features')
Ileum <- AddModuleScore( object = Ileum,  features = goblet_features,  ctrl = 100,  name = 'Goblet_Features')
Ileum <- AddModuleScore( object = Ileum, features = tuft_features, ctrl = 100, name = 'Tuft_Features')
Ileum <- AddModuleScore(object = Ileum, features = stem_features, ctrl = 100, name = 'StemCell_Features')
Ileum <- AddModuleScore(object = Ileum, features = entero_features, ctrl = 100, name = 'Enteroendocrine_Features')
Ileum <- AddModuleScore(object = Ileum, features = entericGlia, ctrl = 100, name = 'EntericGlia_Features')
```

```{r}
theme <-theme(panel.background = element_rect(fill = "white", 
            colour = NA), panel.border = element_rect(fill = NA, 
            colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
            panel.grid.minor = element_line(linewidth = rel(0.5)), 
            strip.background = element_rect(fill = "grey85", 
                colour = "grey20"), axis.title.x = element_blank(),  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
 VlnPlot(Ileum, features = c("Olfm4"), pt.size = 0, group.by = "WNN_clusters") + NoLegend() +  theme
 FeaturePlot_scCustom(Ileum, "StemCell_Features1", reduction = "wnn.umap_cc", na_cutoff = NA)
  VlnPlot(Ileum, features = c("Chgb"), pt.size = 0, group.by = "WNN_clusters") + NoLegend() +  theme
 FeaturePlot_scCustom(Ileum, "Enteroendocrine_Features1", reduction = "wnn.umap_cc", na_cutoff = NA)
  VlnPlot(Ileum, features = c("Muc2"), pt.size = 0, group.by = "WNN_clusters") + NoLegend() +  theme
 FeaturePlot_scCustom(Ileum, "Goblet_Features1", reduction = "wnn.umap_cc", na_cutoff = NA)
 ggsave( "FeaturePlot_Goblet_Features1.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
  VlnPlot(Ileum, features = c("Tslp"), pt.size = 0, group.by = "WNN_clusters") + NoLegend() +  theme
 FeaturePlot_scCustom(Ileum, "Tuft_Features1", reduction = "wnn.umap_cc", na_cutoff = NA)
  VlnPlot(Ileum, features = c("Lyz1"), pt.size = 0, group.by = "WNN_clusters") + NoLegend() +  theme
 FeaturePlot_scCustom(Ileum, "Paneth_Features1", reduction = "wnn.umap_cc", na_cutoff = NA)
 FeaturePlot_scCustom(Ileum, "EntericGlia_Features1", reduction = "wnn.umap_cc", na_cutoff = NA)
```

#Transfer Labels from Prior object
```{r}
#Make the the cell names between objects is compatible. The below code resolved a prior duplication issue 
old_names <- colnames(Ileum_old)
new_cell_names <- sub("^[^_]+_[^_]+_", "", old_names)
colnames(Ileum_old) <- new_cell_names

matching_cells <- match(colnames(Ileum), colnames(Ileum_old))
Ileum$TransferredCellType <- "Unknown"
Ileum$TransferredCellType[!is.na(matching_cells)] <- Ileum_old$FineCellType[matching_cells[!is.na(matching_cells)]]

sum(Ileum$TransferredCellType == "Unknown") #The cells labeled as unknown are the cells from the new samples that are being added

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = c( "TransferredCellType"), combine = TRUE, label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = c(  "Monocle_clusters_cc"), combine = TRUE, label.size = 2, label = TRUE, raster = F) + NoLegend()

old_names <- colnames(MLN_old)
new_cell_names <- sub("^[^_]+_[^_]+_", "", old_names)
colnames(MLN_old) <- new_cell_names

matching_cells <- match(colnames(MLN), colnames(MLN_old))
MLN$TransferredCellType <- "Unknown"
MLN$TransferredCellType[!is.na(matching_cells)] <- MLN_old$FineCellType[matching_cells[!is.na(matching_cells)]]

sum(MLN$TransferredCellType == "Unknown") #The cells labeled as unknown are the cells from the new samples that are being added

DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = c( "TransferredCellType"), combine = TRUE, label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = c(  "Monocle_clusters_cc"), combine = TRUE, label.size = 2, label = TRUE, raster = F) + NoLegend()


```

```{r}
library(BiocParallel)
Idents(MLN_old) <- "CoarseCellType"
unique(MLN_old$`Level 1 Cell Type Annotation`)

scRNA_counts <- GetAssayData(MLN_old, assay = "RNA", layer = "data")
FineCellType_vector <- as.character(MLN_old$`Level 1 Cell Type Annotation`)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

MLN_counts <- GetAssayData(MLN, assay = "RNA", layer = "data")
singleR_results_MLN <- SingleR(test = MLN_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(10))

MLN$Coarse_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_MLN, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_MLN$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_MLN$pruned.labels[is.na(singleR_results_MLN$pruned.labels)] <- "Unknown"
  MLN$Coarse_SingleR<-singleR_results_MLN$labels

#Add the Score for visualization
MLN$Coarse_SingleR_score <- NA
MLN$Coarse_SingleR_score <-singleR_results_MLN$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
MLN$Coarse_SingleR_deltaNext <- NA
MLN$Coarse_SingleR_deltaNext<-singleR_results_MLN$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_MLN , ncol = 3)
plotScoreHeatmap(singleR_results_MLN,  show.pruned = T)


FeaturePlot_scCustom(MLN, reduction = "wnn.umap_cc", features = "Coarse_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(MLN, reduction = "wnn.umap_cc", features = "Coarse_SingleR_deltaNext", na_cutoff = 0) 
DimPlot(MLN, reduction = "wnn.umap_cc", group.by = "Coarse_SingleR" , label = T , repel = T, raster = T) 
table(MLN$Coarse_SingleR) 

library(reticulate)
use_condaenv("/home/hartandrew/.conda/envs/scvi-env", required = TRUE)
Sys.setenv(MPLBACKEND = "Agg")
MLN <- FindClusters(MLN, resolution = 1, cluster.name = "Low_res_clusters", algorithm = 4,  group.singletons = T, graph.name = "harmony_rna_cc", random.seed = 20)
DimPlot(MLN, reduction = "wnn.umap_cc", group.by = "Low_res_clusters" , label = T , repel = T, raster = T) 
DimPlot(MLN, reduction = "wnn.umap_cc", group.by = "Coarse_SingleR" , label = T , repel = T, raster = T) 
MLN <- FindSubCluster( MLN,  "28",  graph.name = "wsnn" , subcluster.name = "Low_res_clusters", resolution = 0.03,  algorithm = 4)
DimPlot_scCustom(MLN, reduction = "wnn.umap_cc", group.by = "Low_res_clusters" , label = T , repel = T, raster = T) 
Idents(MLN) <- "Low_res_clusters"
MLN <- FindSubCluster( MLN,  "17",  graph.name = "wsnn" , subcluster.name = "Low_res_clusters", resolution = 0.03,  algorithm = 4)
DimPlot_scCustom(MLN, reduction = "wnn.umap_cc", group.by = "Low_res_clusters" , label = T , repel = T, raster = T) 
FeaturePlot_scCustom(MLN, reduction = "wnn.umap_cc", features = "Rorc", na_cutoff = 0) 
table(MLN$Low_res_clusters)
```

```{r}
library(BiocParallel)
Idents(Ileum_old) <- "CoarseCellType"
unique(Ileum_old$CoarseCellType)

scRNA_counts <- GetAssayData(Ileum_old, assay = "RNA", layer = "data")
FineCellType_vector <- as.character(Ileum_old$CoarseCellType)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

Ileum_counts <- GetAssayData(Ileum, assay = "RNA", layer = "data")
singleR_results_Ileum <- SingleR(test = Ileum_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(10))

Ileum$Coarse_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_Ileum, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_Ileum$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_Ileum$pruned.labels[is.na(singleR_results_Ileum$pruned.labels)] <- "Unknown"
  Ileum$Coarse_SingleR<-singleR_results_Ileum$labels

#Add the Score for visualization
Ileum$Coarse_SingleR_score <- NA
Ileum$Coarse_SingleR_score <-singleR_results_Ileum$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
Ileum$Coarse_SingleR_deltaNext <- NA
Ileum$Coarse_SingleR_deltaNext<-singleR_results_Ileum$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_Ileum , ncol = 3)
plotScoreHeatmap(singleR_results_Ileum,  show.pruned = T)


FeaturePlot_scCustom(Ileum, reduction = "wnn.umap_cc", features = "Coarse_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(Ileum, reduction = "wnn.umap_cc", features = "Coarse_SingleR_deltaNext", na_cutoff = 0) 
DimPlot(Ileum, reduction = "wnn.umap_cc", group.by = "Coarse_SingleR" , label = T , repel = T, raster = T) 
table(Ileum$Coarse_SingleR) 

library(reticulate)
use_condaenv("/home/hartandrew/.conda/envs/scvi-env", required = TRUE)
Sys.setenv(MPLBACKEND = "Agg")
Ileum <- FindClusters(Ileum, resolution = 0.5, cluster.name = "Low_res_clusters", algorithm = 4,  group.singletons = T, graph.name = "harmony_rna_cc", random.seed = 20)
DimPlot(Ileum, reduction = "wnn.umap_cc", group.by = "Low_res_clusters" , label = T , repel = T, raster = T) 
```
#Assign Coarse Cell types
```{r}
#Refine Identification to most common among cluster
cluster_labels <- as.factor(MLN$Monocle_clusters_cc) 
predicted_labels <- MLN$Coarse_SingleR 

# Function to find the most common predicted label within each cluster
get_best_label <- function(cluster) {
  cluster_cells <- names(cluster_labels[cluster_labels == cluster])
  cluster_pred_labels <- predicted_labels[cluster_cells]
  best_label <- names(sort(table(cluster_pred_labels), decreasing = TRUE))[1]
  return(best_label)
}

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)

# Map the best labels back to the Seurat object
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(MLN$Monocle_clusters_cc) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
MLN$CoarseCellType <-monocle_clusts
DimPlot(MLN, reduction = "wnn.umap_cc", group.by = "CoarseCellType")

#When we double check,The plasma cells need defined, The stromal cells need separated
MLN$CoarseCellType[MLN$Monocle_clusters_cc == "88"] <- "Plasma cells"
MLN$CoarseCellType[MLN$Low_res_clusters == "28_2" | MLN$Low_res_clusters == "28_3"] <- "Stromal cells"
MLN$CoarseCellType[MLN$CoarseCellType %in% c("CD4 T Cells", "CD8 T Cells", "gd T cells", "Innate Lymphoid Cells") ] <- "T cells and ILCs"
DimPlot(MLN, reduction = "wnn.umap_cc", group.by = "CoarseCellType")

Idents(MLN) <- "Low_res_clusters"
library(presto)
Markers_temp <- FindAllMarkers(MLN, assay = "RNA", test.use = "wilcox", only.pos = T)

```


```{r}
#Refine Identification to most common among cluster
cluster_labels <- as.factor(Ileum$Monocle_clusters_cc)  
predicted_labels <- Ileum$Coarse_SingleR  

# Function to find the most common predicted label within each cluster
get_best_label <- function(cluster) {
  cluster_cells <- names(cluster_labels[cluster_labels == cluster])
  cluster_pred_labels <- predicted_labels[cluster_cells]
  best_label <- names(sort(table(cluster_pred_labels), decreasing = TRUE))[1]
  return(best_label)
}

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)

# Map the best labels back to the Seurat object
monocle_clusts <- as.character(Ileum$Monocle_clusters_cc) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
Ileum$CoarseCellType <-monocle_clusts
DimPlot(Ileum, reduction = "wnn.umap_cc", group.by = "CoarseCellType")

#When we double check, It is clear that the MAST cells have been looped in with the T cells and ILCs
Ileum$CoarseCellType[Ileum$Monocle_clusters_cc == "157"] <- "Myeloid cells"
DimPlot(Ileum, reduction = "wnn.umap_cc", group.by = "CoarseCellType")
```

#Export to H5ad and CellxGene for manual check
```{r}

RNA_only <- MLN
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["WNN_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["WNN_clusters"] <- eh 
eh <-unlist(RNA_only[["Monocle_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_clusters"] <- eh 
eh <-unlist(RNA_only[["Monocle_clusters_cc"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_clusters_cc"] <- eh 
eh <-unlist(RNA_only[["WNN_clusters_cc"]])
eh <- as.character(eh) 
RNA_only@meta.data["WNN_clusters_cc"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")

RNA_only[["RNA"]]$scale.data <- NULL 
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)


RNA_only <- Ileum
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["WNN_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["WNN_clusters"] <- eh 
eh <-unlist(RNA_only[["Monocle_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_clusters"] <- eh 
eh <-unlist(RNA_only[["Monocle_clusters_cc"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_clusters_cc"] <- eh 
eh <-unlist(RNA_only[["WNN_clusters_cc"]])
eh <- as.character(eh) 
RNA_only@meta.data["WNN_clusters_cc"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)

```


#Visualization 

```{r}
#Visualizing the UMAPS
umap_data <- Ileum[["wnn.umap_cc"]]@cell.embeddings
cluster_data <- Ileum$CoarseCellType

# Combine UMAP coordinates and cluster information into a single dataframe
umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cluster_data)
ggplot(umap_data, aes(x = wnnUMAPcc_1, y = wnnUMAPcc_2, color = cluster)) +
  geom_point(alpha = 0.2, size = 0.2) +  
  ylim(-15, 20) + xlim(-15, 20) +
  geom_mark_hull(aes(fill = cluster, alpha = 0.01), alpha = 0.1, color = "NA", concavity = 500, 
                  con.type = "none") +  
  scale_color_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  labs(title = "UMAP with Cluster Hulls") +
     theme(panel.background = element_rect(fill = "white", 
            colour = NA), panel.border = element_rect(fill = NA, 
            colour = "grey20"), panel.grid = element_blank(), 
            panel.grid.minor = element_blank(), 
          legend.position = "bottom",
            strip.background = element_blank())
ggsave( "Ileum_UMAP_HULL_Trial.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


umap_data <- MLN[["wnn.umap_cc"]]@cell.embeddings
cluster_data <- MLN$CoarseCellType


umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cluster_data)
ggplot(umap_data, aes(x = wnnUMAPcc_1, y = wnnUMAPcc_2, color = cluster)) +
  geom_point(alpha = 0.2, size = 0.2) + 
  ylim(-12, 12) + xlim(-10, 12) +
  geom_mark_hull(aes(fill = cluster), alpha = 0.1, color = "NA", concavity = 500) +  
  scale_color_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  labs(title = "UMAP with Cluster Hulls") +
     theme(panel.background = element_rect(fill = "white", 
            colour = NA), panel.border = element_rect(fill = NA, 
            colour = "grey20"), panel.grid = element_blank(), 
            panel.grid.minor = element_blank(), 
          legend.position = "bottom",
            strip.background = element_blank())
ggsave( "MLN_UMAP_HULL_Trial.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
# Increase concavity to make hulls tighter
```



#Sub cluster Ileum T cells
```{r}

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum) <- "CoarseCellType"
Ileum_T <- subset(x = Ileum, idents = c("T cells and ILCs"))

#Repeat analysis Process
Ileum_T[['RNA']] <-  JoinLayers(Ileum_T[['RNA']])
Ileum_T[['RNA']] <-  split(x = Ileum_T[['RNA']], f = Ileum_T$Infection_Organ)
#Ileum_T <- NormalizeData(Ileum_T)
Ileum_T <- FindVariableFeatures(Ileum_T, selection.method = "vst", nfeatures = 3000) 
#all.genes <- rownames(Ileum_T[["RNA"]])
#Ileum_T <- ScaleData(Ileum_T, features = all.genes)
Ileum_T <- RunPCA(Ileum_T, verbose = FALSE, reduction.name = "T_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 5000 * 1024^2)
Ileum_T <- IntegrateLayers( object = Ileum_T, method = HarmonyIntegration,  orig.reduction = "T_RNA_pca", new.reduction = "integrated.harmony.T", verbose = FALSE)
DefaultAssay(Ileum_T) <- "ADT"
Ileum_T[['ADT']] <-  JoinLayers(Ileum_T[['ADT']])
Ileum_T[['ADT']] <-  split(x = Ileum_T[['ADT']], f = Ileum_T$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_T), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_T)[!rownames(Ileum_T) %in% Isotype_labels] #remove isotypes
#Ileum_T <- ScaleData(Ileum_T, assay = "ADT")
Ileum_T <- RunPCA(Ileum_T, features =prots, reduction.name = "T_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_T <- IntegrateLayers( object = Ileum_T, method = HarmonyIntegration, features = prots,  orig.reduction = "T_ADT_pca", new.reduction = "integrated.ADT.harmony.T", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_T) <- "RNA"
Ileum_T = FindMultiModalNeighbors(Ileum_T, reduction.list = list("integrated.harmony.T", "integrated.ADT.harmony.T"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_T",  snn.graph.name = "wsnn_T", weighted.nn.name = "weighted_T.nn", verbose = TRUE)

Ileum_T <- RunUMAP(Ileum_T, nn.name = "weighted_T.nn", reduction.name = "wnn.T.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_T <- FindClusters(Ileum_T, graph.name = "wsnn_T",  resolution = 2, verbose = FALSE, cluster.name = "Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_T, reduction = "wnn.T.umap",   group.by = c( "Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_T <-JoinLayers(Ileum_T) 
 cds <- as.cell_data_set(Ileum_T, reductions = "wnn.T.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_T[["Monocle_T_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_T[["Monocle_T_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_T, reduction = "wnn.T.umap",   group.by = c( "Tcell_clusters", "Monocle_T_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_T, c("Cd4","Cd8a", "Cd3e", "Trdc", "TCR.Bchain"), reduction = "wnn.T.umap")
VlnPlot(Ileum_T, "Malat1", group.by = "Tcell_clusters")
FeatureScatter(Ileum_T, "Malat1", "rb.prop")
```


```{r}
#Remove Low quality outlier cells. In clustering, the cells which have < 0.05 ribosomal gene tend to have very higher Malat1 etc. These outliars are technical artifacts that stand apart from other cells. ~2300 in this data set and they are removed
Low_quality_Tcells <- as.data.frame(WhichCells(subset(x = Ileum_T, subset = rb.prop < 0.05)))
write.csv(Low_quality_Tcells, paste0(CSV, "/Ileum_Tcell_LowQuality.csv"))
Ileum_T <- subset(x = Ileum_T, subset = rb.prop >0.05)

DefaultAssay(Ileum_T) <- "RNA"
Ileum_T = FindMultiModalNeighbors(Ileum_T, reduction.list = list("integrated.hamrony.T", "integrated.ADT.hamrony.T"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_T",  snn.graph.name = "wsnn_T", weighted.nn.name = "weighted_T.nn", verbose = TRUE)

Ileum_T <- RunUMAP(Ileum_T, nn.name = "weighted_T.nn", reduction.name = "wnn.T.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_T <- FindClusters(Ileum_T, graph.name = "wsnn_T",  resolution = 2, verbose = FALSE, cluster.name = "Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_T, reduction = "wnn.T.umap",   group.by = c( "Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_T <-JoinLayers(Ileum_T) 
 cds <- as.cell_data_set(Ileum_T, reductions = "wnn.T.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_T[["Monocle_T_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_T[["Monocle_T_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_T, reduction = "wnn.T.umap",   group.by = c( "Tcell_clusters", "Monocle_T_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)

```

```{r}

library(BiocParallel)
Idents(Ileum_old) <- "CoarseCellType"
unique(Ileum_old$CoarseCellType)
Ileum_T_old <- subset(Ileum_old, idents = "T cells and ILCs")
unique(Ileum_T_old$FineCellType)
scRNA_counts <- GetAssayData(Ileum_T_old, assay = "RNA", layer = "data")
#scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(Ileum_T_old$FineCellType)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

Ileum_T_counts <- GetAssayData(Ileum_T, assay = "RNA", layer = "data")
singleR_results_Ileum_T <- SingleR(test = Ileum_T_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

Ileum_T$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_Ileum_T, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_Ileum_T$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_Ileum_T$pruned.labels[is.na(singleR_results_Ileum_T$pruned.labels)] <- "Unknown"
  Ileum_T$Fine_SingleR<-singleR_results_Ileum_T$labels

#Add the Score for visualization
Ileum_T$Fine_SingleR_score <- NA
Ileum_T$Fine_SingleR_score <-singleR_results_Ileum_T$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
Ileum_T$Fine_SingleR_deltaNext <- NA
Ileum_T$Fine_SingleR_deltaNext<-singleR_results_Ileum_T$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_Ileum_T , ncol = 3)
plotScoreHeatmap(singleR_results_Ileum_T,  show.pruned = T)


FeaturePlot_scCustom(Ileum_T, reduction = "wnn.T.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(Ileum_T, reduction = "wnn.T.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(Ileum_T, reduction = "wnn.T.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(Ileum_T$Monocle_T_clusters)  
predicted_labels <- Ileum_T$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(Ileum_T$Monocle_T_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
Ileum_T$FineCellType <-monocle_clusts
DimPlot_scCustom(Ileum_T, reduction = "wnn.T.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_T, reduction = "wnn.T.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_T, reduction = "wnn.T.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
RNA_only <- Ileum_T
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_T_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_T_clusters"] <- eh 
eh <-unlist(RNA_only[["Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_T_filt, "/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediateFiltered.qs2")

```

```{r}
#Cluster 38 contains a portion of cells that are rb.prop ow. It also contains mixed doublets - markers like pla2g7 and tgm2 and slc7a11 which are exclusively found in other cell populations are found in this tiny group - also doublets with ILC2s
#Several cycling clusters will be renamed as such
#Cluster 346 IL17rb + NKT cells? These cells Express higher IL17rb and Il1rl1 (ST2)  and Gata3. They  have reduced IL12rb, reduced CCr5, and other markers. I still believe they are NKT cells because they generally express Zbtb16. There is some evidence of a distinct IL17rb expressing population  PMID: 19015310 PMID: 37258582 (lower Zbtb16 in NKT2)
Idents(Ileum_T) <- "FineCellType"
VlnPlot( subset(Ileum_T, idents = c("NKT Cells", "CD103 low CD4 Trm")), "Zbtb16", group.by = "Monocle_T_clusters")
VlnPlot( subset(Ileum_T, idents = c("NKT Cells", "CD103 low CD4 Trm")), "Gata3", group.by = "Monocle_T_clusters")
Ileum_T$FineCellType[Ileum_T$Monocle_T_clusters %in% c("21", "81", "159", "158", "284")] <- "Cycling ILCS"
Ileum_T$FineCellType[Ileum_T$Monocle_T_clusters %in% c("278", "306", "127", "33")] <- "Cycling T Cells"
Ileum_T$FineCellType[Ileum_T$Monocle_T_clusters %in% c("164")] <- "IFN Stimulated CD8 T Cells"
Ileum_T$FineCellType[Ileum_T$Monocle_T_clusters %in% c("346")] <- "NKT2 Cells"

DimPlot_scCustom(Ileum_T, reduction = "wnn.T.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```
```{r}
#Save final version
RNA_only <- Ileum_T
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_T_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_T_clusters"] <- eh 
eh <-unlist(RNA_only[["Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_T, "/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.qs2")
```


#Sub cluster Ileum B cells
```{r}

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum) <- "CoarseCellType"
Ileum_Bcell <- subset(x = Ileum, idents = c("B cells"))

#Repeat analysis Process
Ileum_Bcell[['RNA']] <-  JoinLayers(Ileum_Bcell[['RNA']])
Ileum_Bcell[['RNA']] <-  split(x = Ileum_Bcell[['RNA']], f = Ileum_Bcell$Infection_Organ)
Ileum_Bcell <- FindVariableFeatures(Ileum_Bcell, selection.method = "vst", nfeatures = 2000) 
Ileum_Bcell <- RunPCA(Ileum_Bcell, verbose = FALSE, reduction.name = "Bcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_Bcell <- IntegrateLayers( object = Ileum_Bcell, method = HarmonyIntegration,  orig.reduction = "Bcell_RNA_pca", new.reduction = "integrated.harmony.Bcell", verbose = FALSE)
DefaultAssay(Ileum_Bcell) <- "ADT"
Ileum_Bcell[['ADT']] <-  JoinLayers(Ileum_Bcell[['ADT']])
Ileum_Bcell[['ADT']] <-  split(x = Ileum_Bcell[['ADT']], f = Ileum_Bcell$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_Bcell), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_Bcell)[!rownames(Ileum_Bcell) %in% Isotype_labels] #remove isotypes
Ileum_Bcell <- RunPCA(Ileum_Bcell, features =prots, reduction.name = "Bcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_Bcell <- IntegrateLayers( object = Ileum_Bcell, method = HarmonyIntegration, features = prots,  orig.reduction = "Bcell_ADT_pca", new.reduction = "integrated.ADT.harmony.Bcell", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_Bcell) <- "RNA"
Ileum_Bcell = FindMultiModalNeighbors(Ileum_Bcell, reduction.list = list("integrated.harmony.Bcell", "integrated.ADT.harmony.Bcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Bcell",  snn.graph.name = "wsnn_Bcell", weighted.nn.name = "weighted_Bcell.nn", verbose = TRUE)

Ileum_Bcell <- RunUMAP(Ileum_Bcell, nn.name = "weighted_Bcell.nn", reduction.name = "wnn.Bcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Bcell <- FindClusters(Ileum_Bcell, graph.name = "wsnn_Bcell",  resolution = 2, verbose = FALSE, cluster.name = "Bcell_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Bcell <-JoinLayers(Ileum_Bcell) 
 cds <- as.cell_data_set(Ileum_Bcell, reductions = "wnn.Bcell.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Bcell[["Monocle_Bcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Bcell[["Monocle_Bcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "Monocle_Bcell_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_Bcell, c( "Cd3e", "Epcam", "rb.prop"), reduction = "wnn.Bcell.umap")
VlnPlot(Ileum_Bcell, "Malat1", group.by = "Bcell_clusters")
FeatureScatter(Ileum_Bcell, "Malat1", "rb.prop")
```


```{r}
#Remove Low quality outlier cells. In clustering, the cells which have < 0.05 ribosomal gene tend to have very higher Malat1 etc. These outliars are technical artifacts that stand apart from other cells. 
Low_quality_Bcells <- as.data.frame(WhichCells(subset(x = Ileum_Bcell, subset = rb.prop < 0.05)))
write.csv(Low_quality_Bcells, paste0(CSV, "/Ileum_Bcell_LowQuality.csv"))
Ileum_Bcell <- subset(x = Ileum_Bcell, subset = rb.prop >0.05)

DefaultAssay(Ileum_Bcell) <- "RNA"
Ileum_Bcell = FindMultiModalNeighbors(Ileum_Bcell, reduction.list = list("integrated.harmony.Bcell", "integrated.ADT.harmony.Bcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Bcell",  snn.graph.name = "wsnn_T", weighted.nn.name = "weighted_Bcell.nn", verbose = TRUE)

Ileum_Bcell <- RunUMAP(Ileum_Bcell, nn.name = "weighted_Bcell.nn", reduction.name = "wnn.Bcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Bcell <- FindClusters(Ileum_Bcell, graph.name = "wsnn_Bcell",  resolution = 2, verbose = FALSE, cluster.name = "Bcell_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Bcell <-JoinLayers(Ileum_Bcell) 
 cds <- as.cell_data_set(Ileum_Bcell, reductions = "wnn.Bcell.umap", default.reduction = "wnn.Bcell.umap",graph =  "wsnn_Bcell",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Bcell[["Monocle_Bcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Bcell[["Monocle_Bcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "Monocle_Bcell_clusters", "IntermediateCellType", "Phase", "InfectionStatus"), combine = FALSE, label.size = 2, label = TRUE, raster = F)
DimPlot(Ileum_Bcell, reduction = "wnn.Bcell.umap",   group.by = c("Monocle_Bcell_clusters"), combine = T, label.size = 2, label = TRUE, raster = F) + NoLegend()

```

```{r}

library(BiocParallel)
Idents(Ileum_old) <- "CoarseCellType"
unique(Ileum_old$CoarseCellType)
Ileum_Bcell_old <- subset(Ileum_old, idents = "B cells")
unique(Ileum_Bcell_old$FineCellType)
scRNA_counts <- GetAssayData(Ileum_Bcell_old, assay = "RNA", layer = "data")
#scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(Ileum_Bcell_old$FineCellType)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

Ileum_Bcell_counts <- GetAssayData(Ileum_Bcell, assay = "RNA", layer = "data")
singleR_results_Ileum_Bcell <- SingleR(test = Ileum_Bcell_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

Ileum_Bcell$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_Ileum_Bcell, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_Ileum_Bcell$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_Ileum_Bcell$pruned.labels[is.na(singleR_results_Ileum_Bcell$pruned.labels)] <- "Unknown"
  Ileum_Bcell$Fine_SingleR<-singleR_results_Ileum_Bcell$labels

#Add the Score for visualization
Ileum_Bcell$Fine_SingleR_score <- NA
Ileum_Bcell$Fine_SingleR_score <-singleR_results_Ileum_Bcell$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
Ileum_Bcell$Fine_SingleR_deltaNext <- NA
Ileum_Bcell$Fine_SingleR_deltaNext<-singleR_results_Ileum_Bcell$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_Ileum_Bcell , ncol = 3)
plotScoreHeatmap(singleR_results_Ileum_Bcell,  show.pruned = T)


FeaturePlot_scCustom(Ileum_Bcell, reduction = "wnn.Bcell.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(Ileum_Bcell, reduction = "wnn.Bcell.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(Ileum_Bcell, reduction = "wnn.Bcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(Ileum_Bcell$Monocle_Bcell_clusters)  
predicted_labels <- Ileum_Bcell$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(Ileum_Bcell$Monocle_Bcell_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
Ileum_Bcell$FineCellType <-monocle_clusts
DimPlot_scCustom(Ileum_Bcell, reduction = "wnn.Bcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Bcell, reduction = "wnn.Bcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Bcell, reduction = "wnn.Bcell.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
RNA_only <- Ileum_Bcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Bcell_clusters"] <- eh 
eh <-unlist(RNA_only[["Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Bcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Bcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Bcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_Bcell, "/path/to/analysis/directory/Seurat_Files/Ileum_Bcell_intermediate.qs2")

```

```{r}
#Remove technical duplicates and Erythroid cells 
# 36, 54, 67, 80, 182, 233, 85, 241, 180, 12 are Klhl14 High - This seems to be the only gene that separates these cells though they come from all populations and infection types
# 52, 76, 33, 69, 99 Hbb high Erythroid
# 106 T cells genes - Mixed
# 255 Myeloid mixed
# 226 244 Mixed populations
Idents(Ileum_Bcell) <- "Monocle_Bcell_clusters"
Low_quality_Bcells2 <- as.data.frame(WhichCells(subset(x = Ileum_Bcell, idents = c("52", "76", "33", "69", "99", "106", "255", "226", "244"))))
colnames(Low_quality_Bcells) <- "CallNames"
colnames(Low_quality_Bcells2) <- "CallNames"
Low_quality_Bcells2 <- rbind(Low_quality_Bcells, Low_quality_Bcells2)
write.csv(Low_quality_Bcells2, paste0(CSV, "/Ileum_Bcell_LowQuality.csv"))
Ileum_Bcell <- subset(x = Ileum_Bcell, idents = c("52", "76", "33", "69", "99", "106", "255", "226", "244") , invert = T)
#After removal, I went 3 blocks up and recluster, drew new umap, etc. 

```


```{r}

DefaultAssay(Ileum_Bcell) <- "RNA"
Ileum_Bcell = FindMultiModalNeighbors(Ileum_Bcell, reduction.list = list("integrated.harmony.Bcell", "integrated.ADT.harmony.Bcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Bcell",  snn.graph.name = "wsnn_T", weighted.nn.name = "weighted_Bcell.nn", verbose = TRUE)

Ileum_Bcell <- RunUMAP(Ileum_Bcell, nn.name = "weighted_Bcell.nn", reduction.name = "wnn.Bcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Bcell <- FindClusters(Ileum_Bcell, graph.name = "wsnn_Bcell",  resolution = 2, verbose = FALSE, cluster.name = "Bcell_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Bcell <-JoinLayers(Ileum_Bcell) 
 cds <- as.cell_data_set(Ileum_Bcell, reductions = "wnn.Bcell.umap", default.reduction = "wnn.Bcell.umap",graph =  "wsnn_Bcell",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Bcell[["Monocle_Bcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Bcell[["Monocle_Bcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "Monocle_Bcell_clusters", "IntermediateCellType", "Phase", "InfectionStatus"), combine = FALSE, label.size = 2, label = TRUE, raster = F)

```



```{r}

Ileum_Bcell$FineCellType[Ileum_Bcell$Monocle_Bcell_clusters %in% c("99", "14", "198", "36", "249", "9", "62")] <- "Klhl14 B cells"

DimPlot_scCustom(Ileum_Bcell, reduction = "wnn.Bcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}

RNA_only <- Ileum_Bcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Bcell_clusters"] <- eh 
eh <-unlist(RNA_only[["Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Bcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Bcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Bcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_Bcell, "/path/to/analysis/directory/Seurat_Files/Ileum_Bcell_intermediate.qs2")
```



#Ileum Stromal Cell Subclustering 
```{r}

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum) <- "CoarseCellType"
Ileum_Stromal <- subset(x = Ileum, idents = c("Stromal cells"))

#Repeat analysis Process
Ileum_Stromal[['RNA']] <-  JoinLayers(Ileum_Stromal[['RNA']])
Ileum_Stromal[['RNA']] <-  split(x = Ileum_Stromal[['RNA']], f = Ileum_Stromal$Infection_Organ)
Ileum_Stromal <- FindVariableFeatures(Ileum_Stromal, selection.method = "vst", nfeatures = 3000) 
Ileum_Stromal <- RunPCA(Ileum_Stromal, verbose = FALSE, reduction.name = "Stromal_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_Stromal <- IntegrateLayers( object = Ileum_Stromal, method = HarmonyIntegration,  orig.reduction = "Stromal_RNA_pca", new.reduction = "integrated.harmony.Stromal", verbose = FALSE)
DefaultAssay(Ileum_Stromal) <- "ADT"
Ileum_Stromal[['ADT']] <-  JoinLayers(Ileum_Stromal[['ADT']])
Ileum_Stromal[['ADT']] <-  split(x = Ileum_Stromal[['ADT']], f = Ileum_Stromal$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_Stromal), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_Stromal)[!rownames(Ileum_Stromal) %in% Isotype_labels] #remove isotypes
Ileum_Stromal <- RunPCA(Ileum_Stromal, features =prots, reduction.name = "Stromal_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_Stromal <- IntegrateLayers( object = Ileum_Stromal, method = HarmonyIntegration, features = prots,  orig.reduction = "Stromal_ADT_pca", new.reduction = "integrated.ADT.harmony.Stromal", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_Stromal) <- "RNA"
Ileum_Stromal = FindMultiModalNeighbors(Ileum_Stromal, reduction.list = list("integrated.harmony.Stromal", "integrated.ADT.harmony.Stromal"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Stromal",  snn.graph.name = "wsnn_Stromal", weighted.nn.name = "weighted_Stromal.nn", verbose = TRUE)

Ileum_Stromal <- RunUMAP(Ileum_Stromal, nn.name = "weighted_Stromal.nn", reduction.name = "wnn.Stromal.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Stromal <- FindClusters(Ileum_Stromal, graph.name = "wsnn_Stromal",  resolution = 2, verbose = FALSE, cluster.name = "Stromal_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Stromal, reduction = "wnn.Stromal.umap",   group.by = c( "Stromal_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Stromal <-JoinLayers(Ileum_Stromal) 
 cds <- as.cell_data_set(Ileum_Stromal, reductions = "wnn.Stromal.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Stromal[["Monocle_Stromal_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Stromal[["Monocle_Stromal_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Stromal, reduction = "wnn.Stromal.umap",   group.by = c( "Stromal_clusters", "Monocle_Stromal_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_Stromal, c( "Cd3e", "Epcam", "rb.prop"), reduction = "wnn.Stromal.umap")
VlnPlot(Ileum_Stromal, "Malat1", group.by = "Stromal_clusters")
FeatureScatter(Ileum_Stromal, "Malat1", "rb.prop")
```


```{r}

library(BiocParallel)
Idents(Ileum_old) <- "CoarseCellType"
unique(Ileum_old$CoarseCellType)
Ileum_Stromal_old <- subset(Ileum_old, idents = "Stromal cells")
unique(Ileum_Stromal_old$FineCellType)
scRNA_counts <- GetAssayData(Ileum_Stromal_old, assay = "RNA", layer = "data")
#scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(Ileum_Stromal_old$FineCellType)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

Ileum_Stromal_counts <- GetAssayData(Ileum_Stromal, assay = "RNA", layer = "data")
singleR_results_Ileum_Stromal <- SingleR(test = Ileum_Stromal_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

Ileum_Stromal$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_Ileum_Stromal, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_Ileum_Stromal$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_Ileum_Stromal$pruned.labels[is.na(singleR_results_Ileum_Stromal$pruned.labels)] <- "Unknown"
  Ileum_Stromal$Fine_SingleR<-singleR_results_Ileum_Stromal$labels

#Add the Score for visualization
Ileum_Stromal$Fine_SingleR_score <- NA
Ileum_Stromal$Fine_SingleR_score <-singleR_results_Ileum_Stromal$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
Ileum_Stromal$Fine_SingleR_deltaNext <- NA
Ileum_Stromal$Fine_SingleR_deltaNext<-singleR_results_Ileum_Stromal$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_Ileum_Stromal , ncol = 3)
plotScoreHeatmap(singleR_results_Ileum_Stromal,  show.pruned = T)


FeaturePlot_scCustom(Ileum_Stromal, reduction = "wnn.Stromal.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(Ileum_Stromal, reduction = "wnn.Stromal.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(Ileum_Stromal, reduction = "wnn.Stromal.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(Ileum_Stromal$Monocle_Stromal_clusters)  
predicted_labels <- Ileum_Stromal$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(Ileum_Stromal$Monocle_Stromal_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
Ileum_Stromal$FineCellType <-monocle_clusts
DimPlot_scCustom(Ileum_Stromal, reduction = "wnn.Stromal.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Stromal, reduction = "wnn.Stromal.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Stromal, reduction = "wnn.Stromal.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
RNA_only <- Ileum_Stromal
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Stromal_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Stromal_clusters"] <- eh 
eh <-unlist(RNA_only[["Stromal_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Stromal_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Stromal_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Stromal_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_Stromal, "/path/to/analysis/directory/Seurat_Files/Ileum_Stromal_intermediate.qs2")

```


# Ileum Myeloid Cells subclustering ileum
```{r}

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum) <- "CoarseCellType"
Ileum_Myeloid <- subset(x = Ileum, idents = c("Myeloid cells"))

#Repeat analysis Process
Ileum_Myeloid[['RNA']] <-  JoinLayers(Ileum_Myeloid[['RNA']])
Ileum_Myeloid[['RNA']] <-  split(x = Ileum_Myeloid[['RNA']], f = Ileum_Myeloid$Infection_Organ)
Ileum_Myeloid <- FindVariableFeatures(Ileum_Myeloid, selection.method = "vst", nfeatures = 3000) 
Ileum_Myeloid <- RunPCA(Ileum_Myeloid, verbose = FALSE, reduction.name = "Myeloid_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_Myeloid <- IntegrateLayers( object = Ileum_Myeloid, method = HarmonyIntegration,  orig.reduction = "Myeloid_RNA_pca", new.reduction = "integrated.harmony.Myeloid", verbose = FALSE)
DefaultAssay(Ileum_Myeloid) <- "ADT"
Ileum_Myeloid[['ADT']] <-  JoinLayers(Ileum_Myeloid[['ADT']])
Ileum_Myeloid[['ADT']] <-  split(x = Ileum_Myeloid[['ADT']], f = Ileum_Myeloid$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_Myeloid), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_Myeloid)[!rownames(Ileum_Myeloid) %in% Isotype_labels] #remove isotypes
Ileum_Myeloid <- RunPCA(Ileum_Myeloid, features =prots, reduction.name = "Myeloid_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_Myeloid <- IntegrateLayers( object = Ileum_Myeloid, method = HarmonyIntegration, features = prots,  orig.reduction = "Myeloid_ADT_pca", new.reduction = "integrated.ADT.harmony.Myeloid", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_Myeloid) <- "RNA"
Ileum_Myeloid = FindMultiModalNeighbors(Ileum_Myeloid, reduction.list = list("integrated.harmony.Myeloid", "integrated.ADT.harmony.Myeloid"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Myeloid",  snn.graph.name = "wsnn_Myeloid", weighted.nn.name = "weighted_Myeloid.nn", verbose = TRUE)

Ileum_Myeloid <- RunUMAP(Ileum_Myeloid, nn.name = "weighted_Myeloid.nn", reduction.name = "wnn.Myeloid.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Myeloid <- FindClusters(Ileum_Myeloid, graph.name = "wsnn_Myeloid",  resolution = 2, verbose = FALSE, cluster.name = "Myeloid_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Myeloid, reduction = "wnn.Myeloid.umap",   group.by = c( "Myeloid_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Myeloid <-JoinLayers(Ileum_Myeloid) 
 cds <- as.cell_data_set(Ileum_Myeloid, reductions = "wnn.Myeloid.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Myeloid[["Monocle_Myeloid_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Myeloid[["Monocle_Myeloid_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Myeloid, reduction = "wnn.Myeloid.umap",   group.by = c( "Myeloid_clusters", "Monocle_Myeloid_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_Myeloid, c( "Cd3e", "Epcam", "rb.prop"), reduction = "wnn.Myeloid.umap")
VlnPlot(Ileum_Myeloid, "Malat1", group.by = "Myeloid_clusters")
FeatureScatter(Ileum_Myeloid, "Malat1", "rb.prop")
```


```{r}

VlnPlot(Ileum_Myeloid, "Malat1", group.by = "Monocle_Myeloid_clusters") + NoLegend()

#Remove the 2 cluster which are extreme outliers - a lack of malat1, higher ribosomal counts - and separated within umap space.
Idents(Ileum_Myeloid) <- "Monocle_Myeloid_clusters"
Low_quality_Myeloid <- as.data.frame(WhichCells(subset(x = Ileum_Myeloid, idents = c("16", "17"))))
write.csv(Low_quality_Myeloid, paste0(CSV, "/Ileum_Myeloid_LowQuality.csv"))
Ileum_Myeloid <- subset(x = Ileum_Myeloid, idents = c("16", "17"), invert = T)


Ileum_Myeloid[['RNA']] <-  JoinLayers(Ileum_Myeloid[['RNA']])
Ileum_Myeloid[['RNA']] <-  split(x = Ileum_Myeloid[['RNA']], f = Ileum_Myeloid$Infection_Organ)
Ileum_Myeloid <- FindVariableFeatures(Ileum_Myeloid, selection.method = "vst", nfeatures = 3000) 
Ileum_Myeloid <- RunPCA(Ileum_Myeloid, verbose = FALSE, reduction.name = "Myeloid_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_Myeloid <- IntegrateLayers( object = Ileum_Myeloid, method = HarmonyIntegration,  orig.reduction = "Myeloid_RNA_pca", new.reduction = "integrated.harmony.Myeloid", verbose = FALSE)
DefaultAssay(Ileum_Myeloid) <- "ADT"
Ileum_Myeloid[['ADT']] <-  JoinLayers(Ileum_Myeloid[['ADT']])
Ileum_Myeloid[['ADT']] <-  split(x = Ileum_Myeloid[['ADT']], f = Ileum_Myeloid$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_Myeloid), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_Myeloid)[!rownames(Ileum_Myeloid) %in% Isotype_labels] #remove isotypes
Ileum_Myeloid <- RunPCA(Ileum_Myeloid, features =prots, reduction.name = "Myeloid_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_Myeloid <- IntegrateLayers( object = Ileum_Myeloid, method = HarmonyIntegration, features = prots,  orig.reduction = "Myeloid_ADT_pca", new.reduction = "integrated.ADT.harmony.Myeloid", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_Myeloid) <- "RNA"
Ileum_Myeloid = FindMultiModalNeighbors(Ileum_Myeloid, reduction.list = list("integrated.harmony.Myeloid", "integrated.ADT.harmony.Myeloid"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Myeloid",  snn.graph.name = "wsnn_Myeloid", weighted.nn.name = "weighted_Myeloid.nn", verbose = TRUE)

Ileum_Myeloid <- RunUMAP(Ileum_Myeloid, nn.name = "weighted_Myeloid.nn", reduction.name = "wnn.Myeloid.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Myeloid <- FindClusters(Ileum_Myeloid, graph.name = "wsnn_Myeloid",  resolution = 2, verbose = FALSE, cluster.name = "Myeloid_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Myeloid, reduction = "wnn.Myeloid.umap",   group.by = c( "Myeloid_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Myeloid <-JoinLayers(Ileum_Myeloid) 
 cds <- as.cell_data_set(Ileum_Myeloid, reductions = "wnn.Myeloid.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.001, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Myeloid[["Monocle_Myeloid_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Myeloid[["Monocle_Myeloid_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Myeloid, reduction = "wnn.Myeloid.umap",   group.by = c( "Myeloid_clusters", "Monocle_Myeloid_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_Myeloid, c( "Cd3e", "Epcam", "rb.prop"), reduction = "wnn.Myeloid.umap")
VlnPlot(Ileum_Myeloid, "Malat1", group.by = "Myeloid_clusters")
FeatureScatter(Ileum_Myeloid, "Malat1", "rb.prop")


```

```{r}

library(BiocParallel)
Idents(Ileum_old) <- "CoarseCellType"
unique(Ileum_old$CoarseCellType)
Ileum_Myeloid_old <- subset(Ileum_old, idents = "Myeloid cells")
unique(Ileum_Myeloid_old$FineCellType)
scRNA_counts <- GetAssayData(Ileum_Myeloid_old, assay = "RNA", layer = "data")
#scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(Ileum_Myeloid_old$FineCellType)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

Ileum_Myeloid_counts <- GetAssayData(Ileum_Myeloid, assay = "RNA", layer = "data")
singleR_results_Ileum_Myeloid <- SingleR(test = Ileum_Myeloid_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

Ileum_Myeloid$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_Ileum_Myeloid, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_Ileum_Myeloid$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_Ileum_Myeloid$pruned.labels[is.na(singleR_results_Ileum_Myeloid$pruned.labels)] <- "Unknown"
  Ileum_Myeloid$Fine_SingleR<-singleR_results_Ileum_Myeloid$labels

#Add the Score for visualization
Ileum_Myeloid$Fine_SingleR_score <- NA
Ileum_Myeloid$Fine_SingleR_score <-singleR_results_Ileum_Myeloid$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
Ileum_Myeloid$Fine_SingleR_deltaNext <- NA
Ileum_Myeloid$Fine_SingleR_deltaNext<-singleR_results_Ileum_Myeloid$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_Ileum_Myeloid , ncol = 3)
plotScoreHeatmap(singleR_results_Ileum_Myeloid,  show.pruned = T)


FeaturePlot_scCustom(Ileum_Myeloid, reduction = "wnn.Myeloid.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(Ileum_Myeloid, reduction = "wnn.Myeloid.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(Ileum_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(Ileum_Myeloid$Monocle_Myeloid_clusters)  
predicted_labels <- Ileum_Myeloid$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(Ileum_Myeloid$Monocle_Myeloid_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
Ileum_Myeloid$FineCellType <-monocle_clusts
DimPlot_scCustom(Ileum_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
#cluster 27 expresses. higher Vim, CD200r1, ccr2, compared to surrounding Dcs and pDCs - This makes it tDC population  PMID: 37414907
#Macrophages are clearly split between Nippo and Yersinia - Nippo macrophages express Mrc1 and Arg1  = M2. In contrast Yersinia derived macrophages express higher B2m, Stat1, Il1b etc. 
#clusters
#https://www.nature.com/articles/s41590-024-01745-9
#https://pmc.ncbi.nlm.nih.gov/articles/PMC10683792/ 

#21 and 19 Retnla high Cd226 + macrophages - high antigen presentation compared to other "monocytes" - similar to a macrophage. unique expression of fcrls which is supposed to exclusively mark microglia  PMID: 38926085. It should be noted that their unique expression of retnla has been associated and called macrophages in other work PMID: 35361907 (type 2 like macrophages). They lack Ly6c2 expression compared to monocytes
#cluster 16 is mgl2 (CD301b) high also CD209a expressing which is has been called many things but is associated with migratory capacity and CD4 activation PMID: 40011469 PMID: 38492222 PMID: 35933399  PMID: 24076051
#Cluster 8 and Cluster 6 are Itgae (CD103) positive which in the intestine is associated with Notch dependent cDC2s PMID: 38351322 Importantly however a portiong of these cells are Lyz2 high which has been associated with Notch independent cDC2 PMID: 38351322
#cluster 25 is a ccr7 expressing cDC subset that also expresses fscn1, ly75, and other unique markers  - trafficking to LN?
#cluster 12 and CLuster 15 are il12a expressing, lower on ccr2, starting to express Adgre1 (F4/80). Inflammatory monocytes?
#Cluster 14, Is Arg1 high macrophages
#Clusters 9, 24, 17, 22, 10 7 are Cx3cr1 high macrophages which several papers refer to as mature/ resident PMID: 32222060. This coincides with loss of ccr2 and high H2-a and C1q expression

Ileum_Myeloid$FineCellType[Ileum_Myeloid$Monocle_Myeloid_clusters %in% c("9", "24", "17", "22", "10", '7' )] <- "Resident Macrophages"
Ileum_Myeloid$FineCellType[Ileum_Myeloid$Monocle_Myeloid_clusters %in% c("14" )] <- "Arg1 High Macrophages"
Ileum_Myeloid$FineCellType[Ileum_Myeloid$Monocle_Myeloid_clusters %in% c("25" )] <- "Ccr7 High cDC2s"
Ileum_Myeloid$FineCellType[Ileum_Myeloid$Monocle_Myeloid_clusters %in% c("8", "6" )] <- "Itgae High cDC2s"
Ileum_Myeloid$FineCellType[Ileum_Myeloid$Monocle_Myeloid_clusters %in% c("16" )] <- "Cd301b High cDCs"
Ileum_Myeloid$FineCellType[Ileum_Myeloid$Monocle_Myeloid_clusters %in% c("21", "19" )] <- "Retnla High Macrophages"
Ileum_Myeloid$FineCellType[Ileum_Myeloid$Monocle_Myeloid_clusters %in% c("27" )] <- "tDCs"

DimPlot_scCustom(Ileum_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()


```




```{r}
RNA_only <- Ileum_Myeloid
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Myeloid_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Myeloid_clusters"] <- eh 
eh <-unlist(RNA_only[["Myeloid_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Myeloid_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Myeloid_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Myeloid_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_Myeloid, "/path/to/analysis/directory/Seurat_Files/Ileum_Myeloid_intermediate.qs2")

```



#Epithelial Ileum Subcluster 

```{r}

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum) <- "CoarseCellType"
Ileum_Epithelial <- subset(x = Ileum, idents = c("Epithelial cells"))

#Repeat analysis Process
Ileum_Epithelial[['RNA']] <-  JoinLayers(Ileum_Epithelial[['RNA']])
Ileum_Epithelial[['RNA']] <-  split(x = Ileum_Epithelial[['RNA']], f = Ileum_Epithelial$Infection_Organ)
Ileum_Epithelial <- FindVariableFeatures(Ileum_Epithelial, selection.method = "vst", nfeatures = 3000) 
Ileum_Epithelial <- RunPCA(Ileum_Epithelial, verbose = FALSE, reduction.name = "Epithelial_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_Epithelial <- IntegrateLayers( object = Ileum_Epithelial, method = HarmonyIntegration,  orig.reduction = "Epithelial_RNA_pca", new.reduction = "integrated.harmony.Epithelial", verbose = FALSE)
DefaultAssay(Ileum_Epithelial) <- "ADT"
Ileum_Epithelial[['ADT']] <-  JoinLayers(Ileum_Epithelial[['ADT']])
Ileum_Epithelial[['ADT']] <-  split(x = Ileum_Epithelial[['ADT']], f = Ileum_Epithelial$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_Epithelial), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_Epithelial)[!rownames(Ileum_Epithelial) %in% Isotype_labels] #remove isotypes
Ileum_Epithelial <- RunPCA(Ileum_Epithelial, features =prots, reduction.name = "Epithelial_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_Epithelial <- IntegrateLayers( object = Ileum_Epithelial, method = HarmonyIntegration, features = prots,  orig.reduction = "Epithelial_ADT_pca", new.reduction = "integrated.ADT.harmony.Epithelial", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_Epithelial) <- "RNA"
Ileum_Epithelial = FindMultiModalNeighbors(Ileum_Epithelial, reduction.list = list("integrated.harmony.Epithelial", "integrated.ADT.harmony.Epithelial"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Epithelial",  snn.graph.name = "wsnn_Epithelial", weighted.nn.name = "weighted_Epithelial.nn", verbose = TRUE)

Ileum_Epithelial <- RunUMAP(Ileum_Epithelial, nn.name = "weighted_Epithelial.nn", reduction.name = "wnn.Epithelial.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Epithelial <- FindClusters(Ileum_Epithelial, graph.name = "wsnn_Epithelial",  resolution = 2, verbose = FALSE, cluster.name = "Epithelial_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Epithelial, reduction = "wnn.Epithelial.umap",   group.by = c( "Epithelial_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Epithelial <-JoinLayers(Ileum_Epithelial) 
 cds <- as.cell_data_set(Ileum_Epithelial, reductions = "wnn.Epithelial.umap", default.reduction = "wnn.Epithelial.umap",graph =  "wsnn_Epithelial",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Epithelial[["Monocle_Epithelial_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Epithelial[["Monocle_Epithelial_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Epithelial, reduction = "wnn.Epithelial.umap",   group.by = c( "Epithelial_clusters", "Monocle_Epithelial_clusters",  "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_Epithelial, c( "Cd3e", "Epcam", "rb.prop", "Malat1"), reduction = "wnn.Epithelial.umap")
VlnPlot(Ileum_Epithelial, "Malat1", group.by = "Epithelial_clusters") + NoLegend()
FeatureScatter(Ileum_Epithelial, "Malat1", "rb.prop")
```




```{r}

library(BiocParallel)
Idents(Ileum_old) <- "CoarseCellType"
unique(Ileum_old$CoarseCellType)
Ileum_Epithelial_old <- subset(Ileum_old, idents = "Epithelial cells")
unique(Ileum_Epithelial_old$FineCellType)
scRNA_counts <- GetAssayData(Ileum_Epithelial_old, assay = "RNA", layer = "data")
FineCellType_vector <- as.character(Ileum_Epithelial_old$FineCellType)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

Ileum_Epithelial_counts <- GetAssayData(Ileum_Epithelial, assay = "RNA", layer = "data")
singleR_results_Ileum_Epithelial <- SingleR(test = Ileum_Epithelial_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

Ileum_Epithelial$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_Ileum_Epithelial, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_Ileum_Epithelial$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_Ileum_Epithelial$pruned.labels[is.na(singleR_results_Ileum_Epithelial$pruned.labels)] <- "Unknown"
  Ileum_Epithelial$Fine_SingleR<-singleR_results_Ileum_Epithelial$labels

#Add the Score for visualization
Ileum_Epithelial$Fine_SingleR_score <- NA
Ileum_Epithelial$Fine_SingleR_score <-singleR_results_Ileum_Epithelial$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
Ileum_Epithelial$Fine_SingleR_deltaNext <- NA
Ileum_Epithelial$Fine_SingleR_deltaNext<-singleR_results_Ileum_Epithelial$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_Ileum_Epithelial , ncol = 3)
plotScoreHeatmap(singleR_results_Ileum_Epithelial,  show.pruned = T)


FeaturePlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(Ileum_Epithelial$Monocle_Epithelial_clusters)  
predicted_labels <- Ileum_Epithelial$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(Ileum_Epithelial$Monocle_Epithelial_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
Ileum_Epithelial$FineCellType <-monocle_clusts
DimPlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}

Ileum_Epithelial$FineCellType[Ileum_Epithelial$FineCellType %in% c("Early Enterocyte", "Middle Enterocyte", "Late Enterocyte" )] <- "Enterocytes"
DimPlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()

#Monocle Clustering
Ileum_Epithelial <-JoinLayers(Ileum_Epithelial) 
 cds <- as.cell_data_set(Ileum_Epithelial, reductions = "wnn.Epithelial.umap", default.reduction = "wnn.Epithelial.umap",graph =  "wsnn_Epithelial",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

cds <- learn_graph(cds, use_partition = F, verbose = F)

#Add the epithelial pseudotime to the seurat object and visualize 
plot_cells(cds,  color_cells_by = "cluster",  label_groups_by_cluster=FALSE,   label_leaves=FALSE,  label_branch_points=FALSE)

cds <- order_cells(cds)
Ileum_Epithelial <- AddMetaData(object = Ileum_Epithelial, metadata = cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "Epithelial_pseudotime")
FeaturePlot(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", "Epithelial_pseudotime" )

unique(Ileum_Epithelial$FineCellType)


Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Epithelial[["Monocle_Epithelial_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Epithelial[["Monocle_Epithelial_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Epithelial, reduction = "wnn.Epithelial.umap",   group.by = c("Monocle_Epithelial_clusters"),combine = FALSE, label.size = 2, label = TRUE, raster = F)

```

```{r}
#Assess appropriate place to mark early middle and late Enterocytes

GeneList <- c("Slc5a1", "Lgr5", "Ada", "Reg3g", "Reg3b")
pseudotime_scores <- Ileum_Epithelial@meta.data$Epithelial_pseudotime
plot_list <- list()
for (gene in unique(GeneList)) {
  
  expression <- FetchData(Ileum_Epithelial, vars = gene)
  data <- data.frame(pseudotime = pseudotime_scores, expression = expression[,1])
  #Bin the pseduotime
  bin_width <- 0.4  
  data$bin <- cut(data$pseudotime, breaks = seq(min(data$pseudotime), max(data$pseudotime), by = bin_width), include.lowest = TRUE)
  
  #Calculate the average  expression for each  bin
  
  bin_averages <- data %>%
    group_by(bin) %>%
    summarise(avg = mean(expression, na.rm = TRUE))
  
  bin_averages$bin <- row_number(bin_averages)
  #replace bin with the center of the bin for plotting
  start_number <- bin_width/2
  step_size <- bin_width
  length_of_vector <- length(bin_averages$bin)
  my_vector <- seq(from = start_number, by = step_size, length.out = length_of_vector)
  bin_averages$bin <- my_vector
  
  
  plot <-ggplot(bin_averages, aes(x = bin, y = avg)) +
    geom_line(color = "#1A0889", size = 0.4) + geom_hline(yintercept = 0, color = "black", size = 1)+
    labs(title = paste0("Average ", gene,  " Expression Across Pseudotime Bins"), 
         x = "Pseudotime (Bins)", y = paste0(gene, " Exp")) +
    theme_minimal() +
    theme(plot.title = element_blank(), axis.title.x = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),  # Remove the panel background
          plot.background = element_blank() )  
  plot_list[[gene]]<- plot
  
  plot_list[[gene]]

}


plot_list[["Lgr5"]]
plot_list[["Reg3b"]]
plot_list[["Slc5a1"]]
plot_list[["Ada"]]
```



```{r}
#Adjust the Early Middle and Late Fine labels and Create Fier scale label
DimPlot(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "FineCellType")

Ileum_Epithelial@meta.data$FineCellType[Ileum_Epithelial$FineCellType %in% c("Early Enteroctye", "Middle Enterocyte", "Late Enterocyte") & Ileum_Epithelial$Epithelial_pseudotime <13.000] <- "Early Enterocyte"
Ileum_Epithelial@meta.data$FineCellType[Ileum_Epithelial$FineCellType %in% c("Early Enteroctye", "Middle Enterocyte", "Late Enterocyte") & Ileum_Epithelial$Epithelial_pseudotime >=13.000 & Ileum_Epithelial$Epithelial_pseudotime <29] <- "Middle Enterocyte"
Ileum_Epithelial@meta.data$FineCellType[Ileum_Epithelial$FineCellType %in% c("Early Enteroctye", "Middle Enterocyte", "Late Enterocyte") & Ileum_Epithelial$Epithelial_pseudotime >29] <- "Late Enterocyte"

DimPlot(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "FineCellType")


Ileum_Epithelial$FineCellType[Ileum_Epithelial$FineCellType == "Differentiating Goblet Cells"] <- "Goblet Cells"

#Additional Cell Types 
Ileum_Epithelial$FinestCellType <- Ileum_Epithelial$FineCellType
Ileum_Epithelial$FinestCellType[Ileum_Epithelial$Epithelial_clusters %in% c("29", "54")] <- "Pla2g4c inflammatory Enterocytes"
Ileum_Epithelial$FinestCellType[Ileum_Epithelial$Monocle_Epithelial_clusters %in% c("8", "43", "123")] <- "Stat1 inflammatory Enterocytes"
Ileum_Epithelial$FinestCellType[Ileum_Epithelial$Monocle_Epithelial_clusters %in% c("45", "85", "27", "116", "48", "60", "126", "111", "21", "93")] <- "Yps inflammatory Enterocytes"
DimPlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "FinestCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(Ileum_Epithelial, reduction = "wnn.Epithelial.umap", group.by = "Monocle_Epithelial_clusters",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```




```{r}
RNA_only <- Ileum_Epithelial
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Epithelial_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Epithelial_clusters"] <- eh 
eh <-unlist(RNA_only[["Epithelial_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Epithelial_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Epithelial_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Epithelial_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_Epithelial, "/path/to/analysis/directory/Seurat_Files/Ileum_Epithelial_intermediate.qs2")

```















```{r}
#Transfer Ileum Labels 
Ileum$FineCellType <- "NA"
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Epithelial$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Bcell$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_T$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Myeloid$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Stromal$FineCellType , col.name = "FineCellType")

DimPlot(Ileum, reduction = "wnn.umap_cc",  group.by = c("FineCellType"), combine = TRUE, label.size = 3, label = TRUE) + NoLegend()
VlnPlot(Ileum, "Il6st", group.by = "FineCellType",  pt.size = 0) + NoLegend() + theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())
ggsave( "VlnPlot_Ileum_Il6st_AllCellTypes_AllSamples.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 12,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
VlnPlot(Ileum, "Il6ra", group.by = "FineCellType", pt.size = 0) + NoLegend() + theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())
ggsave( "VlnPlot_Ileum_Il6ra_AllCellTypes_AllSamples.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 12,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```







#MLN T cells and ILCs 

```{r}

DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(MLN) <- "CoarseCellType"
MLN_Tcell <- subset(x = MLN, idents = c("T cells and ILCs"))

#Repeat analysis Process
MLN_Tcell[['RNA']] <-  JoinLayers(MLN_Tcell[['RNA']])
MLN_Tcell[['RNA']] <-  split(x = MLN_Tcell[['RNA']], f = MLN_Tcell$InfectionStatus)
MLN_Tcell <- FindVariableFeatures(MLN_Tcell, selection.method = "vst", nfeatures = 3000) 
MLN_Tcell <- RunPCA(MLN_Tcell, verbose = FALSE, reduction.name = "Tcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_Tcell <- IntegrateLayers( object = MLN_Tcell, method = HarmonyIntegration,  orig.reduction = "Tcell_RNA_pca", new.reduction = "integrated.harmony.Tcell", verbose = FALSE)
DefaultAssay(MLN_Tcell) <- "ADT"
MLN_Tcell[['ADT']] <-  JoinLayers(MLN_Tcell[['ADT']])
MLN_Tcell[['ADT']] <-  split(x = MLN_Tcell[['ADT']], f = MLN_Tcell$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_Tcell)[!rownames(MLN_Tcell) %in% Isotype_labels] #remove isotypes
MLN_Tcell <- RunPCA(MLN_Tcell, features =prots, reduction.name = "Tcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_Tcell <- IntegrateLayers( object = MLN_Tcell, method = HarmonyIntegration, features = prots,  orig.reduction = "Tcell_ADT_pca", new.reduction = "integrated.ADT.harmony.Tcell", verbose = TRUE)
#Recluster
DefaultAssay(MLN_Tcell) <- "RNA"
MLN_Tcell = FindMultiModalNeighbors(MLN_Tcell, reduction.list = list("integrated.harmony.Tcell", "integrated.ADT.harmony.Tcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Tcell",  snn.graph.name = "wsnn_Tcell", weighted.nn.name = "weighted_Tcell.nn", verbose = TRUE)

MLN_Tcell <- RunUMAP(MLN_Tcell, nn.name = "weighted_Tcell.nn", reduction.name = "wnn.Tcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_Tcell <- FindClusters(MLN_Tcell, graph.name = "wsnn_Tcell",  resolution = 2, verbose = FALSE, cluster.name = "Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_Tcell, reduction = "wnn.Tcell.umap",   group.by = c( "Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_Tcell <-JoinLayers(MLN_Tcell) 
 cds <- as.cell_data_set(MLN_Tcell, reductions = "wnn.Tcell.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_Tcell[["Monocle_Tcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_Tcell[["Monocle_Tcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(MLN_Tcell, reduction = "wnn.Tcell.umap",   group.by = c( "Tcell_clusters", "Monocle_Tcell_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_Tcell, c( "Cd3e", "Cd4", "Cd8a", "Malat1"), reduction = "wnn.Tcell.umap")
VlnPlot(MLN_Tcell, "Malat1", group.by = "Monocle_Tcell_clusters",  pt.size = 0) + NoLegend()
FeatureScatter(MLN_Tcell, "Malat1", "rb.prop", pt.size = 0)
```


```{r}
Idents(MLN_Tcell) <- "Monocle_Tcell_clusters"
Low_quality_Tcell_MLN <- as.data.frame(WhichCells(subset(x = MLN_Tcell, idents = c("75", "27", "90")))) #Extreme quality outliers as evidenced by Malat1 etc
write.csv(Low_quality_Tcell_MLN, paste0(CSV, "/MLN_Tcell_LowQuality.csv"))
MLN_Tcell <- subset(x = MLN_Tcell, idents = c("75", "27", "90"), invert = T)

#Repeat analysis Process
MLN_Tcell[['RNA']] <-  JoinLayers(MLN_Tcell[['RNA']])
MLN_Tcell[['RNA']] <-  split(x = MLN_Tcell[['RNA']], f = MLN_Tcell$InfectionStatus)
MLN_Tcell <- FindVariableFeatures(MLN_Tcell, selection.method = "vst", nfeatures = 3000) 
MLN_Tcell <- RunPCA(MLN_Tcell, verbose = FALSE, reduction.name = "Tcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_Tcell <- IntegrateLayers( object = MLN_Tcell, method = HarmonyIntegration,  orig.reduction = "Tcell_RNA_pca", new.reduction = "integrated.harmony.Tcell", verbose = FALSE)
DefaultAssay(MLN_Tcell) <- "ADT"
MLN_Tcell[['ADT']] <-  JoinLayers(MLN_Tcell[['ADT']])
MLN_Tcell[['ADT']] <-  split(x = MLN_Tcell[['ADT']], f = MLN_Tcell$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_Tcell)[!rownames(MLN_Tcell) %in% Isotype_labels] #remove isotypes
MLN_Tcell <- RunPCA(MLN_Tcell, features =prots, reduction.name = "Tcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_Tcell <- IntegrateLayers( object = MLN_Tcell, method = HarmonyIntegration, features = prots,  orig.reduction = "Tcell_ADT_pca", new.reduction = "integrated.ADT.harmony.Tcell", verbose = TRUE)
#Recluster
DefaultAssay(MLN_Tcell) <- "RNA"
MLN_Tcell = FindMultiModalNeighbors(MLN_Tcell, reduction.list = list("integrated.harmony.Tcell", "integrated.ADT.harmony.Tcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Tcell",  snn.graph.name = "wsnn_Tcell", weighted.nn.name = "weighted_Tcell.nn", verbose = TRUE)

MLN_Tcell <- RunUMAP(MLN_Tcell, nn.name = "weighted_Tcell.nn", reduction.name = "wnn.Tcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_Tcell <- FindClusters(MLN_Tcell, graph.name = "wsnn_Tcell",  resolution = 2, verbose = FALSE, cluster.name = "Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_Tcell, reduction = "wnn.Tcell.umap",   group.by = c( "Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_Tcell <-JoinLayers(MLN_Tcell) 
 cds <- as.cell_data_set(MLN_Tcell, reductions = "wnn.Tcell.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_Tcell[["Monocle_Tcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_Tcell[["Monocle_Tcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(MLN_Tcell, reduction = "wnn.Tcell.umap",   group.by = c( "Tcell_clusters", "Monocle_Tcell_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_Tcell, c( "Cd3e", "Cd4", "Cd8a", "Malat1"), reduction = "wnn.Tcell.umap")
VlnPlot(MLN_Tcell, "Malat1", group.by = "Monocle_Tcell_clusters",  pt.size = 0) + NoLegend()
FeatureScatter(MLN_Tcell, "Malat1", "rb.prop", pt.size = 0)
```

```{r}

library(BiocParallel)
Idents(MLN_old) <- "Level 1 Cell Type Annotation"
unique(MLN_old$`Level 1 Cell Type Annotation`)
MLN_Tcell_old <- subset(MLN_old, idents = c("CD4 T Cells", "CD8 T Cells", "gd T cells", "Innate Lymphoid Cells"))

scRNA_counts <- GetAssayData(MLN_Tcell_old, assay = "RNA", layer = "data")
#scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(MLN_Tcell_old$`Level 3 Cell Type Annotation`)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

MLN_Tcell_counts <- GetAssayData(MLN_Tcell, assay = "RNA", layer = "data")
singleR_results_MLN_Tcell <- SingleR(test = MLN_Tcell_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

MLN_Tcell$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_MLN_Tcell, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_MLN_Tcell$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_MLN_Tcell$pruned.labels[is.na(singleR_results_MLN_Tcell$pruned.labels)] <- "Unknown"
  MLN_Tcell$Fine_SingleR<-singleR_results_MLN_Tcell$labels

#Add the Score for visualization
MLN_Tcell$Fine_SingleR_score <- NA
MLN_Tcell$Fine_SingleR_score <-singleR_results_MLN_Tcell$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
MLN_Tcell$Fine_SingleR_deltaNext <- NA
MLN_Tcell$Fine_SingleR_deltaNext<-singleR_results_MLN_Tcell$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_MLN_Tcell , ncol = 3)
plotScoreHeatmap(singleR_results_MLN_Tcell,  show.pruned = T)


FeaturePlot_scCustom(MLN_Tcell, reduction = "wnn.Tcell.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(MLN_Tcell, reduction = "wnn.Tcell.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(MLN_Tcell, reduction = "wnn.Tcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(MLN_Tcell$Monocle_Tcell_clusters)  
predicted_labels <- MLN_Tcell$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(MLN_Tcell$Monocle_Tcell_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
MLN_Tcell$FineCellType <-monocle_clusts
DimPlot_scCustom(MLN_Tcell, reduction = "wnn.Tcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Tcell, reduction = "wnn.Tcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Tcell, reduction = "wnn.Tcell.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
RNA_only <- MLN_Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Tcell, "/path/to/analysis/directory/Seurat_Files/MLN_Tcell_intermediate.qs2")

```



```{r}
# Separate out the CD4 T cells and NKt Cells for separate clustering _ Can we further resolve some of the activated populations
unique(MLN_Tcell$FineCellType)
Idents(MLN_Tcell) <- "FineCellType"
MLN_CD4Tcell <- subset(x = MLN_Tcell, idents = c("NKT cells", "Th17 cells", "Tfh cells", "CD103 high regulatory cells", "CD103 low regulatory T cells", 
                                           "CD44 high effector CD4 T cells", "Activated CD4 T cells Early", "Naive CD4 T cells", "Activated CD4 T cells Late", 
                                           "Activated CD4 T cells Intermediate", "IFN stimulated CD4 T cells"))

#Repeat analysis Process
MLN_CD4Tcell[['RNA']] <-  JoinLayers(MLN_CD4Tcell[['RNA']])
MLN_CD4Tcell[['RNA']] <-  split(x = MLN_CD4Tcell[['RNA']], f = MLN_CD4Tcell$InfectionStatus)
MLN_CD4Tcell <- FindVariableFeatures(MLN_CD4Tcell, selection.method = "vst", nfeatures = 2000) 
MLN_CD4Tcell <- RunPCA(MLN_CD4Tcell, verbose = FALSE, reduction.name = "CD4Tcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_CD4Tcell <- IntegrateLayers( object = MLN_CD4Tcell, method = HarmonyIntegration,  orig.reduction = "CD4Tcell_RNA_pca", new.reduction = "integrated.harmony.CD4Tcell", verbose = FALSE)
DefaultAssay(MLN_CD4Tcell) <- "ADT"
MLN_CD4Tcell[['ADT']] <-  JoinLayers(MLN_CD4Tcell[['ADT']])
MLN_CD4Tcell[['ADT']] <-  split(x = MLN_CD4Tcell[['ADT']], f = MLN_CD4Tcell$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_CD4Tcell)[!rownames(MLN_CD4Tcell) %in% Isotype_labels] #remove isotypes
MLN_CD4Tcell <- RunPCA(MLN_CD4Tcell, features =prots, reduction.name = "CD4Tcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_CD4Tcell <- IntegrateLayers( object = MLN_CD4Tcell, method = HarmonyIntegration, features = prots,  orig.reduction = "CD4Tcell_ADT_pca", new.reduction = "integrated.ADT.harmony.CD4Tcell", verbose = TRUE)
#Recluster
DefaultAssay(MLN_CD4Tcell) <- "RNA"
MLN_CD4Tcell = FindMultiModalNeighbors(MLN_CD4Tcell, reduction.list = list("integrated.harmony.CD4Tcell", "integrated.ADT.harmony.CD4Tcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_CD4Tcell",  snn.graph.name = "wsnn_CD4Tcell", weighted.nn.name = "weighted_CD4Tcell.nn", verbose = TRUE)

MLN_CD4Tcell <- RunUMAP(MLN_CD4Tcell, nn.name = "weighted_CD4Tcell.nn", reduction.name = "wnn.CD4Tcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_CD4Tcell <- FindClusters(MLN_CD4Tcell, graph.name = "wsnn_CD4Tcell",  resolution = 2, verbose = FALSE, cluster.name = "CD4Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap",   group.by = c( "CD4Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_CD4Tcell <-JoinLayers(MLN_CD4Tcell) 
 cds <- as.cell_data_set(MLN_CD4Tcell, reductions = "wnn.CD4Tcell.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.001, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_CD4Tcell[["Monocle_CD4Tcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_CD4Tcell[["Monocle_CD4Tcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap",   group.by = c( "CD4Tcell_clusters", "Monocle_CD4Tcell_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_CD4Tcell, c( "Cd3e", "Cd4", "Cd8a", "TCR.Bchain"), reduction = "wnn.CD4Tcell.umap")
VlnPlot(MLN_CD4Tcell, "Malat1", group.by = "Monocle_CD4Tcell_clusters",  pt.size = 0) + NoLegend()
DimPlot(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap",   group.by = c( "FineCellType"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)

FeaturePlot_scCustom(MLN_CD4Tcell, c( "Cd3e", "Cd4", "Cd8a", "TCR.Bchain"), reduction = "wnn.CD4Tcell.umap")
```

```{r}

scRNA_counts <- GetAssayData(MLN_Tcell_old, assay = "RNA", layer = "data")
#scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(MLN_Tcell_old$`Level 3 Cell Type Annotation`)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

MLN_CD4Tcell_counts <- GetAssayData(MLN_CD4Tcell, assay = "RNA", layer = "data")
singleR_results_MLN_CD4Tcell <- SingleR(test = MLN_CD4Tcell_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

MLN_CD4Tcell$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_MLN_CD4Tcell, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_MLN_CD4Tcell$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_MLN_CD4Tcell$pruned.labels[is.na(singleR_results_MLN_CD4Tcell$pruned.labels)] <- "Unknown"
  MLN_CD4Tcell$Fine_SingleR<-singleR_results_MLN_CD4Tcell$labels

#Add the Score for visualization
MLN_CD4Tcell$Fine_SingleR_score <- NA
MLN_CD4Tcell$Fine_SingleR_score <-singleR_results_MLN_CD4Tcell$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
MLN_CD4Tcell$Fine_SingleR_deltaNext <- NA
MLN_CD4Tcell$Fine_SingleR_deltaNext<-singleR_results_MLN_CD4Tcell$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_MLN_CD4Tcell , ncol = 3)
plotScoreHeatmap(singleR_results_MLN_CD4Tcell,  show.pruned = T)


FeaturePlot_scCustom(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(MLN_CD4Tcell$Monocle_CD4Tcell_clusters)  
predicted_labels <- MLN_CD4Tcell$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(MLN_CD4Tcell$Monocle_CD4Tcell_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
MLN_CD4Tcell$FineCellType <-monocle_clusts
DimPlot_scCustom(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()

```






```{r, warning=F}
RNA_only <- MLN_CD4Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD4Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD4Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_CD4Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_CD4Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_CD4Tcell, "/path/to/analysis/directory/Seurat_Files/MLN_CD4Tcell_intermediate.qs2")


```


```{r}
#adjust labels 

#PMID: 32286271 Effector memory populations expressing ccl, il7r etc 
unique(MLN_CD4Tcell$FineCellType)
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$Monocle_CD4Tcell_clusters %in% "88"] <- "Rorc Migratory Tregs"
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$Monocle_CD4Tcell_clusters %in% "22"] <- "Klrg1 Migratory Tregs"
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("59", "15")] <- "Tfr cells"
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("46", "16", "1", "65")] <- "Sell high Tregs"
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$FineCellType == "Th17 cells"] <- "NKT17 cells"
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$FineCellType == "CD44 high effector CD4 T cells"] <- "CD4 Tem"
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("34")] <- "CD4 Tem"
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("87", "44", "41", "51", "50")] <- "Activated CD4 T cells Early" # fos, Jun etc activated
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("64")] <- "Activated CD4 T cells Intermediate" 
MLN_CD4Tcell$FineCellType[MLN_CD4Tcell$FineCellType == "CD103 low regulatory T cells"] <- "Tfr cells"
DimPlot_scCustom(MLN_CD4Tcell, reduction = "wnn.CD4Tcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()

```
```{r}

RNA_only <- MLN_CD4Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD4Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD4Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_CD4Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_CD4Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_CD4Tcell, "/path/to/analysis/directory/Seurat_Files/MLN_CD4Tcell_intermediate.qs2")
```

```{r}
#Transfer labels to the MLn T cell 

MLN_Tcell <- AddMetaData(object = MLN_Tcell, metadata = MLN_CD4Tcell$FineCellType , col.name = "FineCellType")
```

```{r}
#Cd8
unique(MLN_Tcell$FineCellType)
Idents(MLN_Tcell) <- "FineCellType"
MLN_CD8Tcell <- subset(x = MLN_Tcell, idents = c("KIR positive CD8 T cells", "Activated CD8 T cells Early", "Naive CD8 T cells",
                                            "Activated CD8 T cells Late", "Effector CD8 T cells", "IFN stimulated CD8 T cells", 
                                            "gdT cells", "Activated CD8 T cells Intermediate"
                                           ))

#Repeat analysis Process
MLN_CD8Tcell[['RNA']] <-  JoinLayers(MLN_CD8Tcell[['RNA']])
MLN_CD8Tcell[['RNA']] <-  split(x = MLN_CD8Tcell[['RNA']], f = MLN_CD8Tcell$InfectionStatus)
MLN_CD8Tcell <- FindVariableFeatures(MLN_CD8Tcell, selection.method = "vst", nfeatures = 2000) 
MLN_CD8Tcell <- RunPCA(MLN_CD8Tcell, verbose = FALSE, reduction.name = "CD8Tcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_CD8Tcell <- IntegrateLayers( object = MLN_CD8Tcell, method = HarmonyIntegration,  orig.reduction = "CD8Tcell_RNA_pca", new.reduction = "integrated.harmony.CD8Tcell", verbose = FALSE)
DefaultAssay(MLN_CD8Tcell) <- "ADT"
MLN_CD8Tcell[['ADT']] <-  JoinLayers(MLN_CD8Tcell[['ADT']])
MLN_CD8Tcell[['ADT']] <-  split(x = MLN_CD8Tcell[['ADT']], f = MLN_CD8Tcell$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_CD8Tcell)[!rownames(MLN_CD8Tcell) %in% Isotype_labels] #remove isotypes
MLN_CD8Tcell <- RunPCA(MLN_CD8Tcell, features =prots, reduction.name = "CD8Tcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_CD8Tcell <- IntegrateLayers( object = MLN_CD8Tcell, method = HarmonyIntegration, features = prots,  orig.reduction = "CD8Tcell_ADT_pca", new.reduction = "integrated.ADT.harmony.CD8Tcell", verbose = TRUE)
#Recluster
DefaultAssay(MLN_CD8Tcell) <- "RNA"
MLN_CD8Tcell = FindMultiModalNeighbors(MLN_CD8Tcell, reduction.list = list("integrated.harmony.CD8Tcell", "integrated.ADT.harmony.CD8Tcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_CD8Tcell",  snn.graph.name = "wsnn_CD8Tcell", weighted.nn.name = "weighted_CD8Tcell.nn", verbose = TRUE)

MLN_CD8Tcell <- RunUMAP(MLN_CD8Tcell, nn.name = "weighted_CD8Tcell.nn", reduction.name = "wnn.CD8Tcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_CD8Tcell <- FindClusters(MLN_CD8Tcell, graph.name = "wsnn_CD8Tcell",  resolution = 2, verbose = FALSE, cluster.name = "CD8Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap",   group.by = c( "CD8Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_CD8Tcell <-JoinLayers(MLN_CD8Tcell) 
 cds <- as.cell_data_set(MLN_CD8Tcell, reductions = "wnn.CD8Tcell.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.001, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_CD8Tcell[["Monocle_CD8Tcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_CD8Tcell[["Monocle_CD8Tcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap",   group.by = c( "CD8Tcell_clusters", "Monocle_CD8Tcell_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_CD8Tcell, c( "Cd3e", "CD8", "Cd8a", "TCR.Bchain"), reduction = "wnn.CD8Tcell.umap")
VlnPlot(MLN_CD8Tcell, "Malat1", group.by = "Monocle_CD8Tcell_clusters",  pt.size = 0) + NoLegend()
DimPlot(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap",   group.by = c( "FineCellType"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)

FeaturePlot_scCustom(MLN_CD8Tcell, c( "Cd3e", "CD8", "Cd8a", "TCR.Bchain"), reduction = "wnn.CD8Tcell.umap")
```

```{r}

scRNA_counts <- GetAssayData(MLN_Tcell_old, assay = "RNA", layer = "data")
#scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(MLN_Tcell_old$`Level 3 Cell Type Annotation`)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

MLN_CD8Tcell_counts <- GetAssayData(MLN_CD8Tcell, assay = "RNA", layer = "data")
singleR_results_MLN_CD8Tcell <- SingleR(test = MLN_CD8Tcell_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

MLN_CD8Tcell$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_MLN_CD8Tcell, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_MLN_CD8Tcell$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_MLN_CD8Tcell$pruned.labels[is.na(singleR_results_MLN_CD8Tcell$pruned.labels)] <- "Unknown"
  MLN_CD8Tcell$Fine_SingleR<-singleR_results_MLN_CD8Tcell$labels

#Add the Score for visualization
MLN_CD8Tcell$Fine_SingleR_score <- NA
MLN_CD8Tcell$Fine_SingleR_score <-singleR_results_MLN_CD8Tcell$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
MLN_CD8Tcell$Fine_SingleR_deltaNext <- NA
MLN_CD8Tcell$Fine_SingleR_deltaNext<-singleR_results_MLN_CD8Tcell$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_MLN_CD8Tcell , ncol = 3)
plotScoreHeatmap(singleR_results_MLN_CD8Tcell,  show.pruned = T)


FeaturePlot_scCustom(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(MLN_CD8Tcell$Monocle_CD8Tcell_clusters)  
predicted_labels <- MLN_CD8Tcell$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(MLN_CD8Tcell$Monocle_CD8Tcell_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
MLN_CD8Tcell$FineCellType <-monocle_clusts
DimPlot_scCustom(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()

```


```{r, warning=F}
RNA_only <- MLN_CD8Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD8Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD8Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_CD8Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_CD8Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_CD8Tcell, "/path/to/analysis/directory/Seurat_Files/MLN_CD8Tcell_intermediate.qs2")


```


```{r}
#adjust labels 

#PMID: 32286271 Effector memory populations expressing ccl, il7r etc 
unique(MLN_CD8Tcell$FineCellType)
MLN_CD8Tcell$FineCellType[MLN_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("6", "9", "20", "65", "38", "59", "62")] <- "Effector CD8 T cells"
MLN_CD8Tcell$FineCellType[MLN_CD8Tcell$Monocle_CD8Tcell_clusters %in% "70"] <- "gdT cells"
MLN_CD8Tcell$FineCellType[MLN_CD8Tcell$Monocle_CD8Tcell_clusters %in% "4"] <- "Activated CD8 T cells Intermediate"
MLN_CD8Tcell$FineCellType[MLN_CD8Tcell$Monocle_CD8Tcell_clusters %in% "67"] <- "IFN stimulated CD8 T cells"
MLN_CD8Tcell$FineCellType[MLN_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("51", "25", "16")] <- "Activated CD8 T cells Early"
MLN_CD8Tcell$FineCellType[MLN_CD8Tcell$Monocle_CD8Tcell_clusters %in% "2"] <- "CD103 positive CD8 T cells"
DimPlot_scCustom(MLN_CD8Tcell, reduction = "wnn.CD8Tcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()

#Tranfer the labels to the MLn T cell object
MLN_Tcell <- AddMetaData(object = MLN_Tcell, metadata = MLN_CD8Tcell$FineCellType , col.name = "FineCellType")
```

```{r}

RNA_only <- MLN_CD8Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD8Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD8Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_CD8Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_CD8Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_CD8Tcell, "/path/to/analysis/directory/Seurat_Files/MLN_CD8Tcell_intermediate.qs2")
```




```{r}
#Visualize final labels
DimPlot(MLN_Tcell, reduction = "wnn.Tcell.umap",   group.by = c("FineCellType"), label = T) + NoLegend()
```


#MLN Myeloid

```{r}

DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(MLN) <- "CoarseCellType"
MLN_Myeloid <- subset(x = MLN, idents = c("Myeloid Cells"))

#Repeat analysis Process
MLN_Myeloid[['RNA']] <-  JoinLayers(MLN_Myeloid[['RNA']])
MLN_Myeloid[['RNA']] <-  split(x = MLN_Myeloid[['RNA']], f = MLN_Myeloid$InfectionStatus)
MLN_Myeloid <- FindVariableFeatures(MLN_Myeloid, selection.method = "vst", nfeatures = 1500) 
MLN_Myeloid <- RunPCA(MLN_Myeloid, verbose = FALSE, reduction.name = "Myeloid_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_Myeloid <- IntegrateLayers( object = MLN_Myeloid, method = HarmonyIntegration,  orig.reduction = "Myeloid_RNA_pca", new.reduction = "integrated.harmony.Myeloid", verbose = FALSE)
DefaultAssay(MLN_Myeloid) <- "ADT"
MLN_Myeloid[['ADT']] <-  JoinLayers(MLN_Myeloid[['ADT']])
MLN_Myeloid[['ADT']] <-  split(x = MLN_Myeloid[['ADT']], f = MLN_Myeloid$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_Myeloid)[!rownames(MLN_Myeloid) %in% Isotype_labels] #remove isotypes
MLN_Myeloid <- RunPCA(MLN_Myeloid, features =prots, reduction.name = "Myeloid_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_Myeloid <- IntegrateLayers( object = MLN_Myeloid, method = HarmonyIntegration, features = prots,  orig.reduction = "Myeloid_ADT_pca", new.reduction = "integrated.ADT.harmony.Myeloid", verbose = TRUE)
#Recluster
DefaultAssay(MLN_Myeloid) <- "RNA"
MLN_Myeloid = FindMultiModalNeighbors(MLN_Myeloid, reduction.list = list("integrated.harmony.Myeloid", "integrated.ADT.harmony.Myeloid"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Myeloid",  snn.graph.name = "wsnn_Myeloid", weighted.nn.name = "weighted_Myeloid.nn", verbose = TRUE)

MLN_Myeloid <- RunUMAP(MLN_Myeloid, nn.name = "weighted_Myeloid.nn", reduction.name = "wnn.Myeloid.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_Myeloid <- FindClusters(MLN_Myeloid, graph.name = "wsnn_Myeloid",  resolution = 2, verbose = FALSE, cluster.name = "Myeloid_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_Myeloid, reduction = "wnn.Myeloid.umap",   group.by = c( "Myeloid_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_Myeloid <-JoinLayers(MLN_Myeloid) 
 cds <- as.cell_data_set(MLN_Myeloid, reductions = "wnn.Myeloid.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 10, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_Myeloid[["Monocle_Myeloid_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_Myeloid[["Monocle_Myeloid_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap",   group.by = c( "Myeloid_clusters", "Monocle_Myeloid_clusters", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_Myeloid, c(  "nFeature_RNA", "Malat1", "Xcr1", "nCount_RNA"), reduction = "wnn.Myeloid.umap")
VlnPlot(MLN_Myeloid, "Malat1", group.by = "Monocle_Myeloid_clusters",  pt.size = 0) + NoLegend()
FeatureScatter(MLN_Myeloid, "Malat1", "rb.prop", pt.size = 0)
FeatureScatter(MLN_Myeloid, "Malat1", "nFeature_RNA", pt.size = 0)
```


```{r}

library(BiocParallel)
Idents(MLN_old) <- "Level 1 Cell Type Annotation"
unique(MLN_old$`Level 1 Cell Type Annotation`)
MLN_Myeloid_old <- subset(MLN_old, idents = c("Myeloid Cells"))

scRNA_counts <- GetAssayData(MLN_Myeloid_old, assay = "RNA", layer = "data")
FineCellType_vector <- as.character(MLN_Myeloid_old$`Level 3 Cell Type Annotation`)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

MLN_Myeloid_counts <- GetAssayData(MLN_Myeloid, assay = "RNA", layer = "data")
singleR_results_MLN_Myeloid <- SingleR(test = MLN_Myeloid_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

MLN_Myeloid$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_MLN_Myeloid, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_MLN_Myeloid$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_MLN_Myeloid$pruned.labels[is.na(singleR_results_MLN_Myeloid$pruned.labels)] <- "Unknown"
  MLN_Myeloid$Fine_SingleR<-singleR_results_MLN_Myeloid$labels

#Add the Score for visualization
MLN_Myeloid$Fine_SingleR_score <- NA
MLN_Myeloid$Fine_SingleR_score <-singleR_results_MLN_Myeloid$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
MLN_Myeloid$Fine_SingleR_deltaNext <- NA
MLN_Myeloid$Fine_SingleR_deltaNext<-singleR_results_MLN_Myeloid$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_MLN_Myeloid , ncol = 3)
plotScoreHeatmap(singleR_results_MLN_Myeloid,  show.pruned = T)


FeaturePlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(MLN_Myeloid$Monocle_Myeloid_clusters)  
predicted_labels <- MLN_Myeloid$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(MLN_Myeloid$Monocle_Myeloid_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
MLN_Myeloid$FineCellType <-monocle_clusts
DimPlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
RNA_only <- MLN_Myeloid
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Myeloid_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Myeloid_clusters"] <- eh 
eh <-unlist(RNA_only[["Myeloid_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Myeloid_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Myeloid_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Myeloid_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Myeloid, "/path/to/analysis/directory/Seurat_Files/MLN_Myeloid_intermediate.qs2")

```

```{r}
MLN_Myeloid$FineCellType[MLN_Myeloid$Monocle_Myeloid_clusters %in% c("51")] <- "Retnla High Macrophages"
MLN_Myeloid$FineCellType[MLN_Myeloid$Monocle_Myeloid_clusters %in% c("24", "33")] <- "Arg1 High Macrophages"
MLN_Myeloid$FineCellType[MLN_Myeloid$Monocle_Myeloid_clusters %in% c("55", "4", "11")] <- "Monocytes"
MLN_Myeloid$FineCellType[MLN_Myeloid$Monocle_Myeloid_clusters %in% c("38", "40", "34", "9", "6", "20", "7", "30")] <- "Macrophages"
DimPlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()

Idents(MLN_Myeloid) <- "Monocle_Myeloid_clusters"
cluster_50_myeloid<- WhichCells(MLN_Myeloid, idents = "50")
MLN_Myeloid <- FindSubCluster( MLN_Myeloid,  "50",  graph.name = "wsnn_Myeloid" , subcluster.name = "sub.cluster_50", resolution = 0.3,  algorithm = 3)
DimPlot(MLN_Myeloid, reduction = "wnn.Myeloid.umap",   group.by = "sub.cluster_50", label.size = 3, label = TRUE, repel = TRUE, raster = FALSE, cells = cluster_50_myeloid)+ NoLegend()
remove(cluster_50_myeloid)

MLN_Myeloid$FineCellType[MLN_Myeloid$sub.cluster_50 %in% c("50_2")] <- "Basophils"
DimPlot_scCustom(MLN_Myeloid, reduction = "wnn.Myeloid.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()

```


```{r}
RNA_only <- MLN_Myeloid
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Myeloid_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Myeloid_clusters"] <- eh 
eh <-unlist(RNA_only[["Myeloid_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Myeloid_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Myeloid_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Myeloid_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Myeloid, "/path/to/analysis/directory/Seurat_Files/MLN_Myeloid_intermediate.qs2")

```



#MLN Stromal

```{r}

DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(MLN) <- "CoarseCellType"
MLN_Stromal <- subset(x = MLN, idents = c("Stromal cells"))

#Repeat analysis Process
MLN_Stromal[['RNA']] <-  JoinLayers(MLN_Stromal[['RNA']])
MLN_Stromal[['RNA']] <-  split(x = MLN_Stromal[['RNA']], f = MLN_Stromal$InfectionStatus)
MLN_Stromal <- FindVariableFeatures(MLN_Stromal, selection.method = "vst", nfeatures = 1500) 
MLN_Stromal <- RunPCA(MLN_Stromal, verbose = FALSE, reduction.name = "Stromal_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_Stromal <- IntegrateLayers( object = MLN_Stromal, method = HarmonyIntegration,  orig.reduction = "Stromal_RNA_pca", new.reduction = "integrated.harmony.Stromal", verbose = FALSE)
DefaultAssay(MLN_Stromal) <- "ADT"
MLN_Stromal[['ADT']] <-  JoinLayers(MLN_Stromal[['ADT']])
MLN_Stromal[['ADT']] <-  split(x = MLN_Stromal[['ADT']], f = MLN_Stromal$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_Stromal)[!rownames(MLN_Stromal) %in% Isotype_labels] #remove isotypes
MLN_Stromal <- RunPCA(MLN_Stromal, features =prots, reduction.name = "Stromal_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_Stromal <- IntegrateLayers( object = MLN_Stromal, method = HarmonyIntegration, features = prots,  orig.reduction = "Stromal_ADT_pca", new.reduction = "integrated.ADT.harmony.Stromal", verbose = TRUE)
#Recluster
DefaultAssay(MLN_Stromal) <- "RNA"
MLN_Stromal = FindMultiModalNeighbors(MLN_Stromal, reduction.list = list("integrated.harmony.Stromal", "integrated.ADT.harmony.Stromal"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Stromal",  snn.graph.name = "wsnn_Stromal", weighted.nn.name = "weighted_Stromal.nn", verbose = TRUE)

MLN_Stromal <- RunUMAP(MLN_Stromal, nn.name = "weighted_Stromal.nn", reduction.name = "wnn.Stromal.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_Stromal <- FindClusters(MLN_Stromal, graph.name = "wsnn_Stromal",  resolution = 2, verbose = FALSE, cluster.name = "Stromal_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_Stromal, reduction = "wnn.Stromal.umap",   group.by = c( "Stromal_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_Stromal <-JoinLayers(MLN_Stromal) 
 cds <- as.cell_data_set(MLN_Stromal, reductions = "wnn.Stromal.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 10, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_Stromal[["Monocle_Stromal_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_Stromal[["Monocle_Stromal_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap",   group.by = c( "Stromal_clusters", "Monocle_Stromal_clusters", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_Stromal, c(  "nFeature_RNA", "Malat1", "Xcr1", "nCount_RNA"), reduction = "wnn.Stromal.umap")
VlnPlot(MLN_Stromal, "Malat1", group.by = "Monocle_Stromal_clusters",  pt.size = 0) + NoLegend()
FeatureScatter(MLN_Stromal, "Malat1", "rb.prop", pt.size = 0)
FeatureScatter(MLN_Stromal, "Malat1", "nFeature_RNA", pt.size = 0)
```


```{r}

library(BiocParallel)
Idents(MLN_old) <- "Level 1 Cell Type Annotation"
unique(MLN_old$`Level 1 Cell Type Annotation`)
MLN_Stromal_old <- subset(MLN_old, idents = c("Fibroblasts"))

scRNA_counts <- GetAssayData(MLN_Stromal_old, assay = "RNA", layer = "data")
FineCellType_vector <- as.character(MLN_Stromal_old$`Level 3 Cell Type Annotation`)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

MLN_Stromal_counts <- GetAssayData(MLN_Stromal, assay = "RNA", layer = "data")
singleR_results_MLN_Stromal <- SingleR(test = MLN_Stromal_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

MLN_Stromal$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_MLN_Stromal, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_MLN_Stromal$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_MLN_Stromal$pruned.labels[is.na(singleR_results_MLN_Stromal$pruned.labels)] <- "Unknown"
  MLN_Stromal$Fine_SingleR<-singleR_results_MLN_Stromal$labels

#Add the Score for visualization
MLN_Stromal$Fine_SingleR_score <- NA
MLN_Stromal$Fine_SingleR_score <-singleR_results_MLN_Stromal$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
MLN_Stromal$Fine_SingleR_deltaNext <- NA
MLN_Stromal$Fine_SingleR_deltaNext<-singleR_results_MLN_Stromal$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_MLN_Stromal , ncol = 3)
plotScoreHeatmap(singleR_results_MLN_Stromal,  show.pruned = T)


FeaturePlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()


#Reduce Labels
cluster_labels <- as.factor(MLN_Stromal$Monocle_Stromal_clusters)  
predicted_labels <- MLN_Stromal$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(MLN_Stromal$Monocle_Stromal_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
MLN_Stromal$FineCellType <-monocle_clusts
DimPlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
RNA_only <- MLN_Stromal
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Stromal_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Stromal_clusters"] <- eh 
eh <-unlist(RNA_only[["Stromal_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Stromal_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Stromal_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Stromal_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Stromal, "/path/to/analysis/directory/Seurat_Files/MLN_Stromal_intermediate.qs2")

```

```{r}
MLN_Stromal$FineCellType[MLN_Stromal$Stromal_clusters %in% c("0")] <- "Blood endothelial cells"
MLN_Stromal$FineCellType[MLN_Stromal$Stromal_clusters %in% c("7")] <- "Lymphatic endothelial cells"
MLN_Stromal$FineCellType[MLN_Stromal$Stromal_clusters %in% c("5")] <- "Mesothelial cells"
MLN_Stromal$FineCellType[MLN_Stromal$Stromal_clusters %in% c("1", "2", "3", "4", "6", "8", "9", "singleton")] <- "Fibroblasts"

DimPlot_scCustom(MLN_Stromal, reduction = "wnn.Stromal.umap", group.by = "FineCellType",  label = T , combine = T, repel = T) + NoLegend()
```


```{r}
RNA_only <- MLN_Stromal
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Stromal_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Stromal_clusters"] <- eh 
eh <-unlist(RNA_only[["Stromal_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Stromal_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Stromal_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Stromal_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Stromal, "/path/to/analysis/directory/Seurat_Files/MLN_Stromal_intermediate.qs2")

```


#MLN Bcell

```{r}

DimPlot(MLN, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(MLN) <- "CoarseCellType"
MLN_Bcell <- subset(x = MLN, idents = c("B Cells"))

#Repeat analysis Process
MLN_Bcell[['RNA']] <-  JoinLayers(MLN_Bcell[['RNA']])
MLN_Bcell[['RNA']] <-  split(x = MLN_Bcell[['RNA']], f = MLN_Bcell$InfectionStatus)
MLN_Bcell <- FindVariableFeatures(MLN_Bcell, selection.method = "vst", nfeatures = 2000) 
MLN_Bcell <- RunPCA(MLN_Bcell, verbose = FALSE, reduction.name = "Bcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_Bcell <- IntegrateLayers( object = MLN_Bcell, method = HarmonyIntegration,  orig.reduction = "Bcell_RNA_pca", new.reduction = "integrated.harmony.Bcell", verbose = FALSE)
DefaultAssay(MLN_Bcell) <- "ADT"
MLN_Bcell[['ADT']] <-  JoinLayers(MLN_Bcell[['ADT']])
MLN_Bcell[['ADT']] <-  split(x = MLN_Bcell[['ADT']], f = MLN_Bcell$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_Bcell)[!rownames(MLN_Bcell) %in% Isotype_labels] #remove isotypes
MLN_Bcell <- RunPCA(MLN_Bcell, features =prots, reduction.name = "Bcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_Bcell <- IntegrateLayers( object = MLN_Bcell, method = HarmonyIntegration, features = prots,  orig.reduction = "Bcell_ADT_pca", new.reduction = "integrated.ADT.harmony.Bcell", verbose = TRUE)
#Recluster
DefaultAssay(MLN_Bcell) <- "RNA"
MLN_Bcell = FindMultiModalNeighbors(MLN_Bcell, reduction.list = list("integrated.harmony.Bcell", "integrated.ADT.harmony.Bcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Bcell",  snn.graph.name = "wsnn_Bcell", weighted.nn.name = "weighted_Bcell.nn", verbose = TRUE)

MLN_Bcell <- RunUMAP(MLN_Bcell, nn.name = "weighted_Bcell.nn", reduction.name = "wnn.Bcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_Bcell <- FindClusters(MLN_Bcell, graph.name = "wsnn_Bcell",  resolution = 2, verbose = FALSE, cluster.name = "Bcell_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_Bcell <-JoinLayers(MLN_Bcell) 
 cds <- as.cell_data_set(MLN_Bcell, reductions = "wnn.Bcell.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 10, cluster_method =  "leiden", resolution = 0.001, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_Bcell[["Monocle_Bcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_Bcell[["Monocle_Bcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "Monocle_Bcell_clusters", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_Bcell, c(  "nFeature_RNA", "Malat1", "Xcr1", "nCount_RNA"), reduction = "wnn.Bcell.umap")
VlnPlot(MLN_Bcell, "Malat1", group.by = "Monocle_Bcell_clusters",  pt.size = 0) + NoLegend()
FeatureScatter(MLN_Bcell, "Malat1", "rb.prop", pt.size = 0)
FeatureScatter(MLN_Bcell, "Malat1", "nFeature_RNA", pt.size = 0)
```



```{r}
RNA_only <- MLN_Bcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Bcell_clusters"] <- eh 
eh <-unlist(RNA_only[["Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Bcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Bcell, "/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.qs2")

```

```{r}

#Remove Doublets and other popualtions
# 201 , 90, 1  - These are Lyz2+ clearly myleoid doublets
# 191 - T cells - doublets CD3e 
# 157, 139 - Malat1 high cells with few genes  - poor quality

Idents(MLN_Bcell) <- "Monocle_Bcell_clusters"
cluster_41_Epithelial<- WhichCells(MLN_Bcell, idents = "41")
MLN_Bcell <- FindSubCluster( MLN_Bcell,  "41",  graph.name = "wsnn_Bcell" , subcluster.name = "sub.cluster_41", resolution = 0.08,  algorithm = 3)
DimPlot(MLN_Bcell, reduction = "wnn.Bcell.umap",   group.by = "sub.cluster_41", label.size = 3, label = TRUE, repel = TRUE, raster = FALSE, cells= cluster_41_Epithelial)+ NoLegend()

# 41_1 T cells 

Idents(MLN_Bcell) <- "sub.cluster_41"
Low_quality_Bcell_MLN <- as.data.frame(WhichCells(subset(x = MLN_Bcell, idents = c("201" , "90", "1", "191", "157", "139", "41_1")))) #Extreme quality outliers as evidenced by Malat1 etc
write.csv(Low_quality_Bcell_MLN, paste0(CSV, "/MLN_Bcell_LowQuality.csv"))
MLN_Bcell <- subset(x = MLN_Bcell, idents = c("201" , "90", "1", "191", "157", "139", "41_1"), invert = T)

MLN_Bcell[['RNA']] <-  JoinLayers(MLN_Bcell[['RNA']])
MLN_Bcell[['RNA']] <-  split(x = MLN_Bcell[['RNA']], f = MLN_Bcell$InfectionStatus)
MLN_Bcell <- FindVariableFeatures(MLN_Bcell, selection.method = "vst", nfeatures = 2000) 
MLN_Bcell <- RunPCA(MLN_Bcell, verbose = FALSE, reduction.name = "Bcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
MLN_Bcell <- IntegrateLayers( object = MLN_Bcell, method = HarmonyIntegration,  orig.reduction = "Bcell_RNA_pca", new.reduction = "integrated.harmony.Bcell", verbose = FALSE)
DefaultAssay(MLN_Bcell) <- "ADT"
MLN_Bcell[['ADT']] <-  JoinLayers(MLN_Bcell[['ADT']])
MLN_Bcell[['ADT']] <-  split(x = MLN_Bcell[['ADT']], f = MLN_Bcell$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN[["ADT"]]), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_Bcell)[!rownames(MLN_Bcell) %in% Isotype_labels] #remove isotypes
MLN_Bcell <- RunPCA(MLN_Bcell, features =prots, reduction.name = "Bcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_Bcell <- IntegrateLayers( object = MLN_Bcell, method = HarmonyIntegration, features = prots,  orig.reduction = "Bcell_ADT_pca", new.reduction = "integrated.ADT.harmony.Bcell", verbose = TRUE)
#Recluster
DefaultAssay(MLN_Bcell) <- "RNA"
MLN_Bcell = FindMultiModalNeighbors(MLN_Bcell, reduction.list = list("integrated.harmony.Bcell", "integrated.ADT.harmony.Bcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Bcell",  snn.graph.name = "wsnn_Bcell", weighted.nn.name = "weighted_Bcell.nn", verbose = TRUE)

MLN_Bcell <- RunUMAP(MLN_Bcell, nn.name = "weighted_Bcell.nn", reduction.name = "wnn.Bcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_Bcell <- FindClusters(MLN_Bcell, graph.name = "wsnn_Bcell",  resolution = 2, verbose = FALSE, cluster.name = "Bcell_clusters", n.start = 20, group.singletons = F)
DimPlot(MLN_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
MLN_Bcell <-JoinLayers(MLN_Bcell) 
 cds <- as.cell_data_set(MLN_Bcell, reductions = "wnn.Bcell.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 10, cluster_method =  "leiden", resolution = 0.001, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_Bcell[["Monocle_Bcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_Bcell[["Monocle_Bcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap",   group.by = c( "Bcell_clusters", "Monocle_Bcell_clusters", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(MLN_Bcell, c(  "nFeature_RNA", "Malat1", "Xcr1", "nCount_RNA"), reduction = "wnn.Bcell.umap")
VlnPlot(MLN_Bcell, "Malat1", group.by = "Monocle_Bcell_clusters",  pt.size = 0) + NoLegend()
FeatureScatter(MLN_Bcell, "Malat1", "rb.prop", pt.size = 0)
FeatureScatter(MLN_Bcell, "Malat1", "nFeature_RNA", pt.size = 0)
```

```{r}

library(BiocParallel)
Idents(MLN_old) <- "Level 1 Cell Type Annotation"
unique(MLN_old$`Level 1 Cell Type Annotation`)
MLN_Bcell_old <- subset(MLN_old, idents = c("B Cells"))

scRNA_counts <- GetAssayData(MLN_Bcell_old, assay = "RNA", layer = "data")
FineCellType_vector <- as.character(MLN_Bcell_old$`Level 3 Cell Type Annotation`)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

MLN_Bcell_counts <- GetAssayData(MLN_Bcell, assay = "RNA", layer = "data")
singleR_results_MLN_Bcell <- SingleR(test = MLN_Bcell_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(20))

MLN_Bcell$Fine_SingleR <- NA
singleR_results_pruned <-pruneScores(singleR_results_MLN_Bcell, nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
print(sum(singleR_results_pruned == TRUE))

#Transfer labels
singleR_results_MLN_Bcell$pruned.labels[singleR_results_pruned] <- "Unknown"
  singleR_results_MLN_Bcell$pruned.labels[is.na(singleR_results_MLN_Bcell$pruned.labels)] <- "Unknown"
  MLN_Bcell$Fine_SingleR<-singleR_results_MLN_Bcell$labels

#Add the Score for visualization
MLN_Bcell$Fine_SingleR_score <- NA
MLN_Bcell$Fine_SingleR_score <-singleR_results_MLN_Bcell$scores

#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
MLN_Bcell$Fine_SingleR_deltaNext <- NA
MLN_Bcell$Fine_SingleR_deltaNext<-singleR_results_MLN_Bcell$delta.next

#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_MLN_Bcell , ncol = 3)
plotScoreHeatmap(singleR_results_MLN_Bcell,  show.pruned = T)
FeaturePlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap", features = "Fine_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap", features = "Fine_SingleR_deltaNext", na_cutoff = 0) 
DimPlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()

#Reduce Labels
cluster_labels <- as.factor(MLN_Bcell$Monocle_Bcell_clusters)  
predicted_labels <- MLN_Bcell$Fine_SingleR  

# Apply this function to each cluster
best_labels <- sapply(levels(cluster_labels), get_best_label)
# Ensure 'best_labels' is a named vector, where names correspond to cluster ids
monocle_clusts <- as.character(MLN_Bcell$Monocle_Bcell_clusters) 
monocle_clusts <- recode(
  as.character(monocle_clusts), 
  !!!best_labels)
MLN_Bcell$FineCellType <-monocle_clusts
DimPlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap", group.by = "Fine_SingleR" , label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()
DimPlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap", group.by = "InfectionStatus",  label = T , combine = T, repel = T, raster = T) + NoLegend()
```

```{r}
RNA_only <- MLN_Bcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Bcell_clusters"] <- eh 
eh <-unlist(RNA_only[["Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Bcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Bcell, "/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.qs2")

```

```{r}
#Fix B cells Labels and Save 
#cluster 117, 75, 110 B cells 
MLN_Bcell$FineCellType[MLN_Bcell$FineCellType %in% c("Naive B cells", "Activated B cells")] <- "B cells"
MLN_Bcell$FineCellType[MLN_Bcell$Monocle_Bcell_clusters %in% c("17", "75", "110", "109", "168")] <- "B cells"
DimPlot_scCustom(MLN_Bcell, reduction = "wnn.Bcell.umap", group.by = "FineCellType",  label = T , combine = T, repel = T, raster = T) + NoLegend()

```


```{r}
RNA_only <- MLN_Bcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Bcell_clusters"] <- eh 
eh <-unlist(RNA_only[["Bcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Bcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(MLN_Bcell, "/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.qs2")

```


#Add the Fine Labels (Level3)

```{r}
unique(Ileum$CoarseCellType)
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Stromal$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_T$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Myeloid$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Bcell$FineCellType , col.name = "FineCellType")
Ileum <- AddMetaData(object = Ileum, metadata = Ileum_Epithelial$FineCellType , col.name = "FineCellType")
Ileum@meta.data$FineCellType[Ileum$CoarseCellType == "Enteric Nervous System" ] <- "Enteric Glia"
Ileum@meta.data$FineCellType[Ileum$CoarseCellType == "Plasma cells" ] <- "Plasma cells"
DimPlot(Ileum, reduction = "wnn.umap_cc", group.by = "FineCellType", label = T, repel = T, label.size = 1) + NoLegend()
table(Ileum$FineCellType) # There are 6351 that are NA - 2883 removed from B cell, 2770 removed from T cells, 696 removed from myeloid cells  = 6349 cells (+ singletons not clustered into these major groups)
Idents(Ileum) <- "FineCellType"
Low_quality_Ileum <- WhichCells(Ileum, idents = "NA")
write.csv(Low_quality_Ileum, paste0(CSV, "/Total_LowQuality_Ileum_Removed_during_Annotation.csv"))

Ileum_filt <- subset(Ileum, idents = c( "NA"), invert = TRUE )
DimPlot(Ileum_filt, reduction = "wnn.umap_cc", group.by = "FineCellType", label = T, repel = T, label.size = 1) + NoLegend()


unique(MLN$CoarseCellType)

MLN <- AddMetaData(object = MLN, metadata = MLN_Bcell$FineCellType , col.name = "FineCellType")
MLN <- AddMetaData(object = MLN, metadata = MLN_T$FineCellType , col.name = "FineCellType")
MLN <- AddMetaData(object = MLN, metadata = MLN_Tcell$FineCellType , col.name = "FineCellType")
MLN <- AddMetaData(object = MLN, metadata = MLN_CD4Tcell$FineCellType , col.name = "FineCellType")
MLN <- AddMetaData(object = MLN, metadata = MLN_CD8Tcell$FineCellType , col.name = "FineCellType")
MLN <- AddMetaData(object = MLN, metadata = MLN_Myeloid$FineCellType , col.name = "FineCellType")
MLN <- AddMetaData(object = MLN, metadata = MLN_Stromal$FineCellType , col.name = "FineCellType")
MLN@meta.data$FineCellType[MLN$CoarseCellType == "Plasma cells" ] <- "Plasma cells"
MLN@meta.data$FineCellType[is.na(MLN$FineCellType)] <- "Junk"
table(MLN$FineCellType) # There are 4588 that are Junk - 2215 removed from B cell, 2373 removed from T cells  = 4588 cells
Idents(MLN) <- "FineCellType"
Low_quality_MLN <- WhichCells(MLN, idents = "Junk")
write.csv(Low_quality_MLN, paste0(CSV, "/Total_LowQuality_MLN_Removed_during_Annotation.csv"))

MLN_filt <- subset(MLN, idents = c( "Junk"), invert = TRUE )
DimPlot(MLN, reduction = "wnn.umap_cc", group.by = "FineCellType", label = T) + NoLegend()
```

#Add MajorCellTypes (an intermediate level of labeling used for graphing)

```{r}
#Ileum 
unique(Ileum_filt$FineCellType)

Gdcells <- c("Gzmb gd T Cells",  "Cd160 gd T Cells" , "gd T17 Cells",  "Ly6c gd T Cells" ,  "Sell high gd T Cells")
CD4abcells <- c("Tfh" ,  "CD103 low CD4 Trm", "Sell high CD4 T Cells",   "Treg",  "CD103 high CD4 Trm",  "Sell high Treg" ,  "IFN stimulated CD4 T Cells" )
CD8abcells <- c("Sell high CD8 T Cells" , "CD103 high CD8 Trm",  "Cycling T Cells",  "Ly6c high CD8 Trm" ,  "IFN Stimulated CD8 T Cells"   )
ILCs <- c("ILC2s" , "LTIs", "ILC3s","ILC1s" ,  "NK Cells",  "Cycling ILCS")
NKT <- c("NKT Cells" , "NKT2 Cells" )
Bcells <- c( "Dark Zone GC B Cells" ,  "Activated B Cells" ,  "Klhl14 B cells" , "B Cells",  "Light Zone GC B Cells", "Sell negative proliferative B Cells",  "IFN stimulated B Cells" )
PlasmaCells <- c( "Plasma cells" )
DendriticCells <- c("Itgae High cDC2s" , "Cd301b High cDCs", "cDC1s",  "tDCs", "Ccr7 High cDC2s")
pDCs<- c("Plasmacytoid DCs")
Macrophages <- c("Resident Macrophages" ,"Arg1 High Macrophages" ,"Retnla High Macrophages"   )
Monocytes<- c( "Monocytes" , "Inflammatory Monocytes"  )
Neutrophils<- c("LP Neutrophils" )
MastCells <- c("Mast Cells" )
Enterocytes <- c("Enterocytes")
GobletCells <- c("Goblet Cells", "Differentiating Goblet Cells")
PanethCells <- c( "Paneth Cells" )
Enteroendocrinecells <- c( "Enteroendocrine Cells", "Enterochromaffin Cells" )
TuftCells <- c("Tuft Cells")
Stemcells<- c("Transit Amplifying Cells",  "Stem Cells" )
EntericGlia<- c( "Enteric Glia")
Mesenchymal <-  c("Trophocytes" , "Interstitial nonvascular smooth muscle", "Visceral nonvascular smooth muscle", "Cd81 negative fibroblasts" , "Ackr4 high trophocytes",  "Pdgfra and Cd81 high fibroblasts" ,  "Inflammatory trophocytes" ,   "Pln high vascular smooth muscle","Chrm2 high vascular smooth muscle",  "Interstitial cells of cajal"  ,"Pericytes" ,  "Telocytes", "Mesothelial cells" )    # (Fibroblasts, ICC, pericytes, telocytes, smooth muscle cells, trophocytes)
Endothelial<-  c("Lymphatic endothelial cells", "Artery BEC",  "Capillary BEC", "Venule BEC")    # (Blood Endothelial cells Lymphatic endothelial cells )



Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Gdcells ] <- "gd T Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% CD4abcells ] <- "CD4 ab T Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% CD8abcells ] <- "CD8 ab T Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% ILCs ] <- "Innate Lymphoid Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% NKT ] <- "NKT Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Bcells ] <- "B Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% PlasmaCells ] <- "Plasma Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% DendriticCells ] <- "Dendritic Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% pDCs ] <- "Plasmacytoid DCs"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Macrophages ] <- "Macrophages"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Monocytes ] <- "Monocytes"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Neutrophils ] <- "Neutrophils"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% MastCells ] <- "Mast Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Enterocytes ] <- "Enterocytes"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% GobletCells ] <- "Goblet Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% PanethCells ] <- "Paneth Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Enteroendocrinecells ] <- "Enteroendocrine Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% TuftCells ] <- "Tuft Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Stemcells ] <- "Stem / Transit Amplifying Cells"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% EntericGlia ] <- "Enteric Nervous System"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Mesenchymal ] <- "Mesenchymal"
Ileum_filt@meta.data$MajorCellType[Ileum_filt$FineCellType %in% Endothelial ] <- "Endothelial"

DimPlot(Ileum_filt, reduction = "wnn.umap_cc", group.by = "MajorCellType", label = T, label.size = 2) + NoLegend()

```


```{r}
#Make the Level 2 filter  Ileum
unique(Ileum_filt$FineCellType)

Enterocytes <- c("Enterocytes"  )          
ILC2 <- c("ILC2s" )     
GobletCells <- c( "Goblet Cells",  "Differentiating Goblet Cells"  )
GCBcells <- c("Dark Zone GC B Cells",  "Light Zone GC B Cells")
TA <- c("Transit Amplifying Cells" )
ILC3 <- c("LTIs" , "ILC3s" ) 
CyclingILC <- c("Cycling ILCS")
CyclingT <- c( "Cycling T Cells")
LEC <- c("Lymphatic endothelial cells")
gdTcells <- c("Gzmb gd T Cells",  "Cd160 gd T Cells" , "gd T17 Cells",  "Ly6c gd T Cells" ,  "Sell high gd T Cells")
EffectorCD4 <- c("Tfh" ,  "CD103 low CD4 Trm",   "Treg",  "CD103 high CD4 Trm",  "Sell high Treg"  )
SmoothMuscle <- c("Interstitial nonvascular smooth muscle", "Visceral nonvascular smooth muscle",  "Pln high vascular smooth muscle","Chrm2 high vascular smooth muscle" )
Bcells <- c(  "Activated B Cells" ,  "Klhl14 B cells" , "B Cells",  "Sell negative proliferative B Cells",  "IFN stimulated B Cells" )
CD8abT <- c( "Sell high CD8 T Cells" ,   "IFN Stimulated CD8 T Cells"  )
Mono <- c( "Monocytes" , "Inflammatory Monocytes"  )
Glia <- c("Enteric Glia"  )
Fibroblasts <- c("Trophocytes" ,  "Cd81 negative fibroblasts" , "Ackr4 high trophocytes",  "Pdgfra and Cd81 high fibroblasts" ,  "Inflammatory trophocytes" )
Stem <- c("Stem Cells")
NKT <- c("NKT Cells" , "NKT2 Cells" )
CD4ab <- c( "Sell high CD4 T Cells", "IFN stimulated CD4 T Cells"  )
BEC <- c(  "Artery BEC",  "Capillary BEC", "Venule BEC")
cDC2 <- c("Itgae High cDC2s" , "Cd301b High cDCs",   "tDCs", "Ccr7 High cDC2s")
cDC1 <-  c("cDC1s" )
plasma <- c( "Plasma cells"  )
Macrophages <- c( "Resident Macrophages" ,"Arg1 High Macrophages" ,"Retnla High Macrophages"   )
bMem <- c("Memory B Cells"  )
cyclingT <- c("Cycling T Cells" )
Paneth <- c( "Paneth Cells" )
Tuft <- c( "Tuft Cells"   )
effectorCD8 <- c( "CD103 high CD8 Trm",   "Ly6c high CD8 Trm"  )
pDC <- c("Plasmacytoid DCs")      
ILC1 <- c("ILC1s"  )
EEC <- c("Enteroendocrine Cells"  ,  "Enterochromaffin Cells"  )
NK <- c(   "NK Cells"  )
Mast <- c( "Mast Cells" )
Neutrophils <- c("LP Neutrophils"    )
Mesothelial <- c( "Mesothelial cells"  )  
Stromal <-c(  "Interstitial cells of cajal"  ,"Pericytes" ,  "Telocytes" )

Ileum_filt$secondlevel <- Ileum_filt$FineCellType 

Ileum_filt$secondlevel[Ileum_filt$FineCellType %in%Enterocytes ] <- "Enterocytes"
  Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% ILC2] <- "ILC2s"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% GobletCells] <- "Goblet Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in%GCBcells ] <- "GC B Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in%TA ] <- "Transit Amplifying Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% ILC3] <- "ILC3s"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% LEC] <- "Lymphatic Endothelial Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in%gdTcells] <- "gd T Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% EffectorCD4] <- "Effector CD4 ab T Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in%SmoothMuscle ] <- "Smooth Muscle Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Bcells] <- "B Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% CD8abT] <- "CD8 ab T Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Mono] <- "Monocytes"
  Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Glia] <- "Enteric Glia"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Fibroblasts] <- "Fibroblasts"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% NKT] <- "NKT Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% CD4ab] <- "CD4 ab T Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% BEC] <- "Blood Endothelial Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% CyclingILC] <- "Cycling ILCs"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% CyclingT] <- "Cycling T Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% cDC1] <- "cDC1s"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% cDC2] <- "cDC2s"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% plasma ] <- "Plasma Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Macrophages] <- "Macrophages"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% bMem] <- "Memory B cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in%Paneth ] <- "Paneth Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Tuft] <- "Tuft Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% effectorCD8] <- "Effector CD8 ab T Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% pDC] <- "Plasmacytoid DCs"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% ILC1 ] <- "ILC1s"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% EEC] <- "Enteroendocrine Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% NK] <- "NK Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Mast] <- "Mast Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Neutrophils] <- "Neutrophils"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Mesothelial] <- "Mesothelial Cells"
Ileum_filt$secondlevel[Ileum_filt$FineCellType %in% Stromal] <- "Stromal Cells"


DimPlot(Ileum_filt, reduction = "wnn.umap_cc", group.by = "secondlevel", label = T, label.size = 2) + NoLegend()
```
```{r}


```





```{r}
#MLN 
unique(MLN_filt$FineCellType)

Gdcells <- c("gdT cells"  )
CD4abcells <- c( "Activated CD4 T cells Early",  "Tfh cells",  "Sell high Tregs", "CD4 Tem" ,  "Naive CD4 T cells", "Activated CD4 T cells Late", "Klrg1 Migratory Tregs", "Tfr cells" ,  "Rorc Migratory Tregs", "IFN stimulated CD4 T cells"  ,   "Activated CD4 T cells Intermediate" )
CD8abcells <- c(  "Naive CD8 T cells",   "Activated CD8 T cells Early",    "Effector CD8 T cells" , "Activated CD8 T cells Late" , "CD103 positive CD8 T cells" , "KIR positive CD8 T cells",  "IFN stimulated CD8 T cells", "Activated CD8 T cells Intermediate"    )
ILCs <- c( "NK cells", "ILC1s",  "ILC2s" ,  "ILC3s")
NKT <- c( "NKT cells",  "NKT17 cells"   )
Bcells <- c("B cells" , "GC Committed B cells", "Memory B cells" , "LZ Germinal Center",  "DZ Germinal Center" , "IFN Activated B cells"  )
PlasmaCells <- c("Plasma cells" )
DendriticCells <- c("Ccr7 cDC2s" ,  "Xcr1 cDC1s" , "Sipra cDC2s"  )
pDCs<- c( "Plasmacytoid DCs")
Macrophages <- c(  "Macrophages" , "Cycling Myeloid" ,     "CD209b Macrophages" ,   "Arg1 High Macrophages" , "Retnla High Macrophages"   )
Monocytes<- c(  "Monocytes" )
Neutrophils<- c("Neutrophils"  )
Basophil<- c( "Basophils"  )
Mesenchymal <-  c(   "Fibroblasts"  , "Mesothelial cells" )    # (Fibroblasts, ICC, pericytes, telocytes, smooth muscle cells, trophocytes)
Endothelial<-  c( "Blood endothelial cells", "Lymphatic endothelial cells"     )    # (Blood Endothelial cells Lymphatic endothelial cells )

                                                                 


MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Gdcells ] <- "gd T Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% CD4abcells ] <- "CD4 ab T Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% CD8abcells ] <- "CD8 ab T Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% ILCs ] <- "Innate Lymphoid Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% NKT ] <- "NKT Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Bcells ] <- "B Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% PlasmaCells ] <- "Plasma Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% DendriticCells ] <- "Dendritic Cells"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% pDCs ] <- "Plasmacytoid DCs"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Macrophages ] <- "Macrophages"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Monocytes ] <- "Monocytes"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Neutrophils ] <- "Neutrophils"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Basophil ] <- "Basophils"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Mesenchymal ] <- "Mesenchymal"
MLN_filt@meta.data$MajorCellType[MLN_filt$FineCellType %in% Endothelial ] <- "Endothelial"

DimPlot(MLN_filt, reduction = "wnn.umap_cc", group.by = "MajorCellType", label = T, label.size = 2) + NoLegend()

```


```{r}
#Make the Level 2 filter  MLN
unique(MLN_filt$FineCellType)


ILC2 <- c( "ILC2s")     
GCBcells <- c("LZ Germinal Center",  "DZ Germinal Center")
ILC3 <- c( "ILC3s") 
LEC <- c("Lymphatic endothelial cells")
gdTcells <- c("gdT cells" )
EffectorCD4 <- c("Tfh cells",  "Klrg1 Migratory Tregs", "Tfr cells",  "Rorc Migratory Tregs" ,  "Sell high Tregs")
CD4Tmem <- c("CD4 Tem" )
Bcells <- c("B cells" , "GC Committed B cells"  , "IFN Activated B cells"  )
CD8abT <- c(  "Naive CD8 T cells",   "Activated CD8 T cells Early",    "Activated CD8 T cells Late" , "IFN stimulated CD8 T cells", "Activated CD8 T cells Intermediate" )
Mono <- c( "Monocytes"  )
Fibroblasts <- c(  "Fibroblasts" )
NKT <- c("NKT cells",  "NKT17 cells"  )
CD4ab <- c(  "Naive CD4 T cells", "IFN stimulated CD4 T cells" ,"Activated CD4 T cells Early",   "Activated CD4 T cells Late",  "Activated CD4 T cells Intermediate" )
BEC <- c("Blood endothelial cells")
cDC2 <- c("Sipra cDC2s" , "Ccr7 cDC2s")
cDC1 <-  c(  "Xcr1 cDC1s")
plasma <- c( "Plasma cells")
Macrophages <- c( "Macrophages" , "Cycling Myeloid" ,     "CD209b Macrophages" ,   "Arg1 High Macrophages" , "Retnla High Macrophages"  )
bMem <- c(  "Memory B cells" )
effectorCD8 <- c(  "Effector CD8 T cells",  "CD103 positive CD8 T cells" , "KIR positive CD8 T cells" )
pDC <- c( "Plasmacytoid DCs")      
ILC1 <- c( "ILC1s" )
NK <- c( "NK cells"  )
Baso <- c( "Basophils" )
Neutrophils <- c( "Neutrophils"  )
Mesothelial <- c("Mesothelial cells" )  


MLN_filt$secondlevel <- MLN_filt$FineCellType 


MLN_filt$secondlevel[MLN_filt$FineCellType %in% ILC2] <- "ILC2s"
MLN_filt$secondlevel[MLN_filt$FineCellType %in%GCBcells ] <- "GC B Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% ILC3] <- "ILC3s"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% LEC] <- "Lymphatic Endothelial Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in%gdTcells] <- "gd T Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% EffectorCD4] <- "Effector CD4 ab T Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% CD4Tmem] <- "Memory CD4 ab T Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% Bcells] <- "B Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% CD8abT] <- "CD8 ab T Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% Mono] <- "Monocytes"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% Fibroblasts] <- "Fibroblasts"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% NKT] <- "NKT Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% CD4ab] <- "CD4 ab T Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% BEC] <- "Blood Endothelial Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% CyclingILC] <- "Cycling ILCs"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% cDC1] <- "cDC1s"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% cDC2] <- "cDC2s"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% plasma ] <- "Plasma Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% Macrophages] <- "Macrophages"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% bMem] <- "Memory B cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% effectorCD8] <- "Effector CD8 ab T Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% pDC] <- "Plasmacytoid DCs"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% ILC1 ] <- "ILC1s"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% NK] <- "NK Cells"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% Baso] <- "Basophils"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% Neutrophils] <- "Neutrophils"
MLN_filt$secondlevel[MLN_filt$FineCellType %in% Mesothelial] <- "Mesothelial Cells"


DimPlot(MLN_filt, reduction = "wnn.umap_cc", group.by = "secondlevel", label = T, label.size = 2, repel = T) + NoLegend()
table()

```


#Rerun UMAP for Ileum after Filter

```{r}
#This object is the loaded Ileum object - which has already been scaled to regress out the Cell cycle genes.Because I eliminated ~7000 cells in filtering, a final umap will be made here

DefaultAssay(Ileum_filt) <- "RNA"
options(future.globals.maxSize = 500000 * 1024^2)

Ileum_filt[["RNA"]] <- as(Ileum_filt[["RNA"]], "Assay5")
Ileum_filt <- RunPCA(Ileum_filt, verbose = FALSE, reduction.name = "RNA_pca_cc", features = VariableFeatures(Ileum_filt, nfeatures = 3000))
Ileum_filt[["RNA"]] <- JoinLayers(Ileum_filt[['RNA']])
Ileum_filt[['RNA']] <-  split(x = Ileum_filt[['RNA']], f = Ileum_filt$Infection_Organ)
Ileum_filt <- IntegrateLayers(object = Ileum_filt, method = HarmonyIntegration, assay = "RNA", features = VariableFeatures(Ileum_filt, nfeatures = 3000),
  orig.reduction = "RNA_pca_cc", new.reduction = "reintegrated_RNA_harmony_cc", verbose = T)

Ileum_filt <- FindNeighbors(Ileum_filt, reduction = "reintegrated_RNA_harmony_cc", dims = 1:35, graph.name = "harmony_rna_cc")
Ileum_filt <- FindClusters(Ileum_filt, resolution = 2, cluster.name = "harmony_integration_clusters_cc",  graph.name = "harmony_rna_cc", group.singletons = F, n.start = 20)
Ileum_filt <- RunUMAP(Ileum_filt, reduction = "reintegrated_RNA_harmony_cc", dims = 1:35, reduction.name = "umap.harmony_cc",  seed.use = 20)


Ileum_filt = FindMultiModalNeighbors(
  Ileum_filt, reduction.list = list("reintegrated_RNA_harmony_cc", "integrated.ADT.harmony"), 
  dims.list = list(1:35, 1:20), 
  modality.weight.name = c("RNA.weight", "ADT.weight"), 
  knn.graph.name = "wknn_cc",
  snn.graph.name = "wsnn_cc",
  weighted.nn.name = "weighted.nn_cc",
  verbose = FALSE)

#compare clusters and Umaps before and after CC regression
DefaultAssay(Ileum_filt) <- "RNA"
Ileum_filt <- RunUMAP(Ileum_filt, nn.name = "weighted.nn_cc", reduction.name = "wnn.umap_cc", reduction.key = "wnnUMAPcc_", seed.use = 20)
Ileum_filt <- FindClusters(Ileum_filt, graph.name = "wsnn_cc",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters_cc", n.start = 25, group.singletons = F)


#Monocle Clusters
library(SeuratWrappers)
Ileum_filt <- JoinLayers(Ileum_filt)

 cds <- as.cell_data_set(
  Ileum_filt,
  reductions = "wnn.umap_cc",
  default.reduction = "wnn.umap_cc",
  graph =  "wsnn",
  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells( cds, reduction_method = "UMAP", k = 15, cluster_method =  "leiden", resolution = 0.00015,
  num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_filt[["Monocle_clusters_cc"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_filt[["Monocle_partitions_cc"]] <- Seurat_monocle@meta.data$monocle3_partitions

#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)
#Save finalized Ileum Object 
qs_save(Ileum_filt, "/path/to/analysis/directory/Seurat_Files/Analysis2_clustered_annotated_Ileum.qs2")
```


```{r}
DimPlot(Ileum_filt, reduction = "wnn.umap_cc", group.by = "CoarseCellType", label = T, label.size = 3, combine = T) + NoLegend()

```




```{r}
DefaultAssay(MLN_filt) <- "RNA"

MLN_filt[["RNA"]] <- as(MLN_filt[["RNA"]], "Assay5")
#MLN_filt[["ADT"]] <- as(MLN_filt[["ADT"]], "Assay5")
MLN_filt <- RunPCA(MLN_filt, verbose = FALSE, reduction.name = "RNA_pca_cc", features = VariableFeatures(MLN_filt, nfeatures = 3000))
DimPlot(MLN_filt, reduction = "RNA_pca_cc" )
DimPlot(MLN_filt, reduction = "unintegrated_RNA_pca")

MLN_filt[['RNA']] <-  split(x = MLN_filt[['RNA']], f = MLN_filt$InfectionStatus)
DefaultAssay(MLN_filt)
MLN_filt <- IntegrateLayers(object = MLN_filt, method = RPCAIntegration, assay = "RNA", features = VariableFeatures(MLN_filt, nfeatures = 3000),
  orig.reduction = "RNA_pca_cc", new.reduction = "reintegrated_RNA_rpca_cc", verbose = T)
#Saved here today
MLN_filt <- FindNeighbors(MLN_filt, reduction = "reintegrated_RNA_rpca_cc", dims = 1:30)
MLN_filt <- FindClusters(MLN_filt, resolution = 2, cluster.name = "rpca_integration_clusters_cc", group.singletons = F)
MLN_filt <- RunUMAP(MLN_filt, reduction = "reintegrated_RNA_rpca_cc", dims = 1:30, reduction.name = "umap.rpca_cc")
#Transfer the previous phase ID 

q1<- DimPlot(MLN_filt, reduction = "umap.rpca_cc", group.by = c( "Phase"), combine = TRUE, label.size = 4, label = TRUE)+ NoLegend()
q2 <- DimPlot(MLN_filt, reduction = "umap.rpca", group.by = c("Phase"), combine = TRUE, label.size = 4, label = TRUE) + NoLegend()
q1 | q2

MLN_filt = FindMultiModalNeighbors(
  MLN_filt, reduction.list = list("reintegrated_RNA_rpca_cc", "integrated.ADT.rpca"), 
  dims.list = list(1:30, 1:20), 
  modality.weight.name = "RNA.weight", 
  knn.graph.name = "wknn_cc",
  snn.graph.name = "wsnn_cc",
  weighted.nn.name = "weighted.nn_cc",
  verbose = FALSE)

#compare clusters and Umaps before and after CC regression
DefaultAssay(MLN_filt) <- "RNA"
MLN_filt <- RunUMAP(MLN_filt, nn.name = "weighted.nn_cc", reduction.name = "wnn.umap_cc", reduction.key = "wnnUMAPcc_")
MLN_filt <- FindClusters(MLN_filt, graph.name = "wsnn_cc",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters_cc", n.start = 25, group.singletons = F)

q1 <- DimPlot(MLN_filt, reduction = "wnn.umap_cc",   group.by = c( "WNN_clusters_cc"),
  combine = T, label.size = 2, label = TRUE) 
q2 <- DimPlot(MLN_filt, reduction = "wnn.umap",   group.by = c("WNN_clusters"),
  combine = T, label.size = 2, label = TRUE) + NoLegend()
q1 | q2

FeaturePlot(MLN_filt, reduction = "wnn.umap_cc", features = "CD200" )
remove(RNA_only)

#Monocle Clusters
library(SeuratWrappers)
MLN_filt <- JoinLayers(MLN_filt)

 cds <- as.cell_data_set(
  MLN_filt,
  reductions = "wnn.umap_cc",
  default.reduction = "wnn.umap_cc",
  graph =  "wsnn",
  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells( cds, reduction_method = "UMAP", k = 15, cluster_method =  "leiden", resolution = 0.00015,
  num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )
#Visualize clusters Monocle - I don't find the partitions useful and don't believe they make sense  
 plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_filt[["Monocle_clusters_cc"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_filt[["Monocle_partitions_cc"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(MLN_filt, reduction = "wnn.umap_cc",   group.by = c( "WNN_clusters_cc", "Monocle_clusters_cc"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)

library(qs2)

qs_save(MLN_filt, "/path/to/analysis/directory/Seurat_Files/Analysis2_MLN_clustered_annotated.qs2")
```

```{r}

DimPlot(Ileum_filt, reduction = "wnn.umap_cc", group.by = "CoarseCellType", label = T, label.size = 3, combine = T) + NoLegend()
DimPlot(Ileum_filt, reduction = "wnn.umap_cc", group.by = "FineCellType", label = T, label.size = 1, combine = T) + NoLegend()
DimPlot(MLN_filt, reduction = "wnn.umap_cc", group.by = "CoarseCellType", label = T, label.size = 3, combine = T) + NoLegend()
DimPlot(MLN_filt, reduction = "wnn.umap_cc", group.by = "FineCellType", label = T, label.size = 1, combine = T) + NoLegend()

```



#Save the Objects 
```{r}
Ileum_filt[["RNA"]] <-  as(Ileum_filt[["RNA"]], "Assay")
Ileum_filt[["ADT"]] <-  as(Ileum_filt[["ADT"]], "Assay")
SaveH5Seurat(Ileum_filt, filename = "/path/to/analysis/directory/Seurat_Files/Analysis_2_Ileum_Annotated_Filtered.h5Seurat", overwrite = T)

MLN_filt[["RNA"]] <-  as(MLN_filt[["RNA"]], "Assay")
MLN_filt[["ADT"]] <-  as(MLN_filt[["ADT"]], "Assay")
SaveH5Seurat(MLN_filt, filename = "/path/to/analysis/directory/Seurat_Files/Analysis_2_MLN_Annotated_Filtered.h5Seurat", overwrite = T)
```






# Plasma Cell and NVS subsetting 
```{r}

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum) <- "CoarseCellType"
Ileum_Plasma <- subset(x = Ileum, idents = c("Plasma cells"))


Ileum_Plasma$Infection_Organ2 <- paste0(Ileum_Plasma$InfectionStatus, "_", Ileum_Plasma$SampleType)
#Repeat analysis Process
Ileum_Plasma[['RNA']] <-  JoinLayers(Ileum_Plasma[['RNA']])
Ileum_Plasma[['RNA']] <-  split(x = Ileum_Plasma[['RNA']], f = Ileum_Plasma$InfectionStatus) # Not split by organ because basically 0 in IEC
Ileum_Plasma <- FindVariableFeatures(Ileum_Plasma, selection.method = "vst", nfeatures = 500) 
Ileum_Plasma <- RunPCA(Ileum_Plasma, verbose = FALSE, reduction.name = "Plasma_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_Plasma <- IntegrateLayers( object = Ileum_Plasma, method = HarmonyIntegration,  orig.reduction = "Plasma_RNA_pca", new.reduction = "integrated.harmony.Plasma", verbose = FALSE)
DefaultAssay(Ileum_Plasma) <- "ADT"
Ileum_Plasma[['ADT']] <-  JoinLayers(Ileum_Plasma[['ADT']])
Ileum_Plasma[['ADT']] <-  split(x = Ileum_Plasma[['ADT']], f = Ileum_Plasma$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(Ileum_Plasma), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_Plasma)[!rownames(Ileum_Plasma) %in% Isotype_labels] #remove isotypes
Ileum_Plasma <- RunPCA(Ileum_Plasma, features =prots, reduction.name = "Plasma_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_Plasma <- IntegrateLayers( object = Ileum_Plasma, method = HarmonyIntegration, features = prots,  orig.reduction = "Plasma_ADT_pca", new.reduction = "integrated.ADT.harmony.Plasma", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_Plasma) <- "RNA"
Ileum_Plasma = FindMultiModalNeighbors(Ileum_Plasma, reduction.list = list("integrated.harmony.Plasma", "integrated.ADT.harmony.Plasma"),  dims.list = list(1:10, 1:5),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Plasma",  snn.graph.name = "wsnn_Plasma", weighted.nn.name = "weighted_Plasma.nn", verbose = TRUE)

Ileum_Plasma <- RunUMAP(Ileum_Plasma, nn.name = "weighted_Plasma.nn", reduction.name = "wnn.Plasma.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Plasma <- FindClusters(Ileum_Plasma, graph.name = "wsnn_Plasma",  resolution = 2, verbose = FALSE, cluster.name = "Plasma_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Plasma, reduction = "wnn.Plasma.umap",   group.by = c( "Plasma_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Plasma <-JoinLayers(Ileum_Plasma) 
 cds <- as.cell_data_set(Ileum_Plasma, reductions = "wnn.Plasma.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Plasma[["Monocle_Plasma_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Plasma[["Monocle_Plasma_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Plasma, reduction = "wnn.Plasma.umap",   group.by = c( "Plasma_clusters", "Monocle_Plasma_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_Plasma, c( "Cd3e", "Epcam", "rb.prop"), reduction = "wnn.Plasma.umap")
VlnPlot(Ileum_Plasma, "Malat1", group.by = "Plasma_clusters")
FeatureScatter(Ileum_Plasma, "Malat1", "rb.prop")
```


```{r}
RNA_only <- Ileum_Plasma
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Plasma_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Plasma_clusters"] <- eh 
eh <-unlist(RNA_only[["Plasma_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Plasma_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Plasma_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Plasma_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_Plasma, "/path/to/analysis/directory/Seurat_Files/Ileum_Plasma_intermediate.qs2")

```

# Nervous Cell and NVS subsetting 
```{r}

DimPlot(Ileum, reduction = "wnn.umap_cc",   group.by = "CoarseCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum) <- "CoarseCellType"
Ileum_Nervous <- subset(x = Ileum, idents = c("Enteric Nervous System"))

#Repeat analysis Process
Ileum_Nervous[['RNA']] <-  JoinLayers(Ileum_Nervous[['RNA']])
Ileum_Nervous[['RNA']] <-  split(x = Ileum_Nervous[['RNA']], f = Ileum_Nervous$InfectionStatus) # Not split by organ because basically 0 in IEC
Ileum_Nervous <- FindVariableFeatures(Ileum_Nervous, selection.method = "vst", nfeatures = 500) 
Ileum_Nervous <- RunPCA(Ileum_Nervous, verbose = FALSE, reduction.name = "Nervous_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_Nervous <- IntegrateLayers( object = Ileum_Nervous, method = HarmonyIntegration,  orig.reduction = "Nervous_RNA_pca", new.reduction = "integrated.harmony.Nervous", verbose = FALSE)
DefaultAssay(Ileum_Nervous) <- "ADT"
Ileum_Nervous[['ADT']] <-  JoinLayers(Ileum_Nervous[['ADT']])
Ileum_Nervous[['ADT']] <-  split(x = Ileum_Nervous[['ADT']], f = Ileum_Nervous$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(Ileum_Nervous), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_Nervous)[!rownames(Ileum_Nervous) %in% Isotype_labels] #remove isotypes
Ileum_Nervous <- RunPCA(Ileum_Nervous, features =prots, reduction.name = "Nervous_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_Nervous <- IntegrateLayers( object = Ileum_Nervous, method = HarmonyIntegration, features = prots,  orig.reduction = "Nervous_ADT_pca", new.reduction = "integrated.ADT.harmony.Nervous", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_Nervous) <- "RNA"
Ileum_Nervous = FindMultiModalNeighbors(Ileum_Nervous, reduction.list = list("integrated.harmony.Nervous", "integrated.ADT.harmony.Nervous"),  dims.list = list(1:10, 1:5),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_Nervous",  snn.graph.name = "wsnn_Nervous", weighted.nn.name = "weighted_Nervous.nn", verbose = TRUE)

Ileum_Nervous <- RunUMAP(Ileum_Nervous, nn.name = "weighted_Nervous.nn", reduction.name = "wnn.Nervous.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_Nervous <- FindClusters(Ileum_Nervous, graph.name = "wsnn_Nervous",  resolution = 2, verbose = FALSE, cluster.name = "Nervous_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_Nervous, reduction = "wnn.Nervous.umap",   group.by = c( "Nervous_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_Nervous <-JoinLayers(Ileum_Nervous) 
 cds <- as.cell_data_set(Ileum_Nervous, reductions = "wnn.Nervous.umap", default.reduction = "wnn.T.umap",graph =  "wsnn_T",  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.0004, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_Nervous[["Monocle_Nervous_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_Nervous[["Monocle_Nervous_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_Nervous, reduction = "wnn.Nervous.umap",   group.by = c( "Nervous_clusters", "Monocle_Nervous_clusters", "IntermediateCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot_scCustom(Ileum_Nervous, c( "Cd3e", "Epcam", "rb.prop"), reduction = "wnn.Nervous.umap")
VlnPlot(Ileum_Nervous, "Malat1", group.by = "Nervous_clusters")
FeatureScatter(Ileum_Nervous, "Malat1", "rb.prop")
```


```{r}
RNA_only <- Ileum_Nervous
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_Nervous_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_Nervous_clusters"] <- eh 
eh <-unlist(RNA_only[["Nervous_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Nervous_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_Nervous_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_Nervous_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_Nervous, "/path/to/analysis/directory/Seurat_Files/Ileum_Nervous_intermediate.qs2")

```

#Repeat T cells split with scale split #CD4 T cells
```{r}
Ileum_T <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.qs2")
DimPlot(Ileum_T, reduction = "wnn.T.umap",   group.by = "FineCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum_T) <- "FineCellType"
unique(Ileum_T$FineCellType)
#Subset out the correct cells 
Ileum_CD4Tcell <- subset(x = Ileum_T, idents = c("Tfh", "NKT Cells", "CD103 low CD4 Trm", "Sell high CD4 T Cells", 
                                               "Treg", "CD103 high CD4 Trm", "NKT2 Cells", "Sell high Treg", "IFN stimulated CD4 T Cells", "Cycling T Cells"))
Ileum_CD4Tcell[['RNA']] <-  JoinLayers(Ileum_CD4Tcell[['RNA']])
expr <- GetAssayData(Ileum_CD4Tcell, layer = "data")
cycling <- Ileum_CD4Tcell$FineCellType == "Cycling T Cells" 
Cd4_expr <- expr["Cd4", ] >= 0.5 
#logical for which cells to keep
cyclingCD4 <- Cd4_expr  & cycling

noncycle <- Ileum_CD4Tcell$FineCellType != "Cycling T Cells" 
cells_to_keep <- cyclingCD4 | noncycle
Ileum_CD4Tcell <- subset(Ileum_CD4Tcell, cells = colnames(Ileum_CD4Tcell)[cells_to_keep])
#Repeat analysis Process
Ileum_CD4Tcell[['RNA']] <-  JoinLayers(Ileum_CD4Tcell[['RNA']])
Ileum_CD4Tcell <- FindVariableFeatures(Ileum_CD4Tcell, selection.method = "vst", nfeatures = 3000) 
Ileum_CD4Tcell[['RNA']] <-  split(x = Ileum_CD4Tcell[['RNA']], f = Ileum_CD4Tcell$Infection_Organ)
Ileum_CD4Tcell <- RunPCA(Ileum_CD4Tcell, verbose = FALSE, reduction.name = "CD4Tcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_CD4Tcell <- IntegrateLayers( object = Ileum_CD4Tcell, method = HarmonyIntegration,  orig.reduction = "CD4Tcell_RNA_pca", new.reduction = "integrated.harmony.CD4Tcell", verbose = FALSE)
DefaultAssay(Ileum_CD4Tcell) <- "ADT"
Ileum_CD4Tcell[['ADT']] <-  JoinLayers(Ileum_CD4Tcell[['ADT']])
Ileum_CD4Tcell[['ADT']] <-  split(x = Ileum_CD4Tcell[['ADT']], f = Ileum_CD4Tcell$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_CD4Tcell), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_CD4Tcell)[!rownames(Ileum_CD4Tcell) %in% Isotype_labels] #remove isotypes
Ileum_CD4Tcell <- RunPCA(Ileum_CD4Tcell, features =prots, reduction.name = "CD4Tcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_CD4Tcell <- IntegrateLayers( object = Ileum_CD4Tcell, method = HarmonyIntegration, features = prots,  orig.reduction = "CD4Tcell_ADT_pca", new.reduction = "integrated.ADT.harmony.CD4Tcell", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_CD4Tcell) <- "RNA"
Ileum_CD4Tcell = FindMultiModalNeighbors(Ileum_CD4Tcell, reduction.list = list("integrated.harmony.CD4Tcell", "integrated.ADT.harmony.CD4Tcell"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_CD4Tcell",  snn.graph.name = "wsnn_CD4Tcell", weighted.nn.name = "weighted_CD4Tcell.nn", verbose = TRUE)

Ileum_CD4Tcell <- RunUMAP(Ileum_CD4Tcell, nn.name = "weighted_CD4Tcell.nn", reduction.name = "wnn.CD4Tcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_CD4Tcell <- FindClusters(Ileum_CD4Tcell, graph.name = "wsnn_CD4Tcell",  resolution = 2, verbose = FALSE, cluster.name = "CD4Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_CD4Tcell, reduction = "wnn.CD4Tcell.umap",   group.by = c( "CD4Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_CD4Tcell[["RNA"]] <-JoinLayers(Ileum_CD4Tcell[["RNA"]]) 
cds <- SeuratWrappers::as.cell_data_set(Ileum_CD4Tcell, reductions = "wnn.CD4Tcell.umap", 
                                        default.reduction = "wnn.CD4Tcell.umap", graph = "wsnn_CD4Tcell",
                                        group.by = NULL, assay = "RNA")
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_CD4Tcell[["Monocle_CD4Tcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_CD4Tcell[["Monocle_CD4Tcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_CD4Tcell, reduction = "wnn.CD4Tcell.umap",   group.by = c( "CD4Tcell_clusters", "Monocle_CD4Tcell_clusters", "FineCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
```
```{r}
RNA_only <- Ileum_CD4Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD4Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD4Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_CD4Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_CD4Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_CD4Tcell, "/path/to/analysis/directory/Seurat_Files/Ileum_CD4Tcell_intermediate.qs2")

```

```{r}
#Finest Cell Type 
Ileum_CD4Tcell$FinestCellType <- "NA"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("36")] <- "Follicular Treg"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("52", "50")] <- "Ccr2 high Treg"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("30", "7")] <- "Klrg1 high Treg"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("4", "40")] <- "Il12rb2 high CD4 Effector"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("21")] <- "Il17rb high CD4 Effector"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("25")] <- "Rorc high CD4 Effector"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("11", "63", "20", "27", "56")] <- "Itgae high CD4 Effector"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("47", "13", "53", "54")] <- "S100a6 high CD4 Effector"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("46", "33")] <- "NKT"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("66")] <- "Il17rb high NKT"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("3")] <- "Il12rb2 high NKT"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("34", "62", "10", "37", "67")] <- "Tfh"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("42")] <- "Cycling CD4 T cells"
Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$Monocle_CD4Tcell_clusters %in% c("12")] <- "Sell high Treg"
DimPlot(Ileum_CD4Tcell, reduction = "wnn.CD4Tcell.umap",   group.by = c( "FinestCellType"),
  combine = TRUE, label.size = 2, label = TRUE, raster = F) +NoLegend()


Ileum_CD4Tcell$FinestCellType[Ileum_CD4Tcell$FinestCellType  == "NA"] <- "Sell high CD4 T cell"
DimPlot(Ileum_CD4Tcell, reduction = "wnn.CD4Tcell.umap",   group.by = c( "FinestCellType"),
  combine = TRUE, label.size = 2, label = TRUE, raster = F) +NoLegend()
```

```{r}
RNA_only <- Ileum_CD4Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD4Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD4Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD4Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_CD4Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_CD4Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_CD4Tcell, "/path/to/analysis/directory/Seurat_Files/Ileum_CD4Tcell_intermediate.qs2")

```


#CD8 T cell subclustering Ileum 
```{r}
#Ileum_T <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.qs2")
DimPlot(Ileum_T, reduction = "wnn.T.umap",   group.by = "FineCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum_T) <- "FineCellType"
unique(Ileum_T$FineCellType)
#Subset out the correct cells 
Ileum_CD8Tcell <- subset(x = Ileum_T, idents = c( "Cycling T Cells", "Sell high CD8 T Cells", "Gzmb gd T Cells", "Cd160 gd T Cells", "gd T17 Cells", "CD103 high CD8 Trm", "Sell high gd T Cells", "Ly6c gd T Cells", "Ly6c high CD8 Trm", "IFN Stimulated CD8 T Cells"))
Ileum_CD8Tcell[['RNA']] <-  JoinLayers(Ileum_CD8Tcell[['RNA']])
expr <- GetAssayData(Ileum_CD8Tcell, layer = "data")
cycling <- Ileum_CD8Tcell$FineCellType == "Cycling T Cells" 
Cd4_expr <- expr["Cd4", ] < 0.5 
#logical for which cells to keep
cyclingCD4 <- Cd4_expr  & cycling

noncycle <- Ileum_CD8Tcell$FineCellType != "Cycling T Cells" 
cells_to_keep <- cyclingCD4 | noncycle
Ileum_CD8Tcell <- subset(Ileum_CD8Tcell, cells = colnames(Ileum_CD8Tcell)[cells_to_keep])
#Repeat analysis Process
Ileum_CD8Tcell[['RNA']] <-  JoinLayers(Ileum_CD8Tcell[['RNA']])
Ileum_CD8Tcell <- FindVariableFeatures(Ileum_CD8Tcell, selection.method = "vst", nfeatures = 3000) 
Ileum_CD8Tcell[['RNA']] <-  split(x = Ileum_CD8Tcell[['RNA']], f = Ileum_CD8Tcell$Infection_Organ)
Ileum_CD8Tcell <- RunPCA(Ileum_CD8Tcell, verbose = FALSE, reduction.name = "CD8Tcell_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_CD8Tcell <- IntegrateLayers( object = Ileum_CD8Tcell, method = HarmonyIntegration,  orig.reduction = "CD8Tcell_RNA_pca", new.reduction = "integrated.harmony.CD8Tcell", verbose = FALSE)
DefaultAssay(Ileum_CD8Tcell) <- "ADT"
Ileum_CD8Tcell[['ADT']] <-  JoinLayers(Ileum_CD8Tcell[['ADT']])
Ileum_CD8Tcell[['ADT']] <-  split(x = Ileum_CD8Tcell[['ADT']], f = Ileum_CD8Tcell$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_CD8Tcell), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_CD8Tcell)[!rownames(Ileum_CD8Tcell) %in% Isotype_labels] #remove isotypes
Ileum_CD8Tcell <- RunPCA(Ileum_CD8Tcell, features =prots, reduction.name = "CD8Tcell_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_CD8Tcell <- IntegrateLayers( object = Ileum_CD8Tcell, method = HarmonyIntegration, features = prots,  orig.reduction = "CD8Tcell_ADT_pca", new.reduction = "integrated.ADT.harmony.CD8Tcell", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_CD8Tcell) <- "RNA"
Ileum_CD8Tcell = FindMultiModalNeighbors(Ileum_CD8Tcell, reduction.list = list("integrated.harmony.CD8Tcell", "integrated.ADT.harmony.CD8Tcell"),  dims.list = list(1:25, 1:15),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_CD8Tcell",  snn.graph.name = "wsnn_CD8Tcell", weighted.nn.name = "weighted_CD8Tcell.nn", verbose = TRUE)

Ileum_CD8Tcell <- RunUMAP(Ileum_CD8Tcell, nn.name = "weighted_CD8Tcell.nn", reduction.name = "wnn.CD8Tcell.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_CD8Tcell <- FindClusters(Ileum_CD8Tcell, graph.name = "wsnn_CD8Tcell",  resolution = 2, verbose = FALSE, cluster.name = "CD8Tcell_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_CD8Tcell, reduction = "wnn.CD8Tcell.umap",   group.by = c( "CD8Tcell_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_CD8Tcell[["RNA"]] <-JoinLayers(Ileum_CD8Tcell[["RNA"]]) 
cds <- SeuratWrappers::as.cell_data_set(Ileum_CD8Tcell, reductions = "wnn.CD8Tcell.umap", 
                                        default.reduction = "wnn.CD8Tcell.umap", graph = "wsnn_CD8Tcell",
                                        group.by = NULL, assay = "RNA")
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.001, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_CD8Tcell[["Monocle_CD8Tcell_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_CD8Tcell[["Monocle_CD8Tcell_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_CD8Tcell, reduction = "wnn.CD8Tcell.umap",   group.by = c( "CD8Tcell_clusters", "Monocle_CD8Tcell_clusters", "FineCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
```
```{r}
RNA_only <- Ileum_CD8Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD8Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD8Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_CD8Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_CD8Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_CD8Tcell, "/path/to/analysis/directory/Seurat_Files/Ileum_CD8Tcell_intermediate.qs2")

```

```{r}
#Finest Cell Type 
# Potentially 2 defined popultions of CD8aa TCRab T cells  - PMID: 32687575 PMID: 28530714. These papers tells us that there is a PD1+ groups which is self reactive. They express CD5, Nr4a1 and other genes. These stand in contrast to the other group which were defined as PD1 negative but found to express Itgae, Cxcr3, 
#gd T cells are complex Frequenly, they are clustered by vg type. for example vg5 are renowned as DETC in the skin. In my data, there is some association of trdv and trgv regions with clustering. A large portion of mostly gzm+ gd t cells appear to be trvg7+ which matches some recent literature PMID: 27641500. These author found the remaining to be mostly vg1 or vg4. In our data, the trgv7- gd T cells are usually trgv1 as well. They report that trvg7+ cells are high for il2rb (CD122) which is also found in my data. Notably, trgv6 is most common among the gd T17 cells in my data. 
#The gd T cells also separate very nicely by cd160 and Gzma/Gzmb. Interesting, a subpopulation of both of these popylations is Themis +. In the Gzm+ gd T cells, this Themis + group is also gzmk+ and other cells are not. What is granzyme k? PMID: 38826230. The Themis+ Cd160+ gd T cells are more of a mystery but are clearly dominated by trvg1 trvg2 and are trvg7-.  Both Themis + groups are thy1+. I wonder if they stem from developmentally distinct populations. 
# the results from the above study tate - Here, we identify a time window early in the development of young mice in which Btnl1 expressed by post-mitotic, small intestinal villus epithelial cells critically and selectively promotes the maturation and expansion of Vg7+ T cells, thereby shaping the IEL compartment. Requiring neither microbial nor food antigens, this process evokes Skint1-mediated DETC selection and ab T cell selection by the MHC



Ileum_CD8Tcell$FinestCellType <- "NA"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("13", "22", "25")] <- "Themis positive Gzma gd T cells"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("42", "33", "43", "50", "7", "14", "46", "16", "51", 
                                                                              "8", "56", "9", "60", "61", "18", "54", "45", "19", "26", "58", "35")] <- "Gzma gd T cells"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("4", "10", "39", "31", "32", "40", "36")] <- "Themis positive Cd160 gd T cells"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("30", "15", "53", "44",  "38", "48", "29", "28", "37", "59", "41")] <- "Cd160 gd T cells"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("3")] <- "Cycling CD8 T cells"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("2")] <- "gd T17"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("20")] <- "Sox4 positive gd T cells"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("23", "34", "57")] <- "Itgae positives CD8aa ab T cells"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("52")] <- "DN ab T cells"

Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("47")] <- "Ikzf2 high CD8 Effector"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("12")] <- "S100a6 high CD8 Effector"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("5", "27")] <- "CD8 Effector"
Ileum_CD8Tcell$FinestCellType[Ileum_CD8Tcell$Monocle_CD8Tcell_clusters %in% c("1", "11", "17", "21", "6", "24", "55", "49")] <- "Sell high CD8 T cells"

DimPlot(Ileum_CD8Tcell, reduction = "wnn.CD8Tcell.umap",   group.by = c( "FinestCellType"),
  combine = TRUE, label.size = 2, label = TRUE, raster = F) +NoLegend()
```

```{r}
RNA_only <- Ileum_CD8Tcell
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_CD8Tcell_clusters"] <- eh 
eh <-unlist(RNA_only[["CD8Tcell_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["CD8Tcell_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_CD8Tcell_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_CD8Tcell_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_CD8Tcell, "/path/to/analysis/directory/Seurat_Files/Ileum_CD8Tcell_intermediate.qs2")

```

#ILCS

```{r}
#Ileum_T <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.qs2")
DimPlot(Ileum_T, reduction = "wnn.T.umap",   group.by = "FineCellType", label.size = 3, label = TRUE, repel = F, raster = FALSE)
Idents(Ileum_T) <- "FineCellType"
unique(Ileum_T$FineCellType)
#Subset out the correct cells 
Ileum_ILC <- subset(x = Ileum_T, idents = c("ILC2s", "LTIs", "ILC3s", "ILC1s", "Cycling ILCS", "NK Cells"))

#Repeat analysis Process
Ileum_ILC[['RNA']] <-  JoinLayers(Ileum_ILC[['RNA']])
Ileum_ILC <- FindVariableFeatures(Ileum_ILC, selection.method = "vst", nfeatures = 3000) 
Ileum_ILC[['RNA']] <-  split(x = Ileum_ILC[['RNA']], f = Ileum_ILC$Infection_Organ)
Ileum_ILC <- RunPCA(Ileum_ILC, verbose = FALSE, reduction.name = "ILC_RNA_pca", assay = "RNA" )
options(future.globals.maxSize = 20000 * 1024^2)
Ileum_ILC <- IntegrateLayers( object = Ileum_ILC, method = HarmonyIntegration,  orig.reduction = "ILC_RNA_pca", new.reduction = "integrated.harmony.ILC", verbose = FALSE)
DefaultAssay(Ileum_ILC) <- "ADT"
Ileum_ILC[['ADT']] <-  JoinLayers(Ileum_ILC[['ADT']])
Ileum_ILC[['ADT']] <-  split(x = Ileum_ILC[['ADT']], f = Ileum_ILC$Infection_Organ)
Isotype_labels <- grep("isotype", rownames(Ileum_ILC), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_ILC)[!rownames(Ileum_ILC) %in% Isotype_labels] #remove isotypes
Ileum_ILC <- RunPCA(Ileum_ILC, features =prots, reduction.name = "ILC_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
Ileum_ILC <- IntegrateLayers( object = Ileum_ILC, method = HarmonyIntegration, features = prots,  orig.reduction = "ILC_ADT_pca", new.reduction = "integrated.ADT.harmony.ILC", verbose = TRUE)
#Recluster
DefaultAssay(Ileum_ILC) <- "RNA"
Ileum_ILC = FindMultiModalNeighbors(Ileum_ILC, reduction.list = list("integrated.harmony.ILC", "integrated.ADT.harmony.ILC"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, knn.graph.name = "wknn_ILC",  snn.graph.name = "wsnn_ILC", weighted.nn.name = "weighted_ILC.nn", verbose = TRUE)

Ileum_ILC <- RunUMAP(Ileum_ILC, nn.name = "weighted_ILC.nn", reduction.name = "wnn.ILC.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
Ileum_ILC <- FindClusters(Ileum_ILC, graph.name = "wsnn_ILC",  resolution = 2, verbose = FALSE, cluster.name = "ILC_clusters", n.start = 20, group.singletons = F)
DimPlot(Ileum_ILC, reduction = "wnn.ILC.umap",   group.by = c( "ILC_clusters", "InfectionStatus", "Phase", "IntermediateCellType"), combine = FALSE, label.size = 2, label = TRUE)

#Monocle Clustering
Ileum_ILC[["RNA"]] <-JoinLayers(Ileum_ILC[["RNA"]]) 
cds <- SeuratWrappers::as.cell_data_set(Ileum_ILC, reductions = "wnn.ILC.umap", 
                                        default.reduction = "wnn.ILC.umap", graph = "wsnn_ILC",
                                        group.by = NULL, assay = "RNA")
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(  cds, reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.005, num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_ILC[["Monocle_ILC_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_ILC[["Monocle_ILC_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_ILC, reduction = "wnn.ILC.umap",   group.by = c( "ILC_clusters", "Monocle_ILC_clusters", "FineCellType", "Phase", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
```
```{r}
RNA_only <- Ileum_ILC
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_ILC_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_ILC_clusters"] <- eh 
eh <-unlist(RNA_only[["ILC_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["ILC_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_ILC_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_ILC_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_ILC, "/path/to/analysis/directory/Seurat_Files/Ileum_ILC_intermediate.qs2")

```

```{r}
#Finest Cell Type 
# Information about ILC2 - PMID: 35354980
Ileum_ILC$FinestCellType <- "NA"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("86", "5", "48")] <- "ILC1s"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("87")] <- "T cells"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("11", "10")] <- "NK cells"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("99")] <- "Cycling NK cells"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("12", "65", "85", "60", "42", "96", "81", "2", "7", "44", "19", "56")] <- "LTIs"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("14", "93", "79", "45", "31", "51", "58", "38", "21", "82", "23", "46", "84")] <- "ILC3s"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("101", "79")] <- "Cycling ILC3s"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("32", "78", "9", "47", "16", "73", "90", "71", "77", "8", "50", "94", "98", "27", "74", "17", "20", "25")] <- "Ccr9 negative resident ILC2s"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("100" , "34", "89", "92", "22", "70", "66", "41", "30")] <- "Cycling ILC2s"
Ileum_ILC$FinestCellType[Ileum_ILC$Monocle_ILC_clusters %in% c("68")] <- "Il1rl1 positive resident ILC2s"
DimPlot(Ileum_ILC, reduction = "wnn.ILC.umap",   group.by = c( "FinestCellType"),
  combine = TRUE, label.size = 2, label = TRUE, raster = F) +NoLegend()

Ileum_ILC$FinestCellType[Ileum_ILC$FinestCellType == "NA"] <- "Ccr9 positive migratory ILC2s"

DimPlot(Ileum_ILC, reduction = "wnn.ILC.umap",   group.by = c( "FinestCellType"),
  combine = TRUE, label.size = 2, label = TRUE, raster = F) +NoLegend()

```

```{r}
RNA_only <- Ileum_ILC
RNA_only[["ADT"]] <- NULL
eh <-unlist(RNA_only[["Monocle_ILC_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["Monocle_ILC_clusters"] <- eh 
eh <-unlist(RNA_only[["ILC_clusters"]])
eh <- as.character(eh) 
RNA_only@meta.data["ILC_clusters"] <- eh 
RNA_only[["RNA"]] <- as(RNA_only[["RNA"]], "Assay")
RNA_only[["RNA"]]$scale.data <- NULL
SaveH5Seurat(RNA_only, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_ILC_intermediate.h5Seurat", overwrite = T)
Convert("/path/to/analysis/directory/Seurat_Files/Ileum_ILC_intermediate.h5Seurat",   assay = "RNA", dest = "h5ad", overwrite = TRUE)
remove(RNA_only)
library(qs2)
qs_save(Ileum_ILC, "/path/to/analysis/directory/Seurat_Files/Ileum_ILC_intermediate.qs2")

```

#Transfer Finest Label to Ileum
```{r}
Ileum_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2_clustered_annotated_Ileum.qs2")

Ileum_filt$FinestCellType <- Ileum_filt$FineCellType
Ileum_filt <- AddMetaData(object = Ileum_filt, metadata = Ileum_CD4Tcell_filt$FinestCellType , col.name = "FinestCellType")
Ileum_filt <- AddMetaData(object = Ileum_filt, metadata = Ileum_CD8Tcell_filt$FinestCellType , col.name = "FinestCellType")
Ileum_filt <- AddMetaData(object = Ileum_filt, metadata = Ileum_ILCs_filt$FinestCellType , col.name = "FinestCellType")
Ileum_filt <- AddMetaData(object = Ileum_filt, metadata = Ileum_Epithelial_filt$FinestCellType , col.name = "FinestCellType")

DimPlot(Ileum_filt, reduction = "wnn.umap_cc",   group.by = c( "FinestCellType"),
  combine = TRUE, label.size = 2, label = TRUE, raster = F) +NoLegend()
unique(Ileum_filt$FinestCellType)

qs_save(Ileum_filt, "/path/to/analysis/directory/Seurat_Files/Analysis2_clustered_annotated_Ileum.qs2")
```

