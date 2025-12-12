#Description: Here I describe how the Xenium data was processed including resegmentation of the cells via Proseg and then analyzed. 
#Figures produced Figure 7A, Figure S7C, Figure S7D, Figure 7B, Figure 7C, Figure 7D, figure 7E, #Figure 7F, Figure 7G, Figure 7H, Figure 7I Figure 7J

#After Xenium machine was finalized, the raw data was then used through Proseg (PMID: 40404994). Proseg is a probabilistic cell segmentation algorithm for spatial data. The Proseg to Xenium Ranger results were utilized here

# Load Libraries----

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
library(RColorBrewer)
library(monocle3)
library(SingleR)
library(gridExtra)
library(gt)
library(arrow) #This is a package used in the custom ReadXenium function needed after10X changed their output nomenclature
library(remotes)
options(future.globals.maxSize = 20000 * 1024^2)
library(FNN)
library(ClusterR)   # For MiniBatchKmeans
library(tibble)
library(qs2)
library(ggrepel)

library(sf) #Used by functions in Seurat to process some types of Xenium data

#Create Directories----
getwd() #/home/hartandrew
setwd("/path/to/data/MIST/Xenium_Yp_2025/")
out <- "/path/to/data/MIST/Xenium_Yp_2025/Output/"
seurat <- "/path/to/data/MIST/Xenium_Yp_2025/Seurat/"
images <- "/path/to/data/MIST/Xenium_Yp_2025/Images/"
CSV <- "/path/to/data/MIST/Xenium_Yp_2025/CSV/"

# Redefine ReadXenium() and LoadXenium functions  - 10X redefined their output directory names and it is no longer compatible with ----
ReadXenium_New <- function (data.dir, outs = c("matrix", "microns"), type = "centroids", 
                            mols.qv.threshold = 20) 
{
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location, 
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}

LoadXenium_New <- function(data.dir, fov = 'fov', assay = 'Xenium') {
  data <- ReadXenium_New(
    data.dir = data.dir,
    type = c("centroids", "segmentations"),
  )
  
  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )
  
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  if("Blank Codeword" %in% names(data$matrix))
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  else
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}



#Load the Objects----
naive1_path <- "/data/beitinglab/Proseg/Naive1/Naive1_Xenium/outs/"
naive2_path <- "/data/beitinglab/Proseg/Naive2/Naive2_Xenium/outs/"
yersinia1_path <- "/data/beitinglab/Proseg/Yersinia1/Yersinia1_Xenium/outs/"
yersinia2_path <- "/data/beitinglab/Proseg/Yersinia2/Yersinia2_Xenium/outs/"

naive1_path
naive2_path
yersinia1_path 
yersinia2_path 



#Remove cells with 0 transcripts
naive1.obj <- LoadXenium(naive1_path, fov = "naive1")
naive1.obj <- subset(naive1.obj, subset = nCount_Xenium > 0)

naive2.obj <- LoadXenium(naive2_path, fov = "naive2")
naive2.obj <- subset(naive2.obj, subset = nCount_Xenium > 0)

yersinia1.obj <- LoadXenium(yersinia1_path, fov = "yersinia1")
yersinia1.obj <- subset(yersinia1.obj, subset = nCount_Xenium > 0)

yersinia2.obj <- LoadXenium(yersinia2_path, fov = "yersinia2")
yersinia2.obj <- subset(yersinia2.obj, subset = nCount_Xenium > 0)

# Add Metadata the Data ----

naive1.obj@meta.data["orig.ident"] <- "naive1"
naive1.obj@meta.data["InfectionStatus"] <- "Naive"
naive1.obj@meta.data["Method"] <- "Proseg"
naive2.obj@meta.data["orig.ident"] <- "naive2"
naive2.obj@meta.data["InfectionStatus"] <- "Naive"
naive2.obj@meta.data["Method"] <- "Proseg"

yersinia1.obj@meta.data["orig.ident"] <- "yersinia1"
yersinia1.obj@meta.data["InfectionStatus"] <- "Yersinia"
yersinia1.obj@meta.data["Method"] <- "Proseg"
yersinia2.obj@meta.data["orig.ident"] <- "yersinia2"
yersinia2.obj@meta.data["InfectionStatus"] <- "Yersinia"
yersinia2.obj@meta.data["Method"] <- "Proseg"


#Qc Before Filtering ----
seurat_objects <- list(
  naive1 = naive1.obj,
  naive2 = naive2.obj,
  yersinia1 = yersinia1.obj,
  yersinia2 = yersinia2.obj
)

# Get number of cells per object
cell_counts <- data.frame(
  Object = names(seurat_objects),
  CellCount = sapply(seurat_objects, ncol)
)

# Bar plot of cell counts
cell_count_barplot <- ggplot(cell_counts, aes(x = Object, y = CellCount, fill = Object)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = CellCount), vjust = -0.3) +
  labs(
    title = "Number of Cells per Sample",
    y = "Number of Cells",
    fill = "Sample"
  ) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), legend.position = "none")
cell_count_barplot
# Save barplot
ggsave(filename = "Unfiltered_cell_counts_barplot.png", plot = cell_count_barplot, path = images, width = 6, height = 5, dpi = 600)
ggsave(filename = "Unfiltered_cell_counts_barplot.svg", plot = cell_count_barplot, path = images, width = 6, height = 5, dpi = 600)

# Combine all metadata for joint violin plots
combined_meta <- do.call(rbind, lapply(names(seurat_objects), function(name) {
  obj <- seurat_objects[[name]]
  meta <- obj@meta.data
  meta$sample <- name
  meta
}))

# Calculate medians for labeling
median_counts <- aggregate(nCount_Xenium ~ sample, data = combined_meta, median)
median_features <- aggregate(nFeature_Xenium ~ sample, data = combined_meta, median)

# Plot nCount per cell
p_ncount <- ggplot(combined_meta, aes(x = sample, y = nCount_Xenium)) +
  geom_violin(fill = "grey", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_text(data = median_counts, aes(x = sample, y = nCount_Xenium, label = paste("Median:", round(nCount_Xenium, 1))),
            inherit.aes = FALSE, vjust = -0.7, size = 4) +
  labs(title = "Transcripts per Cell", y = "nCount_Xenium", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))
p_ncount
ggsave("Unfiltered_Transcripts_perCell_Violin.png", plot = p_ncount, path = images, width = 8, height = 5, dpi = 600)
ggsave("Unfiltered_Transcripts_perCell_Violin.svg", plot = p_ncount, path = images, width = 8, height = 5, dpi = 600)

# Plot nFeature per cell
p_nfeature <- ggplot(combined_meta, aes(x = sample, y = nFeature_Xenium)) +
  geom_violin(fill = "grey", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_text(data = median_features, aes(x = sample, y = nFeature_Xenium, label = paste("Median:", round(nFeature_Xenium, 1))),
            inherit.aes = FALSE, vjust = -0.7, size = 4) +
  labs(title = "Features per Cell", y = "nFeature_Xenium", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))
p_nfeature
ggsave("Unfiltered_Feature_perCell_Violin.png", plot = p_nfeature, path = images, width = 8, height = 5, dpi = 600)
ggsave("Unfiltered_Feature_perCell_Violin.svg", plot = p_nfeature, path = images, width = 8, height = 5, dpi = 600)


seurat_objects_unfiltered <- list(
  naive1 = naive1.obj,
  naive2 = naive2.obj,
  yersinia1 = yersinia1.obj,
  yersinia2 = yersinia2.obj
)


# Filter Each object by MAD ----
filter_by_mad <- function(seurat_object, feature, n_mad = 2.5) {
  # Calculate median and MAD
  feature_values <- seurat_object[[feature]][[1]]
  feature_median <- median(feature_values)
  feature_mad <- mad(feature_values)
  # Define lower and upper bounds
  lower_bound <- feature_median - n_mad * feature_mad
  upper_bound <- feature_median + n_mad * feature_mad
  # Filter cells
  cells_to_keep <- Cells(seurat_object)[feature_values >= lower_bound & feature_values <= upper_bound]
  return(subset(seurat_object, cells = cells_to_keep))
}

# Apply filtering to each object

naive1.obj <- filter_by_mad(naive1.obj, 'nCount_Xenium')
naive1.obj <- filter_by_mad(naive1.obj, 'nFeature_Xenium')
naive2.obj <- filter_by_mad(naive2.obj, 'nCount_Xenium')
naive2.obj <- filter_by_mad(naive2.obj, 'nFeature_Xenium')

yersinia1.obj <- filter_by_mad(yersinia1.obj, 'nCount_Xenium')
yersinia1.obj <- filter_by_mad(yersinia1.obj, 'nFeature_Xenium')
yersinia2.obj <- filter_by_mad(yersinia2.obj, 'nCount_Xenium')
yersinia2.obj <- filter_by_mad(yersinia2.obj, 'nFeature_Xenium')


#Repeat QC After filtering 
# List of Seurat objects (assumed already loaded)
seurat_objects <- list(
  naive1 = naive1.obj,
  naive2 = naive2.obj,
  yersinia1 = yersinia1.obj,
  yersinia2 = yersinia2.obj
)

# Get number of cells per object
cell_counts <- data.frame(
  Object = names(seurat_objects),
  CellCount = sapply(seurat_objects, ncol)
)

# Bar plot of cell counts
cell_count_barplot <- ggplot(cell_counts, aes(x = Object, y = CellCount, fill = Object)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = CellCount), vjust = -0.3) +
  labs(
    title = "Number of Cells per Sample",
    y = "Number of Cells",
    fill = "Sample"
  ) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), legend.position = "none")
cell_count_barplot
# Save barplot
ggsave(filename = "Filtered_cell_counts_barplot.png", plot = cell_count_barplot, path = images, width = 6, height = 5, dpi = 600)
ggsave(filename = "Filtered_cell_counts_barplot.svg", plot = cell_count_barplot, path = images, width = 6, height = 5, dpi = 600)

# Combine all metadata for joint violin plots
combined_meta <- do.call(rbind, lapply(names(seurat_objects), function(name) {
  obj <- seurat_objects[[name]]
  meta <- obj@meta.data
  meta$sample <- name
  meta
}))

# Calculate medians for labeling
median_counts <- aggregate(nCount_Xenium ~ sample, data = combined_meta, median)
median_features <- aggregate(nFeature_Xenium ~ sample, data = combined_meta, median)

# Plot nCount per cell
p_ncount <- ggplot(combined_meta, aes(x = sample, y = nCount_Xenium)) +
  geom_violin(fill = "grey", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_text(data = median_counts, aes(x = sample, y = nCount_Xenium, label = paste("Median:", round(nCount_Xenium, 1))),
            inherit.aes = FALSE, vjust = -0.7, size = 4) +
  #stat_summary(fun = median, geom = "hline", aes(yintercept = ..y..), linetype = "dashed", color = "black") +
  labs(title = "Transcripts per Cell", y = "nCount_Xenium", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))
p_ncount
ggsave("Filtered_Transcripts_perCell_Violin.png", plot = p_ncount, path = images, width = 8, height = 5, dpi = 600)
ggsave("Filtered_Transcripts_perCell_Violin.svg", plot = p_ncount, path = images, width = 8, height = 5, dpi = 600)

# Plot nFeature per cell
p_nfeature <- ggplot(combined_meta, aes(x = sample, y = nFeature_Xenium)) +
  geom_violin(fill = "grey", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_text(data = median_features, aes(x = sample, y = nFeature_Xenium, label = paste("Median:", round(nFeature_Xenium, 1))),
            inherit.aes = FALSE, vjust = -0.7, size = 4) +
  #stat_summary(fun = median, geom = "hline", aes(yintercept = ..y..), linetype = "dashed", color = "black") +
  labs(title = "Features per Cell", y = "nFeature_Xenium", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))
p_nfeature
ggsave("Filtered_Feature_perCell_Violin.png", plot = p_nfeature, path = images, width = 8, height = 5, dpi = 600)
ggsave("Filtered_Feature_perCell_Violin.svg", plot = p_nfeature, path = images, width = 8, height = 5, dpi = 600)

#Merge----
all_objects_filtered <- merge(naive1.obj, y = c(naive2.obj,  yersinia1.obj, yersinia2.obj),
                     add.cell.ids = c("naive1", "naive2", "yersinia1", "yersinia2"), project = "Combined", merge.data = TRUE)
all_objects_filtered[["Xenium"]] <- JoinLayers(all_objects_filtered[["Xenium"]])
all_objects_filtered[["Xenium"]] <- split(x = all_objects_filtered[["Xenium"]],
                                         f = all_objects_filtered$orig.ident)

# Save the complete object with all fov and remove from environment
library(qs2)
qs_save(all_objects_filtered, "/path/to/data/MIST/Xenium_Yp_2025/Seurat/Seurat_all_proseg_fov.qs2")
#Filter by Gene and Transcript Counts

ncol(all_objects_filtered) 
all_objects_filtered <- subset(all_objects_filtered, subset = nFeature_Xenium > 50 &
                                nCount_Xenium > 100)
ncol(all_objects_filtered) 
# SCT Transformation and Clustering ----
#SCT Transform takes ~10 minutes per sample to run 

all_objects_filtered <- SCTransform(all_objects_filtered, assay = "Xenium",  variable.features.n = 1500 ) # Chooses 3000 Variable Features by default
all_objects_filtered
length(VariableFeatures(all_objects_filtered)) 
#VariableFeaturePlot(all_objects_filtered, assay = "SCT")
z <- all_objects_filtered@assays$SCT@SCTModel.list[[1]]@feature.attributes
qplot(z$gmean, z$residual_variance) +
  scale_y_log10() +
  scale_x_log10()

# RunPCA
all_objects_filtered <- RunPCA(all_objects_filtered, npcs = 30, features = rownames(all_objects_filtered)) 
ElbowPlot(all_objects_filtered, ndims = 30, reduction = "pca")
all_objects_filtered <- RunUMAP(all_objects_filtered, dims = 1:30, n.neighbors = 100L, min.dist = 0.2) 
all_objects_filtered <- FindNeighbors(all_objects_filtered, reduction = "pca", dims = 1:30)
all_objects_filtered <- FindClusters(all_objects_filtered, resolution = 1.5,  n.start = 30) 

save.image("/path/to/data/MIST/Xenium_Yp_2025/Yersinia_analysis.RData")


Idents(all_objects_filtered) <- "seurat_clusters"

#Monocle Clustering ----
library(SeuratWrappers)
cds <- as.cell_data_set(
  all_objects_filtered,
  reductions = "umap",
  default.reduction = "umap",
  graph =  "SCT_snn",
  group.by = NULL)

reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(
  cds,
  reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", 
  resolution = 0.00001,
  num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )


#Visualize clusters Monocle - I don't find the partitions useful and don't believe they make sense  
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1|p2 

Seurat_monocle <- as.Seurat(cds, assay = NULL)
all_objects_filtered[["Monocle_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
all_objects_filtered[["Monocle_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(all_objects_filtered, reduction = "umap",   group.by = c(  "Monocle_clusters"),
        combine = T, label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(all_objects_filtered, reduction = "umap",   group.by = c(  "seurat_clusters"),
        combine = T, label.size = 2, label = TRUE, raster = F) + NoLegend()
ggsave( "Umap_SeuratClusters_AllData.png",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)


#Coarse Annotations SingleR----
all_objects_filtered <- qs2::qs_read("/path/to/data/MIST/Xenium_Yp_2025/Seurat/Seurat_all_proseg_fov.qs2")
library(BiocParallel)
Ileum_annotated <-qs2::qs_read( "/path/to/analysis/directory/Seurat_Files/Analysis2.5_Ileum.qs2")
Idents(Ileum_annotated) <- "CoarseCellType"
unique(Ileum_annotated$CoarseCellType)

cell_types <- unique(Ileum_annotated$CoarseCellType)

# Find markers for each cell type vs. all others
markers <- FindAllMarkers(
  Ileum_annotated,
  only.pos = TRUE,         
  min.pct = 0.25,          
  logfc.threshold = 0.25   
)

#Filter out genes that are also expression in some f the 
markers <- markers %>%
  group_by(cluster) %>%
  filter(pct.1 >= 0.25 & pct.2 <= 0.05) %>%  # expressed in ≥25% of group cells, ≤5% in others
  arrange(cluster, desc(avg_log2FC))



# Save to file if needed
write.csv(markers, paste0(out, "/CoarseCellType_SpecificGenes.csv"), row.names = FALSE)

shared_markers <- markers[markers$gene %in% rownames(all_objects_filtered@assays$SCT),]
#Test the SingleR approach
scRNA_counts <- GetAssayData(Ileum_annotated, assay = "RNA", layer = "data")
scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(all_objects_filtered@assays$SCT), ]
FineCellType_vector <- as.character(Ileum_annotated$CoarseCellType)
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)
head(reference$counts)
head(reference$labels)
all_objects_filtered <- JoinLayers(all_objects_filtered, assay = "Xenium")

#Loop over samples to assign SingleR labels
originals <-c("naive1", "naive2", "yersinia1", "yersinia2")
singleR_results_list <- list()
xenium_counts <- GetAssayData(all_objects_filtered, assay = "SCT", layer = "data")[, Cells(all_objects_filtered[["naive1"]])]

library(SingleR)
reference_se <- SummarizedExperiment(
  assays = list(logcounts = as.matrix(scRNA_counts)),
  colData = DataFrame(label = FineCellType_vector)
)
for (i in originals) {
  xenium_counts <- GetAssayData(all_objects_filtered, assay = "SCT", layer = "data")[, Cells(all_objects_filtered[[i]])]
  singleR_results <- SingleR(test = xenium_counts, ref = reference$counts, labels = reference$labels, 
                             de.method="wilcox",   BPPARAM=MulticoreParam(10))
  
  singleR_results_list[[i]] <- singleR_results
}

singleR_results_pruned <- list()
Idents(all_objects_filtered) <- "orig.ident"
all_objects_filtered$Coarse_SingleR <- NA

for (i in originals){
  singleR_results_pruned[[i]] <-pruneScores(singleR_results_list[[i]], nmads = 2, min.diff.med = -Inf,  min.diff.next = 0.05, get.thresholds = FALSE)
  print(sum(singleR_results_pruned[[i]] == TRUE))
}
#The above filters out 8769, 5175, 5966, 5105, 
for (i in originals){
  singleR_results_list[[i]]$pruned.labels[singleR_results_pruned[[i]]] <- "Unknown"
  singleR_results_list[[i]]$pruned.labels[is.na(singleR_results_list[[i]]$pruned.labels)] <- "Unknown"
  all_objects_filtered$Coarse_SingleR[all_objects_filtered$orig.ident == i ]<-singleR_results_list[[i]]$pruned.labels
}
#Add the Score for visualization
all_objects_filtered$Coarse_SingleR_score <- NA
for (i in originals){
  all_objects_filtered$Coarse_SingleR_score[all_objects_filtered$orig.ident == i ]<-singleR_results_list[[i]]$scores
}
#Add the delta next metric for visualization (better than raw score and will tell us which populations struggle in annotation)
all_objects_filtered$Coarse_SingleR_deltaNext <- NA
for (i in originals){
  all_objects_filtered$Coarse_SingleR_deltaNext[all_objects_filtered$orig.ident == i ]<-singleR_results_list[[i]]$delta.next
}
#Visualize results and quality control of singleR assignments
plotDeltaDistribution(singleR_results_list[["naive1"]] , ncol = 3)
plotDeltaDistribution(singleR_results_list[["yersinia1"]] , ncol = 3)


plotScoreHeatmap(singleR_results_list[["yersinia2"]],   show.pruned = T)
plotScoreDistribution(singleR_results_list[["yersinia2"]])

FeaturePlot_scCustom(all_objects_filtered, reduction = "umap", features = "Coarse_SingleR_score", na_cutoff = NULL) 
FeaturePlot_scCustom(all_objects_filtered, reduction = "umap", features = "Coarse_SingleR_deltaNext", na_cutoff = 0) 
DimPlot(all_objects_filtered, reduction = "umap", group.by = "Coarse_SingleR" , label = T , repel = T, raster = T) 

FeaturePlot_scCustom(all_objects_filtered, reduction = "umap", features = "nCount_Xenium", na_cutoff = 500) 
table(all_objects_filtered$Coarse_SingleR) 

#Create Filtered Barcharts for Supplement----
combined_files2 <- all_objects_filtered@meta.data

combined_files <-combined_files2 %>%
  group_by(orig.ident) %>%
  summarise(
    n_cells = n(),
    median_counts = median(nCount_Xenium, na.rm = TRUE),
    median_genes = median(nFeature_Xenium, na.rm = TRUE),
    InfectionStatus = InfectionStatus[1]
  ) %>%
  ungroup()

infection_levels <- c("Naive",  "Yersinia")

ReadNumber <- ggplot(combined_files, aes(x = orig.ident, y = median_counts, fill = InfectionStatus)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") +
  scale_fill_manual(
    values = colors8,
    breaks = infection_levels
  ) +
  scale_x_discrete(labels = sample_labels) +  # map sample names → infection labels
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12),
    strip.background = element_rect(fill = "grey85", colour = "grey20"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ylab("median UMI per cell")

ReadNumber
#Figure S7C
ggsave(paste0(images, "/filtered_median_nCount_Xenium.svg"), plot = last_plot(), width = 11, height = 4)




ReadNumber <- ggplot(combined_files, aes(x = orig.ident, y = median_genes, fill = InfectionStatus)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") +
  scale_fill_manual(
    values = colors8,
    breaks = infection_levels
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12),
    strip.background = element_rect(fill = "grey85", colour = "grey20"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ylab("Median Genes per Cell")

ReadNumber
#Figure S7C
ggsave(paste0(images, "/filtered_Median_genes_per_cell.svg"), plot = last_plot(), width = 11, height = 4)

ReadNumber <- ggplot(combined_files, aes(x = orig.ident, y = n_cells, fill = InfectionStatus)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") +
  scale_fill_manual(
    values = colors8,
    breaks = infection_levels
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12),
    strip.background = element_rect(fill = "grey85", colour = "grey20"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ylab("Filtered # of Cells")

ReadNumber
#Figure S7C
ggsave(paste0(images, "/filtered_cell_number_xenium.svg"), plot = last_plot(), width = 11, height = 4)

#Save a Umap of the Coarse assignments ----
library(ggrepel)
umap_data <- all_objects_filtered[["umap"]]@cell.embeddings
cluster_data <- all_objects_filtered$Coarse_SingleR

# Combine UMAP coordinates and cluster information into a single dataframe
umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cluster_data)
unique()
centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(umap_1 = mean(umap_1), umap_2 = mean(umap_2))
umap_data$Infection <- all_objects_filtered$InfectionStatus
umap_data$orig.ident <- all_objects_filtered$orig.ident

spatial_umap <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(alpha = 0.1, size = 1.0) +  # Add points for each cell
  ylim(-15, 20) + xlim(-15, 20) +

  scale_color_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "black")) + 
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "black")) + 
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.grid = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
spatial_umap + geom_text_repel(data = centroids, aes(label = cluster), size = 4, fontface = "bold",    color = "black",
                             min.segment.length = 0
) 
ggsave( "spatial_umap_UMAP_CoarseCellTypel.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

ggsave( "Ileum_UMAP_CoarseCellType_NoLabels.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

ggplot(umap_data[!umap_data$cluster == "Unknown",], aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(alpha = 0.4, size = 1.0) +  # Add points for each cell
  ylim(-15, 20) + xlim(-15, 20) +
  scale_color_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.grid = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
#Figure S7D
ggsave( "Ileum_UMAP_CoarseCellType_NoLabels_Unlabeled_removed.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


cluster_colors <- c(
  "#CC79A7", "#E60000", "#56B4E9", "#009E73", 
  "#F0E442", "#0072B2", "#D55E00", "black"
)




#Save the objects 
library(qs2)
qs_save(all_objects_filtered, "/path/to/data/MIST/Xenium_Yp_2025/Seurat/CoarseAnnotatedSeurat.qs2")
write.csv(singleR_results, paste0(out, "/Coarse_SingleR_Results.csv"))

DimPlot(all_objects_filtered, reduction = "umap",   group.by = c(  "Coarse_SingleR"),
        combine = T, label.size = 3, label = TRUE, raster = F) + NoLegend()
ggsave( "Umap_AllData_CoarseSingleR.png",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

#Export the Swiss Roll Coordinates of every cell----

coords <- GetTissueCoordinates(all_objects_filtered, image = "naive1") 
head(coords)
# Merge with meta.data
coords_with_meta <- cbind(
  cell_id = rownames(coords),
  coords)

write.csv(coords_with_meta[, c("cell_id", "x", "y")],
          paste0(out,"/naive1_xenium_cell_coordinates.csv"),
          row.names = FALSE)

coords <- GetTissueCoordinates(all_objects_filtered, image = "naive2") 
head(coords)
# Merge with meta.data
coords_with_meta <- cbind(
  cell_id = rownames(coords),
  coords)

write.csv(coords_with_meta[, c("cell_id", "x", "y")],
          paste0(out,"/naive2_xenium_cell_coordinates.csv"),
          row.names = FALSE)

coords <- GetTissueCoordinates(all_objects_filtered, image = "yersinia1") 
head(coords)
# Merge with meta.data
coords_with_meta <- cbind(
  cell_id = rownames(coords),
  coords)

write.csv(coords_with_meta[, c("cell_id", "x", "y")],
          paste0(out,"/yersinia1_xenium_cell_coordinates.csv"),
          row.names = FALSE)

coords <- GetTissueCoordinates(all_objects_filtered, image = "yersinia2") 
head(coords)
# Merge with meta.data
coords_with_meta <- cbind(
  cell_id = rownames(coords),
  coords)

write.csv(coords_with_meta[, c("cell_id", "x", "y")],
          paste0(out,"/yersinia2_xenium_cell_coordinates.csv"),
          row.names = FALSE)
#Subset and Recluster pt 1----
all_objects_filtered <- qs_read( "/path/to/data/MIST/Xenium_Yp_2025/Seurat/CoarseAnnotatedSeurat.qs2")


# Subset the object by each CoarseCellType and store them in a list
CoarseClusts <- unique(all_objects_filtered$Coarse_SingleR)
subset_objects <- lapply(CoarseClusts, function(cell_type) {subset(all_objects_filtered, Coarse_SingleR == cell_type)})
names(subset_objects) <- CoarseClusts


#Loop over the data and subset and recluster each coarse group
for (cell_type in CoarseClusts) {
  # Get the subset for this cell type
  subset_obj <- subset_objects[[cell_type]]
  
  # Repeat the normalization and identification of variable features with SCT TRansform 
  subset_obj <- SCTransform(subset_obj, assay = "Xenium", variable.features.n = 500) 
  subset_obj <- ScaleData(subset_obj)
  
  # PCA for dimensionality reduction
  subset_obj <- RunPCA(subset_obj, npcs = 30) 
  ElbowPlot(subset_obj, ndims = 30, reduction = "pca")
  subset_obj <- RunUMAP(subset_obj, dims = 1:30, n.neighbors = 50L, min.dist = 0.2) 
  subset_obj <- FindNeighbors(subset_obj, reduction = "pca", dims = 1:30) 
  subset_obj <- FindClusters(subset_obj, resolution = 1.5,  n.start = 30, 
                             cluster.name = paste0(gsub(" ", "_", cell_type), "_clusters")  )
  
  library(SeuratWrappers)
  cds <- as.cell_data_set(
    subset_obj,
    reductions = "umap",
    default.reduction = "umap",
    graph =  "SCT_snn",
    group.by = NULL)
  
  reducedDimNames(cds) <- "UMAP"
  cds <-cluster_cells(
    cds,
    reduction_method = "UMAP", 
    k = 15, cluster_method =  "leiden", 
    #resolution = 0.0001,
    num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )
  
  Seurat_monocle <- as.Seurat(cds, assay = NULL)
  subset_obj[["Monocle_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
  
  #Delete the temporary Seurat monocle object
  remove(Seurat_monocle)
  remove(cds)
  
  # Store the updated Seurat object back into the list
  subset_objects[[cell_type]] <- subset_obj
  
}

save.image("/path/to/data/MIST/Xenium_Yp_2025/Subclustered_Xenium.RData")



 
# Subset and Recluster part 2 Fine labeling  ---- 

#Visualize results
DimPlot(subset_objects[["Epithelial cells"]], reduction = "umap", group.by = "seurat_clusters")
DimPlot(subset_objects[["Myeloid cells"]], reduction = "umap", group.by = "seurat_clusters")
DimPlot(subset_objects[["Stromal cells"]], reduction = "umap", group.by = "seurat_clusters")
DimPlot(subset_objects[["B cells"]], reduction = "umap", group.by = "seurat_clusters")
DimPlot(subset_objects[["Plasma cells"]], reduction = "umap", group.by = "seurat_clusters")
DimPlot(subset_objects[["T cells and ILCs"]], reduction = "umap", group.by = "seurat_clusters")
DimPlot(subset_objects[["Enteric Nervous System"]], reduction = "umap", group.by = "seurat_clusters")
FeaturePlot(subset_objects[["Epithelial cells"]], reduction = "umap", features  = "Saa1")

#Label Cells with SingleR
# Create the reference data for SingleR

#Epithelial Cells ----

Idents(Ileum_annotated) <- "CoarseCellType"
unique(Ileum_annotated$CoarseCellType)
Ileum_Epithelial <- subset(Ileum_annotated, idents = "Epithelial cells")
unique(Ileum_Epithelial$FineCellType)
scRNA_counts <- GetAssayData(Ileum_Epithelial, assay = "RNA", slot = "data")
scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(subset_objects[["Epithelial cells"]]@assays$SCT), ]
unique(Ileum_Epithelial$FineCellType)
FineCellType_vector <- Ileum_Epithelial$FineCellType
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

# Extract Xenium info
#subset_objects[["Epithelial cells"]]@assays$SCT <- JoinLayers(subset_objects[["Epithelial cells"]]@assays$SCT)
xenium_counts <- GetAssayData(subset_objects[["Epithelial cells"]], assay = "SCT", slot = "data") 
DimPlot(subset_objects[["Epithelial cells"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()

# Run SingleR
subset_objects[["Epithelial cells"]] <-  split(x = subset_objects[["Epithelial cells"]], f = subset_objects[["Epithelial cells"]]$orig.ident) 

singleR_results <- SingleR(test = xenium_counts, ref = reference$counts, labels = reference$labels, 
                           de.method="wilcox",
                           BPPARAM=MulticoreParam(25))


subset_objects[["Epithelial cells"]]$SingleR_Fine_Scores <- singleR_results$scores
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "SingleR_Fine_Scores", na_cutoff = 0)

# Extract cell type predictions (e.g., the most likely cell type for each cell)
subset_objects[["Epithelial cells"]]$FineLabelsSingleR <- singleR_results$pruned.labels
DimPlot(subset_objects[["Epithelial cells"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of Epithelial cells")) + NoLegend()
table(subset_objects[["Epithelial cells"]]$InfectionStatus, subset_objects[["Epithelial cells"]]$FineLabelsSingleR)
unique(subset_objects[["Epithelial cells"]]$FineLabelsSingleR)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Myh11", na_cutoff = 0)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Cd3e", na_cutoff = 0)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Lum", na_cutoff = 0)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Pou2f3", na_cutoff = 0)

ImageDimPlot(subset_objects[["Epithelial cells"]], group.by = "FineLabelsSingleR", fov = "Yersinia2")
ImageDimPlot(subset_objects[["Epithelial cells"]], group.by = "FineLabelsSingleR", fov = "Naive1")
DimPlot(subset_objects[["Epithelial cells"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()
remove(Ileum_Epithelial)

#It is clear the there are cells labeled Epithelial which contain contaminating transcripts likely from proximal cells. This can be visualized by the Myh11, Lum, Cd3e etc. These clusters are 103, 89, 11, 121, 122, 125, 126, 157, 46, and 2
#I will remove these cells from the object and they will become "Unknown" in the larger objects 
contaminant <- c("103", "89", "11", "121", "122", "125", "126", "157", "46",  "2", "111")
Idents(subset_objects[["Epithelial cells"]]) <- "Monocle_clusters"
contaminant_ids <- WhichCells(subset_objects[["Epithelial cells"]], idents =  contaminant)
all_objects_filtered$Coarse_SingleR[colnames(all_objects_filtered) %in% contaminant_ids] <- "Unknown"
Idents(subset_objects[["Epithelial cells"]]) <- "Monocle_clusters"
subset_objects[["Epithelial cells filt"]] <- subset(subset_objects[["Epithelial cells"]], idents = contaminant, invert = TRUE)

# Recluster the data and visalize 
DimPlot(subset_objects[["Epithelial cells filt"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of Epithelial cells")) + NoLegend()
DimPlot(subset_objects[["Epithelial cells filt"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()


subset_obj <-subset_objects[["Epithelial cells filt"]]
  
  # Repeat the normalization and identification of variable features with SCT TRansform 
subset_obj <- SCTransform(subset_obj, assay = "Xenium", variable.features.n = 500) 
subset_obj <- ScaleData(subset_obj)
  
  # PCA for dimensionality reduction
  subset_obj <- RunPCA(subset_obj, npcs = 30) # 8 minutes
  ElbowPlot(subset_obj, ndims = 30, reduction = "pca")
  subset_obj <- RunUMAP(subset_obj, dims = 1:30, n.neighbors = 50L, min.dist = 0.2) 
  subset_obj <- FindNeighbors(subset_obj, reduction = "pca", dims = 1:30) 
  subset_obj <- FindClusters(subset_obj, resolution = 1.5,  n.start = 30, 
                             cluster.name = paste0(gsub(" ", "_", cell_type), "_clusters")  )
  
  library(SeuratWrappers)
  cds <- as.cell_data_set(
    subset_obj,
    reductions = "umap",
    default.reduction = "umap",
    graph =  "SCT_snn",
    group.by = NULL)
  
  reducedDimNames(cds) <- "UMAP"
  cds <-cluster_cells(
    cds,
    reduction_method = "UMAP", 
    k = 15, cluster_method =  "leiden", 
    #resolution = 0.0001,
    num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )
  
  Seurat_monocle <- as.Seurat(cds, assay = NULL)
  subset_obj[["Monocle_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
  
  #Delete the temporary Seurat monocle object
  remove(Seurat_monocle)
  remove(cds)
  
  # Store the updated Seurat object back into the list
  subset_objects[["Epithelial cells filt"]] <- subset_obj
  
remove(subset_obj)
DimPlot(subset_objects[["Epithelial cells filt"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of Epithelial cells")) + NoLegend()
DimPlot(subset_objects[["Epithelial cells filt"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()
FeaturePlot_scCustom(subset_objects[["Epithelial cells filt"]], reduction =  "umap", "Tph1", na_cutoff = 0)


# Create Epithelial  Barchart & Umap ----
umap_data <- subset_objects[["Epithelial cells filt"]][["umap"]]@cell.embeddings
cluster_data <- subset_objects[["Epithelial cells filt"]]$FineLabelsSingleR
infection <- subset_objects[["Epithelial cells filt"]]$InfectionStatus
umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cluster_data)
umap_data$infection <- factor(infection)
centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(umap_1 = mean(umap_1), umap_2 = mean(umap_2))

cbf_18_colors <- c( "#1170AA",   "#FCB13F",   "#60BFC1",   "#EF6F6A",   "#937860", "#D17A00",   "#B78CBA",  
  "#B3B3B3",   "#64B200",  "#D45E00",   "#7E8CCF",   "#E6A0C4",   "#568B3F",   "#C44E52", "#5FA2A3",  "#CCB974", 
  "#D0A6BE",  "#4E84C4"  )


spatial_umap <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(alpha = 0.1, size = 1.0) +  # Add points for each cell
  ylim(-15, 20) + xlim(-15, 20) +
  #geom_mark_hull(aes(fill = cluster, alpha = 0.01), alpha = 0.1, color = "NA", concavity = 10,    con.type = "none") +  # Add shaded hulls for each cluster
  scale_color_manual(values = cbf_18_colors) + # Custom cluster colors
  scale_fill_manual(values = cbf_18_colors) + # Fill color for hulls
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.grid = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
spatial_umap + geom_text_repel(data = centroids, aes(label = cluster), size = 4, fontface = "bold",    color = "black",
                               min.segment.length = 0
) 
ggsave( "Epithelial_subcluster_UMAP_FineCell.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
spatial_umap
ggsave( "Epithelial_subcluster_UMAP_FineCel_NoLabels.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



spatial_umap <- ggplot(umap_data[umap_data$infection == "Naive",], aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.0) + 
  ylim(-15, 20) + xlim(-15, 20) +
  scale_color_manual(values = cbf_18_colors) + 
  scale_fill_manual(values = cbf_18_colors) + 
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.grid = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
spatial_umap + geom_text_repel(data = centroids, aes(label = cluster), size = 4, fontface = "bold",    color = "black",
                               min.segment.length = 0
) 
ggsave( "Epithelial_subcluster_Naive_UMAP_FineCell.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
spatial_umap
#Figure 7A
ggsave( "Epithelial_subcluster_Naive_UMAP_FineCel_NoLabels.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
ggsave( "Epithelial_subcluster_Naive_UMAP_FineCel_NoLabels.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
spatial_umap +theme(legend.position = "bottom")
ggsave( "Epithelial_subcluster_Naive_UMAP_FineCel_legend.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

spatial_umap <- ggplot(umap_data[umap_data$infection == "Yersinia",], aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.0) + 
  ylim(-15, 20) + xlim(-15, 20) +
  scale_color_manual(values = cbf_18_colors) + 
  scale_fill_manual(values = cbf_18_colors) + 
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.grid = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
spatial_umap + geom_text_repel(data = centroids, aes(label = cluster), size = 4, fontface = "bold",    color = "black",
                               min.segment.length = 0
) 
ggsave( "Epithelial_subcluster_Yersinia_UMAP_FineCell.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
spatial_umap
#Figure 7A
ggsave( "Epithelial_subcluster_Yersinia_UMAP_FineCel_NoLabels.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
ggsave( "Epithelial_subcluster_Yersinia_UMAP_FineCel_NoLabels.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



cluster_freq <- umap_data %>%
  group_by(infection, cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(infection) %>%
  mutate(percent = count / sum(count) * 100)
cluster_freq <- cluster_freq[!cluster_freq$cluster == "NA",]
ggplot(cluster_freq, aes(x = factor(infection), y = percent, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  #facet_wrap(~ infection) +
  scale_fill_manual(values = cbf_18_colors) +
  labs(x = "Cluster", y = "Percentage of Cells") +
  #scale_y_break(c(25, 60)) +  
  theme_minimal() +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave( "Barchart_spatial_CoarseCellType_Epithelial.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


# What does finest cell type look llike
Idents(Ileum_annotated) <- "CoarseCellType"
unique(Ileum_annotated$CoarseCellType)
Ileum_Epithelial <- subset(Ileum_annotated, idents = "Epithelial cells")
unique(Ileum_Epithelial$FinestCellType)
scRNA_counts <- GetAssayData(Ileum_Epithelial, assay = "RNA", slot = "data")
scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(subset_objects[["Epithelial cells filt"]]@assays$SCT), ]
unique(Ileum_Epithelial$FinestCellType)
FinestCellType_vector <- Ileum_Epithelial$FinestCellType
reference <- list(counts = scRNA_counts, labels = FinestCellType_vector)

# Extract Xenium info
subset_objects[["Epithelial cells filt"]]@assays$SCT <- JoinLayers(subset_objects[["Epithelial cells filt"]]@assays$SCT)
xenium_counts <- GetAssayData(subset_objects[["Epithelial cells filt"]], assay = "SCT", layer = "data") 
DimPlot(subset_objects[["Epithelial cells filt"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()

# Run SingleR
subset_objects[["Epithelial cells filt"]] <-  split(x = subset_objects[["Epithelial cells filt"]], f = subset_objects[["Epithelial cells filt"]]$orig.ident) 

singleR_results <- SingleR(test = xenium_counts, ref = reference$counts, labels = reference$labels, 
                           de.method="wilcox",
                           BPPARAM=MulticoreParam(25))


subset_objects[["Epithelial cells filt"]]$SingleR_finest_Scores <- singleR_results$scores
FeaturePlot_scCustom(subset_objects[["Epithelial cells filt"]], reduction =  "umap", "SingleR_finest_Scores", na_cutoff = 0)

# Extract cell type predictions (e.g., the most likely cell type for each cell)
subset_objects[["Epithelial cells filt"]]$finestLabelsSingleR <- singleR_results$pruned.labels
DimPlot(subset_objects[["Epithelial cells filt"]], reduction = "umap", group.by = "finestLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of Epithelial cells")) + NoLegend()


all_objects_filtered <- AddMetaData(all_objects_filtered, subset_objects[["Epithelial cells filt"]]$finestLabelsSingleR , col.name = "FineCellType")
ImageDimPlot(all_objects_filtered, fov = "yersinia1", group.by = "FineCellType")

obj <- all_objects_filtered

# Create a new column for highlight status
obj$HighlightGroup <- ifelse(
  obj$FineCellType == "Saa1 inflammatory Enterocytes",
  "Highlighted",
  "Other"
)

# Plot
ImageDimPlot(
  obj,
  fov = "yersinia1",
  group.by = "HighlightGroup",
  cols = c("Highlighted" = "red", "Other" = "grey80")
)

table(subset_objects[["Epithelial cells"]]$InfectionStatus, subset_objects[["Epithelial cells"]]$finestLabelsSingleR)
unique(subset_objects[["Epithelial cells"]]$finestLabelsSingleR)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Myh11", na_cutoff = 0)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Cd3e", na_cutoff = 0)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Lum", na_cutoff = 0)
FeaturePlot_scCustom(subset_objects[["Epithelial cells"]], reduction =  "umap", "Pou2f3", na_cutoff = 0)


# Ucell enrichment of Epithelial Cells for Yps+ Enterocytes - The script "Analysis10_Xenium_UCell_Yp.R" was run and produced the csv with cell id and score. (Analysis10_Xenium_UCell_Yp.R)

YP_Enrich <- data.table::fread("/path/to/data/MIST/Xenium_Yp_2025/CSV/UCell_PRE_202508.csv")
cell_id_keep <- YP_Enrich$cell_id
all(colnames(subset_objects[["Epithelial cells filt"]]) %in% rownames(YP_Enrich) )
YP_Enrich <- YP_Enrich[, 2:ncol(YP_Enrich)] # remove cell id column 
rownames(YP_Enrich) <- cell_id_keep
all(colnames(subset_objects[["Epithelial cells filt"]]) %in% rownames(YP_Enrich) )
all_objects_filtered@meta.data$Saa1_Gene_Lists <- NULL
subset_objects[["Epithelial cells filt"]] <- AddMetaData(subset_objects[["Epithelial cells filt"]], metadata = YP_Enrich)
all_objects_filtered <-  AddMetaData(all_objects_filtered, metadata = YP_Enrich)
Yersinia_subset_epi <-  AddMetaData(Yersinia_subset_epi, metadata = YP_Enrich)
#all_objects_filtered$Saa1_Gene_Lists[all_objects_filtered$Saa1_Gene_Lists =="NA" ] <- NA
#range(all_objects_filtered$Saa1_Gene_Lists)


Hallmark_Enrich <- data.table::fread("/path/to/data/MIST/Xenium_Yp_2025/CSV/UCell_Hallmark_202511.csv")
cell_id_keep <- Hallmark_Enrich$cell_id
all(colnames(subset_objects[["Epithelial cells filt"]]) %in% rownames(Hallmark_Enrich) )
Hallmark_Enrich <- Hallmark_Enrich[, 2:ncol(Hallmark_Enrich)] # remove cell id column 
rownames(Hallmark_Enrich) <- cell_id_keep
all(colnames(subset_objects[["Epithelial cells filt"]]) %in% rownames(Hallmark_Enrich) )

subset_objects[["Epithelial cells filt"]] <- AddMetaData(subset_objects[["Epithelial cells filt"]], metadata = Hallmark_Enrich)
all_objects_filtered <-  AddMetaData(all_objects_filtered, metadata = Hallmark_Enrich)
all_objects_filtered$Saa1_Gene_Lists[all_objects_filtered$Saa1_Gene_Lists =="NA" ] <- NA
range(all_objects_filtered$Saa1_Gene_Lists)


ImageFeaturePlot(all_objects_filtered, fov = "yersinia1", features = "Saa1")
ImageFeaturePlot(all_objects_filtered, fov = "naive1", features = "Saa1")

expr_vals <- FetchData(all_objects_filtered, "Saa1")[, 1]
min_val <- min(expr_vals, na.rm = TRUE)
max_val <- max(expr_vals, na.rm = TRUE)

ImageFeaturePlot(
  all_objects_filtered,
  fov = "yersinia1",
  features = "Saa1",
  min.cutoff = min_val,
  max.cutoff = max_val
)+ scale_fill_gradient(low = "grey", high = "#D10000", limits = c(0, 4))
#figure 7J
ggsave( "Saa1_FeatureImage_Yersinia1_Image.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


ImageFeaturePlot(
  all_objects_filtered,
  fov = "naive1",
  features = "Saa1"
) + scale_fill_gradient(low = "grey", high = "#D10000", limits = c(0, 4))
#figure 7J
ggsave( "Saa1_FeatureImage_Naive1_Image.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)




ImageFeaturePlot(all_objects_filtered, fov = "yersinia1", features = "Saa1_Gene_Lists")
ImageFeaturePlot(subset_objects[["Epithelial cells filt"]], fov = "yersinia1", features = "Saa1_Gene_Lists", 
                 cols = c("navy", "lightgrey", "darkred"))


ggsave( "Saa1_geneSet_Enrich_Yersinia1_Image.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



ImageFeaturePlot(subset_objects[["Epithelial cells filt"]], fov = "naive1", features = "Saa1_Gene_Lists", 
                 cols = c("navy", "lightgrey", "darkred"))

ggsave( "Saa1_geneSet_Enrich_Naive1_Image.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


#Utilize cell coordinates and Intestline to unravel the small intestine ----
#This analysis uses the cell cell coordinates previously exported abov and loads them into "Intestline" (PMID: 36943334). In Intestline, I selected a series of close points along the outer edge of a segment of the intestine - beginning inside and working outwork when possible. I then chose a second set of consecutive points for the same segment of tissue that followed along the inside edge. These outer and inner edge points were loaded into R as below where I joined the ends to form a trapezoid containing the region of tissue I would work with. I then calculated  midline from the first pair of points (outer and inner) to the last pair of point following along the shape of the ribbon/trapezoid. I then assigned a relative distance of 0 and 100 to the beginning and end of that line. Each cell was assigned a distance on this axis based where on the centerline it was closest to.  Thus. Cells were roughly ordered on the axis of the intestine from proximal to distal
# Naive1 Coordinate Based Analysis ----



naive1_outer <- fread(paste0(CSV, "naive1_outer_sep10.csv"))
naive1_outer <- naive1_outer[, c(1,2)]
naive1_outer$x <- as.numeric(naive1_outer$x)
naive1_outer$y <- as.numeric(naive1_outer$y)
naive1_inner <- fread(paste0(CSV, "naive1_inner_sep10.csv"))
naive1_inner <- naive1_inner[, c(1,2)]
naive1_inner$x <- as.numeric(naive1_inner$x)
naive1_inner$y <- as.numeric(naive1_inner$y)

library(sf)
# Make sure both are in same order along ribbon

ribbon_coords <- rbind(
  naive1_outer,
  naive1_inner[nrow(naive1_inner):1, ], # reverse inner coords
  naive1_outer[1, ]  # close polygon by repeating first point
)

ribbon_polygon <- st_polygon(list(as.matrix(ribbon_coords))) %>%
  st_sfc(crs = NA_crs_) %>%   
  st_make_valid()

query_df <- GetTissueCoordinates(all_objects_filtered, image = "naive1")
query_points <- st_as_sf(query_df, coords = c("x", "y"), crs = NA_crs_)

# Find points inside polygon
inside_idx <- st_within(query_points, ribbon_polygon, sparse = FALSE)[,1]

# Subset coordinates
inside_cells <- query_df[inside_idx, ]


Naive_subset <- subset(all_objects_filtered, cells = inside_cells$cell)

# Make centerline
n_inner <- nrow(naive1_inner)
n_outer <- nrow(naive1_outer)
n_min <- min(n_inner, n_outer)

# Match up shorter spiral to longer one
inner_resampled <- naive1_inner[round(seq(1, n_inner, length.out = n_min)), ]
outer_resampled <- naive1_outer[round(seq(1, n_outer, length.out = n_min)), ]

centerline <- data.frame(
  x = (inner_resampled$x + outer_resampled$x) / 2,
  y = (inner_resampled$y + outer_resampled$y) / 2
)

# Cumulative distances along centerline
distances <- sqrt(diff(centerline$x)^2 + diff(centerline$y)^2)
cumdist <- c(0, cumsum(distances))
total_dist <- max(cumdist)

# Interpolate to get, say, 10x more points along the centerline
interp_factor <- 50
interp_idx <- seq(0, total_dist, length.out = n_min * interp_factor)

# Cumulative distance already computed: cumdist
centerline_interp <- data.frame(
  x = approx(cumdist, centerline$x, xout = interp_idx)$y,
  y = approx(cumdist, centerline$y, xout = interp_idx)$y
)

# Update distances
distances_interp <- sqrt(diff(centerline_interp$x)^2 + diff(centerline_interp$y)^2)
cumdist_interp <- c(0, cumsum(distances_interp))
total_dist_interp <- max(cumdist_interp)

assign_progress <- function(px, py) {
  dists_to_center <- sqrt((centerline_interp$x - px)^2 + (centerline_interp$y - py)^2)
  nearest_idx <- which.min(dists_to_center)
  progress <- (cumdist_interp[nearest_idx] / total_dist_interp) * 100
  return(progress)
}

inside_cells$progress <- mapply(assign_progress, inside_cells$x, inside_cells$y)

#Save the new distance to the seurat object
Naive_subset$progress <- inside_cells$progress[match(Cells(Naive_subset), inside_cells$cell)]


ImageFeaturePlot(Naive_subset, "progress")
#Figure S7F
ggsave( "Naive1_Spiral_Progression_Visualization.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


naive1_subset_epi <- subset(subset_objects[["Epithelial cells filt"]], cells = inside_cells$cell)
naive1_subset_epi$progress <- inside_cells$progress[match(Cells(naive1_subset_epi), inside_cells$cell)]
#ImageFeaturePlot(naive1_subset_epi, "progress")
#ImageFeaturePlot(naive1_subset_epi, fov = "naive1", features = "Saa1_Gene_Lists",   cols = c("navy", "lightgrey", "yellow"))


#ggsave( "naive1_Spiral_Saa1_enrichment_visualization.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
j <- as.matrix(all_objects_filtered@assays$SCT@scale.data)

#line plot
# Extract data from Seurat object
plot_data <- data.frame(
  progress = naive1_subset_epi$progress,
  Saa1_Gene_Lists = FetchData(naive1_subset_epi, "Saa1_Gene_Lists")[,1]
)



# Remove rows with NA values
plot_data <- na.omit(plot_data)

# Bin "progress" into intervals of 0.2
avg_data <- plot_data %>%
  mutate(progress_bin = floor(progress / 0.2) * 0.2) %>%
  group_by(progress_bin) %>%
  summarise(
    mean_Saa1 = mean(Saa1_Gene_Lists, na.rm = TRUE),  # or mean(...)
    .groups = "drop"
  )

library(zoo)  # for rollapply

# Sort and smooth using a rolling median with a window of ~6 points
avg_data_smooth <- avg_data %>%
  arrange(progress_bin) %>%
  mutate(mean_Saa1_smooth = rollapply(mean_Saa1, width = 6, FUN = mean, fill = NA, align = "center"))


ggplot(avg_data_smooth, aes(x = progress_bin, y = mean_Saa1_smooth)) +
  geom_line(color = "#B78CBA", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean Yp Saa1+ GeneSet Enrichment",
    title = "Smoothed mean Yp Saa1+ GeneSet Enrichment"
  ) +
  ylim(0.25, 0.4)+
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )

ggsave( "LineGraph_Saa1_enichment_along_Distance_naive1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



#TNF Expression Across Distance
# Extract data for Tnf from Seurat
plot_data_tnf <- data.frame(
  progress = Naive_subset$progress,
  Tnf = FetchData(Naive_subset, "Tnf")[,1],
  Il1a = FetchData(Naive_subset, "Il1a")[,1]
) %>%
  na.omit()

# Bin "progress" into intervals of 0.2 and compute average
avg_data_tnf <- plot_data_tnf %>%
  mutate(progress_bin = floor(progress / 0.2) * 0.2) %>%
  group_by(progress_bin) %>%
  summarise(
    mean_Tnf = mean(Tnf, na.rm = TRUE),
    mean_Il1a =  mean(Il1a, na.rm = TRUE),
    .groups = "drop"
  )

# Smooth using rolling median 
avg_data_tnf_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_Tnf_smooth = rollapply(mean_Tnf, width = 6, FUN = median, fill = NA, align = "center"))


ggplot(avg_data_tnf_smooth, aes(x = progress_bin, y = mean_Tnf_smooth)) +
  geom_line(color = "#B78CBA", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Median Yp Tnf Expression",
    title = "Smoothed Median Yp Tnf Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_Tnf_exp_along_Distance_naive1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

avg_data_il1a_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_il1a_smooth = rollapply(mean_Il1a, width = 6, FUN = median, fill = NA, align = "center"))

ggplot(avg_data_il1a_smooth, aes(x = progress_bin, y = mean_il1a_smooth)) +
  geom_line(color = "#B78CBA", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean IL1a Expression",
    title = "Smoothed  IL1a Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_IL1a_exp_along_Distance_naive1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


# Extract cell data from Seurat
plot_data_myeloid <- data.frame(
  progress = Naive_subset$progress,
  cell_type = Naive_subset$Coarse_SingleR  
) %>%
  na.omit()

unique(Naive_subset$Coarse_SingleR)
# Bin "progress" into intervals of 0.2 and compute percent myeloid cells
freq_data_myeloid <- plot_data_myeloid %>%
  mutate(progress_bin = floor(progress / 0.2) * 0.2) %>%
  group_by(progress_bin) %>%
  summarise(
    myeloid_percent = mean(cell_type %in% "Myeloid cells") * 100,
    .groups = "drop"
  )

# Smooth using rolling median
freq_data_myeloid_smooth <- freq_data_myeloid %>%
  arrange(progress_bin) %>%
  mutate(myeloid_percent_smooth = rollapply(myeloid_percent, width = 6, FUN = median, fill = NA, align = "center"))


ggplot(freq_data_myeloid_smooth, aes(x = progress_bin, y = myeloid_percent_smooth)) +
  geom_line(color = "#B78CBA", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Percent Myeloid Cells",
    title = "Smoothed Percent Myeloid Cells Across Progress"
  ) +
  theme_minimal() +
  ylim(0,40) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_Percent_Myeloid_along_Distance_naive1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

# remove what was created
remove(ribbon_coords)
remove(ribbon_polygon)
remove(naive1_inner)
remove(naive1_outer)
#remove(Naive_subset)
remove(query_df)
remove(query_points)
remove(inside_cells)
remove(inside_idx)
remove(n_inner)
remove(n_outer)
#remove(naive1_subset_epi)
remove(inner_resampled)
remove(outer_resampled)
remove(centerline)
remove(distances)
remove(cumdist)
remove(distances_interp)
remove(total_dist)
remove(interp_factor)
remove(plot_data)
remove(avg_data)
remove(avg_data_smooth)
remove(interp_factor)
remove(plot_data)
remove(plot_data_myeloid)
remove(freq_data_myeloid)
remove(freq_data_myeloid_smooth)



























#Yersinia1 Coordinate based analysis----
yp1_outer <- fread(paste0(CSV, "/yp1_outer_sep2.csv"))
yp1_outer <- yp1_outer[, c(1,2)]
yp1_outer$x <- as.numeric(yp1_outer$x)
yp1_outer$y <- as.numeric(yp1_outer$y)
yp1_inner <- fread(paste0(CSV, "/yp1_inner_sep2.csv"))
yp1_inner <- yp1_inner[, c(1,2)]
yp1_inner$x <- as.numeric(yp1_inner$x)
yp1_inner$y <- as.numeric(yp1_inner$y)

library(sf)
# Make sure both are in same order along ribbon

ribbon_coords <- rbind(
  yp1_outer,
  yp1_inner[nrow(yp1_inner):1, ], # reverse inner coords
  yp1_outer[1, ]  # close polygon by repeating first point
)

ribbon_polygon <- st_polygon(list(as.matrix(ribbon_coords))) %>%
  st_sfc(crs = NA_crs_) %>%   # use NA_crs_ for undefined Cartesian CRS
  st_make_valid()

query_df <- GetTissueCoordinates(all_objects_filtered, image = "yersinia1")
query_points <- st_as_sf(query_df, coords = c("x", "y"), crs = NA_crs_)

# Find points inside polygon
inside_idx <- st_within(query_points, ribbon_polygon, sparse = FALSE)[,1]

# Subset coordinates
inside_cells <- query_df[inside_idx, ]


Yersinia_subset <- subset(all_objects_filtered, cells = inside_cells$cell)

# Make centerline
n_inner <- nrow(yp1_inner)
n_outer <- nrow(yp1_outer)
n_min <- min(n_inner, n_outer)

# Match up shorter spiral to longer one
inner_resampled <- yp1_inner[round(seq(1, n_inner, length.out = n_min)), ]
outer_resampled <- yp1_outer[round(seq(1, n_outer, length.out = n_min)), ]

centerline <- data.frame(
  x = (inner_resampled$x + outer_resampled$x) / 2,
  y = (inner_resampled$y + outer_resampled$y) / 2
)

# Cumulative distances along centerline
distances <- sqrt(diff(centerline$x)^2 + diff(centerline$y)^2)
cumdist <- c(0, cumsum(distances))
total_dist <- max(cumdist)

# Interpolate to get, say, 10x more points along the centerline
interp_factor <- 50
interp_idx <- seq(0, total_dist, length.out = n_min * interp_factor)

# Cumulative distance already computed cumdist
centerline_interp <- data.frame(
  x = approx(cumdist, centerline$x, xout = interp_idx)$y,
  y = approx(cumdist, centerline$y, xout = interp_idx)$y
)

# Update distances
distances_interp <- sqrt(diff(centerline_interp$x)^2 + diff(centerline_interp$y)^2)
cumdist_interp <- c(0, cumsum(distances_interp))
total_dist_interp <- max(cumdist_interp)

assign_progress <- function(px, py) {
  dists_to_center <- sqrt((centerline_interp$x - px)^2 + (centerline_interp$y - py)^2)
  nearest_idx <- which.min(dists_to_center)
  progress <- (cumdist_interp[nearest_idx] / total_dist_interp) * 100
  return(progress)
}

inside_cells$progress <- mapply(assign_progress, inside_cells$x, inside_cells$y)

#Save the new distance to the seurat object
Yersinia_subset$progress <- inside_cells$progress[match(Cells(Yersinia_subset), inside_cells$cell)]


ImageFeaturePlot(Yersinia_subset, "progress")
#Figure S7F
ggsave( "Yersinia1_Spiral_Progression_Visualization.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

Idents(all_objects_filtered) <- "Coarse_SingleR"
Yersinia_subset_epi<- subset(all_objects_filtered, idents = "Epithelial cells")
Yersinia_subset_epi <- subset(Yersinia_subset_epi, cells = inside_cells$cell)

Yersinia_subset_epi$progress <- inside_cells$progress[match(Cells(Yersinia_subset_epi), inside_cells$cell)]
ImageFeaturePlot(Yersinia_subset_epi, "progress")
ImageFeaturePlot(Yersinia_subset_epi, fov = "yersinia1", features = "Saa1_Gene_Lists", 
                 cols = c("navy", "lightgrey", "yellow"))


ggsave( "Yersinia1_Spiral_Saa1_enrichment_visualization.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


#line plot
# Extract data from Seurat object
plot_data <- data.frame(
  progress = Yersinia_subset_epi$progress,
  Saa1_Gene_Lists = FetchData(Yersinia_subset_epi, "Saa1_Gene_Lists")[,1]
)

# Remove rows with NA values
plot_data <- na.omit(plot_data)

# Bin "progress" into intervals of 0.2
avg_data <- plot_data %>%
  mutate(progress_bin = floor(progress / 0.2) * 0.2) %>%
  group_by(progress_bin) %>%
  summarise(
    mean_Saa1 = mean(Saa1_Gene_Lists, na.rm = TRUE),  # or mean(...)
    .groups = "drop"
  )

library(zoo)  # for rollapply

# Sort and smooth using a rolling median with a window of ~6 points
avg_data_smooth <- avg_data %>%
  arrange(progress_bin) %>%
  mutate(mean_Saa1_smooth = rollapply(mean_Saa1, width = 6, FUN = mean, fill = NA, align = "center"))

highlight_regions <- data.frame(
  xmin = c(12, 28, 51),
  xmax = c(16, 32, 55),
  ymin = -Inf,
  ymax = Inf
)

ggplot(avg_data_smooth, aes(x = progress_bin, y = mean_Saa1_smooth)) +
  # Add grey background for highlighted regions
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color = "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean Yp Saa1+ GeneSet Enrichment",
    title = "Smoothed mean Yp Saa1+ GeneSet Enrichment"
  ) +
  theme_minimal() +
  ylim(0.25, 0.40)+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )

ggsave( "LineGraph_Saa1_enichment_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



Yersinia_subset_epi$Progress_Highlight <- ifelse(
    (Yersinia_subset_epi$progress >= 12 & Yersinia_subset_epi$progress <= 16 |
       Yersinia_subset_epi$progress >= 28 & Yersinia_subset_epi$progress <= 32 |
      Yersinia_subset_epi$progress >= 51 & Yersinia_subset_epi$progress <= 55),
  "Highlighted", "Other"
)

# Check the counts
table(Yersinia_subset_epi$Progress_Highlight)

# Now you can use it for coloring in ImageDimPlot
ImageDimPlot(
  Yersinia_subset_epi,
  fov = "yersinia1",
  group.by = "Progress_Highlight",
  cols = c("grey80", "red") 
)

ggsave( "Yersinia1_Highlighted_image_of peakSaa1_enrich.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

unique(Yersinia_subset$secondlevel)


#TNF Expression Across Distance
# Extract data for Tnf from Seurat
plot_data_tnf <- data.frame(
  progress = Yersinia_subset$progress,
  Tnf = FetchData(Yersinia_subset, "Tnf")[,1],
  Il1a = FetchData(Yersinia_subset, "Il1a")[,1],
  Il18 = FetchData(Yersinia_subset, "Il18")[,1],
  Il18r1 = FetchData(Yersinia_subset, "Il18r1")[,1],
  Il18rap = FetchData(Yersinia_subset, "Il18rap")[,1],
  Il22 = FetchData(Yersinia_subset, "Il22")[,1]
) %>%
  na.omit()

# Bin "progress" into intervals of 0.2 and compute average
avg_data_tnf <- plot_data_tnf %>%
  mutate(progress_bin = floor(progress / 0.2) * 0.2) %>%
  group_by(progress_bin) %>%
  summarise(
    mean_Tnf = mean(Tnf, na.rm = TRUE),
    mean_Il1a =  mean(Il1a, na.rm = TRUE),
    mean_Il18 =  mean(Il18, na.rm = TRUE),
    mean_Il18r1 =  mean(Il18r1, na.rm = TRUE),
    mean_Il18rap =  mean(Il18rap, na.rm = TRUE),
    mean_Il22 = mean(Il22, na.rm = TRUE),
    .groups = "drop"
  )

# Smooth using rolling median 
avg_data_tnf_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_Tnf_smooth = rollapply(mean_Tnf, width = 6, FUN = median, fill = NA, align = "center"))


# Plot
ggplot(avg_data_tnf_smooth, aes(x = progress_bin, y = mean_Tnf_smooth)) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color =  "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Median Yp Tnf Expression",
    title = "Smoothed Median Yp Tnf Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_Tnf_exp_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

avg_data_il1a_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_il1a_smooth = rollapply(mean_Il1a, width = 6, FUN = median, fill = NA, align = "center"))

ggplot(avg_data_il1a_smooth, aes(x = progress_bin, y = mean_il1a_smooth)) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color =  "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean IL1a Expression",
    title = "Smoothed  IL1a Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_IL1a_exp_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


avg_data_Il18_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_Il18_smooth = rollapply(mean_Il18, width = 6, FUN = median, fill = NA, align = "center"))

ggplot(avg_data_Il18_smooth, aes(x = progress_bin, y = mean_Il18_smooth)) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color =  "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean Il18 Expression",
    title = "Smoothed  Il18 Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_Il18_exp_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


avg_data_Il18r1_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_Il18r1_smooth = rollapply(mean_Il18r1, width = 6, FUN = median, fill = NA, align = "center"))

ggplot(avg_data_Il18r1_smooth, aes(x = progress_bin, y = mean_Il18r1_smooth)) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color = "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean Il18r1 Expression",
    title = "Smoothed  Il18r1 Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_Il18r1_exp_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



avg_data_Il18rap_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_Il18rap_smooth = rollapply(mean_Il18rap, width = 6, FUN = median, fill = NA, align = "center"))

ggplot(avg_data_Il18rap_smooth, aes(x = progress_bin, y = mean_Il18rap_smooth)) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color = "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean Il18rap Expression",
    title = "Smoothed  Il18rap Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_Il18rap_exp_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


avg_data_il22_smooth <- avg_data_tnf %>%
  arrange(progress_bin) %>%
  mutate(mean_Il22_smooth = rollapply(mean_Il22, width = 6, FUN = median, fill = NA, align = "center"))

ggplot(avg_data_il22_smooth, aes(x = progress_bin, y = mean_Il22_smooth)) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color = "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Mean IL22 Expression",
    title = "Smoothed  IL22 Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )


ggsave( "LineGraph_IL22_exp_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



# Plot the myeloid cells line graph
# Extract cell data from Seurat
plot_data_myeloid <- data.frame(
  progress = Yersinia_subset$progress,
  cell_type = Yersinia_subset$Coarse_SingleR  
) %>%
  na.omit()

unique(Yersinia_subset$Coarse_SingleR)
# Bin "progress" into intervals of 0.2 and compute percent myeloid cells
freq_data_myeloid <- plot_data_myeloid %>%
  mutate(progress_bin = floor(progress / 0.2) * 0.2) %>%
  group_by(progress_bin) %>%
  summarise(
    myeloid_percent = mean(cell_type %in% "Myeloid cells") * 100,
    .groups = "drop"
  )

# Smooth using rolling median
freq_data_myeloid_smooth <- freq_data_myeloid %>%
  arrange(progress_bin) %>%
  mutate(myeloid_percent_smooth = rollapply(myeloid_percent, width = 6, FUN = median, fill = NA, align = "center"))

ggplot(freq_data_myeloid_smooth, aes(x = progress_bin, y = myeloid_percent_smooth)) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.3) +
  geom_line(color = "#D17A00", linewidth = 1.5) +
  labs(
    x = "Progress",
    y = "Percent Myeloid Cells",
    title = "Smoothed Percent Myeloid Cells Across Progress"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray90")
  )
ggsave( "LineGraph_Percent_Myeloid_along_Distance_Yersinia1.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)




# Naive 2 Coordinate based analysis ----

# This sample is split into two ploygons for my purposes because of the quality of the swiss roll


# Combine them (Yersinia and Naive) into plots ----
#Sample 1 for both
library(zoo)
library(dplyr)
bin_and_smooth <- function(seu, feature, bin_width = 0.2, k = 6) {
  df <- data.frame(
    progress = seu$progress,
    value = FetchData(seu, feature)[,1]
  ) %>% na.omit()
  
  df %>%
    mutate(progress_bin = floor(progress / bin_width) * bin_width) %>%
    group_by(progress_bin) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(progress_bin) %>%
    mutate(smooth = rollapply(mean_val, width = k, FUN = mean, fill = NA, align = "center"))
}


#Set Highlight
highlight_regions <- data.frame(
  xmin = c(12, 28, 51),
  xmax = c(16, 32, 55),
  ymin = -Inf,
  ymax = Inf
)

unique(colnames(all_objects_filtered@meta.data))
#Saa1
dat_saa1_y <- bin_and_smooth(Yersinia_subset_epi, "Saa1_Gene_Lists") %>%
  mutate(Condition = "Yersinia")
dat_saa1_n <- bin_and_smooth(naive1_subset_epi, "Saa1_Gene_Lists") %>%
  mutate(Condition = "Naive")
dat_saa1 <- bind_rows(dat_saa1_y, dat_saa1_n )

summary(dat_saa1_y)
range(dat_saa1_y$mean_val, na.rm = TRUE)
range(dat_saa1_y$smooth, na.rm = TRUE)
saa1_values <- as.data.frame(all_objects_filtered$Saa1_Gene_Lists)

saa1_values %>% na.omit() %>% range()


#TNf_genelist
dat_tnf_enrich_y <- bin_and_smooth(Yersinia_subset, "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  mutate(Condition = "Yersinia")
dat_tnf_enrich_n <- bin_and_smooth(Naive_subset, "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  mutate(Condition = "Naive")
dat_tnf_enrich <- bind_rows(dat_tnf_enrich_y, dat_tnf_enrich_n )

#Ifng_Genelist
dat_ifng_enrich_y <- bin_and_smooth(Yersinia_subset, "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  mutate(Condition = "Yersinia")
dat_ifng_enrich_n <- bin_and_smooth(Naive_subset, "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  mutate(Condition = "Naive")
dat_ifng_enrich <- bind_rows(dat_ifng_enrich_y, dat_ifng_enrich_n )




#Tnf
dat_tnf_y <- bin_and_smooth(Yersinia_subset, "Tnf") %>%
  mutate(Condition = "Yersinia")
dat_tnf_n <- bin_and_smooth(Naive_subset, "Tnf") %>%
  mutate(Condition = "Naive")
dat_tnf <- bind_rows(dat_tnf_y, dat_tnf_n)

# Il1a
dat_il1a_y <- bin_and_smooth(Yersinia_subset, "Il1a") %>%
  mutate(Condition = "Yersinia")
dat_il1a_n <- bin_and_smooth(Naive_subset, "Il1a") %>%
  mutate(Condition = "Naive")
dat_il1a <- bind_rows(dat_il1a_y, dat_il1a_n)


# Il22
dat_il22_y <- bin_and_smooth(Yersinia_subset, "Il22") %>%
  mutate(Condition = "Yersinia")
dat_il22_n <- bin_and_smooth(Naive_subset, "Il22") %>%
  mutate(Condition = "Naive")
dat_il22 <- bind_rows(dat_il22_y, dat_il22_n)


# Il18r1
dat_il18r1_y <- bin_and_smooth(Yersinia_subset, "Il18r1") %>%
  mutate(Condition = "Yersinia")
dat_il18r1_n <- bin_and_smooth(Naive_subset, "Il18r1") %>%
  mutate(Condition = "Naive")
dat_il18r1 <- bind_rows(dat_il18r1_y, dat_il18r1_n)


# Il18rap
dat_Il18rap_y <- bin_and_smooth(Yersinia_subset, "Il18rap") %>%
  mutate(Condition = "Yersinia")
dat_Il18rap_n <- bin_and_smooth(Naive_subset, "Il18rap") %>%
  mutate(Condition = "Naive")
dat_Il18rap <- bind_rows(dat_Il18rap_y, dat_Il18rap_n)



# Il18
dat_il18_y <- bin_and_smooth(Yersinia_subset_epi, "Il18") %>%
  mutate(Condition = "Yersinia")
dat_il18_n <- bin_and_smooth(naive1_subset_epi, "Il18") %>%
  mutate(Condition = "Naive")
dat_il18 <- bind_rows(dat_il18_y, dat_il18_n)

# saa1
dat_saa1_y <- bin_and_smooth(Yersinia_subset_epi, "Saa1") %>%
  mutate(Condition = "Yersinia")
dat_saa1_n <- bin_and_smooth(naive1_subset_epi, "Saa1") %>%
  mutate(Condition = "Naive")
dat_saa1b <- bind_rows(dat_saa1_y, dat_saa1_n)


# Car4
dat_Car4_y <- bin_and_smooth(Yersinia_subset_epi, "Car4") %>%
  mutate(Condition = "Yersinia")
dat_Car4_n <- bin_and_smooth(naive1_subset_epi, "Car4") %>%
  mutate(Condition = "Naive")
dat_Car4 <- bind_rows(dat_Car4_y, dat_Car4_n)

# Nos2
dat_Hif1a1_y <- bin_and_smooth(Yersinia_subset_epi, "Hif1a") %>%
  mutate(Condition = "Yersinia")
dat_Hif1a1_n <- bin_and_smooth(naive1_subset_epi, "Hif1a") %>%
  mutate(Condition = "Naive")
dat_Hif1a <- bind_rows(dat_Hif1a1_y, dat_Hif1a1_n)

#percent myeloid 
compute_myeloid <- function(seu, bin_width = 0.2, k = 6) {
  df <- data.frame(
    progress = seu$progress,
    cell_type = seu$Coarse_SingleR
  ) %>% na.omit()
  
  df %>%
    mutate(progress_bin = floor(progress / bin_width) * bin_width) %>%
    group_by(progress_bin) %>%
    summarise(pct = mean(cell_type %in% "Myeloid cells") * 100, .groups = "drop") %>%
    arrange(progress_bin) %>%
    mutate(smooth = rollapply(pct, width = k, FUN = median, fill = NA, align = "center"))
}

dat_my_y <- compute_myeloid(Yersinia_subset) %>% mutate(Condition = "Yersinia")
dat_my_n <- compute_myeloid(Naive_subset) %>% mutate(Condition = "Naive")

dat_my <- bind_rows(dat_my_y, dat_my_n)

compute_t <- function(seu, bin_width = 0.2, k = 6) {
  df <- data.frame(
    progress = seu$progress,
    cell_type = seu$Coarse_SingleR
  ) %>% na.omit()
  
  df %>%
    mutate(progress_bin = floor(progress / bin_width) * bin_width) %>%
    group_by(progress_bin) %>%
    summarise(pct = mean(cell_type %in% "T cells and ILCs") * 100, .groups = "drop") %>%
    arrange(progress_bin) %>%
    mutate(smooth = rollapply(pct, width = k, FUN = median, fill = NA, align = "center"))
}
dat_t_y <- compute_t(Yersinia_subset) %>% mutate(Condition = "Yersinia")
dat_t_n <- compute_t(Naive_subset) %>% mutate(Condition = "Naive")

dat_t <- bind_rows(dat_t_y, dat_t_n)


#plot 
#Saa1
ggplot(dat_saa1, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +
  ylim(0.25,0.4)+ 
  geom_vline(xintercept = c(43, 62, 65, 95),
                              linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean Saa1 Enrichment") +
  theme(axis.text = element_text(color = "black"))
#Figure 7C
ggsave( "LineGraph_Yps_enrichment_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



#Enrichment IFNG 
ggplot(dat_ifng_enrich, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +
  geom_vline(xintercept = c(43, 62, 65, 95),
             linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Ifng Enrichmentt") +
  theme(axis.text = element_text(color = "black"))

ggsave( "LineGraph_Ifng_enrichment_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


#Enrichment TNFa 
ggplot(dat_tnf_enrich, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +
  geom_vline(xintercept = c(43, 62, 65, 95),
             linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed TNFa Enrichmentt") +
  theme(axis.text = element_text(color = "black"))

ggsave( "LineGraph_TNFa_enrichment_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)






# Il1a
ggplot(dat_il1a, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +  geom_vline(xintercept = c(43, 62, 65, 95),
                                linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean Il1a Expression") +
  theme(axis.text = element_text(color = "black"))
  #Figure 7B
ggsave( "LineGraph_Il1a_exp_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


ggplot(dat_il22, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +  geom_vline(xintercept = c(43, 62, 65, 95),
                                linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean Il22 Expression") +
  theme(axis.text = element_text(color = "black"))
  #Figure 7G
ggsave( "LineGraph_Il22_exp_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)




# il18r components
ggplot(dat_il18r1, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +  geom_vline(xintercept = c(43, 62, 65, 95),
                                linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean Il18r1 Expression") +
  theme(axis.text = element_text(color = "black"))
ggsave( "LineGraph_Il18r1_exp_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

ggplot(dat_Il18rap, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean Il18rap Expression") +
  theme(axis.text = element_text(color = "black"))
  #Figure 7E
ggsave( "LineGraph_Il18rap_exp_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



#tnf
ggplot(dat_tnf, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +  geom_vline(xintercept = c(43, 62, 65, 95),
                                linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean TNfa Expression") +
  theme(axis.text = element_text(color = "black"))
  #Figure 7B
ggsave( "LineGraph_TNFa_exp_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

#il18
ggplot(dat_il18, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +  geom_vline(xintercept = c(43, 62, 65, 95),
                                linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean Il18 Expression Epithelial") +
  theme(axis.text = element_text(color = "black"))
  #Figure 7E
ggsave( "LineGraph_Il18_exp_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


#saa1
ggplot(dat_saa1b, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +  geom_vline(xintercept = c(43, 62, 65, 95),
                                linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed Mean Saa1 Expression") +
  theme(axis.text = element_text(color = "black"))
  #figure 7I
ggsave( "LineGraph_Saa1_exp_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



#Myeloid
ggplot(dat_my, aes(x = progress_bin, y = smooth, color = Condition)) +
  geom_line(linewidth = 1.4) +
  geom_rect(data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  theme_minimal() +  geom_vline(xintercept = c(43, 62, 65, 95),
                                linetype = "dotted", linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("#B78CBA", "#D17A00"), breaks = c("Naive", "Yersinia"))+
  labs(x = "Progress",
       y = "Smoothed % Myeloid") +
  theme(axis.text = element_text(color = "black"))
#Figure 7B
ggsave( "LineGraph_Freq_Myeloid_along_Distance_Sample1_combined.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 1.5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)




#Examine the Spread of individual genes in the Saa1 gene list ----
#Are the genes that seem to be overlapping the PGS all uniformly distributed? If not, what is the local expression distribution of enterocyte-related Yps expression. 
library(GSVA)
library(GSEABase) 
Saa1_genes <- getGmt("/path/to/data/Genesets/Saa1_geneList_Aug14_2025_Xenium.gmt")
Saa1_genes <- Saa1_genes[["Saa1_Gene_Lists"]]
Saa1_genes <- Saa1_genes@geneIds

# Extract the data
expr_data <- FetchData(
  Yersinia_subset_epi,
  vars = c(Saa1_genes, "progress")
)

# Filter for the range containing the PG of interest
expr_data <- expr_data %>%
  filter(progress >= 45 & progress <= 61)

library(tidyr)
spread_stats <- expr_data %>%
  pivot_longer(cols = all_of(Saa1_genes), names_to = "gene", values_to = "expression") %>%
  group_by(gene) %>%
  summarise(
    center = sum(progress * expression) / sum(expression), # weighted mean
    spread = sqrt(sum(expression * (progress - center)^2) / sum(expression)), # weighted SD
    .groups = "drop"
  )


#Create a MAtrix
bin_size <- 0.2
heatmap_data <- expr_data %>%
  mutate(progress_bin = floor(progress / bin_size) * bin_size) %>%
  pivot_longer(cols = all_of(Saa1_genes), names_to = "gene", values_to = "expression") %>%
  group_by(gene, progress_bin) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = progress_bin, values_from = mean_expr, values_fill = 0)



library(pheatmap)

# Convert to matrix for pheatmap
mat <- as.matrix(heatmap_data[,-1])
rownames(mat) <- heatmap_data$gene

pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(50),
  main = "Expression along progress (45-61)"
)


heatmap_data_ordered <- heatmap_data %>%
  arrange(match(gene, spread_stats$gene[order(spread_stats$spread)]))

mat_ordered <- as.matrix(heatmap_data_ordered[,-1])
rownames(mat_ordered) <- heatmap_data_ordered$gene

pheatmap(
  mat_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(50),
  main = "Genes ordered by spread from peak", 
  border_color = NA
)

#Scaled
heatmap_data_ordered <- heatmap_data %>%
  arrange(match(gene, spread_stats$gene[order(spread_stats$spread)]))

# Convert to matrix and scale each row (gene)
mat_ordered <- as.matrix(heatmap_data_ordered[,-1])
rownames(mat_ordered) <- heatmap_data_ordered$gene

# Row-wise scaling: subtract mean, divide by SD for each gene
mat_scaled <- t(scale(t(mat_ordered)))  # t() tricks scale() into working per row

# Heatmap
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy",  "yellow"))(50),
  main = "Genes ordered by spread from peak (row scaled)", 
  border_color = NA
)
#Figure 7H
svg(paste0(images, "/Saa1_genes_heatmap_dispersion_along_PG.svg"), width = 4, height = 6)

# Plot heatmap
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy",  "yellow"))(50),
  main = "Genes ordered by spread from peak (row scaled)", 
  border_color = NA
)


dev.off()

# Examine the Spread of all genes - new genes which are locally defined to the PG within Epithelial cells? ----
# Extract the data
genes <- rownames(all_objects_filtered)
expr_data <- FetchData(
  Yersinia_subset_epi,
  vars = c(genes, "progress")
)

expr_data <- expr_data %>%
  filter(progress >= 45 & progress <= 61) # Filter for the range containing the PG of interest

spread_stats <- expr_data %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  group_by(gene) %>%
  summarise(
    center = sum(progress * expression) / sum(expression), # weighted mean
    spread = sqrt(sum(expression * (progress - center)^2) / sum(expression)), # weighted SD
    .groups = "drop"
  )


#Create a MAtrix
bin_size <- 0.2
heatmap_data <- expr_data %>%
  mutate(progress_bin = floor(progress / bin_size) * bin_size) %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  group_by(gene, progress_bin) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = progress_bin, values_from = mean_expr, values_fill = 0)

heatmap_data <- heatmap_data[rowMeans(heatmap_data[,2:ncol(heatmap_data)]) > 0.2, ]

# Convert to matrix for pheatmap
mat <- as.matrix(heatmap_data[,-1])
rownames(mat) <- heatmap_data$gene

pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(50),
  main = "Expression along progress (45-61)"
)


heatmap_data_ordered <- heatmap_data %>%
  arrange(match(gene, spread_stats$gene[order(spread_stats$spread)]))

mat_ordered <- as.matrix(heatmap_data_ordered[,-1])
rownames(mat_ordered) <- heatmap_data_ordered$gene

pheatmap(
  mat_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(50),
  main = "Genes ordered by spread from peak", 
  border_color = NA
)

#Scaled
heatmap_data_ordered <- heatmap_data %>%
  arrange(match(gene, spread_stats$gene[order(spread_stats$spread)]))

# Convert to matrix and scale each row (gene)
mat_ordered <- as.matrix(heatmap_data_ordered[1:200,-1])
rownames(mat_ordered) <- heatmap_data_ordered[1:200,]$gene

# Row-wise scaling: subtract mean, divide by SD for each gene
mat_scaled <- t(scale(t(mat_ordered)))  

# Heatmap
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy",  "yellow"))(50),
  main = "Genes ordered by spread from peak (row scaled)", 
  border_color = NA
)

svg(paste0(images, "/Saa1_genes_heatmap_dispersion_along_PG.svg"), width = 4, height = 6)
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy",  "yellow"))(50),
  main = "Genes ordered by spread from peak (row scaled)", 
  border_color = NA
)
dev.off()


# Myeloid Cell Further Annotation ----
Idents(Ileum_annotated) <- "CoarseCellType"
unique(Ileum_annotated$CoarseCellType)
Ileum_Myeloid <- subset(Ileum_annotated, idents = "Myeloid cells")
unique(Ileum_Myeloid$FineCellType)
scRNA_counts <- GetAssayData(Ileum_Myeloid, assay = "RNA", slot = "data")
scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(subset_objects[["Myeloid cells"]]@assays$SCT), ]
FineCellType_vector <- Ileum_Myeloid$secondlevel
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)
unique(Ileum_Myeloid$FineCellType)
# Extract Xenium info
subset_objects[["Myeloid cells"]]@assays$SCT <- JoinLayers(subset_objects[["Myeloid cells"]]@assays$SCT)
MyeloidMarkers <- FindAllMarkers(subset_objects[["Myeloid cells"]], group.by = "Monocle_clusters", test.use = "wilcox", 
                                 only.pos = TRUE, min.pct = 0.10, assay = "SCT")
xenium_counts <- GetAssayData(subset_objects[["Myeloid cells"]], assay = "SCT", slot = "data") 
DimPlot(subset_objects[["Myeloid cells"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()
DimPlot(subset_objects[["Myeloid cells"]], reduction = "umap", group.by = "Myeloid_cells_clusters" , label = T)  + NoLegend()
FeaturePlot(subset_objects[["Myeloid cells"]], reduction = "umap", c("Itgam", "Itgax", "Sirpa", "Xcr1", "H2-DMb2"))
FeaturePlot(subset_objects[["Myeloid cells"]], reduction = "umap", c("Epcam", "Cd3e", "Cd19", "Lum", "Tagln"))
FeaturePlot(subset_objects[["Myeloid cells"]], reduction = "umap", c("Plvap", "Eng", "Lyve1", "Lum", "Mki67"))
FeaturePlot(subset_objects[["Myeloid cells"]], reduction = "umap", c("Jchain", "Igha", "Igkc", "Fgfr2" ))

# Missegmentation is a problem. In the Myeloid cells, there are clusters formed by contaminating transcript leakage (misassigned transcripts)  - 
# Cluster 16 is tagln hi - muscle cells
# Cluster 5 is high in Lum - a fibroblast marker
# cluster 14 is Cd3e high for T cell links
# cluster 9 Eng and lvap high endothelial cells contaminant
# Cluster 11 - Clearly Cdc1 - Xcr1 positive
# Cluster 17 lyve1, tagn lum - soe mix of stromal contaminating transcripts
# cluster 12 is cycling by Mki67
# cluster 15 is plasma cells 

# Run SingleR
library(BiocParallel)
subset_objects[["Myeloid cells"]] <-  split(x = subset_objects[["Myeloid cells"]], f = subset_objects[["Myeloid cells"]]$orig.ident) 

singleR_results <- SingleR(test = xenium_counts, ref = reference$counts, labels = reference$labels, 
                           de.method="wilcox", 
                           BPPARAM=MulticoreParam(25))


subset_objects[["Myeloid cells"]]$SingleR_Fine_Scores <- singleR_results$scores
FeaturePlot_scCustom(subset_objects[["Myeloid cells"]], reduction =  "umap", "SingleR_Fine_Scores", na_cutoff = 0)

# Extract cell type predictions (e.g., the most likely cell type for each cell)
subset_objects[["Myeloid cells"]]$FineLabelsSingleR <- singleR_results$pruned.labels
DimPlot(subset_objects[["Myeloid cells"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of Myeloid cells")) + NoLegend()

FeaturePlot(subset_objects[["Myeloid cells"]], reduction = "umap", c("Ly6g", "Cxcr2", "Xcr1", "Ccr2", "Adgre1", "Itgax"))
table(subset_objects[["Myeloid cells"]]$InfectionStatus, subset_objects[["Myeloid cells"]]$FineLabelsSingleR)
unique(subset_objects[["Myeloid cells"]]$FineLabelsSingleR)



#Doesnt work 
#ImageDimPlot(subset_objects[["Myeloid cells"]], group.by = "FineLabelsSingleR", fov = "Naive1")

remove(Ileum_Myeloid)

umap_data <- subset_objects[["Myeloid cells"]][["umap"]]@cell.embeddings
cluster_data <- subset_objects[["Myeloid cells"]]$FineLabelsSingleR
infection <- subset_objects[["Myeloid cells"]]$InfectionStatus
umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cluster_data)
umap_data$infection <- factor(infection)
centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(umap_1 = mean(umap_1), umap_2 = mean(umap_2))

cluster_freq <- umap_data %>%
  group_by(infection, cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(infection) %>%
  mutate(percent = count / sum(count) * 100)

ggplot(cluster_freq, aes(x = cluster, y = count, fill = infection)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), color = "black") +
  #facet_wrap(~ infection) +
  scale_fill_manual(values = cbf_18_colors) +
  labs(x = "Cluster", y = "Percentage of Cells") +
  #scale_y_break(c(25, 60)) +  
  theme_minimal() +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background = element_blank(), axis.text.x = element_text(angle = 90),
    legend.position = "none"
  )

ggsave( "Barchart_spatial_Myeloid_CoarseCellType.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)






#T cells and ILCs further annotation ----
library(BiocParallel)
Idents(Ileum_annotated) <- "CoarseCellType"

unique(Ileum_annotated$FineCellType)
Ileum_T <- subset(Ileum_annotated, idents = "T cells and ILCs")
unique(Ileum_T$secondlevel)
scRNA_counts <- GetAssayData(Ileum_T, assay = "RNA", slot = "data")
scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(subset_objects[["T cells and ILCs"]]@assays$SCT), ]
FineCellType_vector <- Ileum_T$secondlevel
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

# Extract Xenium info
subset_objects[["T cells and ILCs"]]@assays$SCT <- JoinLayers(subset_objects[["T cells and ILCs"]]@assays$SCT)
xenium_counts <- GetAssayData(subset_objects[["T cells and ILCs"]], assay = "SCT", slot = "data") 
DimPlot(subset_objects[["T cells and ILCs"]], reduction = "umap", group.by = "T_cells_and_ILCs_clusters" , label = T)  + NoLegend()

# Run SingleR
#subset_objects[["T cells and ILCs"]] <- JoinLayers(subset_objects[["T cells and ILCs"]] )
subset_objects[["T cells and ILCs"]] <-  split(x = subset_objects[["T cells and ILCs"]], f = subset_objects[["T cells and ILCs"]]$orig.ident) 

singleR_results <- SingleR(test = xenium_counts, ref = reference$counts, labels = reference$labels, 
                           de.method="wilcox",
                           BPPARAM=MulticoreParam(25))


subset_objects[["T cells and ILCs"]]$SingleR_Fine_Scores <- singleR_results$scores
FeaturePlot_scCustom(subset_objects[["T cells and ILCs"]], reduction =  "umap", "SingleR_Fine_Scores", na_cutoff = 0)

# Extract cell type predictions 
subset_objects[["T cells and ILCs"]]$FineLabelsSingleR <- singleR_results$pruned.labels
DimPlot(subset_objects[["T cells and ILCs"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T) + ggtitle(paste("UMAP of T cells and ILCs")) + NoLegend()

FeaturePlot(subset_objects[["T cells and ILCs"]], reduction = "umap", features  = "Mki67")

unique(subset_objects[["T cells and ILCs"]]$FineLabelsSingleR)
remove(Ileum_T)

DimPlot(subset_objects[["T cells and ILCs"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of T cells and ILCs")) + NoLegend()
FeaturePlot(subset_objects[["T cells and ILCs"]], reduction = "umap", c("Rorc", "Ncam1", "Ccr6", "Ncr1", "Gata3", "Sell"))
table(subset_objects[["T cells and ILCs"]]$FineLabelsSingleR )




#Stromal Cells further annotation ----
Idents(Ileum_annotated) <- "CoarseCellType"
unique(Ileum_annotated$CoarseCellType)
Ileum_Stromal <- subset(Ileum_annotated, idents = "Stromal cells")
unique(Ileum_Stromal$FineCellType)
scRNA_counts <- GetAssayData(Ileum_Stromal, assay = "RNA", slot = "data")
scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(subset_objects[["Stromal cells"]]@assays$SCT), ]
FineCellType_vector <- Ileum_Stromal$secondlevel
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

# Extract Xenium info
#subset_objects[["Stromal cells"]]@assays$SCT <- JoinLayers(subset_objects[["Stromal cells"]]@assays$SCT)
xenium_counts <- GetAssayData(subset_objects[["Stromal cells"]], assay = "SCT", slot = "data") 
DimPlot(subset_objects[["Stromal cells"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()

# Run SingleR
subset_objects[["Stromal cells"]] <-  split(x = subset_objects[["Stromal cells"]], f = subset_objects[["Stromal cells"]]$orig.ident) 

singleR_results <- SingleR(test = xenium_counts, ref = reference$counts, labels = reference$labels, 
                           de.method="wilcox",
                           BPPARAM=MulticoreParam(25))


subset_objects[["Stromal cells"]]$SingleR_Fine_Scores <- singleR_results$scores
FeaturePlot_scCustom(subset_objects[["Stromal cells"]], reduction =  "umap", "SingleR_Fine_Scores", na_cutoff = 0)

# Extract cell type predictions (e.g., the most likely cell type for each cell)
subset_objects[["Stromal cells"]]$FineLabelsSingleR <- singleR_results$pruned.labels
DimPlot(subset_objects[["Stromal cells"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of Stromal cells")) + NoLegend()
table(subset_objects[["Stromal cells"]]$InfectionStatus, subset_objects[["Stromal cells"]]$FineLabelsSingleR)
unique(subset_objects[["Stromal cells"]]$FineLabelsSingleR)
remove(Ileum_Stromal)





#B cells annotation ----

Idents(Ileum_annotated) <- "CoarseCellType"
unique(Ileum_annotated$CoarseCellType)
Ileum_B <- subset(Ileum_annotated, idents = "B cells")
unique(Ileum_B$secondlevel)
scRNA_counts <- GetAssayData(Ileum_B, assay = "RNA", slot = "data")
scRNA_counts <- scRNA_counts[rownames(scRNA_counts) %in% rownames(subset_objects[["B cells"]]@assays$SCT), ]
FineCellType_vector <- Ileum_B$secondlevel
reference <- list(counts = scRNA_counts, labels = FineCellType_vector)

# Extract Xenium info
#subset_objects[["B cells"]]@assays$SCT <- JoinLayers(subset_objects[["B cells"]]@assays$SCT)
xenium_counts <- GetAssayData(subset_objects[["B cells"]], assay = "SCT", slot = "data") 
DimPlot(subset_objects[["B cells"]], reduction = "umap", group.by = "Monocle_clusters" , label = T)  + NoLegend()

# Run SingleR
subset_objects[["B cells"]] <-  split(x = subset_objects[["B cells"]], f = subset_objects[["B cells"]]$orig.ident) 

singleR_results <- SingleR(test = xenium_counts, ref = reference$counts, labels = reference$labels, 
                           de.method="wilcox",
                           BPPARAM=MulticoreParam(25))


subset_objects[["B cells"]]$SingleR_Fine_Scores <- singleR_results$scores
FeaturePlot_scCustom(subset_objects[["B cells"]], reduction =  "umap", "SingleR_Fine_Scores", na_cutoff = 0)

# Extract cell type predictions (e.g., the most likely cell type for each cell)
subset_objects[["B cells"]]$FineLabelsSingleR <- singleR_results$pruned.labels
DimPlot(subset_objects[["B cells"]], reduction = "umap", group.by = "FineLabelsSingleR" , label = T, repel = T) + ggtitle(paste("UMAP of B cells")) + NoLegend()
FeaturePlot(subset_objects[["B cells"]], reduction = "umap", c("Cd3e", "Trac", "Epcam", "Lum", "Itgax"))


ImageDimPlot(subset_objects[["B cells"]], group.by = "FineLabelsSingleR", fov = "Macp2")
table(subset_objects[["B cells"]]$InfectionStatus, subset_objects[["B cells"]]$FineLabelsSingleR)
unique(subset_objects[["B cells"]]$FineLabelsSingleR)
remove(Ileum_B)




#Enteric Nervous System----

#These cells were not numerous enought to adequately quantify in the prior data sets. Here I will begin by defining major Gene expression changes
DimPlot(subset_objects[["Enteric Nervous System"]], reduction = "umap", group.by = "Monocle_clusters", label = T, repel = T ) + ggtitle(paste("UMAP of Enteric Nervous System")) + NoLegend()
Idents(subset_objects[["Enteric Nervous System"]]) <- "Monocle_clusters"
markers_all <- FindAllMarkers(subset_objects[["Enteric Nervous System"]], only.pos = TRUE, logfc.threshold = 0.25, assay = "SCT")

top_genes <- markers_all %>%
  group_by(cluster) %>%  # Group by cluster
  arrange(cluster, p_val_adj) %>%  # Arrange by cluster and p-value
  slice_head(n = 5) %>%  # Select the top 5 genes per cluster
  pull(gene)

top_genes
top_genes <- unique(top_genes)

dotplot <- DotPlot(subset_objects[["Enteric Nervous System"]], features = top_genes, ) +
  theme_minimal() +  # Use a minimal theme
  scale_color_viridis_c() +  # Attractive color scale
  ggtitle("Top 5 Differential Genes Across Clusters") +  # Title for the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability



FeaturePlot(subset_objects[["Enteric Nervous System"]], reduction = "umap", c("Nefl",
                                                                              "Tubb3",
                                                                              "Syn1",
                                                                              "Sox10",
                                                                              "Chat",
                                                                              "Pdgfra", 
                                                                              "Scn7a", "Plp1"))



FeaturePlot(subset_objects[["Enteric Nervous System"]], reduction = "umap", c("Lum",
                                                                              "Pecam1",
                                                                              "Actg2",
                                                                              "Tagln",
                                                                              "Myh11",
                                                                              "Epcam"))

# Display the plot
dotplot


# Attach the Preliminary Assignments to each Cell Type in the whole data set ----
all_objects_filtered <- AddMetaData(all_objects_filtered, subset_objects[["Enteric Nervous System"]]$Coarse_SingleR , col.name = "secondlevel")
all_objects_filtered <- AddMetaData(all_objects_filtered, subset_objects[["B cells"]]$FineLabelsSingleR , col.name = "secondlevel")
all_objects_filtered <- AddMetaData(all_objects_filtered, subset_objects[["Stromal cells"]]$FineLabelsSingleR, col.name = "secondlevel" )
all_objects_filtered <- AddMetaData(all_objects_filtered, subset_objects[["T cells and ILCs"]]$FineLabelsSingleR , col.name = "secondlevel")
all_objects_filtered <- AddMetaData(all_objects_filtered, subset_objects[["Myeloid cells"]]$FineLabelsSingleR , col.name = "secondlevel")
all_objects_filtered <- AddMetaData(all_objects_filtered, subset_objects[["Epithelial cells"]]$FineLabelsSingleR , col.name = "secondlevel")

DimPlot(all_objects_filtered, reduction = "umap", group.by = "secondlevel", label = T, repel = T) + NoLegend()



# Prep objects for CellChat ----
# I will approach the analysis by comparing separate regions. one from naive sample and two from Yersinia. 

# Naive Coordinate Based Analysis 
naive1_outer <- fread(paste0(CSV, "naive1_outer_sep10.csv"))
naive1_outer <- naive1_outer[, c(1,2)]
naive1_outer$x <- as.numeric(naive1_outer$x)
naive1_outer$y <- as.numeric(naive1_outer$y)
naive1_inner <- fread(paste0(CSV, "naive1_inner_sep10.csv"))
naive1_inner <- naive1_inner[, c(1,2)]
naive1_inner$x <- as.numeric(naive1_inner$x)
naive1_inner$y <- as.numeric(naive1_inner$y)

library(sf)
# Make sure both are in same order along ribbon

ribbon_coords <- rbind(
  naive1_outer,
  naive1_inner[nrow(naive1_inner):1, ], # reverse inner coords
  naive1_outer[1, ]  # close polygon by repeating first point
)

ribbon_polygon <- st_polygon(list(as.matrix(ribbon_coords))) %>%
  st_sfc(crs = NA_crs_) %>%   # use NA_crs_ for undefined Cartesian CRS
  st_make_valid()

query_df <- GetTissueCoordinates(all_objects_filtered, image = "naive1")
query_points <- st_as_sf(query_df, coords = c("x", "y"), crs = NA_crs_)

# Find points inside polygon
inside_idx <- st_within(query_points, ribbon_polygon, sparse = FALSE)[,1]

# Subset coordinates
inside_cells <- query_df[inside_idx, ]


Naive_subset <- subset(all_objects_filtered, cells = inside_cells$cell)

# Make centerline
n_inner <- nrow(naive1_inner)
n_outer <- nrow(naive1_outer)
n_min <- min(n_inner, n_outer)

# Match up shorter spiral to longer one
inner_resampled <- naive1_inner[round(seq(1, n_inner, length.out = n_min)), ]
outer_resampled <- naive1_outer[round(seq(1, n_outer, length.out = n_min)), ]

centerline <- data.frame(
  x = (inner_resampled$x + outer_resampled$x) / 2,
  y = (inner_resampled$y + outer_resampled$y) / 2
)

# Cumulative distances along centerline
distances <- sqrt(diff(centerline$x)^2 + diff(centerline$y)^2)
cumdist <- c(0, cumsum(distances))
total_dist <- max(cumdist)

# Interpolate to get, say, 10x more points along the centerline
interp_factor <- 50
interp_idx <- seq(0, total_dist, length.out = n_min * interp_factor)

# Cumulative distance already computed: cumdist
centerline_interp <- data.frame(
  x = approx(cumdist, centerline$x, xout = interp_idx)$y,
  y = approx(cumdist, centerline$y, xout = interp_idx)$y
)

# Update distances
distances_interp <- sqrt(diff(centerline_interp$x)^2 + diff(centerline_interp$y)^2)
cumdist_interp <- c(0, cumsum(distances_interp))
total_dist_interp <- max(cumdist_interp)

assign_progress <- function(px, py) {
  dists_to_center <- sqrt((centerline_interp$x - px)^2 + (centerline_interp$y - py)^2)
  nearest_idx <- which.min(dists_to_center)
  progress <- (cumdist_interp[nearest_idx] / total_dist_interp) * 100
  return(progress)
}

inside_cells$progress <- mapply(assign_progress, inside_cells$x, inside_cells$y)

#Save the new distance to the seurat object
Naive_subset$progress <- inside_cells$progress[match(Cells(Naive_subset), inside_cells$cell)]


ImageFeaturePlot(Naive_subset, "progress")
ImageDimPlot(Naive_subset,  group.by = "Coarse_SingleR" )

#Separate out a portion of the swiss roll for CellCellCommunication analysis
Naive_cc <- subset(Naive_subset, progress > 30 & progress < 60)
ImageFeaturePlot(Naive_cc, "progress")
table(Naive_cc$secondlevel)
remove(Naive_subset)



#Yersinia Coordinate based analysis
yp1_outer <- fread(paste0(CSV, "yp1_outer_sep2.csv"))
yp1_outer <- yp1_outer[, c(1,2)]
yp1_outer$x <- as.numeric(yp1_outer$x)
yp1_outer$y <- as.numeric(yp1_outer$y)
yp1_inner <- fread(paste0(CSV, "yp1_inner_sep2.csv"))
yp1_inner <- yp1_inner[, c(1,2)]
yp1_inner$x <- as.numeric(yp1_inner$x)
yp1_inner$y <- as.numeric(yp1_inner$y)

library(sf)
# Make sure both are in same order along ribbon

ribbon_coords <- rbind(
  yp1_outer,
  yp1_inner[nrow(yp1_inner):1, ], # reverse inner coords
  yp1_outer[1, ]  # close polygon by repeating first point
)

ribbon_polygon <- st_polygon(list(as.matrix(ribbon_coords))) %>%
  st_sfc(crs = NA_crs_) %>%   # use NA_crs_ for undefined Cartesian CRS
  st_make_valid()

query_df <- GetTissueCoordinates(all_objects_filtered, image = "yersinia1")
query_points <- st_as_sf(query_df, coords = c("x", "y"), crs = NA_crs_)

# Find points inside polygon
inside_idx <- st_within(query_points, ribbon_polygon, sparse = FALSE)[,1]

# Subset coordinates
inside_cells <- query_df[inside_idx, ]


Yersinia_subset <- subset(all_objects_filtered, cells = inside_cells$cell)

# Make centerline
n_inner <- nrow(yp1_inner)
n_outer <- nrow(yp1_outer)
n_min <- min(n_inner, n_outer)

# Match up shorter spiral to longer one
inner_resampled <- yp1_inner[round(seq(1, n_inner, length.out = n_min)), ]
outer_resampled <- yp1_outer[round(seq(1, n_outer, length.out = n_min)), ]

centerline <- data.frame(
  x = (inner_resampled$x + outer_resampled$x) / 2,
  y = (inner_resampled$y + outer_resampled$y) / 2
)

# Cumulative distances along centerline
distances <- sqrt(diff(centerline$x)^2 + diff(centerline$y)^2)
cumdist <- c(0, cumsum(distances))
total_dist <- max(cumdist)


interp_factor <- 50
interp_idx <- seq(0, total_dist, length.out = n_min * interp_factor)

# Cumulative distance already computed: cumdist
centerline_interp <- data.frame(
  x = approx(cumdist, centerline$x, xout = interp_idx)$y,
  y = approx(cumdist, centerline$y, xout = interp_idx)$y
)

# Update distances
distances_interp <- sqrt(diff(centerline_interp$x)^2 + diff(centerline_interp$y)^2)
cumdist_interp <- c(0, cumsum(distances_interp))
total_dist_interp <- max(cumdist_interp)

assign_progress <- function(px, py) {
  dists_to_center <- sqrt((centerline_interp$x - px)^2 + (centerline_interp$y - py)^2)
  nearest_idx <- which.min(dists_to_center)
  progress <- (cumdist_interp[nearest_idx] / total_dist_interp) * 100
  return(progress)
}

inside_cells$progress <- mapply(assign_progress, inside_cells$x, inside_cells$y)

#Save the new distance to the seurat object
Yersinia_subset$progress <- inside_cells$progress[match(Cells(Yersinia_subset), inside_cells$cell)]


ImageFeaturePlot(Yersinia_subset, "progress")
Yersinia_cc <- subset(Yersinia_subset, progress > 65 & progress < 95)
ImageFeaturePlot(Yersinia_cc, "progress")
table(Yersinia_cc$secondlevel)

#PGcells <- WhichCells(Yersinia_subset, progress > 43 & progress < 62)
Yersinia_pg <- subset(Yersinia_subset, progress > 43 & progress < 62)
ImageFeaturePlot(Yersinia_pg, "progress")
table(Yersinia_pg$secondlevel)


Yersinia_subset$Progress_Class <- ifelse(
  Yersinia_subset$progress >= 43 & Yersinia_subset$progress <= 62, "red",
  ifelse(Yersinia_subset$progress >= 65 & Yersinia_subset$progress <= 95, "blue", "grey")
)

# Check the counts
table(Yersinia_subset$Progress_Class)

# Now you can use it for coloring in ImageDimPlot
ImageDimPlot(
  Yersinia_subset,
  fov = "yersinia1",
  group.by = "Progress_Class",
  cols = c("#D17A00", "navy", "grey") 
)

ggsave( "Yersinia1_Highlighted_image_of_CCC_areas_PG_nonPG.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 6,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)




#Perform CEllChat analysis -----

#all_objects_filtered <- qs_read( "/path/to/data/MIST/Xenium_Yp_2025/Seurat/Seurat_all_proseg_fov.qs2")
#Idents(all_objects_filtered) <- "orig.ident"
#naive1 <- subset(all_objects_filtered, idents = c("naive1"))
#yersinia1 <- subset(all_objects_filtered, idents = c("yersinia1"))

table(Yersinia_pg$secondlevel)
table(Yersinia_cc$secondlevel)
table(Naive_cc$secondlevel)

anyNA(Yersinia_pg@meta.data$secondlevel)
anyNA(rownames(Yersinia_pg))
anyNA(Naive_cc@meta.data$secondlevel)
anyNA(rownames(Naive_cc))
anyNA(Yersinia_cc@meta.data$secondlevel)
table(Yersinia_pg@meta.data$secondlevel)
Yersinia_pg$secondlevel[is.na(Yersinia_pg@meta.data$secondlevel)] <- "Unknown"
Yersinia_cc$secondlevel[is.na(Yersinia_cc@meta.data$secondlevel)] <- "Unknown"
Naive_cc$secondlevel[is.na(Naive_cc@meta.data$secondlevel)] <- "Unknown"

yp_pg_celltypes <- unique(Yersinia_pg$secondlevel)
yp_cc_celltypes <- unique(Yersinia_cc$secondlevel)
na_cc_celltypes <- unique(Naive_cc$secondlevel)
celltypes <- Reduce(intersect, list(yp_pg_celltypes,yp_cc_celltypes, na_cc_celltypes ))
celltypes <- setdiff(celltypes, "Unknown")

Idents(Yersinia_pg) <- "secondlevel"
Yersinia_pg <- subset(Yersinia_pg, idents = celltypes)
Idents(Naive_cc) <- "secondlevel"
Naive_cc <- subset(Naive_cc, idents = celltypes)
Idents(Yersinia_cc) <- "secondlevel"
Yersinia_cc <- subset(Yersinia_cc, idents = celltypes)


Yersinia_pg$sample <- "Yersinia_pg"
Yersinia_cc$sample <- "Yersinia_cc"
data.input1 = Seurat::GetAssayData(Yersinia_pg, slot = "data", assay = "SCT" ) 
data.input2 = Seurat::GetAssayData(Yersinia_cc, slot = "data", assay = "SCT") 
data.input <- cbind(data.input1, data.input2)

genes.common <- intersect(rownames(data.input1), rownames(data.input2))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ])

meta1 = data.frame(Yersinia_pg@meta.data) # manually create a dataframe consisting of the cell labels
meta2 = data.frame(Yersinia_cc@meta.data) 
meta <- rbind(meta1, meta2)
rownames(meta) <- colnames(data.input)

# a factor level should be defined for the `meta$labels` and `meta$samples`
meta$labels <- factor(meta$secondlevel)
meta$samples <- factor(meta$sample, levels = c("Yersinia_pg", "Yersinia_cc"))

spatial.locs1 <-  GetTissueCoordinates(Yersinia_pg) 
spatial.locs2<-  GetTissueCoordinates(Yersinia_cc) 
spatial.locs <- rbind(spatial.locs1, spatial.locs2)
spatial.locs <- spatial.locs[,1:2]
# If spatial.locs is a data.frame, convert it:
spatial.locs <- as.matrix(spatial.locs)

rownames(spatial.locs) <- colnames(data.input)


spatial_factors <- data.frame(
  ratio = c(0.2125, 0.2125),  
  tol = c(32.5, 32.5) 
  #ratio = c(0.566, 0.653),  
)
rownames(spatial_factors) <- c("Yersinia_pg", "Yersinia_cc")

library(CellChat)
print(table(meta$labels))
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors =spatial_factors )

CellChatDB <- CellChatDB.mouse 
unique(CellChatDB$interaction$annotation)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact"), key = "annotation") # use Secreted Signaling
# set the used database in the object


library(CellChat)
db <- CellChatDB.use
# Create a new interaction row mimicking CellChat structure
custom_IL22 <- data.frame(
  interaction_name = "IL22_IL22RA1",
  pathway_name = "IL10",
  ligand = "Il22",
  receptor = "Il22ra1",
  agonist = "",
  antagonist = "",
  co_A_receptor = "",
  co_I_receptor = "",
  evidence = "UserDefined",
  annotation = "Secreted Signaling",
  interaction_name_2 = "Il22 - Il22ra1",
  is_neurotransmitter = FALSE,
  ligand.symbol = "Il22",
  ligand.family = "IL-10",
  ligand.location = "Secreted",
  ligand.keyword = "Signal, Reference proteome, Disulfide bond",
  ligand.secreted_type = "cytokine",
  ligand.transmembrane = FALSE,
  receptor.symbol = "Il22ra1",
  receptor.family = "Type II cytokine receptor",
  receptor.location = "Membrane",
  receptor.keyword = "Membrane, Signal, Reference proteome",
  receptor.surfaceome_main = "Receptors",
  receptor.surfaceome_sub = "CytokineR;IG;Type2",
  receptor.adhesome = "",
  receptor.secreted_type = "",
  receptor.transmembrane = TRUE,
  version = "Custom_v1",
  stringsAsFactors = FALSE
)
# Bind to interaction table
CellChatDB.use$interaction <- rbind(db$interaction, custom_IL22)
number_rownames <- length(rownames(CellChatDB.use$interaction)) - 1
rownames(CellChatDB.use$interaction) <- c(rownames(CellChatDB.use$interaction)[1:number_rownames], 
                                          "IL22_IL22RA1")


cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Part 2 inference of Cell Cell Communication network  - Importantly, spatial restriction was not applied because I had separated the objects into regions of interest already as distinct "samples"
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 150, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100
)

cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)


#Calculate aggregated cell chat network 
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

pathways.show <- c( "IL1" ) 
# Circle plot
unique(cellchat@netP$pathways)
par(mfrow=c(1,1), xpd=TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "Yersinia_pg", layout = "spatial", edge.width.max = 2, vertex.size.max = 0.5, alpha.image = 0.01, vertex.label.cex = 0)

netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "Yersinia_cc", layout = "spatial", edge.width.max = 2, vertex.size.max = 0.5, alpha.image = 0.01, vertex.label.cex = 0)

# Extract all weights for the datasets
weights_list <- lapply(c("Yersinia_pg", "Yersinia_cc"), function(sample) {
  w <- cellchat@netP$weight[[sample]]
  w[is.na(w)] <- 0   # replace NA with 0
  as.vector(w)
})

# Find the global maximum weight
global_max <- max(unlist(weights_list), na.rm = TRUE)

# Scale each dataset by the global maximum (0-1)
cellchat@netP$weight[["Yersinia_pg"]] <- cellchat@netP$weight[["Yersinia_pg"]] / global_max
cellchat@netP$weight[["Yersinia_cc"]] <- cellchat@netP$weight[["Yersinia_cc"]] / global_max

# Plot
svg(paste0(images, "/TNF_Spatial_signals_plot_Yersinia_PG.svg"), 
    width = 12, height = 19)

netVisual_aggregate(cellchat, signaling = "TNF", sample.use = "Yersinia_pg",
                    layout = "spatial", edge.width.max = 2,  # max width in plot
                    vertex.size.max = 0.5, alpha.image = 0.01,
                    vertex.label.cex = 0)

dev.off()

svg(paste0(images, "/TNF_Spatial_signals_plot_Yersinia_nonPG.svg"), 
    width = 12, height = 19)
netVisual_aggregate(cellchat, signaling = "TNF", sample.use = "Yersinia_cc",
                    layout = "spatial", edge.width.max = 2,
                    vertex.size.max = 0.5, alpha.image = 0.01,
                    vertex.label.cex = 0)
dev.off()



#Run cell chat individually on each objects 


spatial.locs1 <- spatial.locs1[,1:2]
spatial.locs1 <- as.matrix(spatial.locs1)# If spatial.locs is a data.frame, convert it:
rownames(spatial.locs1) <- colnames(data.input1)
cellchat_pg <- createCellChat(object = data.input1, meta = meta1, group.by = "secondlevel",
                           datatype = "spatial", coordinates = spatial.locs1, spatial.factors =spatial_factors[1,] )

cellchat_pg@DB <- CellChatDB.use
cellchat_pg <- subsetData(cellchat_pg) 
cellchat_pg <- identifyOverExpressedGenes(cellchat_pg)
cellchat_pg <- identifyOverExpressedInteractions(cellchat_pg)

cellchat_pg <- computeCommunProb(cellchat_pg, type = "triMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL, contact.dependent = TRUE, contact.range = 100)

cellchat_pg <- filterCommunication(cellchat_pg, min.cells = 1)
cellchat_pg.df <- subsetCommunication(cellchat_pg) 
cellchat_pg <- computeCommunProbPathway(cellchat_pg)
cellchat_pg <- aggregateNet(cellchat_pg)
cellchat_pg <- netAnalysis_computeCentrality(cellchat_pg, slot.name = "netP")

spatial.locs2 <- spatial.locs2[,1:2]
spatial.locs2 <- as.matrix(spatial.locs2)
rownames(spatial.locs2) <- colnames(data.input2)
cellchat_cc <- createCellChat(object = data.input2, meta = meta2, group.by = "secondlevel",
                              datatype = "spatial", coordinates = spatial.locs2, spatial.factors =spatial_factors[2,] )

cellchat_cc@DB <- CellChatDB.use
cellchat_cc <- subsetData(cellchat_cc) 
cellchat_cc <- identifyOverExpressedGenes(cellchat_cc)
cellchat_cc <- identifyOverExpressedInteractions(cellchat_cc)

cellchat_cc <- computeCommunProb(cellchat_cc, type = "triMean", trim = 0.1, 
                                 distance.use = FALSE, interaction.range = 250, scale.distance = NULL, contact.dependent = TRUE, contact.range = 100)

cellchat_cc <- filterCommunication(cellchat_cc, min.cells = 1)
cellchat_cc.df <- subsetCommunication(cellchat_cc) 
cellchat_cc <- computeCommunProbPathway(cellchat_cc)
cellchat_cc <- aggregateNet(cellchat_cc)
cellchat_cc <- netAnalysis_computeCentrality(cellchat_cc, slot.name = "netP")


# Naive sample
unique(Naive_cc$secondlevel)
unique(Yersinia_pg$secondlevel)
unique(Yersinia_cc$secondlevel)

Naive_cc$secondlevel[is.na(Naive_cc@meta.data$secondlevel)] <- "Unknown"
Naive_cc$sample <- "Naive_cc"
data.input3 = Seurat::GetAssayData(Naive_cc, slot = "data", assay = "SCT" ) # normalized data matrix
meta3 = data.frame(Naive_cc@meta.data) 
spatial.locs3 <-  GetTissueCoordinates(Naive_cc) 
spatial.locs3 <- spatial.locs3[,1:2]
spatial.locs3 <- as.matrix(spatial.locs3)# If spatial.locs is a data.frame, convert it:
rownames(spatial.locs3) <- colnames(data.input3)
cellchat_na <- createCellChat(object = data.input3, meta = meta3, group.by = "secondlevel",
                              datatype = "spatial", coordinates = spatial.locs3, spatial.factors =spatial_factors[2,] )

cellchat_na@DB <- CellChatDB.use
cellchat_na <- subsetData(cellchat_na) 
cellchat_na <- identifyOverExpressedGenes(cellchat_na)
cellchat_na <- identifyOverExpressedInteractions(cellchat_na)

cellchat_na <- computeCommunProb(cellchat_na, type = "triMean", trim = 0.1, 
                                 distance.use = FALSE, interaction.range = 250, scale.distance = NULL, contact.dependent = TRUE, contact.range = 100)

cellchat_na <- filterCommunication(cellchat_na, min.cells = 1)
cellchat_na.df <- subsetCommunication(cellchat_na) 
cellchat_na <- computeCommunProbPathway(cellchat_na)
cellchat_na <- aggregateNet(cellchat_na)
cellchat_na <- netAnalysis_computeCentrality(cellchat_na, slot.name = "netP")




#cellchat object list
object.list <- list(Yersinia_pg = cellchat_pg, Yersinia_cc = cellchat_cc, Naive_cc = cellchat_na)
save(object.list, file = paste0(out,"cellchat_object.list_unfiltered_secondlevel.RData"))
#I will create the object list that is not restrained by overlapping cell types
cellchat3 <- mergeCellChat(object.list, add.names = names(object.list))

#cellchat<- mergeCellChat(object.list, add.names = names(object.list))
object.list2 <- list(Yersinia_pg = cellchat_pg, Yersinia_cc = cellchat_cc, Naive_cc = cellchat_na)
#Identify CellTypes to keep - If you include cell types which are not present in different groups, you cannot perform the comparisons
a<- as.character(sort(unique(cellchat_pg.df$source))) 
b <-as.character(sort(unique(cellchat_cc.df$source)))
c <-as.character(sort(unique(cellchat_na.df$source)))

#Identify the cell types that all  samples share
overlap <- Reduce(intersect, list(a, b, c))
overlap <- setdiff(overlap, "Unknown")

object.list2[[1]] <- subsetCellChat(object.list2[[1]], idents.use = overlap)
object.list2[[2]] <- subsetCellChat(object.list2[[2]], idents.use = overlap)
object.list2[[3]] <- subsetCellChat(object.list2[[3]], idents.use = overlap)


common_levels <- overlap

for (i in seq_along(object.list2)) {
  # Standardize idents
  object.list2[[i]]@idents <- factor(as.character(object.list2[[i]]@idents),
                                     levels = common_levels)
  
  object.list2[[i]]@idents <- droplevels(object.list2[[i]]@idents)
  
  # Tag dataset in metadata (important for merge)
  object.list2[[i]]@meta$dataset <- names(object.list2)[i]
}

cellchat2 <- mergeCellChat(object.list2, add.names = names(object.list2), cell.prefix = TRUE)
cellchat2 <- liftCellChat(cellchat2, group.new = common_levels)

lapply(object.list2, function(x) table(x@idents))



gg1 <- compareInteractions(cellchat2, show.legend = F, group = c(1,2, 3))
gg2 <- compareInteractions(cellchat2, show.legend = F, group = c(1,2, 3), measure = "weight")
gg1 + gg2

#Visualize differential interactions

par(mfrow = c(1,2) , xpd=T)
netVisual_diffInteraction(cellchat2, weight.scale = T,  comparison = c("Yersinia_pg","Yersinia_cc"),  top = 0.05)
netVisual_diffInteraction(cellchat2, weight.scale = T, comparison = c("Yersinia_pg","Yersinia_cc"), measure = "weight",   top = 0.05)


par(mfrow = c(1,2) , xpd=T)
netVisual_diffInteraction(cellchat2, weight.scale = T,targets.use = c("Early Enterocyte", "Middle Enterocyte", "Late Enterocyte"),  comparison = c("Yersinia_pg","Yersinia_cc"),  top = 0.20)
netVisual_diffInteraction(cellchat2, weight.scale = T,targets.use = c("Early Enterocyte", "Middle Enterocyte", "Late Enterocyte"),  comparison = c("Yersinia_pg","Naive_cc"),  top = 0.20)

gg1 <- netVisual_heatmap(cellchat2, comparison = c("Yersinia_pg", "Yersinia_cc"), title.name = "Yersinia vs Naive Ileum")
gg2 <- netVisual_heatmap(cellchat2,comparison = c("Yersinia_pg","Yersinia_cc"), title.name = "Yersinia vs Naive Ileum", measure = "weight")
gg1 + gg2
png(paste0(images,"/Total_interactions_compare_yp_PG_vs_nonPG_heatmap.png"), width = 8, height = 8, units = "in", res = 600)
gg1
dev.off()


gg1 <- netVisual_heatmap(cellchat2, comparison = c("Yersinia_pg", "Naive_cc"), title.name = "Yersinia PG vs Naive Ileum")
gg2 <- netVisual_heatmap(cellchat2,comparison = c("Yersinia_pg","Naive_cc"), title.name = "Yersinia Pg vs Naive Ileum", measure = "weight")
gg1 + gg2
png(paste0(images,"/Total_interactions_compare_yp_PG_vs_naive_heatmap.png"), width = 8, height = 8, units = "in", res = 600)
gg1
dev.off()




weight.max <- getMaxWeight(object.list2, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_circle(object.list2[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list2)[i]))
}

# Aggregate Scatter 
# in the below code you can separately graph for each cell type its ingoing and outgoing relative number of interactions and compare across plots
object.list[[1]]@netP$pathways # access pathways
levels(object.list[[1]]@idents) #Access your cell types
num.link <- lapply(object.list2[1:3], function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

# Combine all values into a single vector
num.link.vec <- unlist(num.link)

# Now compute min and max safely
weight.MinMax <- c(min(num.link.vec), max(num.link.vec))
gg <- list()
for (i in 1:3) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list2[[i]], title = names(object.list2)[i], weight.MinMax = weight.MinMax) + scale_y_continuous(limits = c(0,5)) + scale_x_continuous(limits = c(0,5))
}
patchwork::wrap_plots(plots = gg) 


png(paste0(images,"/Net_Analysis_Scatter_Yp_PG.png"), width = 7, height = 7, units = "in", res = 600)
gg[[1]]
dev.off()
png(paste0(images,"/Net_Analysis_Scatter_yp_nonPG.png"), width = 7, height = 7, units = "in", res = 600)
gg[[2]]
dev.off()
png(paste0(images,"/Net_Analysis_Scatter_Naive.png"), width = 7, height = 7, units = "in", res = 600)
gg[[3]]
dev.off()




unique(cellchat@idents)
anyNA(cellchat@idents)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat2, 
                                            idents.use = "Early Enterocyte", 
                                            comparison = c(1,2),
                                            color.use = c("grey10", "#F8766D", "#00BFC4"))
 gg1
 
png(paste0(images,"/Net_Analysis_Scatter_Enterocytes_Naive_v_Candida.png"), width = 14, height = 12, units = "in", res = 600)
gg1
dev.off()
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "Stem Cells", comparison = c(1,3),
                                            color.use = c( "grey10","#F8766D", "#00BFC4"))
gg1
png(paste0(images,"/Net_Analysis_Scatter_Enterocytes_Naive_v_Crypto.png"), width = 14, height = 12, units = "in", res = 600)
gg1
dev.off()



library(reticulate)
use_python("/home/hartandrew/.conda/envs/sccellfie_env/bin/python", required = TRUE)
ptm = Sys.time()
set.seed(40)
cellchat2 <- computeNetSimilarityPairwise(cellchat2, type = "functional")
cellchat2 <- netEmbedding(cellchat2, type = "functional")
cellchat2 <- netClustering(cellchat2, type = "functional")








gg1 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  
               cutoff.pvalue = 0.05, slot.name = "net",
               stacked = T, do.stat = TRUE, comparison = c(1,3))
gg1
gg1$data <- subset(gg1$data, pvalues < 0.05)
gg1
svg(paste0(images, "/Ranknet_allsignaling_YersiniaPG_vs_Naive_net.svg"), 
    width = 4, height = 9)
print(gg1 )
dev.off()

gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  
               cutoff.pvalue = 0.10,slot.name = "net",
               stacked = T, do.stat = TRUE, comparison = c(1,2))
gg2
gg2$data <- subset(gg2$data, pvalues <= 0.1)
gg2
svg(paste0(images, "/Ranknet_allsignaling_YersiniaPG_vs_YersiniaCC_net.svg"), 
    width = 4, height = 9)
print(gg2 )
dev.off()

# Investigate 

unique(object.list[[3]]@netP$pathways)
unique(object.list[[1]]@netP$pathways)
unique(object.list[[2]]@netP$pathways)

# IL1 signaling is found within the PG but not other groups - What are the cell types providing IL1 signals
pathways.show <- c("IL1") 
netVisual_aggregate(cellchat2, signaling = pathways.show, sample.use = "Yersinia_pg",
                    layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 0.5, alpha.image = 0.01,
                    vertex.label.cex = 0)

pathway <- "IL1"   # or "IL6", "TGFB", etc.
pairLR.use <- subsetDB(object.list[[1]]@DB, search = "IL18_IL18R1_IL18RAP")
svg(paste0(images, "/Circle_Yersinia_PG_IL1_all_aggregate.svg"), 
    width = 10, height = 10)

netVisual_chord_gene(
  object.list[[1]],
  signaling = c("IL1"),
  slot.name = "netP",
  layout = "circle",
  edge.width.max = 2,
  vertex.size.max = 0.5,
  alpha.image = 0.01,
  vertex.label.cex = 0
)
dev.off()



svg(paste0(images, "/Circle_gene_Yersinia_pg_IL1_all.svg"), 
    width = 10, height = 10)
netVisual_chord_cell(object.list[[1]],
                     signaling = c("IL1"),
                     big.gap = 0, 
                     title.name = paste0("IL1 Pathways - ", names(object.list)[1]))

dev.off()


# Extract all significant communications for IL1
il1_comm <- subsetCommunication(cellchat, signaling = "IL1")
unique(il1_comm$target)
active.groups <- c("Early Enterocyte",  "Late Enterocyte",   "Middle Enterocyte",
                   "Cycling T Cells",   "Effector CD4 ab T Cells", "gd T Cells ",  "ILC1s",         
                   "ILC2s" ,   "ILC3s",      "NK Cells",    "NKT Cells") 



netVisual_aggregate(cellchat, signaling = "IL1",
                    sample.use = "Yersinia_pg",
                    layout = "circle",
                    edge.width.max = 7,
                    vertex.size.max = 0.2,
                    point.size = 0.6,
                    alpha.image = 0.1,
                    vertex.label.cex = 0.5,
                    #top = 0.25,
                    show.legend = FALSE)

#Figure 7F
svg(paste0(images, "/Circle_gene_Yersinia_pg_IL1_all.svg"), 
    width = 10, height = 10)
netVisual_aggregate(cellchat, signaling = "IL1",
                    sample.use = "Yersinia_pg",
                    layout = "circle",
                    edge.width.max = 10,
                    vertex.size.max = 0.2,
                    point.size = 0.6,
                    alpha.image = 0.1,
                    vertex.label.cex = 0.5,
                    #top = 0.25,
                    show.legend = FALSE)

dev.off()

# Evaluate the Ranknet plots more thoroughly
gg1 <- rankNet(cellchat2, mode = "comparison", measure = "weight",  slot.name = "net",
               cutoff.pvalue = 0.05,sources.use = c("Early Enterocyte", "Middle Enterocyte", "Late Enterocyte"),
               stacked = T, do.stat = TRUE, comparison = c(1,3))
gg1

svg(paste0(images, "/Ranknet_Enterocyte_Ligands_YersiniaPG_vs_Naive_net.svg"), 
    width = 4, height = 9)
print(gg1 )
dev.off()

gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight",   slot.name = "net",
               cutoff.pvalue = 0.05,sources.use = c("Early Enterocyte", "Middle Enterocyte", "Late Enterocyte"),
               stacked = T, do.stat = TRUE, comparison = c(1,2))
gg2
# Figure 7D
svg(paste0(images, "/Ranknet_Enterocyte_Ligands_YersiniaPG_vs_YersiniaCC_net.svg"), 
    width = 4, height = 9)
print(gg2 )
dev.off()

unique(Yersinia_pg$secondlevel)

gg1 <- rankNet(cellchat2, mode = "comparison", measure = "weight", 
               #slot.name = "net",
               cutoff.pvalue = 0.05,targets.use = c("Early Enterocyte", "Middle Enterocyte", "Late Enterocyte"),
               stacked = T, do.stat = TRUE, comparison = c(1,3))
gg1

svg(paste0(images, "/Ranknet_Enterocyte_Receptors_YersiniaPG_vs_Naive_net.svg"), 
    width = 4, height = 9)
print(gg1 )
dev.off()

gg2 <- rankNet(cellchat2, mode = "comparison", measure = "weight", 
               #slot.name = "net", 
               cutoff.pvalue = 0.05,targets.use  = c( "Early Enterocyte", "Middle Enterocyte", "Late Enterocyte"),
               stacked = T, do.stat = TRUE, comparison = c(1,2))
gg2

svg(paste0(images, "/Ranknet_Enterocyte_receptors_YersiniaPG_vs_YersiniaCC_net.svg"), 
    width = 4, height = 9)
print(gg2 )
dev.off()




il1_comm <- subsetCommunication(cellchat2, signaling = "IL-10")
unique(il1_comm$target)
active.groups <- c("Early Enterocyte",  "Late Enterocyte",   "Middle Enterocyte",
                   "Cycling T Cells",   "Effector CD4 ab T Cells", "gd T Cells ",  "ILC1s",         
                   "ILC2s" ,   "ILC3s",      "NK Cells",    "NKT Cells") 

unique(cellchat@netP$pathways)

netVisual_aggregate(cellchat, signaling = "IL1",
                    sample.use = "Yersinia_pg",
                    layout = "circle",
                    edge.width.max = 7,
                    vertex.size.max = 0.2,
                    point.size = 0.6,
                    alpha.image = 0.1,
                    vertex.label.cex = 0.5,
                    #top = 0.25,
                    show.legend = FALSE)


svg(paste0(images, "/Circle_gene_Yersinia_pg_IL1_all.svg"), 
    width = 10, height = 10)
netVisual_aggregate(cellchat, signaling = "IL1",
                    sample.use = "Yersinia_pg",
                    layout = "circle",
                    edge.width.max = 10,
                    vertex.size.max = 0.2,
                    point.size = 0.6,
                    alpha.image = 0.1,
                    vertex.label.cex = 0.5,
                    #top = 0.25,
                    show.legend = FALSE)

dev.off()


