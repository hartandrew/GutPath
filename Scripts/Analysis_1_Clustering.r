# This code takes the Seurat object that was made and normalizes and scales data, performs PCA, clustering, integreation, cell cycle scoring, doublet removal, and WNN dual modality analyses

---
title: "SingleCell_Analysis"
author: "AndrewHart"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Load the Libraries----
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
library(RColorBrewer)
library(monocle3)
library(gt)
options(future.globals.maxSize = 200000 * 1024^2)
```

# Directory set up and data loading
```{r}
getwd() #/home/hartandrew
setwd("/path/to/analysis/directory")
out <- "/path/to/analysis/directory/Output"
version <- "/path/to/analysis/directory/Version"
seurat <- "/path/to/analysis/directory/Seurat_Files"
images <- "/path/to/analysis/directory/Images"
CSV <- "/path/to/analysis/directory/CSV"

#Load Data from Processing01.rmd script - The R environment for the merged objects prior to downstream analysis was loaded - Ileum_MLN_combined_prior_to_anlysis.RData in Seurat_Files for this project
```

# Color Palettes
```{r}
colors6 <- c("#941339", "#379FFB", "#D6A206", "#029102", "#6313B7", "#AFAFAF")
colors7 <- c("#941339", "#379FFB", "#D6A206", "#029102", "#6313B7", "#AFAFAF","#542600" )
colors8 <- c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
colors9 <- c("#332288", "#88CCEE", "#117733","#44AA99", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
colors5 <- c("#941339", "#379FFB", "#D6A206", "#029102", "#6313B7")
```

# QC Averages across samples
```{r}

#Theme for Plots 
my_theme <-   theme(panel.background = element_rect(fill = "white", 
            colour = NA), panel.border = element_rect(fill = NA, 
            colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
            panel.grid.minor = element_line(linewidth = rel(0.5)), 
          axis.title.x = element_blank(),
            strip.background = element_rect(fill = "grey85", 
                colour = "grey20"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

```




```{r}
#Add Variable for 10x kit chemistry type
#v3 <- c("CHMI2482_MLN", "CHMI2482_LP", "CHMI2482_IEL", "CHMI2489_MLN", "CHMI2489_IEL", "CHMI2489_LP")
MLN_data@meta.data$Chemistry  <- ifelse(MLN_data$orig.ident  %in% c("CHMI2482_MLN", "CHMI2482_LP", "CHMI2482_IEL", "CHMI2489_MLN", "CHMI2489_IEL", "CHMI2489_LP"), "10X_v3", "10X_v4")
unique(MLN_data@meta.data$Chemistry )
Ileum_data@meta.data$Chemistry  <- ifelse(Ileum_data$orig.ident  %in% c("CHMI2482_MLN", "CHMI2482_LP", "CHMI2482_IEL", "CHMI2489_MLN", "CHMI2489_IEL", "CHMI2489_LP"), "10X_v3", "10X_v4")

```

# Remove Samples not used in publication
```{r}
Idents(MLN_data) <- "InfectionStatus"
MLN_data <- subset(MLN_data, idents = c("H.poly", "SFB_CU"), invert = T)
Idents(Ileum_data) <- "InfectionStatus"
Ileum_data <- subset(Ileum_data, idents = c( "SFB_CU"), invert = T)

#Save Cell Numbers
df <- as.data.frame(table(Ileum_data$InfectionStatus, Ileum_data$orig.ident, Ileum_data$Tissue ))
df <- df[df$Freq != 0 ,]
colnames(df) <- c("Infection", "Sample", "Tissue", "Freq")
write.csv(df, paste0(CSV, "/Ileum_Table_cells_per_sample_filtered_Analysis1.csv"))
df <- as.data.frame(table(MLN_data$InfectionStatus, MLN_data$orig.ident, MLN_data$Tissue ))
df <- df[df$Freq != 0 ,]
colnames(df) <- c("Infection", "Sample", "Tissue", "Freq")
write.csv(df, paste0(CSV, "/MLN_Table_cells_per_sample_filtered_Analysis1.csv"))
remove(df)

```


```{r}
#At this point, the RNA assay has only count (raw) layers and the ADT assay has both count (raw) and data (normalized) layers 
head(colnames(MLN_data[["ADT"]]))
head(colnames(MLN_data[["RNA"]]))

#Split the combined and filtered data into the two assay for individual Integration
MLN_RNA <- CreateSeuratObject(
  MLN_data[["RNA"]],
  project = "SeuratProject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = MLN_data@meta.data)
MLN_RNA
#Check original identity
table(MLN_RNA$orig.ident)

#Split the combined and filtered data into the two assay for individual Integration
Ileum_RNA <- CreateSeuratObject(
  Ileum_data[["RNA"]],
  project = "SeuratProject",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = Ileum_data@meta.data)
Ileum_RNA
#Check original identity
table(Ileum_RNA$orig.ident)
```

```{r}
MLN_data[['ADT']] <-  split(x = MLN_data[['ADT']], f = MLN_data$orig.ident)
Layers(MLN_data[['ADT']])

MLN_ADT <- CreateSeuratObject(
  MLN_data[["ADT"]],
  project = "SeuratProject",
  assay = "ADT",
  names.delim = "_",
  meta.data = MLN_data@meta.data)
MLN_ADT

#Check original identity
table(MLN_ADT$orig.ident)
table(MLN_ADT$Tissue)
table(MLN_ADT$InfectionStatus)
#Save number of cells per sample
 write.csv(table(MLN_ADT$orig.ident), paste0(CSV, "/MLN_Table_cells_per_sample_filtered.csv"))
 write.csv(table(MLN_ADT$Tissue), paste0(CSV, "/MLN_Table_cells_per_Tissue_filtered.csv"))
 write.csv(table(MLN_ADT$InfectionStatus), paste0(CSV,  "/MLN_Table_cells_per_InfectionStatus_filtered.csv"))
 
 Ileum_data[['ADT']] <-  split(x = Ileum_data[['ADT']], f = Ileum_data$orig.ident)
Layers(Ileum_data[['ADT']])

Ileum_ADT <- CreateSeuratObject(
  Ileum_data[["ADT"]],
  project = "SeuratProject",
  assay = "ADT",
  names.delim = "_",
  meta.data = Ileum_data@meta.data)
Ileum_ADT

#Check original identity
table(Ileum_ADT$orig.ident)
table(Ileum_ADT$Tissue)
table(Ileum_ADT$InfectionStatus)
#Save number of cells per sample
 write.csv(table(Ileum_ADT$orig.ident), paste0(CSV, "/Ileum_Table_cells_per_sample_filtered.csv"))
 write.csv(table(Ileum_ADT$Tissue), paste0(CSV, "/Ileum_Table_cells_per_Tissue_filtered.csv"))
 write.csv(table(Ileum_ADT$InfectionStatus), paste0(CSV,  "/Ileum_Table_cells_per_InfectionStatus_filtered.csv"))
```



#RNA-BASED CLUSTERING before integration
```{r}
# Run the analysis on the MLN to visualize
DefaultAssay(MLN_data) <- "RNA" 
MLN_data[['RNA']] <-  JoinLayers(MLN_data[['RNA']])
MLN_data[['RNA']] <-  split(x = MLN_data[['RNA']], f = MLN_data$InfectionStatus)
MLN_data <- NormalizeData(MLN_data)
MLN_data <- FindVariableFeatures(MLN_data, selection.method = "vst", nfeatures = 3000) 

VariableFeaturePlot(MLN_data, cols = c("black", "red"),
  pt.size = 0.1, log = NULL, selection.method = NULL, raster = FALSE) + scale_y_log10()

all.genes <- rownames(MLN_data[["RNA"]])
MLN_data <- ScaleData(MLN_data, features = all.genes)

MLN_data <- RunPCA(MLN_data, verbose = FALSE, reduction.name = "unintegrated_RNA_pca")
ElbowPlot(MLN_data, reduction = "unintegrated_RNA_pca", ndims = 50)
# Determine percent of variation associated with each PC
pct <- MLN_data[["unintegrated_RNA_pca"]]@stdev / sum(MLN_data[["unintegrated_RNA_pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1 

#Method 2 to visualize PCs and variance changes
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2 
#Visualize the total percent variation captured by the cumultaive PC chosen
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
  ggplot(plot_df, aes(cumu, pct, label = rank)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

MLN_data <- FindNeighbors(MLN_data, dims = 1:30, reduction = "unintegrated_RNA_pca", graph.name = "RNA_unint")
MLN_data <- FindClusters(MLN_data, resolution = 2, verbose = FALSE, cluster.name = "unintegrated_RNA_clusters", group.singletons = F, graph.name = "RNA_unint")
MLN_data <- RunUMAP(MLN_data, dims = 1:30, reduction = "unintegrated_RNA_pca", reduction.name = "unintegrated_RNA_umap")
DimPlot(MLN_data, reduction = "unintegrated_RNA_umap", group.by = c("unintegrated_RNA_clusters", "Tissue", "InfectionStatus"), combine = FALSE)

remove(co2, pct, cumu, co1, plot_df)
```



```{r}
# Run the analysis on the Ileum before integration to visualize
Ileum_data$SampleType <- Ileum_data$Tissue
Ileum_data$SampleType<- paste0("Ileum ", Ileum_data$Tissue)
Ileum_data$SampleType[Ileum_data$orig.ident %in% c("CHMI25134_LP_D","CHMI25132_LP_D" )] <- "Duodenum Lamina Propria"
Ileum_data$SampleType[Ileum_data$orig.ident %in% c("CHMI25134_IEL_D","CHMI25132_IEL_D" )] <- "Duodenum Intestinal Epithelial Cells"
Ileum_data$Infection_Organ <- paste0(Ileum_data$InfectionStatus, " ", Ileum_data$SampleType)

DefaultAssay(Ileum_data) <- "RNA"
Ileum_data[['RNA']] <-  JoinLayers(Ileum_data[['RNA']])
Ileum_data[['RNA']] <-  split(x = Ileum_data[['RNA']], f = Ileum_data$Infection_Organ)
Ileum_data <- NormalizeData(Ileum_data)
Ileum_data <- FindVariableFeatures(Ileum_data, selection.method = "vst", nfeatures = 3000) 

VariableFeaturePlot(Ileum_data, cols = c("black", "red"),
  pt.size = 0.1, log = NULL, selection.method = NULL, raster = FALSE) + scale_y_log10()

all.genes <- rownames(Ileum_data[["RNA"]])
Ileum_data <- ScaleData(Ileum_data, features = all.genes)
Ileum_data <- RunPCA(Ileum_data, verbose = FALSE, reduction.name = "unintegrated_RNA_pca")
ElbowPlot(Ileum_data, reduction = "unintegrated_RNA_pca", ndims = 50)

pct <- Ileum_data[["unintegrated_RNA_pca"]]@stdev / sum(Ileum_data[["unintegrated_RNA_pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1

#Method 2 to determine the PC at which point the percentage of change in PC drops to <0.1 
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
  ggplot(plot_df, aes(cumu, pct, label = rank)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

  
Ileum_data <- FindNeighbors(Ileum_data, dims = 1:35, reduction = "unintegrated_RNA_pca", graph.name = "RNA_unint" )
Ileum_data <- FindClusters(Ileum_data, resolution = 2, verbose = FALSE, cluster.name = "unintegrated_RNA_clusters", group.singletons = F, graph.name = "RNA_unint")
Ileum_data <- RunUMAP(Ileum_data, dims = 1:35, reduction = "unintegrated_RNA_pca", reduction.name = "unintegrated_RNA_umap")
DimPlot(Ileum_data, reduction = "unintegrated_RNA_umap", group.by = c("unintegrated_RNA_clusters", "Tissue", "InfectionStatus"), combine = FALSE)
DimPlot(Ileum_data, reduction = "unintegrated_RNA_umap", group.by = c("SampleType"), combine = FALSE)

remove(co2, pct, cumu, co1, plot_df)

```




#INTEGRATE RNA 
```{r}
#chunk can take a long time to run depending on operating system - Memory requirements are high
# I examined the results of multiple integration methods to ascertain what might be the best fit without losing biological signal
options(future.globals.maxSize = 200000 * 1024^2)

MLN_data <- IntegrateLayers(
  object = MLN_data, method = HarmonyIntegration, 
  orig.reduction = "unintegrated_RNA_pca", new.reduction = "integrated.harmony",
  verbose = FALSE)

DefaultAssay(MLN_data) <- "RNA"
MLN_data <- FindNeighbors(MLN_data, reduction = "integrated.harmony", dims = 1:30, graph.name = "harmony_rna")
MLN_data <- FindClusters(MLN_data, resolution = 2, cluster.name = "harmony_integration_clusters", group.singletons = F, graph.name = "harmony_rna")
MLN_data <- RunUMAP(MLN_data, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(MLN_data, reduction = "umap.harmony", group.by = c("harmony_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)


Ileum_data <- IntegrateLayers(
  object = Ileum_data, method = HarmonyIntegration, 
  orig.reduction = "unintegrated_RNA_pca", new.reduction = "integrated.harmony",
  verbose = FALSE)

DefaultAssay(Ileum_data) <- "RNA"
Ileum_data <- FindNeighbors(Ileum_data, reduction = "integrated.harmony", dims = 1:35, graph.name = "harmony_rna")
Ileum_data <- FindClusters(Ileum_data, resolution = 2, cluster.name = "harmony_integration_clusters", group.singletons = F, graph.name = "harmony_rna")
Ileum_data <- RunUMAP(Ileum_data, reduction = "integrated.harmony", dims = 1:35, reduction.name = "umap.harmony")
DimPlot(Ileum_data, reduction = "umap.harmony", group.by = c("harmony_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
FeaturePlot(Ileum_data, reduction = "umap.harmony", "Epcam")
FeaturePlot(Ileum_data, reduction = "umap.harmony", "Gata3")


options(future.globals.maxSize = 200000 * 1024^2)
MLN_data <- IntegrateLayers(
  object = MLN_data, method = RPCAIntegration, 
  orig.reduction = "unintegrated_RNA_pca", new.reduction = "integrated.rpca",
  verbose = FALSE)

DefaultAssay(MLN_data) <- "RNA"
MLN_data <- FindNeighbors(MLN_data, reduction = "integrated.rpca", dims = 1:30, graph.name = "rpca_rna")
MLN_data <- FindClusters(MLN_data, resolution = 2, cluster.name = "rpca_integration_clusters", group.singletons = F, graph.name = "rpca_rna")
MLN_data <- RunUMAP(MLN_data, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(MLN_data, reduction = "umap.rpca", group.by = c("rpca_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)


Ileum_data <- IntegrateLayers(
  object = Ileum_data, method = RPCAIntegration, 
  orig.reduction = "unintegrated_RNA_pca", new.reduction = "integrated.rpca",
  verbose = FALSE)

DefaultAssay(Ileum_data) <- "RNA"
Ileum_data <- FindNeighbors(Ileum_data, reduction = "integrated.rpca", dims = 1:35, graph.name = "rpca_rna")
Ileum_data <- FindClusters(Ileum_data, resolution = 2, cluster.name = "rpca_integration_clusters", group.singletons = F, graph.name = "rpca_rna")
Ileum_data <- RunUMAP(Ileum_data, reduction = "integrated.rpca", dims = 1:35, reduction.name = "umap.rpca")
DimPlot(Ileum_data, reduction = "umap.rpca", group.by = c("rpca_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
FeaturePlot(Ileum_data, reduction = "umap.rpca", "Epcam")
FeaturePlot(Ileum_data, reduction = "umap.rpca", "Gata3")

```

```{r}
#SCVI Integration method - PMID: 34949812 - one of the better and more scaleable approaches
library(reticulate)
use_condaenv("/home/hartandrew/.conda/envs/scvi-env", required = TRUE)
Sys.setenv(MPLBACKEND = "Agg")
library(SeuratWrappers)
options(future.globals.maxSize = 100000 * 1024^2)
MLN_data <- IntegrateLayers(
  object = MLN_data, method = scVIIntegration, 
  conda_env = "/home/hartandrew/.conda/envs/scvi-env", 
  orig.reduction = "unintegrated_RNA_pca", new.reduction = "integrated.scvi",
  verbose = FALSE)

DefaultAssay(MLN_data) <- "RNA"
MLN_data <- FindNeighbors(MLN_data, reduction = "integrated.scvi", dims = 1:30, graph.name = "scvi_rna")
MLN_data <- FindClusters(MLN_data, resolution = 2, cluster.name = "scvi_integration_clusters", group.singletons = F, graph.name = "scvi_rna")
MLN_data <- RunUMAP(MLN_data, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
DimPlot(MLN_data, reduction = "umap.scvi", group.by = c("scvi_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)


Ileum_data <- IntegrateLayers(
  object = Ileum_data, method = scVIIntegration, 
  conda_env = "/home/hartandrew/.conda/envs/scvi-env", 
  orig.reduction = "unintegrated_RNA_pca", new.reduction = "integrated.scvi",
  features = VariableFeatures(Ileum_data),
  verbose = FALSE)

DefaultAssay(Ileum_data) <- "RNA"
Ileum_data <- FindNeighbors(Ileum_data, reduction = "integrated.scvi", dims = 1:30, graph.name = "scvi_rna")
Ileum_data <- FindClusters(Ileum_data, resolution = 2, cluster.name = "scvi_integration_clusters", group.singletons = F, graph.name = "scvi_rna")
Ileum_data <- RunUMAP(Ileum_data, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", graph = "scvi_rna" )
DimPlot(Ileum_data, reduction = "umap.scvi", group.by = c("scvi_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(Ileum_data, reduction = "umap.scvi", group.by = c("Mouse"), combine = FALSE, label.size = 4, label = TRUE)
FeaturePlot(Ileum_data, reduction = "umap.scvi", "Epcam")
FeaturePlot(Ileum_data, reduction = "umap.scvi", "Ly6m")
```

```{r}

#Visualize the Integration methods and clustering
DimPlot(MLN_data, reduction = "umap.harmony", group.by = c("harmony_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(MLN_data, reduction = "umap.rpca", group.by = c("rpca_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(MLN_data, reduction = "umap.scvi", group.by = c("scvi_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(Ileum_data, reduction = "umap.harmony", group.by = c("harmony_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(Ileum_data, reduction = "umap.rpca", group.by = c("rpca_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(Ileum_data, reduction = "umap.scvi", group.by = c("scvi_integration_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
```

# Cell Cycle Scoring 
```{r}
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  They were updated in 2019
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Convert human to mouse gene names
library(biomaRt)
convertHumanGeneList <- function(x){
require("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}
m.s.genes <- convertHumanGeneList(s.genes)
m.g2m.genes <- convertHumanGeneList(g2m.genes)

DefaultAssay(MLN_data) <- "RNA"
MLN_data <- JoinLayers(MLN_data)
MLN_data <- CellCycleScoring(MLN_data, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
DefaultAssay(Ileum_data) <- "RNA"
Ileum_data <- JoinLayers(Ileum_data)
Ileum_data <- CellCycleScoring(Ileum_data, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
```

```{r}
DimPlot(MLN_data, reduction = "umap.harmony",   group.by = c( "Phase"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
DimPlot(MLN_data, reduction = "umap.rpca",   group.by = c( "Phase"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
DimPlot(MLN_data, reduction = "umap.scvi",   group.by = c( "Phase"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
DimPlot(Ileum_data, reduction = "umap.harmony",   group.by = c( "Phase"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
DimPlot(Ileum_data, reduction = "umap.rpca",   group.by = c( "Phase"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
DimPlot(Ileum_data, reduction = "umap.scvi",   group.by = c( "Phase"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)

```

```{r}
#Evaluate RNA-derived Clusters Manually
# below I sought look for epithelial populations and large groups of lineages such as Ptprc+ hematopoietic cells
FeaturePlot(MLN_data, reduction = "umap.rpca",features = c( "Epcam", "Olfm4", "Lgr5", "Dclk1", "Fcgr3", "Ncam1"))
FeaturePlot(MLN_data, reduction = "umap.rpca",features = c("Ighd",  "Cd3e", "Cd14",  "Cd4", "Ptprc", "Cd8a", "Ly6c1", "Klrg1"))
FeaturePlot(Ileum_data, reduction = "umap.scvi",features = c( "Epcam", "Olfm4", "Lgr5", "Dclk1", "Fcgr3", "Ncam1"))
FeaturePlot(Ileum_data, reduction = "umap.scvi",features = c("Ighd",  "Cd3e", "Cd14",  "Cd4", "Ptprc", "Cd8a", "Ly6c1", "Klrg1"))
```

#ADT Clustering before integration
```{r}
DefaultAssay(MLN_data) <- "ADT"
rownames(MLN_data)
MLN_data
Isotype_labels <- grep("isotype", rownames(MLN_data), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_data)[!rownames(MLN_data) %in% Isotype_labels] #remove isotypes

#MLN_data <- FindVariableFeatures(MLN_data) # Use all features and the dsb normalized instead
MLN_data[['ADT']] <-  JoinLayers(MLN_data[['ADT']])
MLN_data[['ADT']] <-  split(x = MLN_data[['ADT']], f = MLN_data$InfectionStatus)
MLN_data <- ScaleData(MLN_data, assay = "ADT")
MLN_data <- RunPCA(MLN_data, features =prots, reduction.name = "unintegrated_pca_adt", reduction.key = "pcaADT_", verbose = FALSE)
DimPlot(MLN_data, reduction = "unintegrated_pca_adt", group.by = "orig.ident")
ElbowPlot(MLN_data, reduction = "unintegrated_pca_adt", ndims = 40)

pct <- MLN_data[["unintegrated_pca_adt"]]@stdev / sum(MLN_data[["unintegrated_pca_adt"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1 

#Method 2 to determine the PC at which point the percentage of change in PC drops to <0.1 
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
  ggplot(plot_df, aes(cumu, pct, label = rank)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


MLN_data <- FindNeighbors(object = MLN_data, dims = 1:20, features = VariableFeatures(object = MLN_data),
                           k.param = 30, verbose = FALSE, reduction = "unintegrated_pca_adt" )

#direct graph clustering 
MLN_data <- FindClusters(MLN_data, resolution = 2,  graph.name = 'ADT_snn', verbose = FALSE, cluster.name = "unintegrated_ADT_clusters", group.singletons = F)

MLN_data <- RunUMAP(object = MLN_data, assay = "ADT", features = prots, 
                      n.neighbors = 30,  verbose = FALSE, reduction.name = "umap.ADT.unintegrated")

DimPlot(MLN_data, reduction = "umap.ADT.unintegrated", group.by = c("unintegrated_ADT_clusters", "Tissue", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE)
```

```{r}
DefaultAssay(Ileum_data) <- "ADT"
rownames(Ileum_data)
Ileum_data
Isotype_labels <- grep("isotype", rownames(Ileum_data), ignore.case = TRUE, value = TRUE)
prots = rownames(Ileum_data)[!rownames(Ileum_data) %in% Isotype_labels] #remove isotypes

#Ileum_data <- FindVariableFeatures(Ileum_data) # Use all features and the dsb normalized instead
Ileum_data[['ADT']] <-  JoinLayers(Ileum_data[['ADT']])
Ileum_data[['ADT']] <-  split(x = Ileum_data[['ADT']], f = Ileum_data$Infection_Organ)
Ileum_data <- ScaleData(Ileum_data, assay = "ADT")
Ileum_data <- RunPCA(Ileum_data, features =prots, reduction.name = "unintegrated_pca_adt", reduction.key = "pcaADT_", verbose = FALSE)
DimPlot(Ileum_data, reduction = "unintegrated_pca_adt", group.by = "orig.ident")
ElbowPlot(Ileum_data, reduction = "unintegrated_pca_adt", ndims = 40)
# Determine percent of variation associated with each PC
pct <- Ileum_data[["unintegrated_pca_adt"]]@stdev / sum(Ileum_data[["unintegrated_pca_adt"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1 

#Method 2 
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank)) +  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + theme_bw()


Ileum_data <- FindNeighbors(object = Ileum_data, dims = 1:20, features = VariableFeatures(object = Ileum_data),
                           k.param = 30, verbose = FALSE, reduction = "unintegrated_pca_adt" )

#direct graph clustering 
Ileum_data <- FindClusters(Ileum_data, resolution = 2,  graph.name = 'ADT_snn', verbose = FALSE, cluster.name = "unintegrated_ADT_clusters", group.singletons = F)
Ileum_data <- RunUMAP(object = Ileum_data, assay = "ADT", features = prots, 
                      n.neighbors = 30,  verbose = FALSE, reduction.name = "umap.ADT.unintegrated")

DimPlot(Ileum_data, reduction = "umap.ADT.unintegrated", group.by = c("unintegrated_ADT_clusters", "Tissue", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE)
```
#Integrate the ADT data -
```{r}
DefaultAssay(MLN_data) <- "ADT"
MLN_data[['ADT']] <-  JoinLayers(MLN_data[['ADT']])
MLN_data[["ADT"]] <- split(MLN_data[["ADT"]], f = MLN_data$InfectionStatus)
MLN_data <- IntegrateLayers(
  object = MLN_data, method = HarmonyIntegration, features = prots, 
  orig.reduction = "unintegrated_pca_adt", new.reduction = "integrated.ADT.harmony", verbose = TRUE)
MLN_data <- FindNeighbors(MLN_data, reduction = "integrated.ADT.harmony", dims = 1:20, graph.name = "harmony_adt")
MLN_data <- FindClusters(MLN_data, resolution = 2, cluster.name = "adt_harmony_clusters", group.singletons = F, graph.name = "harmony_adt")
MLN_data <- RunUMAP(MLN_data, reduction = "integrated.ADT.harmony", dims = 1:20, reduction.name = "umap.harmony.ADT")
DimPlot(MLN_data, reduction = "umap.harmony.ADT", group.by = c("adt_harmony_clusters", "Tissue", "InfectionStatus"),combine = FALSE, label.size = 2)

MLN_data <- IntegrateLayers(
  object = MLN_data, method = RPCAIntegration, features = prots, 
  orig.reduction = "unintegrated_pca_adt", new.reduction = "integrated.ADT.rpca", verbose = TRUE)
MLN_data <- FindNeighbors(MLN_data, reduction = "integrated.ADT.rpca", dims = 1:20, graph.name = "rpca_adt")
MLN_data <- FindClusters(MLN_data, resolution = 2, cluster.name = "adt_rpca_clusters", group.singletons = F, graph.name = "rpca_adt")
MLN_data <- RunUMAP(MLN_data, reduction = "integrated.ADT.rpca", dims = 1:20, reduction.name = "umap.rpca.ADT")
DimPlot(MLN_data, reduction = "umap.rpca.ADT", group.by = c("adt_rpca_clusters", "Tissue", "InfectionStatus"),combine = FALSE, label.size = 2)
DimPlot(MLN_data, reduction = "umap.rpca.ADT", group.by = c("adt_rpca_clusters"),combine = FALSE, label.size = 2)
ggsave("MLN_data_based_integrated_rpca_UMAP_ADT_clusters.png",plot = last_plot(), device = NULL,
  path = images,scale = 1, width = 5, height = 5,  units = c("in"),dpi = 600, limitsize = TRUE, bg = NULL)
```

```{r}
DefaultAssay(Ileum_data) <- "ADT"
Ileum_data[['ADT']] <-  JoinLayers(Ileum_data[['ADT']])
Ileum_data[["ADT"]] <- split(Ileum_data[["ADT"]], f = Ileum_data$Infection_Organ)
Ileum_data <- IntegrateLayers(
  object = Ileum_data, method = RPCAIntegration, features = prots, 
  orig.reduction = "unintegrated_pca_adt", new.reduction = "integrated.ADT.rpca", verbose = TRUE)
Ileum_data <- FindNeighbors(Ileum_data, reduction = "integrated.ADT.rpca", dims = 1:20, graph.name = "rpca_adt" )
Ileum_data <- FindClusters(Ileum_data, resolution = 2, cluster.name = "adt_rpca_clusters", group.singletons = F, graph.name = "rpca_adt" )
Ileum_data <- RunUMAP(Ileum_data, reduction = "integrated.ADT.rpca", dims = 1:20, reduction.name = "umap.rpca.ADT")

Ileum_data <- IntegrateLayers(
  object = Ileum_data, method = HarmonyIntegration, features = prots, 
  orig.reduction = "unintegrated_pca_adt", new.reduction = "integrated.ADT.harmony", verbose = TRUE)
Ileum_data <- FindNeighbors(Ileum_data, reduction = "integrated.ADT.harmony", dims = 1:20, graph.name = "harmony_adt" )
Ileum_data <- FindClusters(Ileum_data, resolution = 2, cluster.name = "adt_harmony_clusters", group.singletons = F, graph.name = "harmony_adt" )
Ileum_data <- RunUMAP(Ileum_data, reduction = "integrated.ADT.harmony", dims = 1:20, reduction.name = "umap.harmony.ADT")

DimPlot(Ileum_data, reduction = "umap.rpca.ADT", group.by = c("adt_rpca_clusters", "Tissue", "InfectionStatus"),combine = FALSE, label.size = 2)
DimPlot(Ileum_data, reduction = "umap.harmony.ADT", group.by = c("adt_harmony_clusters", "Tissue", "InfectionStatus"),combine = FALSE, label.size = 2)

DimPlot(Ileum_ADT, reduction = "umap.harmony.ADT", group.by = c("adt_harmony_clusters"),combine = TRUE, label = T, label.size = 3) + NoLegend()
ggsave("Ileum_UMAP_ADT_clusters.png",plot = last_plot(), device = NULL,
  path = images,scale = 1, width = 5, height = 5,  units = c("in"),dpi = 600, limitsize = TRUE, bg = NULL)
```



```{r}
DimPlot(MLN_data, reduction = "umap.harmony.ADT", group.by = c("adt_harmony_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(MLN_data, reduction = "umap.rpca.ADT", group.by = c("adt_rpca_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)

DimPlot(Ileum_data, reduction = "umap.harmony.ADT", group.by = c("adt_harmony_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)
DimPlot(Ileum_data, reduction = "umap.rpca.ADT", group.by = c("adt_rpca_clusters", "Tissue", "InfectionStatus"), combine = FALSE, label.size = 4, label = TRUE)

```


#Combine Data and perform WNN multimodal neighbor analysis and Umap Visualization
```{r}
DefaultAssay(MLN_data) <- "RNA"

MLN_data = FindMultiModalNeighbors(
  MLN_data, reduction.list = list("integrated.harmony", "integrated.ADT.harmony"), 
  dims.list = list(1:30, 1:20), 
  modality.weight.name = c("RNA.weight", "ADT.weight"), 
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn_harmony",
  verbose = T)

MLN_data <- RunUMAP(MLN_data, nn.name = "weighted.nn_harmony", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
MLN_data <- FindClusters(MLN_data, graph.name = "wsnn",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters", n.start = 20, group.singletons = F)

DimPlot(MLN_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "Tissue", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE)

```

```{r}
#Evaluate Cluster Genes
Idents(MLN_data) <- "WNN_clusters"
FeaturePlot_scCustom(MLN_data, reduction = "wnn.umap",features = c("Ighd",  "Cd3e", "Cd14", "Cd4", "Ptprc", "Cd8a", "Ly6c1", "Epcam"), raster = F)
DefaultAssay(MLN_data) <- "RNA"
FeaturePlot_scCustom(MLN_data, reduction = "wnn.umap",features = c("Ptprc"), raster = F)
FeaturePlot_scCustom(MLN_data, reduction = "wnn.umap",features = c("Epcam"), raster = F) 
ggsave("MLN_Epcam_expression_WNNumap.png",plot = last_plot(), device = NULL, path = images,
  scale = 1, width = 5,height = 5, units = c("in"), dpi = 600, limitsize = TRUE, bg = NULL)
``` 


```{r}
#Cluster Method #2 Using Monocle
library(SeuratWrappers)
 cds <- as.cell_data_set(MLN_data, reductions = "wnn.umap",
  default.reduction = "wnn.umap", graph =  "wsnn", group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(cds,
  reduction_method = "UMAP", resolution =  NULL, k = 15, cluster_method =  "leiden",
  num_iter = 2, partition_qval = 0.25, weight = TRUE, verbose = T )
```


```{r}
#Visualize clusters Monocle - I don't find the partitions useful and don't believe they make sense  
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1|p2 

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_data[["Monocle_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_data[["Monocle_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(MLN_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "Monocle_clusters"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
```

```{r}
#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)
```

```{r}
#WNN with the Ileum data

DefaultAssay(Ileum_data) <- "RNA"
Ileum_data = FindMultiModalNeighbors(
  Ileum_data, reduction.list = list("integrated.harmony", "integrated.ADT.harmony"), 
  dims.list = list(1:30, 1:20), 
  modality.weight.name = c("RNA.weight", "ADT.weight"), 
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn.harmony",
  verbose = FALSE)

Ileum_data <- RunUMAP(Ileum_data, nn.name = "weighted.nn.harmony", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Ileum_data <- FindClusters(Ileum_data, graph.name = "wsnn",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters", n.start = 25, group.singletons = F)

DimPlot(Ileum_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "Tissue", "InfectionStatus"),
  combine = FALSE, label.size = 2, label = TRUE)
```

```{r}
#Evaluate Cluster Genes
Idents(Ileum_data) <- "WNN_clusters"
FeaturePlot_scCustom(Ileum_data, reduction = "wnn.umap",features = c("Ighd",  "Cd3e", "Cd14", "Cd4", "Ptprc", "Cd8a", "Ly6c1", "Epcam"), raster = F)
DefaultAssay(Ileum_data) <- "RNA"
FeaturePlot_scCustom(Ileum_data, reduction = "wnn.umap",features = c("Ptprc"), raster = F)
FeaturePlot_scCustom(Ileum_data, reduction = "wnn.umap",features = c("Epcam"), raster = F) 
ggsave("Ileum_Epcam_expression_WNNumap.png",plot = last_plot(), device = NULL, path = images,
  scale = 1, width = 5,height = 5, units = c("in"), dpi = 600, limitsize = TRUE, bg = NULL)
``` 


```{r}
Ileum_data
#Cluster Method #2 Using Monocle
library(SeuratWrappers)
 cds <- as.cell_data_set(Ileum_data, reductions = "wnn.umap",
  default.reduction = "wnn.umap", graph =  "wsnn", group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(cds,
  reduction_method = "UMAP", resolution =  0.00015, k = 15, cluster_method =  "leiden",
  num_iter = 2, partition_qval = 0.25, weight = TRUE, verbose = T )
```


```{r}
#Visualize clusters Monocle - I don't find the partitions useful and don't believe they make sense  
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1|p2 

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_data[["Monocle_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_data[["Monocle_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters"),
  combine = T, label.size = 2, label = TRUE, raster = F) + NoLegend()
DimPlot(Ileum_data, reduction = "wnn.umap",   group.by = c( "Monocle_clusters"),
  combine = T, label.size = 2, label = TRUE, raster = F) + NoLegend()
```

```{r}
#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)
```


#Advanced Quality Control - Cell Cycle Effects and doublet removal
```{r}
#Visualize the CEll cycle across the UMAP
DimPlot(MLN_data, reduction = "wnn.umap",   group.by = c(  "Phase"),
  combine = T, label.size = 2, label = TRUE, raster = F) + NoLegend()
FeaturePlot(MLN_data, reduction = "wnn.umap",features = c("Ighd",  "Cd3e", "Cd14", "Cd4", 
    "Cd8a", "Ms4a1", ), raster = F)
VlnPlot_scCustom(MLN_data, "Malat1", group.by = "Monocle_clusters", pt.size = 0) + NoLegend()
DimPlot(Ileum_data, reduction = "wnn.umap",   group.by = c( "Phase"),
  combine = T, label.size = 2, label = TRUE, raster = F)+ NoLegend()
FeaturePlot(Ileum_data, reduction = "wnn.umap",features = c( "Cd3e", "Cd14", "Cd4", "Ptprc",
    "Cd8a", "Ms4a1", "Epcam"), raster = F)
VlnPlot_scCustom(Ileum_data, "Malat1", group.by = "Monocle_clusters", pt.size = 0) + NoLegend()
```




```{r}
# Doublet removal
library(SeuratDisk)
library(Seurat)
library(BiocNeighbors) #manually loaded version2.0.1
library(scDblFinder)

#Doublet Identification
MLN_data[["RNA"]] <- JoinLayers(MLN_data[["RNA"]])
MLN_data.sce <- as.SingleCellExperiment(MLN_data, assay = "RNA")

MLN_data.sce <- scDblFinder(MLN_data.sce, samples="orig.ident", clusters="Monocle_clusters")
RDS <- as.Seurat(MLN_data.sce)
MLN_data <- JoinLayers(MLN_data)
MLN_data[["RNA"]] <-  as(MLN_data[["RNA"]], "Assay5")
MLN_data[["ADT"]] <-  as(MLN_data[["ADT"]], "Assay5")

```

```{r}
#Visualize the results of the doublet call
table(MLN_data.sce$scDblFinder.class) #189561 singlets 24388 Doublets
DimPlot(RDS, reduction = "WNN.UMAP", group.by = "scDblFinder.class")
DimPlot(RDS, reduction = "WNN.UMAP", group.by = c("scDblFinder.class", "Monocle_clusters"),
  combine = TRUE, label.size = 3, label = TRUE) + NoLegend()
DimPlot(RDS, reduction = "WNN.UMAP", group.by = c("scDblFinder.class", "Phase"),
  combine = TRUE, label.size = 3, label = TRUE) + NoLegend()
FeaturePlot_scCustom(RDS, reduction = "WNN.UMAP", features = c("Cd19", "Cd4", "Cd8a", "Ly6g") ) #overlapping CD19 and CD8/CD4 demonstrate true doublet call
```

```{r}
#Transfer Doublet labels to Seurat object and remove doublets from the data set
MLN_data[["Doublet"]] <-RDS[["scDblFinder.class"]] 
MLN_data[["scDblFinder.weighted"]] <-RDS[["scDblFinder.weighted"]] 
MLN_data[["scDblFinder.score"]] <-RDS[["scDblFinder.score"]]
MLN_data[["scDblFinder.difficulty"]] <-RDS[["scDblFinder.difficulty"]]
MLN_data[["scDblFinder.cxds_score"]] <-RDS[["scDblFinder.cxds_score"]]
library(qs2)
qs_save(MLN_data, "/path/to/analysis/directory/Seurat_Files/MLN_analysis1_w_doublets.qs2")
#saveRDS(MLN_data, file = "/path/to/analysis/directory/Seurat_Files/MLN_analysis1_w_doublets.rds")
```


```{r}
#Remove Doublets
MLN_data <-qs2::qs_read("/path/to/analysis/directory/Seurat_Files/MLN_analysis1_w_doublets.qs2")
Idents(MLN_data) <- "Doublet"

doublets_mln <- WhichCells(MLN_data, idents = c("doublet"))
write.csv(doublets_mln,  paste0(CSV, "/MLN_doubletCellNames.csv"))
write.csv(table(MLN_data$Doublet), paste0(CSV, "/MLN_Doublet_Calls.csv"))
MLN_data <- subset(x = MLN_data, idents = "singlet")
remove(RDS) 
remove(MLN_data.sce)
remove(RNA_only)
table(MLN_data$Doublet) 
```


```{r}
#Save the Object WITH Doublets

#Doublet Identification
Ileum_data[["RNA"]] <- JoinLayers(Ileum_data[["RNA"]])
Ileum_data.sce <- as.SingleCellExperiment(Ileum_data, assay = "RNA")
library(scDblFinder)
Ileum_data.sce <- scDblFinder(Ileum_data.sce, samples="orig.ident", clusters="Monocle_clusters")
RDS <- as.Seurat(Ileum_data.sce)
#Visualize the results of the doublet call
table(Ileum_data.sce$scDblFinder.class) #327334 singlets 33843 Doublets
DimPlot(RDS, reduction = "WNN.UMAP", group.by = "scDblFinder.class")
DimPlot(RDS, reduction = "WNN.UMAP", group.by = c("scDblFinder.class", "WNN_clusters"),
  combine = TRUE, label.size = 3, label = TRUE) + NoLegend()
DimPlot(RDS, reduction = "WNN.UMAP", group.by = c("scDblFinder.class", "Phase"),
  combine = TRUE, label.size = 3, label = TRUE) + NoLegend()
FeaturePlot(RDS, reduction = "WNN.UMAP", features = c("Cd19", "Cd4", "Cd8a", "Epcam", "Cd3e") ) #overlapping CD19 and CD8/CD4 demonstrate true doublet call
#Ileum_data[["RNA"]] <-  as(Ileum_data[["RNA"]], "Assay5")
#Ileum_data[["ADT"]] <-  as(Ileum_data[["ADT"]], "Assay5")
#SaveH5Seurat(Ileum_data, filename = "/path/to/analysis/directory/Seurat_Files/Ileum_analysis1_w_doublets.h5Seurat", overwrite = T)

```

```{r}
#Transfer Doublet labels to Seurat object and remove doublets from the data set
Ileum_data[["Doublet"]] <-RDS[["scDblFinder.class"]] 
Ileum_data[["scDblFinder.weighted"]] <-RDS[["scDblFinder.weighted"]] 
Ileum_data[["scDblFinder.score"]] <-RDS[["scDblFinder.score"]]
Ileum_data[["scDblFinder.difficulty"]] <-RDS[["scDblFinder.difficulty"]]
Ileum_data[["scDblFinder.cxds_score"]] <-RDS[["scDblFinder.cxds_score"]]

qs_save(Ileum_data, "/path/to/analysis/directory/Seurat_Files/Ileum_analysis1_w_doublets.qs2")
#Remove Doublets
Idents(Ileum_data) <- "Doublet"
table(Ileum_data$Doublet) 
doublets_ileum <- WhichCells(Ileum_data, idents = c("doublet"))
write.csv(doublets_ileum,  paste0(CSV, "/ILEUM_doubletCellNames.csv"))
write.csv(table(Ileum_data$Doublet), paste0(CSV, "/Ileum_Doublet_Calls.csv"))
Ileum_data <- subset(x = Ileum_data, idents = "singlet")
remove(RDS) 
remove(Ileum_data.sce)
remove(RNA_only)
table(Ileum_data$Doublet) 
```

```{r}
#Save the Cell Numbers 
write.csv(table(Ileum_data$Tissue, Ileum_data$orig.ident, Ileum_data$InfectionStatus) , paste0(CSV, "/Ileum_Cell_Breakdown.csv"))
table(Ileum_data$Tissue, Ileum_data$orig.ident, Ileum_data$InfectionStatus) 
write.csv(table(MLN_RNA$Tissue, MLN_RNA$orig.ident, MLN_RNA$InfectionStatus) , paste0(CSV, "/MLN_Cell_Breakdown.csv"))
```

#Recluster following doublet removal
```{r}
DefaultAssay(MLN_data) <- "RNA"

MLN_data = FindMultiModalNeighbors(
  MLN_data, reduction.list = list("integrated.harmony", "integrated.ADT.harmony"), 
  dims.list = list(1:30, 1:20), 
  modality.weight.name = c("RNA.weight", "ADT.weight"), 
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn_harmony",
  verbose = T)

MLN_data <- RunUMAP(MLN_data, nn.name = "weighted.nn_harmony", reduction.name = "wnn.umap", reduction.key = "wnnUMAP", seed.use = 12)
MLN_data <- FindClusters(MLN_data, graph.name = "wsnn",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters", n.start = 20, group.singletons = F)
table(MLN_data$WNN_clusters)
DimPlot(MLN_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "InfectionStatus", "Phase"),
  combine = FALSE, label.size = 2, label = TRUE)
```

```{r}
#Monocle Clustering
library(SeuratWrappers)
 cds <- as.cell_data_set(
  MLN_data,
  reductions = "wnn.umap",
  default.reduction = "wnn.umap",
  graph =  "wsnn",
  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(
  cds,
  reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.00015,
  num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )
```


```{r}  
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1|p2 

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_data[["Monocle_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_data[["Monocle_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(MLN_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "Monocle_clusters"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)
```


```{r}

DefaultAssay(Ileum_data) <- "RNA"
Ileum_data = FindMultiModalNeighbors(
  Ileum_data, reduction.list = list("integrated.harmony", "integrated.ADT.harmony"), 
  dims.list = list(1:30, 1:20), 
  modality.weight.name = c("RNA.weight", "ADT.weight"), 
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn.harmony",
  verbose = FALSE)

Ileum_data <- RunUMAP(Ileum_data, nn.name = "weighted.nn.harmony", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 12, n.neighbors = 15L)
Ileum_data <- FindClusters(Ileum_data, graph.name = "wsnn",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters", n.start = 20, group.singletons = F)

table(Ileum_data$WNN_clusters)
DimPlot(Ileum_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "Tissue", "InfectionStatus", "Phase"),
  combine = FALSE, label.size = 2, label = TRUE)
```

```{r}
library(SeuratWrappers)
 cds <- as.cell_data_set(
  Ileum_data,
  reductions = "wnn.umap",
  default.reduction = "wnn.umap",
  graph =  "wsnn",
  group.by = NULL)
 

reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells(
  cds,
  reduction_method = "UMAP", 
  k = 15, cluster_method =  "leiden", resolution = 0.00015,
  num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )
```


```{r}
#Visualize clusters Monocle - I don't find the partitions useful and don't believe they make sense  
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1|p2 

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_data[["Monocle_clusters"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_data[["Monocle_partitions"]] <- Seurat_monocle@meta.data$monocle3_partitions
DimPlot(Ileum_data, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "Monocle_clusters"),
  combine = FALSE, label.size = 2, label = TRUE, raster = F)
FeaturePlot(Ileum_data, "Saa1", reduction = "wnn.umap",)
#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)
```

#Save Data 
```{r}
library(qs2)
qs_save(MLN_data, "/path/to/analysis/directory/Seurat_Files/Analysis1_MLN_clustered.qs2")
#SaveH5Seurat(Ileum_RNA, filename = "/path/to/analysis/directory/Seurat_Files/Analysis1_Ileum_clustered.h5Seurat", overwrite = T)
qs_save(Ileum_data, "/path/to/analysis/directory/Seurat_Files/Analysis1_Ileum_clustered.qs2")
```

#Regress out Cell cycle MLN
```{r}
#The following was done separately in an R Script 'MLN_ccRegression.R'
# nohup Rscript ./MLN_ccRegression.R > ./MLN_regression.R.log 2>&1 &
MLN_data <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis1_MLN_clustered.qs2")

MLN_data[["RNA"]] <- as(MLN_data[["RNA"]], "Assay")

DefaultAssay(MLN_data) <- "RNA"
options(future.globals.maxSize = 400000 * 1024^2)
plan(strategy =  "multisession", workers = 30 )
all.genes <- rownames(MLN_data[["RNA"]]) 
MLN_data <-  ScaleData(MLN_data, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes, assay = "RNA", verbose = T, block.size = 11000, min.cells.to.block = 5000)

qs_save(MLN_data, "/path/to/analysis/directory/Seurat_Files/Analysis1_MLN_clustered_regressed_int.qs2")

MLN_data[["RNA"]] <- as(MLN_data[["RNA"]], "Assay5")
MLN_data <- RunPCA(MLN_data, verbose = FALSE, reduction.name = "RNA_pca_cc", features = VariableFeatures(MLN_data, nfeatures = 3000))

MLN_data[['RNA']] <-  split(x = MLN_data[['RNA']], f = MLN_data$InfectionStatus)
MLN_data <- IntegrateLayers(object = MLN_data, method = HarmonyIntegration, assay = "RNA", features = VariableFeatures(MLN_data, nfeatures = 3000),
  orig.reduction = "RNA_pca_cc", new.reduction = "reintegrated_RNA_harmony_cc", verbose = T)

MLN_data <- FindNeighbors(MLN_data, reduction = "reintegrated_RNA_harmony_cc", dims = 1:30, graph.name = "harmony_rna_cc")
MLN_data <- FindClusters(MLN_data, resolution = 2, cluster.name = "harmony_integration_clusters_cc",  graph.name = "harmony_rna_cc", group.singletons = F, n.start = 20)
MLN_data <- RunUMAP(MLN_data, reduction = "reintegrated_RNA_harmony_cc", dims = 1:30, reduction.name = "umap.harmony_cc",  seed.use = 12)


MLN_data = FindMultiModalNeighbors(
  MLN_data, reduction.list = list("reintegrated_RNA_harmony_cc", "integrated.ADT.harmony"), 
  dims.list = list(1:30, 1:20), 
  modality.weight.name = c("RNA.weight", "ADT.weight"), 
  knn.graph.name = "wknn_cc",
  snn.graph.name = "wsnn_cc",
  weighted.nn.name = "weighted.nn_cc",
  verbose = FALSE)

#compare clusters and Umaps before and after CC regression
DefaultAssay(MLN_data) <- "RNA"
MLN_data <- RunUMAP(MLN_data, nn.name = "weighted.nn_cc", reduction.name = "wnn.umap_cc", reduction.key = "wnnUMAPcc_")
MLN_data <- FindClusters(MLN_data, graph.name = "wsnn_cc",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters_cc", n.start = 25, group.singletons = F)


#Monocle Clusters
library(SeuratWrappers)
MLN_data <- JoinLayers(MLN_data)

 cds <- as.cell_data_set(
  MLN_data,
  reductions = "wnn.umap_cc",
  default.reduction = "wnn.umap_cc",
  graph =  "wsnn",
  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells( cds, reduction_method = "UMAP", k = 15, cluster_method =  "leiden", resolution = 0.00015,
  num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
MLN_data[["Monocle_clusters_cc"]] <- Seurat_monocle@meta.data$monocle3_clusters
MLN_data[["Monocle_partitions_cc"]] <- Seurat_monocle@meta.data$monocle3_partitions

#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)
#Save
qs_save(MLN_data, "/path/to/analysis/directory/Seurat_Files/Analysis1_MLN_clustered_regressed.qs2")

```


#Regress out Cell cycle Ileum
```{r}
#The following was done separately in an R Script 'Ileum_ccRegression.R'
# nohup Rscript ./Ileum_ccRegression.R > ./Ileum_regression.R.log 2>&1 &
Ileum_data <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis1_Ileum_clustered.qs2")

Ileum_data[["RNA"]] <- as(Ileum_data[["RNA"]], "Assay")

DefaultAssay(Ileum_data) <- "RNA"
options(future.globals.maxSize = 500000 * 1024^2)
plan(strategy =  "multisession", workers = 30 )
all.genes <- rownames(Ileum_data[["RNA"]]) 
Ileum_data <-  ScaleData(Ileum_data, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes, assay = "RNA", verbose = T, block.size = 11000, min.cells.to.block = 5000)

qs_save(Ileum_data, "/path/to/analysis/directory/Seurat_Files/Analysis1_Ileum_clustered_regressed_int.qs2")

Ileum_data[["RNA"]] <- as(Ileum_data[["RNA"]], "Assay5")
Ileum_data <- RunPCA(Ileum_data, verbose = FALSE, reduction.name = "RNA_pca_cc", features = VariableFeatures(Ileum_data, nfeatures = 3000))
Ileum_data[["RNA"]] <- JoinLayers(Ileum_data[['RNA']])
Ileum_data[['RNA']] <-  split(x = Ileum_data[['RNA']], f = Ileum_data$Infection_Organ)
Ileum_data <- IntegrateLayers(object = Ileum_data, method = HarmonyIntegration, assay = "RNA", features = VariableFeatures(Ileum_data, nfeatures = 3000),
  orig.reduction = "RNA_pca_cc", new.reduction = "reintegrated_RNA_harmony_cc", verbose = T)

Ileum_data <- FindNeighbors(Ileum_data, reduction = "reintegrated_RNA_harmony_cc", dims = 1:35, graph.name = "harmony_rna_cc")
Ileum_data <- FindClusters(Ileum_data, resolution = 2, cluster.name = "harmony_integration_clusters_cc",  graph.name = "harmony_rna_cc", group.singletons = F, n.start = 20)
Ileum_data <- RunUMAP(Ileum_data, reduction = "reintegrated_RNA_harmony_cc", dims = 1:35, reduction.name = "umap.harmony_cc",  seed.use = 20)


Ileum_data = FindMultiModalNeighbors(
  Ileum_data, reduction.list = list("reintegrated_RNA_harmony_cc", "integrated.ADT.harmony"), 
  dims.list = list(1:35, 1:20), 
  modality.weight.name = c("RNA.weight", "ADT.weight"), 
  knn.graph.name = "wknn_cc",
  snn.graph.name = "wsnn_cc",
  weighted.nn.name = "weighted.nn_cc",
  verbose = FALSE)

#compare clusters and Umaps before and after CC regression
DefaultAssay(Ileum_data) <- "RNA"
Ileum_data <- RunUMAP(Ileum_data, nn.name = "weighted.nn_cc", reduction.name = "wnn.umap_cc", reduction.key = "wnnUMAPcc_", seed.use = 20)
Ileum_data <- FindClusters(Ileum_data, graph.name = "wsnn_cc",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters_cc", n.start = 25, group.singletons = F)


#Monocle Clusters
library(SeuratWrappers)
Ileum_data <- JoinLayers(Ileum_data)

 cds <- as.cell_data_set(
  Ileum_data,
  reductions = "wnn.umap_cc",
  default.reduction = "wnn.umap_cc",
  graph =  "wsnn",
  group.by = NULL)
 
reducedDimNames(cds) <- "UMAP"
cds <-cluster_cells( cds, reduction_method = "UMAP", k = 15, cluster_method =  "leiden", resolution = 0.00015,
  num_iter = 1, partition_qval = 0.25, weight = TRUE, verbose = T )

Seurat_monocle <- as.Seurat(cds, assay = NULL)
Ileum_data[["Monocle_clusters_cc"]] <- Seurat_monocle@meta.data$monocle3_clusters
Ileum_data[["Monocle_partitions_cc"]] <- Seurat_monocle@meta.data$monocle3_partitions

#Delete the temporary Seurat monocle object
remove(Seurat_monocle)
remove(cds)
#Save
qs_save(Ileum_data, "/path/to/analysis/directory/Seurat_Files/Analysis1_Ileum_clustered_regressed.qs2")
```
