##### Load Libraries----
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(scCustomize)
library(forcats)
library(ggbreak)
library(presto)
library(glmGamPoi)
library(future)
library(dsb)
library(SeuratDisk)
library(RColorBrewer)
library(monocle3)
library(gt)
library(pheatmap)
library(circlize)
library(SeuratWrappers)
options(future.globals.maxSize = 50000 * 1024^2)
library(BiocParallel)
library(GSVA)
library(GSEABase)
library(qs2)
library(AUCell)
library(escape)
library(UCell)






#Perform GSEA ----

# Load your Xenium Seurat object (already loaded in your environment)
load("/path/to/data/MIST/Xenium_Yp_2025/Subclustered_Xenium.RData")

# Load Gene Sets
genesets_PRE <- getGmt("/path/to/data/Genesets/Saa1_geneList_Aug14_2025_Xenium.gmt", geneIdType=SymbolIdentifier())


# Subset the Seurat object
Epithelial <- subset_objects[["Epithelial cells filt"]]

# Extract the scaled expression matrix directly
expr_matrix <- as.matrix(Epithelial@assays$SCT@scale.data)

# Load gene sets
genesets_PRE <- getGmt("/path/to/data/Genesets/Saa1_geneList_Aug14_2025_Xenium.gmt", geneIdType = SymbolIdentifier())

# Run ssGSEA
enrichment.scores <- escape.matrix(expr_matrix, 
                                   gene.sets = genesets_PRE, 
                                   groups = 9000, 
                                   min.size = 5,
                                   method = "UCell", 
                                   normalize = FALSE)

# Convert to data.frame and save
enrichment.scores.out <- as.data.frame(enrichment.scores)
enrichment.scores.out$cell_id <- rownames(enrichment.scores.out)
enrichment.scores.out <- enrichment.scores.out[, c("cell_id", setdiff(colnames(enrichment.scores.out), "cell_id"))]

fwrite(enrichment.scores.out, "/path/to/data/MIST/Xenium_Yp_2025/CSV/UCell_PRE_202508.csv")