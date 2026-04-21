#Load Libraries 

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


#Load Gene set 
genesets_Hallmark<- getGmt("/path/to/data/Genesets/Mouse_Hallmark/Hallmark.v2023.Mm.gmt", geneIdType=SymbolIdentifier())
genesets_Th2<- getGmt("/path/to/data/Genesets/Th2_Gata3.Mm.gmt", geneIdType=SymbolIdentifier())

#Load Data 
Ileum <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2_clustered_annotated_Ileum.qs2")

#Run the analysis 
#Ileum[["RNA"]] <- JoinLayers(Ileum[["RNA"]])

Ileum[["RNA"]] <- as(Ileum[["RNA"]], "Assay")
print(paste0("Starting enrichment"))
enrichment.scores <- escape.matrix(Ileum, 
                                   gene.sets = genesets_Hallmark, 
                                   groups = 9000, 
                                   min.size = 5,
                                   method = "UCell", 
                                   normalize = FALSE)

                                   


# Convert rownames to a column
enrichment.scores.out <- as.data.frame(enrichment.scores)
enrichment.scores.out$cell_id <- rownames(enrichment.scores.out)

# Move cell_id to the first column (optional)
enrichment.scores.out <- enrichment.scores.out[, c("cell_id", setdiff(colnames(enrichment.scores.out), "cell_id"))]

# Save to CSV with cell names
fwrite(enrichment.scores.out, "/path/to/analysis/directory/CSV/Hallmark_ssGSEA_UCell.csv")


enrichment.scores <- escape.matrix(Ileum, 
                                   gene.sets = genesets_Th2, 
                                   groups = 9000, 
                                   min.size = 5,
                                   method = "UCell", 
                                   normalize = FALSE)

                                   


# Convert rownames to a column
enrichment.scores.out <- as.data.frame(enrichment.scores)
enrichment.scores.out$cell_id <- rownames(enrichment.scores.out)

# Move cell_id to the first column (optional)
enrichment.scores.out <- enrichment.scores.out[, c("cell_id", setdiff(colnames(enrichment.scores.out), "cell_id"))]

# Save to CSV with cell names
fwrite(enrichment.scores.out, "/path/to/analysis/directory/CSV/Th2_ssGSEA_Ucell.csv")





genesets_Th17<- getGmt("/path/to/data/Genesets/GOBP_T_HELPER_17_CELL_DIFFERENTIATION.v2025.1.Mm.gmt", geneIdType=SymbolIdentifier())
enrichment.scores <- escape.matrix(Ileum, 
                                   gene.sets = genesets_Th17, 
                                   groups = 9000, 
                                   min.size = 5,
                                   method = "UCell", 
                                   normalize = FALSE)

                                   


# Convert rownames to a column
enrichment.scores.out <- as.data.frame(enrichment.scores)
enrichment.scores.out$cell_id <- rownames(enrichment.scores.out)

# Move cell_id to the first column (optional)
enrichment.scores.out <- enrichment.scores.out[, c("cell_id", setdiff(colnames(enrichment.scores.out), "cell_id"))]

# Save to CSV with cell names
fwrite(enrichment.scores.out, "/path/to/analysis/directory/CSV/Th17_ssGSEA_Ucell.csv")



#Load Gene set 
genesets_Hallmark<- getGmt("/path/to/data/Genesets/Mouse_Hallmark/Hallmark.v2023.Mm.gmt", geneIdType=SymbolIdentifier())
genesets_Th2<- getGmt("/path/to/data/Genesets/Th2_Gata3.Mm.gmt", geneIdType=SymbolIdentifier())

#Load Data 
MLN <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2_MLN_clustered_annotated.qs2")

#Run the analysis 
MLN[["RNA"]] <- JoinLayers(MLN[["RNA"]])

MLN[["RNA"]] <- as(MLN[["RNA"]], "Assay")
print(paste0("Starting enrichment"))
enrichment.scores <- escape.matrix(MLN, 
                                   gene.sets = genesets_Hallmark, 
                                   groups = 9000, 
                                   min.size = 5,
                                   method = "UCell", 
                                   normalize = FALSE)

                                   


# Convert rownames to a column
enrichment.scores.out <- as.data.frame(enrichment.scores)
enrichment.scores.out$cell_id <- rownames(enrichment.scores.out)

# Move cell_id to the first column (optional)
enrichment.scores.out <- enrichment.scores.out[, c("cell_id", setdiff(colnames(enrichment.scores.out), "cell_id"))]

# Save to CSV with cell names
fwrite(enrichment.scores.out, "/path/to/analysis/directory/CSV/Hallmark_ssGSEA_UCell_MLN.csv")


enrichment.scores <- escape.matrix(MLN, 
                                   gene.sets = genesets_Th2, 
                                   groups = 9000, 
                                   min.size = 5,
                                   method = "UCell", 
                                   normalize = FALSE)

                                   


# Convert rownames to a column
enrichment.scores.out <- as.data.frame(enrichment.scores)
enrichment.scores.out$cell_id <- rownames(enrichment.scores.out)

# Move cell_id to the first column (optional)
enrichment.scores.out <- enrichment.scores.out[, c("cell_id", setdiff(colnames(enrichment.scores.out), "cell_id"))]

# Save to CSV with cell names
fwrite(enrichment.scores.out, "/path/to/analysis/directory/CSV/Th2_ssGSEA_Ucell_MLN.csv")





genesets_Th17<- getGmt("/path/to/data/Genesets/GOBP_T_HELPER_17_CELL_DIFFERENTIATION.v2025.1.Mm.gmt", geneIdType=SymbolIdentifier())
enrichment.scores <- escape.matrix(MLN, 
                                   gene.sets = genesets_Th17, 
                                   groups = 9000, 
                                   min.size = 5,
                                   method = "UCell", 
                                   normalize = FALSE)

                                   


# Convert rownames to a column
enrichment.scores.out <- as.data.frame(enrichment.scores)
enrichment.scores.out$cell_id <- rownames(enrichment.scores.out)

# Move cell_id to the first column (optional)
enrichment.scores.out <- enrichment.scores.out[, c("cell_id", setdiff(colnames(enrichment.scores.out), "cell_id"))]

# Save to CSV with cell names
fwrite(enrichment.scores.out, "/path/to/analysis/directory/CSV/Th17_ssGSEA_Ucell_MLN.csv")

