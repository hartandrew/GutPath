#Load Libraries 

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
library(glmGamPoi) 
library(future)
library(dsb) 
library(SeuratDisk) 
library(RColorBrewer)
library(monocle3)
library(gt)
library(SeuratWrappers)
library(grid)
library(gridExtra)
library(qs2)
library(tidyr)
library(tibble)
options(future.globals.maxSize = 5000 * 1024^2)

#Establish Directories
getwd() #/home/hartandrew
setwd("/path/to/analysis/directory")
out <- "/path/to/analysis/directory/Output"
version <- "/path/to/analysis/directory/Version"
seurat <- "/path/to/analysis/directory/Seurat_Files"
images <- "/path/to/analysis/directory/Images"
CSV <- "/path/to/analysis/directory/CSV"


#Load Objects

Ileum_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2.5_Ileum.qs2")

# Compare Each Cell Type for Each condition

library(MAST)
Idents(Ileum_filt) <- "SampleType"
Ileum_filt_i <-  subset(Ileum_filt, idents = c("Ileum Intestinal Epithelial Cells", "Ileum Lamina Propria")) #remove the duodenum
#For Every condition, compare each cell type against the naive. Limit this to cell types where there are at least 10 cells in each condition
all_conditions <- setdiff(unique(Ileum_filt_i$InfectionStatus), "Naive")
all_celltypes <- unique(Ileum_filt_i$secondlevel)

all_degs_list <- list()
results_list <- list()

 for (cond in all_conditions) {
  for (ct in all_celltypes) {
    print(paste0("starting DEG for ", cond, " ", ct))
    cell_counts <- table(Ileum_filt_i$secondlevel, Ileum_filt_i$InfectionStatus)
    
    # Check that both condition and Naive have at least one cell
    if (ct %in% rownames(cell_counts) &&
        cond %in% colnames(cell_counts) &&
        "Naive" %in% colnames(cell_counts) &&
        cell_counts[ct, cond] > 0 &&
        cell_counts[ct, "Naive"] > 0) {
      
      # Subset to relevant cell type and condition
      subset_obj <- subset(Ileum_filt_i, subset = secondlevel == ct & InfectionStatus %in% c("Naive", cond))
      counts <- table(subset_obj$InfectionStatus)
      Idents(subset_obj) <- "InfectionStatus"
      subset_obj <- subset(subset_obj, downsample = 5000)
        if (all(c("Naive", cond) %in% names(counts)) && all(counts[c("Naive", cond)] >= 10)) {
        
        Idents(subset_obj) <- "InfectionStatus"
        
        degs <- FindMarkers(
          subset_obj,
          ident.1 = cond,
          ident.2 = "Naive",
          group.by = "InfectionStatus",
          logfc.threshold = 0,
          test.use = "MAST",
          min.pct = 0.1, 
          latent.vars = c("nFeature_RNA", "nCount_RNA", "rb.prop", "mt.prop")
        )
        
        degs <- degs %>%
          rownames_to_column("gene") %>%
          mutate(condition = cond, celltype = ct)
        
        all_degs_list[[paste(cond, ct, sep = "_")]] <- degs  # Store raw DEGs (for global correction)
        
      } else {
        print(paste0("Not enough cells in ", cond, " ", ct))
        results_list[[paste(cond, ct, sep = "_")]] <- tibble(
          condition = cond,
          celltype = ct,
          total_deg = NA,
          up_deg = NA,
          down_deg = NA
        )
      }
    } else {
      print(paste0("Skipping ", cond, " ", ct, " — no cells found"))
      results_list[[paste(cond, ct, sep = "_")]] <- tibble(
        condition = cond,
        celltype = ct,
        total_deg = NA,
        up_deg = NA,
        down_deg = NA
      )
    }
  }
}
# Combine all DEGs and apply global p correction BH
all_degs_df <- bind_rows(all_degs_list)
all_degs_df$global_padj <- p.adjust(all_degs_df$p_val, method = "BH")
write.csv(all_degs_df, "/path/to/analysis/directory/CSV/all_degs_df.csv")
