# Purpose: This script will make marker gene profiles for both the MLN and the Ileum at high levels (coarse cell types) and refined levels (Finest cell types)
# Figure Panels: This script generates Supplemental Data files for marker gene identification pt 1
#Load Libraries ----

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
library(pheatmap)
library(grid)
library(gridExtra)
library(qs2)
library(tidyr)
library(tibble)
options(future.globals.maxSize = 5000 * 1024^2)


#Establish Directories-----
 
getwd() #/home/hartandrew
setwd("/path/to/analysis/directory")
out <- "/path/to/analysis/directory/Output"
version <- "/path/to/analysis/directory/Version"
seurat <- "/path/to/analysis/directory/Seurat_Files"
images <- "/path/to/analysis/directory/Images"
CSV <- "/path/to/analysis/directory/CSV"


#Set colors----

colors7 <- c("#980339", "#379FFB", "#D6A206", "#027102", "#6313B7", "#AFAFAF","#542600" )
colors6 <- c("#980339", "#379FFB", "#D6A206", "#027102", "#6313B7", "#AFAFAF")
colors5 <- c("#980339", "#379FFB", "#D6A206", "#027102", "#6313B7")
colors8 <- c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
colors9 <- c("#332288", "#88CCEE", "#117733","#44AA99", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
cbf_18_colors <- c(
  "#1170AA",   "#FCB13F",  "#60BFC1",  "#EF6F6A",  "#937860",  "#D17A00",  "#B78CBA",   "#B3B3B3",   "#64B200",  "#D45E00",  
  "#7E8CCF",  "#E6A0C4",  "#568B3F",  "#C44E52","#5FA2A3",  "#CCB974",  "#D0A6BE", "#4E84C4"  )

#Load Objects----

MLN_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2_MLN_clustered_annotated.qs2")
Ileum_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2.5_Ileum.qs2")


# Make Coarse and secondlevel marker Lists----
pval = 0.05
Idents(Ileum_filt) <- "CoarseCellType"

Ileum_Coarse_Markers <- FindAllMarkers(
  Ileum_filt,
  only.pos = TRUE,
  min.pct = 0.1,
  test.use = "wilcox",
  return.thresh = pval
)
Idents(Ileum_filt) <- "secondlevel"

Ileum_intermediate_Markers <- FindAllMarkers(
  Ileum_filt,
  only.pos = TRUE,
  min.pct = 0.1,
  test.use = "wilcox",
  return.thresh = pval
)

Idents(MLN_filt) <- "CoarseCellType"

MLN_Coarse_Markers <- FindAllMarkers(
  MLN_filt,
  only.pos = TRUE,
  min.pct = 0.1,
  test.use = "wilcox",
  return.thresh = pval
)
Idents(MLN_filt) <- "secondlevel"

MLN_intermediate_Markers <- FindAllMarkers(
  MLN_filt,
  only.pos = TRUE,
  min.pct = 0.1,
  test.use = "wilcox",
  return.thresh = pval
)

library(openxlsx)
wb <- createWorkbook()

addWorksheet(wb, "Ileum_Coarse")
writeData(wb, "Ileum_Coarse", Ileum_Coarse_Markers)

addWorksheet(wb, "Ileum_intermediate")
writeData(wb, "Ileum_intermediate", Ileum_intermediate_Markers)

addWorksheet(wb, "MLN_Coarse")
writeData(wb, "MLN_Coarse", MLN_Coarse_Markers)

addWorksheet(wb, "MLN_intermediate")
writeData(wb, "MLN_intermediate", MLN_intermediate_Markers)

saveWorkbook(wb, paste0(CSV,"/Coarse_secondlevel_MarkerGenes_FindMarkers.xlsx"), overwrite = TRUE)


# Generate Marker lists for the subclustered objects and the fine cell type definitions. In this case, the labels were defined on clustering the data after each coarse cell type was subset from the data.


get_finest_markers_within_coarse <- function(seurat_obj, coarse_col, finest_col,
                                             logfc, pval) {
  
  coarse_vals <- unique(seurat_obj[[coarse_col, drop = TRUE]])
  
  result_list <- lapply(coarse_vals, function(cc) {
    message("Processing coarse cell type: ", cc)
    
    # subset the coarse cell type
    sub_obj <- subset(seurat_obj, subset = !!as.name(coarse_col) == cc)
    
    # set identity to FinestCellType
    Idents(sub_obj) <- finest_col
    
    if (length(unique(Idents(sub_obj))) < 2) {
      message("Skipping ", cc, " — only one FinestCellType present.")
      return(NULL)
    }
    
    FindAllMarkers(
      sub_obj,
      only.pos = TRUE,
      min.pct = 0.1,
      test.use = "wilcox",
      min.diff.pct = 0.15
    )
  })
  
  names(result_list) <- coarse_vals
  return(result_list)
}

Ileum_FinestMarkers_byCoarse <- get_finest_markers_within_coarse(
  seurat_obj = Ileum_filt,
  coarse_col = "CoarseCellType",
  finest_col = "FinestCellType",
  pval = 0.05
)

MLN_FinestMarkers_byCoarse <- get_finest_markers_within_coarse(
  seurat_obj = MLN_filt,
  coarse_col = "CoarseCellType",
  finest_col = "FineCellType",
  pval = 0.05
)

library(openxlsx)
wb <- openxlsx::createWorkbook()
for (cc in names(Ileum_FinestMarkers_byCoarse)) {
  addWorksheet(wb, cc)
  writeData(wb, cc, Ileum_FinestMarkers_byCoarse[[cc]])
}
saveWorkbook(wb, paste0(CSV, "/Ileum_FinestMarkers_byCoarse.xlsx"), overwrite = TRUE)

wb <- openxlsx::createWorkbook()
for (cc in names(MLN_FinestMarkers_byCoarse)) {
  addWorksheet(wb, cc)
  writeData(wb, cc, MLN_FinestMarkers_byCoarse[[cc]])
}
saveWorkbook(wb, paste0(CSV, "/MLN_FinestMarkers_byCoarse.xlsx"), overwrite = TRUE)

# Check to make sure that each cell type produced significant markers 
check_finest_marker_coverage <- function(marker_list) {
  data.frame(
    CoarseCellType = names(marker_list),
    n_markers = sapply(marker_list, function(x) if(is.null(x)) 0 else nrow(x)),
    has_markers = sapply(marker_list, function(x) !is.null(x) && nrow(x) > 0)
  )
}
ileum_check <- check_finest_marker_coverage(Ileum_FinestMarkers_byCoarse)
mln_check   <- check_finest_marker_coverage(MLN_FinestMarkers_byCoarse)

ileum_check
mln_check


#Save the Seurat files as H5ad files to then use
library(reticulate)
use_python("/home/hartandrew/.conda/envs/sccellfie_env/bin/python", required = TRUE)
library(sceasy)


Ileum_filt[["RNA"]] <- JoinLayers(Ileum_filt[["RNA"]])
Ileum_filt[["ADT"]] <- JoinLayers(Ileum_filt[["ADT"]])
DefaultAssay(Ileum_filt) <- "RNA"
sceasy::convertFormat(Ileum_filt, from = "seurat", to = "anndata", main_layer = "counts",  outFile = "/path/to/analysis/directory/Seurat_Files/Analysis2.5_labeled_Ileum.h5ad")

MLN_filt[["RNA"]] <- JoinLayers(MLN_filt[["RNA"]])
MLN_filt[["RNA"]] <- as(MLN_filt[["RNA"]], "Assay")
MLN_filt[["ADT"]] <- JoinLayers(MLN_filt[["ADT"]])
DefaultAssay(MLN_filt) <- "RNA"
sceasy::convertFormat(MLN_filt, from = "seurat", to = "anndata", main_layer = "counts",  outFile = "/path/to/analysis/directory/Seurat_Files/Analysis2.5_labeled_MLN.h5ad")

