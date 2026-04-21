# This standalon3 script eprforms cell type deconvolution on the Stereoseq data

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(forcats)
library(ggbreak)
library(presto)
library(ggforce)
library(SingleR)
library(glmGamPoi)
library(future)
library(org.Mm.eg.db)
library(Matrix)
library(jsonlite)
library(patchwork)
library(tidyr)
library(spacexr)
library(qs2)

# Perform a deconvolution and cell type annotation 

load("/path/to/directory/STOmics/Intestine/Naive_yersinia.RData")
Ileum <- qs2::qs_read("/path/to/directory/Seurat_Files/Analysis2.5_Ileum.qs2")

counts_ref <- GetAssayData(Ileum, assay = "RNA", slot = "counts")

cluster_labels <- factor(Ileum$secondlevel)
unique(cluster_labels)
# RCTD automatically handles the signature matrix generation
reference <- spacexr::Reference(counts_ref, cluster_labels)

# raw counts from the object - currently using cellbinned data
counts_spatial <- GetAssayData(combined_stereo_filt, assay = "Spatial", slot = "counts")

# pull out the cell bin centers
# In SAW-generated Seurat objects, coordinates are in the @images slot
image_names <- names(combined_stereo_filt@images) 

# Pull coordinates for each image 
all_coords <- do.call(rbind, lapply(image_names, function(x) {
  df <- GetTissueCoordinates(combined_stereo_filt, image = x)
  return(df)
}))

print(nrow(all_coords))
table(combined_stereo_filt$InfectionStatus)
# format coordinates, rctd needs x and y
colnames(all_coords) <- c("x", "y") 

# Create spatial RNA objects
puck <- spacexr::SpatialRNA(all_coords, counts_spatial)



# RCTD process
myRCTD <- spacexr::create.RCTD(puck, reference)
myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = 'doublet')# doublet_mode = 'doublet' to detect cell-state mixing

save.image("/path/to/directory/STOmics/Intestine/Naive_yersinia_rctd.RData")
# pull out the weights/ cell types probabilities
weights <- myRCTD@results$weights
norm_weights <- weights / rowSums(weights)
# Add to seurat object
# This adds a column for every cell type in your reference
combined_stereo_filt <- AddMetaData(combined_stereo_filt, metadata = as.data.frame(as.matrix(norm_weights)))

# New metadata of the best predicted cell type
cell_type_names <- colnames(norm_weights)
combined_stereo_filt$RCTD_prediction <- cell_type_names[max.col(norm_weights)]
 
saveRDS(myRCTD, "RCTD_CellBin_Results.rds")
saveRDS(combined_stereo_filt, "CellBin_Deconvolved.rds")
save.image("/path/to/directory/STOmics/Intestine/Naive_yersinia_rctd.RData")
