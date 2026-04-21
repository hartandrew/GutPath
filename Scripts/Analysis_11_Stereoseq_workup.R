# Description -----
# We have performed probe-based Xenium analysis on intestinal tissue and we have more recently partnered with Complete Genomics to perform Stereoseq analyses of intestinal tissue as well as brain tissue. The SAW pipeline produced by Complete Genomics is used for analysis of STEREOSeq samples. I will attempt to read in the data and explore it
# Figures S9G
# Load libraries ----
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
library(harmony)
library(glmGamPoi)
library(future)
library(RColorBrewer)
library(org.Mm.eg.db)
library(Matrix)
library(jsonlite)
library(patchwork)
library(tidyr)

# Select color palettes
c8_a <- c( "#000000",
  "#E69F00", 
  "#56B4E9",  
  "#009E73",  
  "#F0E442",  
  "#0072B2",  
  "#D55E00", 
  "#CC79A7"  )
c8_b <- c( "#1B9E77",
  "#D95F02",
  "#7570B3",
  "#E7298A",
  "#66A61E",
  "#E6AB02",
  "#A6761D",
  "#666666")

# Create an output directory
images <- "/path/to/directory/Projects/STOmics/Intestine/Images"

# Load the Yersinia and Naive data. The data were generated using the adjusted cell bins and specific notes can be found in the txt file "Notes_2_Stereoseq_methods_pt1_manuscript.txt" for how seurat compatible objects were made ----


# The RDS data is generated following a generous pipeline [https://github.com/limin321/addImg2annData/]. The data creates a Visium-like object 
Yp_SI_stereo <- readRDS("/path/to/directory/SAW_pipeline/runs/Yp_SI_2025_custom/format_conversion/addimage/Y01679L2.addimg.rds")
Naive_stereo <- readRDS("/path/to/directory/SAW_pipeline/runs/B6_SI_2025/format_conversion/addimage/A05035E1.addimg.rds")



# Map the gene Symbols ----
# Create a Mapping function to get gene names where available

convert_and_rename <- function(obj, status_label) {
  ensembl_ids <- rownames(obj)
  gene_map <- select(org.Mm.eg.db, keys = ensembl_ids, 
                     columns = c("SYMBOL"), keytype = "ENSEMBL")
  
  # Map symbols, keep Ensembl if Symbol is missing
  mapped_symbols <- gene_map$SYMBOL[match(ensembl_ids, gene_map$ENSEMBL)]
  mapped_symbols[is.na(mapped_symbols)] <- ensembl_ids[is.na(mapped_symbols)]
  
  # Collapse counts for duplicate symbols
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  row_names_factor <- factor(mapped_symbols)
  map_matrix <- sparse.model.matrix(~ 0 + row_names_factor)
  colnames(map_matrix) <- levels(row_names_factor)
  collapsed_counts <- t(map_matrix) %*% counts
  
  new_obj <- CreateSeuratObject(counts = collapsed_counts, assay = "Spatial")
  
  # Transfer metadata and images
  new_obj@images <- obj@images
  new_obj@meta.data <- obj@meta.data
  new_obj$InfectionStatus <- status_label
  
  return(new_obj)
}


#Convert the gene symbols 
Yp_SI_symbols <- convert_and_rename(Yp_SI_stereo, "Yersinia")
Naive_symbols <- convert_and_rename(Naive_stereo, "Naive")

Yp_SI_symbols <- UpdateSeuratObject(Yp_SI_symbols)
Naive_symbols <- UpdateSeuratObject(Naive_symbols)

#Merge samples
combined_stereo <- merge(Yp_SI_symbols, y = Naive_symbols, 
                         add.cell.ids = c("Yersinia", "Naive"), 
                         project = "STOmics")

table(combined_stereo$InfectionStatus)# Naive 177378 Yersinia 214806


# Calculate QC metrics and filter ----
combined_stereo[["Spatial"]] <- JoinLayers(combined_stereo[["Spatial"]])
combined_stereo$Transcript_Sums <- colSums(GetAssayData(combined_stereo, assay = "Spatial", layer = "counts"))
combined_stereo$nFeature_RNA <- colSums(GetAssayData(combined_stereo, assay = "Spatial", layer = "counts") > 0)

# Temporarily clear images to avoid subset errors with Seurat
unfilt1 <- VlnPlot(combined_stereo, features = "Transcript_Sums", group.by = "InfectionStatus", pt.size = 0)
unfilt2 <- VlnPlot(combined_stereo, features = "nFeature_RNA", group.by = "InfectionStatus", pt.size = 0)
img_backup <- combined_stereo@images
combined_stereo_filt <- subset(
  combined_stereo,
  subset = nFeature_RNA > 40 &  nFeature_RNA < 450 &
    Transcript_Sums > 150 & Transcript_Sums < 2000 &
    percent.mito < 8
)
combined_stereo_filt@images <- img_backup
# Remove filtered cells from images 
remaining_cells <- colnames(combined_stereo_filt)
for (img_name in names(combined_stereo_filt@images)) {
  coords <- GetTissueCoordinates(combined_stereo_filt@images[[img_name]])
  valid_coords <- coords[rownames(coords) %in% remaining_cells, ]
  combined_stereo_filt@images[[img_name]]@coordinates <- valid_coords
}

names(combined_stereo_filt@images) <- c("Yersinia", "Naive")
names(combined_stereo@images) <- c("Yersinia", "Naive")

SpatialFeaturePlot(combined_stereo, features = "Transcript_Sums", images = "Yersinia")
SpatialFeaturePlot(combined_stereo_filt, features = "Transcript_Sums", images = "Yersinia")

# Compare counts by condition
filt1 <- VlnPlot(combined_stereo_filt, features = "Transcript_Sums", group.by = "InfectionStatus", pt.size = 0)
filt2 <- VlnPlot(combined_stereo_filt, features = "nFeature_RNA", group.by = "InfectionStatus", pt.size = 0)

unfilt1 +filt1
unfilt2 +filt2

SpatialFeaturePlot(combined_stereo_filt, features = "Transcript_Sums", images = "Yersinia", pt.size = 10, min.cutoff = 0.001)
SpatialFeaturePlot(combined_stereo_filt, features = "Reg3b", images = "Yersinia", slot = "counts", 
                   pt.size = 2000,  max.cutoff = 10)
# Change the names to actual conditions
names(combined_stereo_filt@images) 
names(combined_stereo_filt@images) <- c("Yersinia", "Naive")
combined_stereo_filt@images$Yersinia@assay <- "Spatial"
combined_stereo_filt@images$Naive@assay <- "Spatial"
names(combined_stereo_filt@images)

# Spatial plot, specify the image
SpatialDimPlot(combined_stereo_filt, group.by = "InfectionStatus", pt.size = 2000)

#Normalization and Integration ----

# Remove filtered cells from images 
remaining_cells <- colnames(combined_stereo_filt)
# Subset the image coordinates for each image 
for (img_name in names(combined_stereo_filt@images)) {
  # Get coordinates table
  coords <- GetTissueCoordinates(combined_stereo_filt@images[[img_name]])
  valid_coords <- coords[rownames(coords) %in% remaining_cells, ]
  combined_stereo_filt@images[[img_name]]@coordinates <- valid_coords
}

validObject(combined_stereo_filt)

# Run SCTransform on the combined object, splitting by sample to handle batch effects
combined_stereo_filt <- SCTransform(combined_stereo_filt, assay = "Spatial", 
                                    vars.to.regress = "InfectionStatus", 
                                    verbose = TRUE)


DefaultAssay(combined_stereo_filt) <- "SCT"
combined_stereo_filt <- FindVariableFeatures(
  combined_stereo_filt,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  mean.function = "FastExpMean",
  dispersion.function = "FastLogVMR",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE)
combined_stereo_filt <- RunPCA(combined_stereo_filt, verbose = FALSE)
combined_stereo_filt <- RunUMAP(combined_stereo_filt, reduction = "pca", dims = 1:25, reduction.name = "Dim25_Umap")
combined_stereo_filt <- FindNeighbors(combined_stereo_filt, reduction = "pca", dims = 1:25)
combined_stereo_filt <- FindClusters(combined_stereo_filt, resolution = 0.5,  cluster.name = "Clusters_0.5")
combined_stereo_filt <- FindClusters(combined_stereo_filt, resolution = 0.75,  cluster.name = "Clusters_0.75")
DimPlot(combined_stereo_filt, reduction = "Dim25_Umap", group.by = "Clusters_0.75") + NoLegend()
FeaturePlot(combined_stereo_filt, "Vil1", reduction = "Dim25_Umap")
SpatialDimPlot(combined_stereo_filt, group.by = "Clusters_0.75", images = "Yersinia", pt.size =2500)

# Save the CellBin based analysis at this stage ----
save.image("/path/to/directory/Projects/STOmics/Intestine/Naive_yersinia.RData")


# Perform a deconvolution and cell type annotation ----
# A separate script, based on the Seurat guide for RCTD deconvolution was written and executed in base R on the objects here. It provided the RCTD_prediction for the object 
# Scripts/RCTD_secondlevel.r

# Load the environment and Examine the Cell type annotations provided 
SpatialDimPlot(combined_stereo_filt, group.by = "RCTD_prediction", images = c("Naive"), pt.size =2000,  combine= TRUE) 
SpatialDimPlot(combined_stereo_filt, group.by = "RCTD_prediction", images = c("Yersinia"), pt.size =2000,  combine= TRUE) 


combined_stereo_filt$CoarseCellPrediction <- combined_stereo_filt$RCTD_prediction
unique(combined_stereo_filt$CoarseCellPrediction)
Myeloid <- c("cDC1s" ,"Macrophages","Monocytes",   "Plasmacytoid DCs" ,  "Neutrophils" ,   "Mast Cells" , "cDC2s"  )
epithelial <- c(  "Goblet Cells" , "Paneth Cells", "Stem Cells" , "Transit Amplifying Cells", "Enterocytes" , "Enteroendocrine Cells", "Tuft Cells"    )
Stromal <- c( "Smooth Muscle Cells"   , "Stromal Cells" , "Fibroblasts"  ,  "Blood Endothelial Cells"  ,   "Mesothelial Cells"  , 
              "Lymphatic Endothelial Cells")
Tcell <- c( "CD8 ab T Cells","Effector CD4 ab T Cells" , "gd T Cells" , "NKT Cells" , 
            "CD4 ab T Cells", "ILC1s" ,  "Cycling ILCs" , "Cycling T Cells" , "ILC3s"  , "Effector CD8 ab T Cells" , "NK Cells"  , "ILC2s"   )
Nrvous <- c( "Enteric Glia" )
Bcell <- c("B Cells" ,  "GC B Cells" , "Plasma Cells"  )

combined_stereo_filt$CoarseCellPrediction[combined_stereo_filt$CoarseCellPrediction %in%Myeloid] <- "Myeloid"
combined_stereo_filt$CoarseCellPrediction[combined_stereo_filt$CoarseCellPrediction %in%epithelial] <- "Epithelial"
combined_stereo_filt$CoarseCellPrediction[combined_stereo_filt$CoarseCellPrediction %in%Stromal] <- "Stromal"
combined_stereo_filt$CoarseCellPrediction[combined_stereo_filt$CoarseCellPrediction %in%Tcell] <- "T Cells and ILCs"
combined_stereo_filt$CoarseCellPrediction[combined_stereo_filt$CoarseCellPrediction %in%Nrvous] <- "ENS"
combined_stereo_filt$CoarseCellPrediction[combined_stereo_filt$CoarseCellPrediction %in%Bcell] <- "B Cells and Plasma Cells"
combined_stereo_filt$CoarseCellPrediction[combined_stereo_filt$CoarseCellPrediction %in%Nrvous] <- "ENS"


library(ggplot2)
SpatialDimPlot(combined_stereo_filt, group.by = "CoarseCellPrediction", images = c("Yersinia"), pt.size =2000,  combine= TRUE) + scale_fill_manual(values = c8_b)

ggsave( "CoarseCellType_Yersinia_SpatialDim.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 10,  height = 8,  units = c("in"),  dpi = 600)

SpatialDimPlot(combined_stereo_filt, group.by = "CoarseCellPrediction", images = c("Naive"), pt.size =2000,  combine= TRUE) + scale_fill_manual(values = c8_b)

ggsave( "CoarseCellType_Naive_SpatialDim.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 10,  height = 8,  units = c("in"),  dpi = 600)



# Add the Microbial data to the objects
library(data.table)
library(Matrix)
Yersinia_microbes <- data.table::fread("/path/to/directory/SAW_pipeline/runs/Yp_SI_2025_custom/Yp_SI_2025_custom/outs/feature_expression/microorganism/Y01679L2_microorganism.genus.gem")


# Create a plot overlaying the microbial data 
# filter to Yersinia genus reads

range_x <- GetTissueCoordinates(combined_stereo_filt, image = "Yersinia")

cell_map <- fread("/path/to/directory/SAW_pipeline/runs/Yp_SI_2025_custom/format_conversion/Y01679L2.adjusted.cellbin.gem.gz")

head(cell_map)
head(colnames(combined_stereo_filt))

yersinia_only <- Yersinia_microbes[geneID == "g__Yersinia"]
# Match  Yersinia coordinates to the CellIDs 
yersinia_with_ids <- merge(yersinia_only, cell_map, by = c("x", "y"))

# Format the CellID to match 
yersinia_with_ids[, Seurat_Barcode := paste0("Yersinia_", CellID)]

# Filter for only the cells in Seurat object
target_barcodes <- intersect(yersinia_with_ids$Seurat_Barcode, colnames(combined_stereo_filt))

all_coords <- GetTissueCoordinates(combined_stereo_filt)
star_coords <- all_coords[target_barcodes, ]


p_base <- SpatialDimPlot(combined_stereo_filt, 
                         group.by = "CoarseCellPrediction", 
                         images = "Yersinia", 
                         pt.size = 2000) + 
  scale_fill_manual(values = c8_b)

p_final <- p_base + 
  geom_point(data = as.data.frame(star_coords), 
             aes(x = imagecol, y = imagerow), 
             color = "yellow", 
             shape = 8, 
             size = 3, 
             inherit.aes = FALSE)

print(p_final) # Notice the Yersinia reads are rotated compared to the image - there is an error with the way the coordinates are input in seurat vs loaded directly for the microbial reads. The below code checks for the correct orientation (which is known because it is available in the STEREOMAP software for the outputs of the data).

library(patchwork)

# The Yersinia reads are rotated relative to the image and this needs to be corrected
base_coords <- as.data.frame(star_coords)
x_range <- range(all_coords$imagecol)
y_range <- range(all_coords$imagerow)

# generate the plot for each test
plot_test <- function(df, title) {
  SpatialDimPlot(combined_stereo_filt, group.by = "CoarseCellPrediction", images = "Yersinia", pt.size = 2000) +
    geom_point(data = df, aes(x = imagecol, y = imagerow), 
               color = "yellow", shape = 8, size = 2, inherit.aes = FALSE) +
    scale_fill_manual(values = c8_b) +
    labs(title = title) +
    theme(legend.position = "none") 
}

# Current coordinates
df_0 <- base_coords

# Fix type a
df_A <- base_coords
df_A$imagerow <- y_range[2] - (df_A$imagerow - y_range[1])

# fix bX-Y Swap 90 Degree Rotation/Transposition
df_B <- data.frame(imagecol = base_coords$imagerow, imagerow = base_coords$imagecol)

# fix c
df_C <- base_coords
df_C$imagecol <- x_range[2] - (df_C$imagecol - x_range[1])
df_C$imagerow <- y_range[2] - (df_C$imagerow - y_range[1])

# plot all fixes
p0 <- plot_test(df_0, "0: Original")
pA <- plot_test(df_A, "A: Y-Mirror (Flip)")
pB <- plot_test(df_B, "B: X-Y Swap (Transpose)")
pC <- plot_test(df_C, "C: 180 Rotation")
diagnostic_grid <- (p0 + pA) / (pB + pC)
print(diagnostic_grid)


# Final flipped image
y_range <- range(all_coords$imagerow)

# reflect
star_coords_final <- as.data.frame(star_coords)
star_coords_final$imagerow <- y_range[2] - (star_coords_final$imagerow - y_range[1])

p_final <- SpatialDimPlot(combined_stereo_filt, 
                          group.by = "CoarseCellPrediction", 
                          images = "Yersinia", 
                          pt.size = 2000) + 
  scale_fill_manual(values = c8_b) +
  geom_point(data = star_coords_final, 
             aes(x = imagecol, y = imagerow), 
             color = "yellow", 
             shape = 8, 
             size = 0.5, 
             inherit.aes = FALSE) +
  labs(title = "Yersinia Localization (Y-Flip Corrected)")

print(p_final)
ggsave( "CoarseCellType_Yersinia_SpatialDim_With_Yersinia_transcripts.pdf",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 15 ,  height = 13,  units = c("in"),  dpi = 600)



SpatialFeaturePlot(combined_stereo_filt, "Nos2", images = "Yersinia", pt.size = 4000, max.cutoff = 0.5)



# Plot the Nos2+ cells as Stars 

nos2_counts <- GetAssayData(combined_stereo_filt, assay = "SCT", slot = "counts")["Nos2", ]
yersinia_image_cells <- colnames(combined_stereo_filt)[combined_stereo_filt$InfectionStatus == "Yersinia" & combined_stereo_filt$CoarseCellPrediction == "Epithelial"]
nos2_pos_cells <- intersect(yersinia_image_cells, names(nos2_counts[nos2_counts > 0]))

coordinates_nos2 <-GetTissueCoordinates(combined_stereo_filt)
coordinates_nos2 <- coordinates_nos2[nos2_pos_cells, ] # only Epithelial Cell Expression is relevant 



p_final_nos2 <- SpatialDimPlot(combined_stereo_filt, 
                          group.by = "CoarseCellPrediction", 
                          images = "Yersinia", 
                          pt.size = 2000) + 
  scale_fill_manual(values = c8_b) +
  geom_point(data = star_coords_final, 
             aes(x = imagecol, y = imagerow), 
             color = "yellow", 
             shape = 8, 
             size = 0.5, 
             inherit.aes = FALSE) +
  geom_point(data = coordinates_nos2, 
             aes(x = imagecol, y = imagerow), 
             color = "red", 
             shape = 8, 
             size = 0.5, 
             inherit.aes = FALSE) +
  labs(title = "Yersinia Localization and Nos2 Expressio")

print(p_final_nos2)



# Correct the coordinate positions for plotting 
coords_all <- GetTissueCoordinates(combined_stereo_filt, image = "Yersinia")
max_row_tissue <- max(coords_all$imagerow)
min_row_tissue <- min(coords_all$imagerow)


p_final_nos2 <- SpatialDimPlot(combined_stereo_filt, 
                               group.by = "CoarseCellPrediction", 
                               images = "Yersinia", 
                               pt.size = 2000) + 
  scale_fill_manual(values = c8_b) +
  geom_point(data = star_coords_final, 
             aes(x = imagecol, y = imagerow), 
             color = "yellow", 
             shape = 8, 
             size = 0.5, 
             inherit.aes = FALSE) +
  geom_point(data = coordinates_nos2, 
             aes(x = imagecol, 
                 y = (max_row_tissue + min_row_tissue) - imagerow), 
             color = "red", 
             shape = 8, 
             size = 0.8, 
             inherit.aes = FALSE) +
  labs(title = "Yersinia Localization (Nos2 Vertically Reflected & Re-Anchored)")

print(p_final_nos2)

ggsave( "CoarseCellType_Yersinia_SpatialFeature_Yps_enrichment_Yersinia_transcripts_Nos2.pdf",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 15 ,  height = 13,  units = c("in"),  dpi = 600)


