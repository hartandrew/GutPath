# Description: This script will take as input, the Cell Ranger pipeline output for each sample. I will load the samples in R using Seurat and it will filter samples and perform normalization of the ADT reads using dsb.  Metadata about each experiment is added and the objects are merged for further analysis
#Figure Panels produced : Figure S1C Model of Analysis Approach
# Load the Libraries 
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

# Define needed variables 
# Define Global Variables and isotype names 
isotype.controls <- c("Isotype_RTK4530", "Isotype_HTK888", "Isotype_RTK2071", 
                      "Isotype_RTK2758", "Isotype_MOPC.173", "Isotype_MOPC.21", "Isotype_MPC.11")

images_dir <- "/path/to/analysis/directory/Images"



Sample_directory <- data.frame(directory = c("/path/CellRanger/CHMI_183/", "/path/CellRanger/CHMI_183/", "/path/CellRanger/CHMI_183/",
                                            "/path/CellRanger/CHMI_183/", "/path/CellRanger/CHMI_183/", "/path/CellRanger/CHMI_183/", 
                                            "/path/CellRanger/CHMI_190/", "/path/CellRanger/CHMI_190/", "/path/CellRanger/CHMI_190/",
                                            "/path/CellRanger/CHMI_190/", "/path/CellRanger/CHMI_190/",
                                            "/path/CellRanger/CHMI_192/", "/path/CellRanger/CHMI_192/", "/path/CellRanger/CHMI_192/",
                                            "/path/CellRanger/CHMI_197/", "/path/CellRanger/CHMI_197/", "/path/CellRanger/CHMI_197/",
                                            "/path/CellRanger/CHMI_204/", 
                                            "/path/CellRanger/CHMI_210/", "/path/CellRanger/CHMI_210/", "/path/CellRanger/CHMI_210/", 
                                            "/path/CellRanger/CHMI_217/", "/path/CellRanger/CHMI_217/", "/path/CellRanger/CHMI_217/", "/path/CellRanger/CHMI_217/", 
                                           "/path/CellRanger/CHMI_229/", "/path/CellRanger/CHMI_229/", "/path/CellRanger/CHMI_229/" ,"/path/CellRanger/CHMI_229/" ,
                                           "/path/CellRanger/CHMI_233/", "/path/CellRanger/CHMI_233/", "/path/CellRanger/CHMI_233/",
                                           "/path/CellRanger/CHMI_239/", "/path/CellRanger/CHMI_239/", "/path/CellRanger/CHMI_239/", 
                                           "/path/CellRanger/CHMI_245/", "/path/CellRanger/CHMI_245/", "/path/CellRanger/CHMI_245/",
                                           "/path/CellRanger/CHMI_250/", "/path/CellRanger/CHMI_250/", "/path/CellRanger/CHMI_250/", 
                                           "/path/CellRanger/CHMI_260/", "/path/CellRanger/CHMI_260/", "/path/CellRanger/CHMI_260/", 
                                           "/path/CellRanger/CHMI_262/", "/path/CellRanger/CHMI_262/", "/path/CellRanger/CHMI_262/", 
                                           "/path/CellRanger/CHMI_266/", "/path/CellRanger/CHMI_266/", "/path/CellRanger/CHMI_266/", 
                                           "/path/CellRanger/CHMI_266/", "/path/CellRanger/CHMI_266/",
                                           "/path/CellRanger/CHMI_267/", "/path/CellRanger/CHMI_267/", "/path/CellRanger/CHMI_267/",
                                           "/path/CellRanger/CHMI_267/", "/path/CellRanger/CHMI_267/"),

                                FolderID = c("MIST2482_MLN", "MIST2482_LP", "MIST2482_IEL", "MIST2489_MLN", "MIST2489_LP", "MIST2489_IEL", 
                                            "MIST24109_MLN", "MIST24109_LP", "MIST24109_IEL", "MIST2493_MLN", "MIST2493_LP", 
                                            "MIST24134_LP", "MIST24134_MLN", "MIST24134_IEL", 
                                            "MIST24142_MLN", "MIST24142_IEL", "MIST24142_LP",
                                            "MIST24180_MLN", 
                                            "MIST24218_MLN", "MIST24218_IEL", "MIST24218_LP", 
                                            "MIST24240_MLN", "MIST24240_IEL", "MIST24240_LP", "MIST24235_MLN", 
                                            "MIST24235_LP_reseq", "MIST24268_MLN", "MIST24268_IEL", "MIST24268_LP", 
                                            "MIST24282_MLN", "MIST24282_IEL", "MIST24282_LP", 
                                            "MIST24311_MLN", "MIST24311_IEL", "MIST24311_LP",
                                            "MIST24353_MLN", "MIST24353_IEL", "MIST24353_LP", 
                                            "MIST25029_MLN", "MIST25029_IEL", "MIST25029_LP", 
                                            "MIST25100_MLN", "MIST25100_IEL", "MIST25100_LP", 
                                            "MIST25114_MLN", "MIST25114_IEL", "MIST25114_LP", 
                                            "MIST25132_MLN" , "MIST25132_IEL", "MIST25132_LP", "MIST25132_IEL_D", "MIST25132_LP_D",
                                            "MIST25134_MLN", "MIST25134_IEL", "MIST25134_LP", "MIST25134_IEL_D", "MIST25134_LP_D"
                                            ), 
                                SampleID = c("CHMI2482_MLN", "CHMI2482_LP", "CHMI2482_IEL", "CHMI2489_MLN", "CHMI2489_LP", "CHMI2489_IEL", 
                                            "CHMI24109_MLN", "CHMI24109_LP", "CHMI24109_IEL", "CHMI2493_MLN", "CHMI2493_LP", 
                                            "CHMI24134_LP", "CHMI24134_MLN", "CHMI24134_IEL", 
                                            "CHMI24142_MLN", "CHMI24142_IEL", "CHMI24142_LP",
                                            "CHMI24180_MLN", 
                                            "CHMI24218_MLN", "CHMI24218_IEL", "CHMI24218_LP", 
                                            "CHMI24240_MLN", "CHMI24240_IEL", "CHMI24240_LP", "CHMI24235_MLN", 
                                            "CHMI24235_LP", "CHMI24268_MLN", "CHMI24268_IEL", "CHMI24268_LP", 
                                            "CHMI24282_MLN", "CHMI24282_IEL", "CHMI24282_LP", 
                                            "CHMI24311_MLN", "CHMI24311_IEL", "CHMI24311_LP",
                                            "CHMI24353_MLN", "CHMI24353_IEL", "CHMI24353_LP", 
                                            "CHMI25029_MLN", "CHMI25029_IEL", "CHMI25029_LP", 
                                            "CHMI25100_MLN", "CHMI25100_IEL", "CHMI25100_LP", 
                                            "CHMI25114_MLN", "CHMI25114_IEL", "CHMI25114_LP", 
                                            "CHMI25132_MLN" , "CHMI25132_IEL", "CHMI25132_LP", "CHMI25132_IEL_D", "CHMI25132_LP_D",
                                            "CHMI25134_MLN", "CHMI25134_IEL", "CHMI25134_LP", "CHMI25134_IEL_D", "CHMI25134_LP_D"), 
                                Mouse = c("CHMI2482", "CHMI2482", "CHMI2482", "CHMI2489", "CHMI2489", "CHMI2489", "CHMI24109", "CHMI24109", "CHMI24109", "CHMI2493", "CHMI2493", "CHMI24134", "CHMI24134", "CHMI24134", "CHMI24142", "CHMI24142", "CHMI24142", "CHMI24180", "CHMI24218", "CHMI24218", "CHMI24218", "CHMI24240", "CHMI24240", "CHMI24240", "CHMI24235", "CHMI24235", "CHMI24268", "CHMI24268", "CHMI24268", "CHMI24282", "CHMI24282", "CHMI24282", "CHMI24311", "CHMI24311", "CHMI24311", "CHMI24353", "CHMI24353", "CHMI24353", "CHMI25029", "CHMI25029", "CHMI25029", "CHMI25100", "CHMI25100", "CHMI25100", "CHMI25114", "CHMI25114", "CHMI25114", "CHMI25132", "CHMI25132", "CHMI25132", "CHMI25132", "CHMI25132", "CHMI25134", "CHMI25134", "CHMI25134", "CHMI25134", "CHMI25134"), 
                                Tissue = c("Mesenteric Lymph Node", "Lamina Propria", "Intestinal Epithelial Cells", "Mesenteric Lymph Node", "Lamina Propria", "Intestinal Epithelial Cells", "Mesenteric Lymph Node", "Lamina Propria", "Intestinal Epithelial Cells", "Mesenteric Lymph Node", "Lamina Propria", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Intestinal Epithelial Cells", "Lamina Propria", "Mesenteric Lymph Node", "Intestinal Epithelial Cells", "Lamina Propria", "Intestinal Epithelial Cells", "Lamina Propria"), 
                                InfectionStatus = c("Naive", "Naive", "Naive", "Cryptosporidium", "Cryptosporidium", "Cryptosporidium", "Naive", "Naive", "Naive", "Cryptosporidium", "Cryptosporidium", "Yersinia", "Yersinia", "Yersinia", "Yersinia", "Yersinia", "Yersinia", "H.poly", "Candida", "Candida", "Candida", "Cryptosporidium", "Cryptosporidium", "Cryptosporidium", "Candida", "Candida", "MNV", "MNV", "MNV", "MNV", "MNV", "MNV", "Yersinia", "Yersinia", "Yersinia", "SFB_CU", "SFB_CU", "SFB_CU", "SFB_CU", "SFB_CU", "SFB_CU", "SFB_YA", "SFB_YA", "SFB_YA", "SFB_YA", "SFB_YA", "SFB_YA", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus", "Nippostrongylus"))

# Assign the filtering - 0.05 for MLN and 0.2 for other tissues
Sample_directory$MtCutoff <- ifelse(
  Sample_directory$Tissue == "Mesenteric Lymph Node", 
  0.05, 
  0.20
)


# Create a function that will load the data, filter, and prep for QC evaluation
process_sample_mist_meta <- function(sample_id,        # Name of the variable 
                                     folder_id,        # Name of the folder 
                                     mouse_id,         
                                     tissue_type,      
                                     infection_status, 
                                     mt_cutoff,        # 0.05 or 0.20
                                     project_dir , 
                                     isotype_controls,
                                     image_output_dir = NULL) {
  
  message(paste("Processing:", sample_id, "| Mouse:", mouse_id, "| Tissue:", tissue_type))
  
  
  raw_path <- file.path(project_dir, folder_id, "outs/raw_feature_bc_matrix/")
  filt_path <- file.path(project_dir, folder_id, "outs/filtered_feature_bc_matrix/")
  
  raw.data <- Read10X(data.dir = raw_path)
  filtered.data <- Read10X(data.dir = filt_path)
  stained_cells <- colnames(filtered.data$`Gene Expression`)
  
  prot <- raw.data$`Antibody Capture`
  rna <- raw.data$`Gene Expression`
  
  # Calculate mitochondrial and ribosomal reads
  mtgene <- grep(pattern = "^mt-", rownames(rna), value = TRUE)
  rb.genes <- c(grep(pattern = "Rps", rownames(rna), value = TRUE), grep(pattern = "Rpl", rownames(rna), value = TRUE))
  
  md <- data.frame(
    rna.size = log10(Matrix::colSums(rna)), 
    prot.size = log10(Matrix::colSums(prot)), 
    n.gene = Matrix::colSums(rna > 0), 
    mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna), 
    rb.prop = Matrix::colSums(rna[rb.genes, ]) / Matrix::colSums(rna)
  )
  
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  cellmd <- md[md$drop.class == 'cell', ]
  # Calculate 2.5 absolute median deviations and filter
  rna.mult <- (2.5 * mad(cellmd$rna.size))
  prot.mult <- (2.5 * mad(cellmd$prot.size))
  rna.lower <- median(cellmd$rna.size) - rna.mult
  rna.upper <- median(cellmd$rna.size) + rna.mult
  prot.lower <- median(cellmd$prot.size) - prot.mult
  prot.upper <- median(cellmd$prot.size) + prot.mult
  
  qc_cells <- rownames(cellmd[cellmd$prot.size > prot.lower & cellmd$n.gene > 750 & 
                              cellmd$prot.size < prot.upper & cellmd$rna.size > rna.lower & 
                              cellmd$rna.size < rna.upper & cellmd$mt.prop < mt_cutoff, ])
  
  # Normalize ADT data
  background_drops <- rownames(md[md$prot.size > 1.5 & md$prot.size < 3 & md$rna.size < 2.5, ])
  background.adt.mtx <- as.matrix(prot[, background_drops])
  cell.adt.raw <- as.matrix(prot[, qc_cells])
  cell.rna.raw <- rna[, qc_cells]
  cellmd_filt <- cellmd[qc_cells, ]
  
  cells.dsb.norm <- DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.raw, empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, use.isotype.control = TRUE, isotype.control.name.vec = isotype_controls
  )
  
  # Create the Seurat Object and add the metadata
  seurat_obj <- Seurat::CreateSeuratObject(counts = cell.rna.raw, meta.data = cellmd_filt, assay = "RNA", min.cells = 20)
  seurat_obj[["ADT"]] <- Seurat::CreateAssayObject(data = cells.dsb.norm)
  seurat_obj$orig.ident <- sample_id # Sets the Identity
  
  # Add the Custom Metadata
  seurat_obj$Mouse <- mouse_id
  seurat_obj$Tissue <- tissue_type
  seurat_obj$InfectionStatus <- infection_status
  
  # Rename Cells
  seurat_obj@meta.data["CellNames"] <- paste0(sample_id, "_", colnames(seurat_obj))

  # Plotting
  if (!is.null(image_output_dir)) {
    plot_qc_metrics(md, seurat_obj, sample_id, image_output_dir)
  }
  
  return(seurat_obj)
}

# Function for Plotting
plot_qc_metrics <- function(full_md, seurat_obj, sample_id, image_dir) {
  
  # Drop class plot (Protein vs nGene)
  p1 <- ggplot(full_md, aes(x = log10(n.gene), y = prot.size)) +
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~drop.class) +
    ggtitle(paste(sample_id, "Protein vs nGene"))
  
  # Drop class plot (RNA vs nGene)
  p2 <- ggplot(full_md, aes(x = log10(n.gene), y = rna.size)) +
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~drop.class) +
    ggtitle(paste(sample_id, "RNA vs nGene"))
  
  # Violin Plots
  feats <- c("nFeature_RNA", "nCount_RNA", "mt.prop")
  p3 <- VlnPlot(seurat_obj, group.by = "orig.ident", features = feats, pt.size = 0.1, combine = T, ncol = 3, cols = "#FB20FE") + NoLegend()
  
  # QC Scatter
  p4 <- ggplot(as.data.frame(seurat_obj@meta.data), aes(x = nFeature_RNA, y = nCount_RNA)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + theme_bw() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    ggtitle(paste(sample_id, "Post-Filter QC"))

  # Save Violin 
  ggsave(paste0(sample_id, "_RNA_QC_violins.png"), plot = p3, device = NULL,
         path = images_dir, scale = 1, width = 7, height = 5,
         units = c("in"), dpi = 600, limitsize = TRUE, bg = NULL)
  
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
}


for(i in 1:nrow(Sample_directory)) {
  
  # Extract variables for readability
  row <- Sample_directory[i, ]
  
  # Run the function
  processed_obj <- process_sample_mist_meta(
    sample_id = row$SampleID,
    folder_id = row$FolderID,
    mouse_id = row$Mouse,
    tissue_type = row$Tissue,
    infection_status = row$InfectionStatus,
    mt_cutoff = row$MtCutoff,
    project_dir = row$directory,
    isotype_controls = isotype.controls,
    image_output_dir = images_dir
  )
  
  assign(row$SampleID, processed_obj)
}

#Merge the data sets - based on combined
MLN_data <- merge(CHMI2482_MLN, y = c(CHMI2489_MLN, CHMI24109_MLN, CHMI2493_MLN,  CHMI24134_MLN,  CHMI24142_MLN, CHMI24180_MLN, CHMI24218_MLN, CHMI24240_MLN, CHMI24235_MLN, CHMI24268_MLN, CHMI24282_MLN, CHMI24311_MLN , CHMI24353_MLN, CHMI25029_MLN, CHMI25100_MLN, CHMI25114_MLN, CHMI25132_MLN, CHMI25134_MLN), add.cell.ids = c("CHMI2482_MLN", "CHMI2489_MLN", "CHMI24109_MLN", "CHMI2493_MLN", "CHMI24134_MLN","CHMI24142_MLN",  "CHMI24180_MLN", "CHMI24218_MLN", "CHMI24240_MLN", "CHMI24235_MLN", "CHMI24268_MLN", "CHMI24282_MLN", "CHMI24311_MLN", "CHMI24353_MLN", "CHMI25029_MLN", "CHMI25100_MLN", "CHMI25114_MLN", "CHMI25132_MLN", "CHMI25134_MLN"), project = "Combined", merge.data = TRUE)

Ileum_data <- merge(CHMI2482_LP, y = c(  CHMI2482_IEL,  CHMI2489_IEL, CHMI2489_LP, CHMI24109_IEL, CHMI24109_LP,  CHMI2493_LP,  CHMI24134_LP, CHMI24134_IEL,  CHMI24142_IEL, CHMI24142_LP,  CHMI24218_LP, CHMI24218_IEL, CHMI24240_LP, CHMI24240_IEL, CHMI24235_LP,CHMI24268_LP, CHMI24268_IEL, CHMI24282_LP, CHMI24282_IEL, CHMI24311_LP, CHMI24311_IEL, CHMI24353_LP, CHMI24353_IEL, CHMI25029_LP, CHMI25029_IEL,  CHMI25100_LP, CHMI25100_IEL, CHMI25114_LP, CHMI25114_IEL, CHMI25132_LP, CHMI25132_IEL, CHMI25132_LP_D, CHMI25132_IEL_D,  CHMI25134_LP, CHMI25134_IEL, CHMI25134_LP_D, CHMI25134_IEL_D), add.cell.ids = c( "CHMI2482_LP", "CHMI2482_IEL",  "CHMI2489_IEL", "CHMI2489_LP", "CHMI24109_IEL", "CHMI24109_LP",  "CHMI2493_LP",  "CHMI24134_LP", "CHMI24134_IEL",  "CHMI24142_IEL", "CHMI24142_LP",  "CHMI24218_LP", "CHMI24218_IEL", "CHMI24240_LP", "CHMI24240_IEL", "CHMI24235_LP","CHMI24268_LP", "CHMI24268_IEL", "CHMI24282_LP", "CHMI24282_IEL", "CHMI24311_LP", "CHMI24311_IEL", "CHMI24353_LP", "CHMI24353_IEL", "CHMI25029_LP", "CHMI25029_IEL", "CHMI25100_LP", "CHMI25100_IEL", "CHMI25114_LP", "CHMI25114_IEL", "CHMI25132_LP", "CHMI25132_IEL", "CHMI25132_LP_D", "CHMI25132_IEL_D",  "CHMI25134_LP", "CHMI25134_IEL", "CHMI25134_LP_D", "CHMI25134_IEL_D"), project = "Combined", merge.data = TRUE)


#The R environment for the merged objects prior to downstream analysis was saved as Ileum_MLN_combined_prior_to_anlysis.RData in Seurat_Files for this project
save.image(file = "Ileum_MLN_combined_prior_to_anlysis.RData")