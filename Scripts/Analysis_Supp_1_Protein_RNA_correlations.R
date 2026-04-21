# Description. This script will look at whether the Protein data and RNA data are aligned. Because of sparcity at the single cell level, each cell type will be pseudobulked for each independent mouse for both the normalized RNA and normalized ADT expression. A simple Pearson correlation is performed to determine the relationship between the protein and the RNA. For Proteins composed on more than one gene transcript, a single transript was chosen. This accounts for less than 10% of the proteins and I do not believe influences the results. Overall the results suggest fairly high correlation for proteins with more abundant proteins having beter correlations. The premade cocktail includes many proteins which might not be heavily expressed in the cells of the gut. These proteins show poor correlation. 

# Figures S3A

# Notably these correlations would likely improve if we removed the cell types for each protein that did not express the protein. (noise)

# Load Libraries ----
library(patchwork)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(scCustomize)
library(forcats)
library(ggbreak)
library(ggforce)
library(future)
library(SeuratDisk)
library(RColorBrewer)
library(gt)
library(SeuratWrappers)
library(qs2)
library(Matrix)
library(tidyverse)
library(tidyverse)


# Load the Samples 
Ileum_filt <- qs_read("/path/to/directory/MIST/MIST_Analysis/Seurat_Files/Analysis2.5_Ileum.qs2")

# Directories 
seurat <- "/path/to/directory/MIST/MIST_Analysis/Seurat_Files"
images <- "/path/to/directory/MIST/MIST_Analysis/Images"


# Establish Protein to Gene linkages 
rownames(Ileum_filt[["ADT"]])


protein_names <- c(
  "CD103", "CD106", "CD107a", "CD115", "CD117", "CD11a", "CD11b", "CD11c", "CD122", "CD127", 
  "CD138", "CD150", "CD163", "CD169", "CD170", "CD172a", "CD182", "CD183", "CD185", "CD19", 
  "CD199", "CD1d", "CD20", "CD200", "CD210", "CD22", "CD223", "CD24", "CD25", "CD27", 
  "CD274", "CD278", "CD279", "CD3", "CD300LG", "CD301b", "CD304", "CD31", "CD317", "CD34-SA376A4", 
  "CD357", "CD36", "CD366", "CD371", "CD38", "CD4", "CD40", "CD41", "CD44", "CD45", 
  "CD45.2", "CD45R-B220", "CD48", "CD49a", "CD49b", "CD49d", "CD49f", "CD5", "CD54", "CD55", 
  "CD61", "CD62L", "CD63", "CD68", "CD69", "CD73", "CD79b", "CD86", "CD88", "CD8a", 
  "CD8b", "CD9", "CD90.2", "CD98", "CX3CR1", "F4-80", "FceRIa", "I.A-I.E", "IgD", "integrin.b7", 
  "KLRG1", "Ly.49A", "Ly.6A-E", "Ly.6C", "Ly.6G", "Ly108", "Ly49D", "Ly49H", "MERTK", "NK.1.1", 
  "TCR.Bchain", "TCR.Vb8.1-8.2", "TCR.Vg2", "Tim.4", "XCR1"
)

# Define  corresponding gene symbols
gene_symbols <- c(
  "Itgae", "Vcam1", "Lamp1", "Csf1r", "Kit", "Itgal", "Itgam", "Itgax", "Il2rb", "Il7r", 
  "Sdc1", "Slamf1", "Cd163", "Siglec1", "Siglecf", "Sirpa", "Cxcr2", "Cxcr3", "Cxcr5", "Cd19", 
  "Ccr9", "Cd1d1", "Ms4a1", "Cd200", "Il10ra", "Cd22", "Lag3", "Cd24a", "Il2ra", "Cd27", 
  "Cd274", "Icos", "Pdcd1", "Cd3e", "Cd300lg", "Mgl2", "Nrp1", "Pecam1", "Bst2", "Cd34", 
  "Tnfrsf18", "Cd36", "Havcr2", "Clec12a", "Cd38", "Cd4", "Cd40", "Itga2b", "Cd44", "Ptprc", 
  "Ptprc", "Ptprc", "Cd48", "Itga1", "Itga2", "Itga4", "Itga6", "Cd5", "Icam1", "Cd55", 
  "Itgb3", "Sell", "Cd63", "Cd68", "Cd69", "Nt5e", "Cd79b", "Cd86", "C5ar1", "Cd8a", 
  "Cd8b1", "Cd9", "Thy1", "Slc3a2", "Cx3cr1", "Adgre1", "Fcer1a", "H2-Ab1", "Ighd", "Itgb7", 
  "Klrg1", "Klra1", "Ly6a", "Ly6c1", "Ly6g", "Slamf6", "Klra4", "Klra8", "Mertk", "Klrb1c", 
  "Trbc1", "Trbv13-1", "Trgv2", "Timd4", "Xcr1"
)

setdiff(gene_symbols, rownames(Ileum_filt[["RNA"]]) )

adt_to_gene <- data.frame(Protein = protein_names, Gene = gene_symbols)


Ileum_filt$pb_cite <- paste(Ileum_filt$Mouse, Ileum_filt$secondlevel, sep = "_")



# Extract normalized layers and aggregate by group (sum and average expression per sample_celltype)
Ileum_filt[["ADT"]] <- JoinLayers(Ileum_filt[["ADT"]])
adt_data <- GetAssayData(Ileum_filt, assay = "ADT", layer = "data")
rna_data <- GetAssayData(Ileum_filt, assay = "RNA", layer = "data")


group_factor <- factor(Ileum_filt$pb_cite)
selection_mat <- sparse.model.matrix(~ 0 + group_factor)


adt_sum <- adt_data %*% selection_mat
rna_sum <- rna_data %*% selection_mat

cell_counts <- as.numeric(table(group_factor))
diag_inv <- Diagonal(x = 1/cell_counts)

adt_pb_avg <- as.matrix(adt_sum %*% diag_inv)
rna_pb_avg <- as.matrix(rna_sum %*% diag_inv)

# Add group names 
colnames(adt_pb_avg) <- levels(group_factor)
colnames(rna_pb_avg) <- levels(group_factor)



# For each pair of the Protein- RNA perform pearson correlation across all pseudobulked levels
coherence_results <- adt_to_gene %>%
  mutate(Coherence = map2_dbl(Protein, Gene, function(p, g) {
    # Check if both exist in our matrices
    if (p %in% rownames(adt_pb_avg) && g %in% rownames(rna_pb_avg)) {
      # Calculate Pearson correlation 
      return(cor(adt_pb_avg[p, ], rna_pb_avg[g, ], method = "pearson"))
    } else {
      return(NA)
    }
  })) %>%
  filter(!is.na(Coherence)) %>%
  arrange(desc(Coherence))

# Save correlations
write.csv(coherence_results, "/path/to/directory/MIST/MIST_Analysis/Coherence_Results.csv", row.names = FALSE)

# Create categories for results
coherence_results <- coherence_results %>%
  mutate(Correlation_Level = case_when(
    Coherence > 0.5 ~ "High (>0.5)",
    Coherence >= 0.3 & Coherence <= 0.5 ~ "Moderate (0.3-0.5)",
    Coherence < 0.3 ~ "Low (<0.3)"
  )) %>%
  # legend order
  mutate(Correlation_Level = factor(Correlation_Level, 
                                    levels = c("High (>0.5)", "Moderate (0.3-0.5)", "Low (<0.3)")))


pdf(paste0(images, "/CITEseq_Coherence_Barplot_3Colors.pdf"), width = 8, height = 15)

ggplot(coherence_results, aes(x = reorder(Protein, Coherence), y = Coherence)) +
  geom_col(aes(fill = Correlation_Level)) +
  scale_fill_manual(values = c(
    "High (>0.5)" = "steelblue", 
    "Moderate (0.3-0.5)" = "orange", 
    "Low (<0.3)" = "firebrick"
  )) +
  coord_flip() +
  labs(
    title = "RNA-Protein Coherence (Pearson Correlation)",
    subtitle = "Categorized by Correlation Strength",
    x = "Protein/Gene Pair",
    y = "Pearson Correlation (r)",
    fill = "Coherence Level"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  )

dev.off()
coherence_results$Correlation_Level <- factor(coherence_results$Correlation_Level, levels = c( "High (>0.5)", "Moderate (0.3-0.5)",  "Low (<0.3)"))
pdf(paste0(images, "/CITEseq_Coherence_3Column_Facets.pdf"), width = 12, height = 8)




coherence_results$Correlation_Level <- factor(
  coherence_results$Correlation_Level, 
  levels = c("High (>0.5)", "Moderate (0.3-0.5)", "Low (<0.3)")
)



data_high <- coherence_results %>% filter(Correlation_Level == "High (>0.5)")
data_mod  <- coherence_results %>% filter(Correlation_Level == "Moderate (0.3-0.5)")
data_low  <- coherence_results %>% filter(Correlation_Level == "Low (<0.3)")

create_fixed_bar <- function(df, color, title) {
  ggplot(df, aes(x = reorder(Protein, Coherence), y = Coherence)) +
    geom_col(fill = color, width = 0.7) + # Width controls the bar 'fatness'
    coord_flip() +
    theme_minimal() +
    labs(title = title, x = NULL, y = "Pearson r") +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}


p1 <- create_fixed_bar(data_high, "steelblue", "High (>0.5)")
p2 <- create_fixed_bar(data_mod, "orange", "Moderate (0.3-0.5)")
p3 <- create_fixed_bar(data_low, "firebrick", "Low (<0.3)")

n1 <- nrow(data_high)
n2 <- nrow(data_mod)
n3 <- nrow(data_low)
max_n <- max(n1, n2, n3)

final_plot <- (
  (p1 / plot_spacer()) + plot_layout(heights = c(n1, max_n - n1)) | 
    (p2 / plot_spacer()) + plot_layout(heights = c(n2, max_n - n2)) | 
    (p3 / plot_spacer()) + plot_layout(heights = c(n3, max_n - n3))
) + plot_annotation(title = "Uniform Bar Thickness Coherence Plot")

pdf(paste0(images, "/CITEseq_Coherence_Faceted_columns_intestine_FigureS3.pdf"), width = 12, height = 6)
print(final_plot)
dev.off()

# Examine how the correlation results are tied to measured protein expression
coherence_results <- coherence_results %>%
  mutate(Correlation_Level = case_when(
    Coherence > 0.5 ~ "High (>0.5)",
    Coherence >= 0.3 & Coherence <= 0.5 ~ "Moderate (0.3-0.5)",
    Coherence < 0.3 ~ "Low (<0.3)"
  ))

adt_long_clean <- as.data.frame(t(as.matrix(adt_data))) %>%
  mutate(Cell_ID = rownames(.)) %>%
  pivot_longer(
    cols = -Cell_ID, 
    names_to = "Protein", 
    values_to = "Expression"
  )
# Join data
adt_long_colored <- adt_long_clean %>%
  left_join(coherence_results %>% select(Protein, Correlation_Level, Coherence), by = "Protein")

# Order according to prior plot 
ordered_proteins <- coherence_results %>%
  arrange(Coherence) %>% 
  pull(Protein)
clean_levels <- ordered_proteins[!is.na(ordered_proteins)]
adt_long_colored$Protein <- factor(adt_long_colored$Protein, levels = clean_levels)
adt_long_colored$Correlation_Level <- factor(adt_long_colored$Correlation_Level, 
                                             levels = c("High (>0.5)", "Moderate (0.3-0.5)", "Low (<0.3)"))
adt_long_colored <- adt_long_colored[!is.na(adt_long_colored$Protein),]

pdf(paste0(images, "/CITEseq_Protein_Expression_Clean_Boxplots.pdf"), width = 8, height = 15)

ggplot(adt_long_colored, aes(x = Expression, y = Protein, fill = Correlation_Level)) +
  # We removed geom_jitter() to eliminate all dots/points
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  scale_fill_manual(values = c(
    "High (>0.5)" = "steelblue", 
    "Moderate (0.3-0.5)" = "orange", 
    "Low (<0.3)" = "firebrick"
  )) +
  xlim(0,40)+
  labs(
    title = "Protein Expression Distribution by Coherence Category",
    subtitle = "Boxplots only (dots and outliers removed)",
    x = "Normalized Expression (Mean per Group)",
    y = "Protein",
    fill = "Coherence Level"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

dev.off()


pdf(paste0(images, "/CITEseq_Protein_Expression_3Column_Facets.pdf"), width = 14, height = 10)

ggplot(adt_long_colored, aes(x = Expression, y = Protein, fill = Correlation_Level)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  
  # split into 3 columns by level
  facet_wrap(~Correlation_Level, ncol = 3, scales = "free_y") +
  
  scale_fill_manual(values = c(
    "High (>0.5)" = "steelblue", 
    "Moderate (0.3-0.5)" = "orange", 
    "Low (<0.3)" = "firebrick"
  )) +
  xlim(0, 40) +
  labs(
    title = "Protein Expression Distribution by Coherence Category",
    subtitle = "Columns split by Pearson Correlation Strength",
    x = "Normalized Expression",
    y = "Protein",
    fill = "Coherence Level"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7), 
    legend.position = "none",            
    panel.spacing = unit(1, "lines"),    
    strip.background = element_rect(fill = "grey90", color = NA), 
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

dev.off()




# Add expression data to a new data frame
ct_groups <- factor(Ileum_filt$secondlevel)
ct_selection_mat <- sparse.model.matrix(~ 0 + ct_groups)

adt_ct_sum <- adt_data %*% ct_selection_mat
ct_counts <- as.numeric(table(ct_groups))
ct_diag_inv <- Diagonal(x = 1/ct_counts)

adt_ct_avg <- as.matrix(adt_ct_sum %*% ct_diag_inv)
colnames(adt_ct_avg) <- levels(ct_groups)

adt_ct_metrics <- data.frame(
  Protein = rownames(adt_ct_avg),
  Mean_ADT_Global = rowMeans(adt_ct_avg),
  Median_ADT_Global = apply(adt_ct_avg, 1, median),
  Range_ADT_Global = apply(adt_ct_avg, 1, function(x) max(x) - min(x))
)


coherence_final_grouped <- coherence_results %>%
  select(Protein, Gene, Coherence, Correlation_Level) %>%
  left_join(adt_ct_metrics, by = "Protein")

coherence_final_grouped$Correlation_Level <- factor(
  coherence_final_grouped$Correlation_Level, 
  levels = c("High (>0.5)", "Moderate (0.3-0.5)", "Low (<0.3)")
)

# Boxplot Protein expression for the categories
pdf(paste0(images, "/CITEseq_Coherence_vs_Expression_Summary.pdf"), width = 8, height = 6)

ggplot(coherence_final_grouped, aes(x = Correlation_Level, y = Mean_ADT_Global, fill = Correlation_Level)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c(
    "High (>0.5)" = "steelblue", 
    "Moderate (0.3-0.5)" = "orange", 
    "Low (<0.3)" = "firebrick"
  )) +
  labs(
    title = "Protein Expression Levels by Coherence Category",
    subtitle = "Distributions of Mean Protein Expression",
    x = "RNA-Protein Coherence Category",
    y = "Mean Protein Expression (DSB Normalized)",
    fill = "Coherence Level"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black")
  )

dev.off()


MLN_filt <- qs_read("/path/to/directory/MIST/MIST_Analysis/Seurat_Files/Analysis2_MLN_clustered_annotated.qs2")

# Extract normalized layers and aggregate by group (sum and average expression per sample_celltype)
MLN_filt[["ADT"]] <- JoinLayers(MLN_filt[["ADT"]])
adt_data <- GetAssayData(MLN_filt, assay = "ADT", layer = "data")
rna_data <- GetAssayData(MLN_filt, assay = "RNA", layer = "data")

MLN_filt$pb_cite <- paste(MLN_filt$Mouse, MLN_filt$secondlevel, sep = "_")
group_factor <- factor(MLN_filt$pb_cite)
selection_mat <- sparse.model.matrix(~ 0 + group_factor)


adt_sum <- adt_data %*% selection_mat
rna_sum <- rna_data %*% selection_mat

cell_counts <- as.numeric(table(group_factor))
diag_inv <- Diagonal(x = 1/cell_counts)

adt_pb_avg <- as.matrix(adt_sum %*% diag_inv)
rna_pb_avg <- as.matrix(rna_sum %*% diag_inv)

# Add group names 
colnames(adt_pb_avg) <- levels(group_factor)
colnames(rna_pb_avg) <- levels(group_factor)



# For each pair of the Protein- RNA perform pearson correlation across all pseudobulked levels
coherence_results <- adt_to_gene %>%
  mutate(Coherence = map2_dbl(Protein, Gene, function(p, g) {
    # Check if both exist in our matrices
    if (p %in% rownames(adt_pb_avg) && g %in% rownames(rna_pb_avg)) {
      # Calculate Pearson correlation 
      return(cor(adt_pb_avg[p, ], rna_pb_avg[g, ], method = "pearson"))
    } else {
      return(NA)
    }
  })) %>%
  filter(!is.na(Coherence)) %>%
  arrange(desc(Coherence))

# Save correlations
write.csv(coherence_results, "/path/to/directory/MIST/MIST_Analysis/Coherence_Results_MLN.csv", row.names = FALSE)

# Create categories for results
coherence_results <- coherence_results %>%
  mutate(Correlation_Level = case_when(
    Coherence > 0.5 ~ "High (>0.5)",
    Coherence >= 0.3 & Coherence <= 0.5 ~ "Moderate (0.3-0.5)",
    Coherence < 0.3 ~ "Low (<0.3)"
  )) %>%
  # legend order
  mutate(Correlation_Level = factor(Correlation_Level, 
                                    levels = c("High (>0.5)", "Moderate (0.3-0.5)", "Low (<0.3)")))


pdf(paste0(images, "/CITEseq_Coherence_Barplot_3Colors_MLN.pdf"), width = 8, height = 15)

ggplot(coherence_results, aes(x = reorder(Protein, Coherence), y = Coherence)) +
  geom_col(aes(fill = Correlation_Level)) +
  scale_fill_manual(values = c(
    "High (>0.5)" = "steelblue", 
    "Moderate (0.3-0.5)" = "orange", 
    "Low (<0.3)" = "firebrick"
  )) +
  coord_flip() +
  labs(
    title = "RNA-Protein Coherence (Pearson Correlation)",
    subtitle = "Categorized by Correlation Strength",
    x = "Protein/Gene Pair",
    y = "Pearson Correlation (r)",
    fill = "Coherence Level"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  )

dev.off()

# Examine how the correlation results are tied to measured protein expression
coherence_results <- coherence_results %>%
  mutate(Correlation_Level = case_when(
    Coherence > 0.5 ~ "High (>0.5)",
    Coherence >= 0.3 & Coherence <= 0.5 ~ "Moderate (0.3-0.5)",
    Coherence < 0.3 ~ "Low (<0.3)"
  ))

# Join data
adt_long_colored <- adt_long_clean %>%
  left_join(coherence_results %>% select(Protein, Correlation_Level), by = "Protein")

# Order according to prior plot 
clean_levels <- rev(ordered_proteins[!is.na(ordered_proteins)])
adt_long_colored$Protein <- factor(adt_long_colored$Protein, levels = clean_levels)
adt_long_colored$Correlation_Level <- factor(adt_long_colored$Correlation_Level, 
                                             levels = c("High (>0.5)", "Moderate (0.3-0.5)", "Low (<0.3)"))


pdf(paste0(images, "/CITEseq_Protein_Expression_Clean_Boxplots_MLN.pdf"), width = 8, height = 15)

ggplot(adt_long_colored, aes(x = Expression, y = Protein, fill = Correlation_Level)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  scale_fill_manual(values = c(
    "High (>0.5)" = "steelblue", 
    "Moderate (0.3-0.5)" = "orange", 
    "Low (<0.3)" = "firebrick"
  )) +
  xlim(0,40)+
  labs(
    title = "Protein Expression Distribution by Coherence Category",
    subtitle = "Boxplots only (dots and outliers removed)",
    x = "Normalized Expression (Mean per Group)",
    y = "Protein",
    fill = "Coherence Level"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

dev.off()







# Add expression data to a new data frame
ct_groups <- factor(MLN_filt$secondlevel)
ct_selection_mat <- sparse.model.matrix(~ 0 + ct_groups)

adt_ct_sum <- adt_data %*% ct_selection_mat
ct_counts <- as.numeric(table(ct_groups))
ct_diag_inv <- Diagonal(x = 1/ct_counts)

adt_ct_avg <- as.matrix(adt_ct_sum %*% ct_diag_inv)
colnames(adt_ct_avg) <- levels(ct_groups)

adt_ct_metrics <- data.frame(
  Protein = rownames(adt_ct_avg),
  Mean_ADT_Global = rowMeans(adt_ct_avg),
  Median_ADT_Global = apply(adt_ct_avg, 1, median),
  Range_ADT_Global = apply(adt_ct_avg, 1, function(x) max(x) - min(x))
)


coherence_final_grouped <- coherence_results %>%
  select(Protein, Gene, Coherence, Correlation_Level) %>%
  left_join(adt_ct_metrics, by = "Protein")

coherence_final_grouped$Correlation_Level <- factor(
  coherence_final_grouped$Correlation_Level, 
  levels = c("High (>0.5)", "Moderate (0.3-0.5)", "Low (<0.3)")
)

# Boxplot Protein expression for the categories
pdf(paste0(images, "/CITEseq_Coherence_vs_Expression_Summary_MLN.pdf"), width = 8, height = 6)

ggplot(coherence_final_grouped, aes(x = Correlation_Level, y = Mean_ADT_Global, fill = Correlation_Level)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c(
    "High (>0.5)" = "steelblue", 
    "Moderate (0.3-0.5)" = "orange", 
    "Low (<0.3)" = "firebrick"
  )) +
  labs(
    title = "Protein Expression Levels by Coherence Category",
    subtitle = "Distributions of Mean Protein Expression",
    x = "RNA-Protein Coherence Category",
    y = "Mean Protein Expression (DSB Normalized)",
    fill = "Coherence Level"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black")
  )

dev.off()
