# The Purpose of this Script is to visulize the Genes in space to validate our scRNA pseudotime define expression patterns 
# Figures S8E, Figure S8F, were made using XeniumExplorer
# Figure S8G Produced in this script
# Load Libraries----
# Analysis 
library(Seurat)
library(Signac)
library(arrow)

library(sf) #Used by functions in Seurat to process some types of Xenium data
#visualization
library(ggplot2)
library(ggbreak)
library(ggforce)
library(ggrepel)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(scCustomize)
library(gridExtra)
library(gt)
library(ggpubr)
#data manipulation
library(dplyr)
library(reshape2)
library(data.table)
library(tidyverse)
library(forcats)
library(tibble)


#system
library(future)
library(remotes)
options(future.globals.maxSize = 20000 * 1024^2)
library(qs2)

# Loaded the all_objects_filtered Data 
Seurat_obj <-qs2::qs_read("/path/to/directoryXenium_Yp_2025/Seurat/Seurat_all_proseg_Niches.qs2")

#  slc genes from Figure 4G w


slc_genes <- c(
  "Slc39a10", "Slc23a3", "Slc41a1", "Slc25a25", "Slc1a2", "Slc24a5", "Slc27a2", 
  "Slc24a3", "Slc22a15", "Slc44a3", "Slc39a8", "Slc2a1", "Slc25a33", "Slc7a1", 
  "Slc25a26", "Slc1a5", "Slco3a1", "Slc7a2", "Slc7a5", "Slc16a7", "Slc16a6", 
  "Slc39a11", "Slc4a7", "Slc25a30", "Slc1a3", "Slc45a4", "Slc2a13", "Slc38a1", 
  "Slc12a8", "Slc5a3", "Slc22a3", "Slc29a1", "Slc25a23", "Slc39a6", "Slc12a2", 
  "Slc14a1", "Slc66a2", "Slc25a39", "Slc25a17", "Slc22a1", "Slc38a2", "Slc39a7", 
  "Slc25a3", "Slco2b1", "Slc40a1", "Slc20a1", "Slc25a11", "Slc25a5", "Slc28a3", 
  "Slc6a6", "Slc1a1", "Slc35b1", "Slc25a24", "Slc4a4", "Slc25a20", "Slc39a5", 
  "Slc16a5", "Slc35g1", "Slc23a2", "Slc26a3", "Slc30a1", "Slc25a4", "Slc51b", 
  "Slc7a8", "Slc51a", "Slc6a20a", "Slc9a3", "Slc25a45", "Slc10a2", "Slc13a1", 
  "Slc6a20b", "Slc2a5", "Slc5a11", "Slc2a2", "Slc28a1", "Slc5a4b", "Slc5a4a", 
  "Slc25a48", "Slc16a10", "Slc5a12", "Slc28a2", "Slc23a4"
)
slc_genes <- c(slc_genes, "Lgr5", "Reg3b", "Reg3a", "Slc5a1", "Ada")
slc_genes <- slc_genes[slc_genes %in% rownames(Seurat_obj[["Xenium"]])]
slc_genes <- unique(slc_genes)
#Seurat_obj <- UpdateSeuratObject(Seurat_obj)




ImageFeaturePlot(
  Seurat_obj, 
  features = "Lgr5", 
  fov = "yersinia1", axes = T, molecules = "Lgr5"
)

ImageFeaturePlot(
  Seurat_obj,
  features = "Lgr5",
  fov = "yersinia1", 
  axes = TRUE, 
  molecules = "Lgr5", size = 2, mols.size = 5,
) + coord_cartesian(xlim = c(5300, 5900), ylim = c(3500, 4500))



# Fli
ImageFeaturePlot(
  Seurat_obj,
  features = "Lgr5",
  fov = "yersinia1", 
  axes = TRUE, 
  #molecules = "Lgr5", 
  size = 2, 
  #mols.size = 5, 
max.cutoff = 1
) + aes(size = Lgr5) + 
  # This controls the range of the dot sizes (e.g., from size 0.5 to 3)
  scale_size_continuous(range = c(0.5, 3)) +
  coord_flip(xlim = c(5300, 5900), ylim = c(3500, 4500)) +
  scale_x_reverse()


ImageFeaturePlot(
  Seurat_obj,
  features = "Lgr5",
  fov = "yersinia1", 
  axes = TRUE, 
  max.cutoff = 1
) + 
  aes(size = .data[["Lgr5"]]) + 
  scale_size_continuous(range = c(1, 4)) +
  coord_flip(xlim = c(5300, 5900), ylim = c(3500, 4500)) +
  scale_x_reverse()

# Generate the plot object
p <- ImageFeaturePlot(Seurat_obj, features = "Lgr5", fov = "yersinia1")

head(p$data)

plot_data <- p$data

ggplot(plot_data, aes(x = x, y = y, color = Lgr5, size = Lgr5)) +
  geom_point() +
  # Apply your max cutoff of 1 to the color scale
  scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0, 6), oob = scales::squish) +
  # Control the dot size range here
  scale_size_continuous(range = c(2, 6)) +
  # Apply your specific crop and rotation
  coord_flip(xlim = c(3500, 4500), ylim = c(5300, 5900)) +
  theme_minimal() +
  labs(title = "Lgr5 Expression (Manual Scale)")

unique(all_objects_filtered$secondlevel)
VlnPlot(all_objects_filtered, "Slc40a1", pt.size = 0)

all_objects_filtered$secondlevel[all_objects_filtered$secondlevel %in% c("Early Enterocyte", "Middle Enterocyte", "Late Enterocyte")] <- "Enterocyte"


# Fetch expression and metadata into a data frame
plot_df <- FetchData(
  all_objects_filtered, 
  vars = c("Slc40a1", "orig.ident", "secondlevel")
)

# filter
plot_df_filtered <- plot_df %>%
  filter(secondlevel == "Enterocyte") %>%
  filter(orig.ident %in% c("yersinia1", "naive1"))

ggplot(plot_df_filtered, aes(x = orig.ident, y = Slc40a1, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) + 
  theme_classic() +
  labs(
    title = "Slc40a1 Expression in Enterocytes",
    subtitle = "Comparison: Naive1 vs Yersinia1",
    x = "Sample Condition",
    y = "Expression Level"
  ) +
  scale_fill_manual(values = c("naive1" = "#66c2a5", "yersinia1" = "#fc8d62"))
+ stat_compare_means(method = "wilcox.test", label = "p.signif")

ggplot(plot_df_filtered, aes(x = orig.ident, y = Slc40a1, fill = orig.ident)) +
  # jitter
  geom_jitter(aes(color = orig.ident), alpha = 0.2, width = 0.2) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  coord_cartesian(ylim = c(0.01, max(plot_df_filtered$Slc40a1))) + 
  theme_minimal()

plot_df_filtered %>%
  group_by(orig.ident) %>%
  summarize(
    percent_expressing = sum(Slc40a1 > 0) / n() * 100,
    mean_expression = mean(Slc40a1)
  )



ggplot(plot_df_filtered, aes(x = orig.ident, y = Slc40a1, fill = orig.ident)) +
  geom_jitter(aes(color = orig.ident), size = 0.5, alpha = 0.3, 
              position = position_jitter(width = 0.15)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white")+
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
  theme_classic() +
  labs(title = "Slc40a1 Expression: Naive vs Yersinia",
       subtitle = "Enterocyte Population Only",
       y = "Log-Normalized Expression",
       x = "Condition") +
  scale_fill_manual(values = c("naive1" = "#999999", "yersinia1" = "#D55E00")) +
  scale_color_manual(values = c("naive1" = "#999999", "yersinia1" = "#D55E00"))

# Create the summary text string
summary_text <- paste0(
  "Naive: 8.5% expr, Mean 0.066\n",
  "Yersinia: 25.5% expr, Mean 0.244"
)

ggplot(plot_df_filtered, aes(x = orig.ident, y = Slc40a1, fill = orig.ident)) +
  geom_jitter(aes(color = orig.ident), size = 0.5, alpha = 0.3, 
              position = position_jitter(width = 0.15)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
  
  # Add the summary box in the upper right
  annotate("label", x = Inf, y = Inf, label = summary_text, 
           hjust = 1.1, vjust = 1.1, fill = "white", 
           size = 3.5, fontface = "italic", alpha = 0.8) +
  
  theme_classic() +
  labs(title = "Slc40a1 Expression: Naive vs Yersinia",
       subtitle = "Enterocyte Population Only",
       y = "Log-Normalized Expression",
       x = "Condition") +
  scale_fill_manual(values = c("naive1" = "#999999", "yersinia1" = "#D55E00")) +
  scale_color_manual(values = c("naive1" = "#999999", "yersinia1" = "#D55E00"))

ggsave(filename = "Slc_visualizatino_violin_Enterocytes_Slc40a1.svg", plot = last_plot(), path = images, width = 6, height = 5, dpi = 600)



