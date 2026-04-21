---
title: "Analysis_2.5_Figures"
author: "AndrewHart"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Description. Using the annotated data for the MLn and Ileu. This script will explore that data and create major visualization. 
# Figure Panels made: Figure 2A, Figure 2B, Figure S2A, Figure S2B, Figure 2D, Figure S2G, Figure 2C, Figure S2F, Figure S2D, Figure S2E, Figure S2C, Figure S1D, Figure S3A, Figure 2E, Figure 3A, Figure 3B
#Load Libraries 
```{r}
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
options(future.globals.maxSize = 5000 * 1024^2)
```

#Establish Directories
```{r} 
getwd() #/home/hartandrew
setwd("/path/to/analysis/directory")
out <- "/path/to/analysis/directory/Output"
version <- "/path/to/analysis/directory/Version"
seurat <- "/path/to/analysis/directory/Seurat_Files"
images <- "/path/to/analysis/directory/Images"
CSV <- "/path/to/analysis/directory/CSV"

```
#Set colors
```{r}
colors7 <- c("#980339", "#379FFB", "#D6A206", "#027102", "#6313B7", "#AFAFAF","#542600" )
colors6 <- c("#980339", "#379FFB", "#D6A206", "#027102", "#6313B7", "#AFAFAF")
colors5 <- c("#980339", "#379FFB", "#D6A206", "#027102", "#6313B7")
colors8 <- c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
colors9 <- c("#332288", "#88CCEE", "#117733","#44AA99", "#A27DF1", "#DDCC77",  "#AA4499", "#CC6677", "#882255")
infectioncolors <-c(  "#B78CBA", 
  "#1170AA",  
  "#FCB13F",  
  "#60BFC1",  
  "#EF6F6A",  
    "#937860", 
  "#D17A00") 
 
```

#Load Objects
```{r}

MLN_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2_MLN_clustered_annotated.qs2")
Ileum_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Analysis2_clustered_annotated_Ileum.qs2")

```


```{r}
#Graphing UMAPS
# Figure2B
library(ggrepel)
umap_data <- Ileum_filt[["wnn.umap_cc"]]@cell.embeddings
cluster_data <- Ileum_filt$CoarseCellType

# Combine UMAP coordinates and cluster information into a single dataframe
umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cluster_data)
centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnUMAPcc_1 = mean(wnnUMAPcc_1), wnnUMAPcc_2 = mean(wnnUMAPcc_2))

ileum_umap <- ggplot(umap_data, aes(x = wnnUMAPcc_1, y = wnnUMAPcc_2, color = cluster)) +
  geom_point(alpha = 0.1, size = 1.0) +  # Add points for each cell
  ylim(-15, 20) + xlim(-15, 20) +
  #geom_mark_hull(aes(fill = cluster, alpha = 0.01), alpha = 0.1, color = "NA", concavity = 10,    con.type = "none") +  
  scale_color_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
     theme(panel.background = element_rect(fill = "white", 
            colour = NA), panel.grid = element_blank(),
            panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
           axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          legend.position = "none",
            strip.background = element_blank())
ileum_umap + geom_text_repel(data = centroids, aes(label = cluster), size = 4, fontface = "bold", 
                  color = "black",
                 #box.padding = 0.1,  
                #  point.padding = 0.5,  
                   min.segment.length = 0
                  ) 
ggsave( "Ileum_UMAP_HULL_Trial.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
ileum_umap
ggsave( "Ileum_UMAP_HULL_Trial_NoLabels.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


ggplot(umap_data, aes(x = wnnUMAPcc_1, y = wnnUMAPcc_2, color = cluster)) +
  geom_point(alpha = 0.2, size = 1.0) +  
  ylim(-15, 20) + xlim(-15, 20) +
  #geom_mark_hull(aes(fill = cluster), alpha = 0.1, color = "NA", concavity = 500) +  
  scale_color_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
     theme(panel.background = element_rect(fill = "black", 
            colour = NA), panel.grid = element_blank(),
            panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
           axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          legend.position = "right",
            strip.background = element_blank())
ggsave( "Ileum_UMAP_HULL_Trial_blackBackground.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

ggsave( "Ileum_UMAP_HULL_Trial_blackBackground.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- MLN_filt[["wnn.umap_cc"]]@cell.embeddings
cluster_data <- MLN_filt$CoarseCellType

# Combine UMAP coordinates and cluster information into a single dataframe
umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cluster_data)
centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnUMAPcc_1 = mean(wnnUMAPcc_1), wnnUMAPcc_2 = mean(wnnUMAPcc_2))

mln_umap <- ggplot(umap_data, aes(x = wnnUMAPcc_1, y = wnnUMAPcc_2, color = cluster)) +
  geom_point(alpha = 0.1, size = 1.0) +  
  ylim(-15, 20) + xlim(-15, 20) +
  #geom_mark_hull(aes(fill = cluster), alpha = 0.1, color = "NA", concavity = 500) + 
  scale_color_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  scale_fill_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
     theme(panel.background = element_rect(fill = "white", 
            colour = NA), panel.grid = element_blank(),
            panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
           axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          legend.position = "none",
            strip.background = element_blank())
mln_umap + geom_text_repel(data = centroids, aes(label = cluster), size = 4, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) 
ggsave( "MLN_UMAP_HULL_Trial.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
mln_umap
ggsave( "MLN_UMAP_HULL_Trial_noLabels.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


ggplot(umap_data, aes(x = wnnUMAPcc_1, y = wnnUMAPcc_2, color = cluster)) +
  geom_point(alpha = 0.2, size = 1.0) +  
  ylim(-15, 20) + xlim(-15, 20) +
  #geom_mark_hull(aes(fill = cluster), alpha = 0.1, color = "NA", concavity = 500) + 
  scale_color_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  scale_fill_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + 
  labs( x = "UMAP 1", y = "UMAP 2") + theme_minimal() +
     theme(panel.background = element_rect(fill = "black", 
            colour = NA), panel.grid = element_blank(),
            panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
           axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          legend.position = "right",
            strip.background = element_blank())
ggsave( "MLN_UMAP_HULL_Trial_blackBackground.png",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

ggsave( "MLN_UMAP_HULL_Trial_blackBackground.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```

```{r}

#barcharts - Total Cells per InfectionStatus
data <- Ileum_filt@meta.data %>%
  group_by(InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Yersinia", "SFB_YA", "Nippostrongylus", "MNV", "Cryptosporidium",  "Candida", "Naive"   ))
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = (count / total)*100)

stacked_bar <- ggplot(data, aes(x = count, y = InfectionStatus, color = "black")) +
  geom_bar(aes(fill = "lightgrey"), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("black")) + 
  scale_fill_manual(values = c("lightgrey")) + 
  labs(x = "Cell Number") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.y = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
             text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
         legend.position = "none")
stacked_bar
ggsave( "TotalCellCounts_by_infection_Ileum.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 4,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)



# Ensure factor order
data <- Ileum_filt@meta.data
data$InfectionStatus <- factor(data$InfectionStatus,levels = c("Yersinia", "SFB_YA", "Nippostrongylus",
    "MNV", "Cryptosporidium", "Candida", "Naive"))

# Calculate mean nFeature_RNA per infection
mean_features <- data %>%
  group_by(InfectionStatus) %>%
  summarise(mean_nFeature = mean(nFeature_RNA, na.rm = TRUE))

ggplot(mean_features, aes(x = InfectionStatus, y = mean_nFeature)) +
  geom_bar(stat = "identity", color = "black", fill = "lightgrey") +
  labs(x = "Infection Status", y = "Mean nFeature_RNA") +

  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )+ coord_flip()

#figure S2A
ggsave( "Mean_nFeature_RNA_by_Infection_Ileum.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 4,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)







data2 <- MLN_filt@meta.data %>%
  group_by(InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
data2$InfectionStatus <- factor(data2$InfectionStatus, levels = c("Yersinia", "SFB_YA", "Nippostrongylus", "MNV", "Cryptosporidium",  "Candida", "Naive"   ))
# Calculate total number of cells per condition for normalization
total_counts <- data2 %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data2 <- data2 %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = (count / total)*100)

stacked_bar <- ggplot(data2, aes(x = count, y = InfectionStatus, color = "black")) +
  geom_bar(aes(fill = "lightgrey"), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("black")) + 
  scale_fill_manual(values = c("lightgrey")) + 
  labs(x = "Cell Number") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.y = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
             text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
         legend.position = "none") 
stacked_bar
ggsave( "TotalCellCounts_by_infection_MLN.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 4,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


# Ensure factor order
data <- MLN_filt@meta.data
data$InfectionStatus <- factor(data$InfectionStatus,levels = c("Yersinia", "SFB_YA", "Nippostrongylus",
    "MNV", "Cryptosporidium", "Candida", "Naive"))

# Calculate mean nFeature_RNA per infection
mean_features <- data %>%
  group_by(InfectionStatus) %>%
  summarise(mean_nFeature = mean(nFeature_RNA, na.rm = TRUE))

ggplot(mean_features, aes(x = InfectionStatus, y = mean_nFeature)) +
  geom_bar(stat = "identity", color = "black", fill = "lightgrey") +
  labs(x = "Infection Status", y = "Mean nFeature_RNA") +

  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )+ coord_flip()

#figure S2A
ggsave( "Mean_nFeature_RNA_by_Infection_MLN.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 4,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```

```{r}
fulldata <- left_join(data, data2, by = "InfectionStatus")
fulldata$count <- fulldata$count.x + fulldata$count.y

fulldata$InfectionStatus <- factor(fulldata$InfectionStatus, levels = c("Yersinia", "SFB_YA", "Nippostrongylus", "MNV", "Cryptosporidium",  "Candida", "Naive"   ))
# Calculate total number of cells per condition for normalization
total_counts <- fulldata %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data2 <- data2 %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = (count / total)*100)

stacked_bar <- ggplot(fulldata, aes(x = count, y = InfectionStatus, color = "black")) +
  geom_bar(aes(fill = "lightgrey"), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("black")) + 
  scale_fill_manual(values = c("lightgrey")) + 
  labs(x = "Cell Number") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.y = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
             text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
         legend.position = "none")
stacked_bar
# Figure S2B
ggsave( "TotalCellCounts_byinfecion_MLN_and_Ileum.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 4,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```

#Graphing Major Cell Types
```{r}
unique(Ileum_filt$MajorCellType)
unique(MLN_filt$MajorCellType)
colors <- c("#CE864D", "#D07C3A","#D27227", "#D36813",  "#D55E00",   "#CC79A7", "#F0E442","#8DD4C0", "#71C9B1", "#55BEA1", "#38B492", "#1CA982", "#009E73","#ABD7F1", "#9AD0EF", "#89C9EE", "#78C2EC", "#67BBEB",  "#56B4E9","#E60000", "#0072B2", "#014971")
MajorCellTypes <- c("gd T Cells", "CD4 ab T Cells", "CD8 ab T Cells", "Innate Lymphoid Cells", "NKT Cells", "B Cells", "Plasma Cells", "Dendritic Cells", "Plasmacytoid DCs", "Macrophages", "Monocytes", "Neutrophils", "Mast Cells", "Enterocytes", "Goblet Cells", "Paneth Cells", "Enteroendocrine Cells", "Tuft Cells", "Stem / Transit Amplifying Cells", "Enteric Nervous System", "Mesenchymal", 
  "Endothelial")

```

```{r}
#Total filtered Cell plot 
tbl_ileum <- as.data.frame(table(Ileum_filt$InfectionStatus))
tbl_mln <- as.data.frame(table(MLN_filt$InfectionStatus))

tbl_both <- right_join(tbl_ileum, tbl_mln, by = "Var1")
tbl_both$Freq.x[is.na(tbl_both$Freq.x)] <-0
tbl_both$Combined <- tbl_both$Freq.x + tbl_both$Freq.y

tbl_both$Var1 <-  factor(tbl_both$Var1 , levels = c("Naive", "Candida", "Cryptosporidium",  "MNV", "Nippostrongylus" , "SFB_YA", "Yersinia"))

ggplot(tbl_both, aes(x = Var1, y = Combined), color = "black") +
  geom_col(color = "black" , fill = "black") +
  labs(title = "Cell Count by Condition") + 
  ylab("Total Cell Number MLN and Ileum") + theme_minimal() +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), 
            panel.grid.minor = element_blank(), axis.title.x = element_blank(),
          legend.position = "none",
            strip.background = element_blank())
ggsave( "TotalCell.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 3,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```


```{r}
library(patchwork)

data <- Ileum_filt@meta.data %>%
  group_by(CoarseCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = (count / total)*100)

stacked_bar <- ggplot(data, aes(x = InfectionStatus, y = count, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar



stacked_bar2 <- ggplot(data, aes(x = InfectionStatus, y = frequency, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar2



data <- Ileum_filt@meta.data %>%
  group_by(MajorCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = count / total)

data$MajorCellType <- factor(data$MajorCellType , levels = MajorCellTypes)
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))
dotplot <- ggplot(data, aes(x = InfectionStatus, y = MajorCellType, size = frequency, color = MajorCellType)) +
  geom_point(aes(color = MajorCellType)) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  scale_color_manual(values = colors, breaks =MajorCellTypes) + 
  #scale_fill_manual(values = colors, breaks = MajorCellTypes) +
   guides(alpha = "none", color = "none") +
  labs(size = "Frequency") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))
dotplot 


# Calculate total count by cell type
total_by_type <- data %>%
  group_by(MajorCellType) %>%
  summarise(total_count = sum(count))
total_by_type$MajorCellType <- factor(total_by_type$MajorCellType , levels = MajorCellTypes)
# Horizontal bar chart # Will Need to be assigned new cell types and colors once we establish what we will show 
horizontal_bar <- ggplot(total_by_type, aes(x = total_count, y = MajorCellType)) +
  geom_bar(aes(fill = total_by_type$MajorCellType), stat = "identity") +
  theme_minimal() +scale_fill_manual(values = colors, breaks = MajorCellTypes) +
  labs(title = "Total Cell Type Counts", fill = "Major Cell Type", x = "Cell Number") +
  theme(axis.text.y = element_text(hjust = 1), legend.position = "none", 
        axis.title.y = element_blank())
horizontal_bar

combined_plot <-(stacked_bar ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "Ileum_filt_Total_fig2_Count.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

horizontal_bar
ggsave( "Ileum_filt_fig2_Number.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

combined_plot <-(stacked_bar2 ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "Ileum_filt_Total_fig2_Frequency.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```

```{r}

data <- Ileum_filt@meta.data %>%
  group_by(Tissue, CoarseCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
data <- data[data$Tissue == "Intestinal Epithelial Cells", ]
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = (count / total)*100)

stacked_bar <- ggplot(data, aes(x = InfectionStatus, y = count, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7",  "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7",  "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar



stacked_bar2 <- ggplot(data, aes(x = InfectionStatus, y = frequency, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7",  "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar2



data <- Ileum_filt@meta.data %>%
  group_by(Tissue, MajorCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
data <- data[data$Tissue == "Intestinal Epithelial Cells", ]
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = count / total)

data$MajorCellType <- factor(data$MajorCellType , levels = MajorCellTypes)
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))
dotplot <- ggplot(data, aes(x = InfectionStatus, y = MajorCellType, size = frequency, color = MajorCellType)) +
  geom_point(aes(color = MajorCellType)) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  scale_color_manual(values = colors, breaks =MajorCellTypes) + 
  #scale_fill_manual(values = colors, breaks = MajorCellTypes) +
   guides(alpha = "none", color = "none") +
  labs(size = "Frequency") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))
dotplot 


# Calculate total count by cell type
total_by_type <- data %>%
  group_by(MajorCellType) %>%
  summarise(total_count = sum(count))
total_by_type$MajorCellType <- factor(total_by_type$MajorCellType , levels = MajorCellTypes)
# Horizontal bar chart # Will Need to be assigned new cell types and colors once we establish what we will show 
horizontal_bar <- ggplot(total_by_type, aes(x = total_count, y = MajorCellType)) +
  geom_bar(aes(fill = total_by_type$MajorCellType), stat = "identity") +
  theme_minimal() +scale_fill_manual(values = colors, breaks = MajorCellTypes) +
  labs(title = "Total Cell Type Counts", fill = "Major Cell Type", x = "Cell Number") +
  theme(axis.text.y = element_text(hjust = 1), legend.position = "none", 
        axis.title.y = element_blank())
horizontal_bar

#Figure 2D and Figure S2G
combined_plot <-(stacked_bar ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "Ileum_filt_IEC_fig2_Count.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

horizontal_bar
ggsave( "Ileum_filt_IEC_fig2_Number.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

combined_plot <-(stacked_bar2 ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "Ileum_filt_IEC_fig2_Frequency.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```

```{r}
data <- Ileum_filt@meta.data %>%
  group_by(Tissue, CoarseCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
data <- data[data$Tissue == "Lamina Propria", ]
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = (count / total)*100)

stacked_bar <- ggplot(data, aes(x = InfectionStatus, y = count, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar



stacked_bar2 <- ggplot(data, aes(x = InfectionStatus, y = frequency, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7","#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7", "#E60000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar2



data <- Ileum_filt@meta.data %>%
  group_by(Tissue, MajorCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
data <- data[data$Tissue == "Lamina Propria", ]
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = count / total)

data$MajorCellType <- factor(data$MajorCellType , levels = MajorCellTypes)
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))
dotplot <- ggplot(data, aes(x = InfectionStatus, y = MajorCellType, size = frequency, color = MajorCellType)) +
  geom_point(aes(color = MajorCellType)) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  scale_color_manual(values = colors, breaks =MajorCellTypes) + 
  #scale_fill_manual(values = colors, breaks = MajorCellTypes) +
   guides(alpha = "none", color = "none") +
  labs(size = "Frequency") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))
dotplot 


# Calculate total count by cell type
total_by_type <- data %>%
  group_by(MajorCellType) %>%
  summarise(total_count = sum(count))
total_by_type$MajorCellType <- factor(total_by_type$MajorCellType , levels = MajorCellTypes)
# Horizontal bar chart # Will Need to be assigned new cell types and colors once we establish what we will show 
horizontal_bar <- ggplot(total_by_type, aes(x = total_count, y = MajorCellType)) +
  geom_bar(aes(fill = total_by_type$MajorCellType), stat = "identity") +
  theme_minimal() +scale_fill_manual(values = colors, breaks = MajorCellTypes) +
  labs(title = "Total Cell Type Counts", fill = "Major Cell Type", x = "Cell Number") +
  theme(axis.text.y = element_text(hjust = 1), legend.position = "none", 
        axis.title.y = element_blank())
horizontal_bar

#Figure 2D and Figure S2G
combined_plot <-(stacked_bar ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "Ileum_filt_LP_fig2_Count.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

horizontal_bar
ggsave( "Ileum_filt_LP_fig2_Number.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

combined_plot <-(stacked_bar2 ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "Ileum_filt_LP_fig2_Frequency.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```
```{r}

# Generate the plots so that we can see the frequency of cells across both tissues - to compare between IEC and LP rather than within LP and IEC

all_celltypes <- MajorCellTypes
all_conditions <- c("Naive", "Candida", "Cryptosporidium", "MNV",
                    "Nippostrongylus", "SFB_YA", "Yersinia")


data_all <- Ileum_filt@meta.data %>%
  group_by(Tissue, MajorCellType, InfectionStatus) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(Tissue,
           MajorCellType = all_celltypes,
           InfectionStatus = all_conditions,
           fill = list(count = 0))

#  Total cells across both tissues
total_per_condition <- data_all %>%
  group_by(InfectionStatus) %>%
  summarise(global_total = sum(count), .groups = "drop")

data_all <- data_all %>%
  left_join(total_per_condition, by = "InfectionStatus") %>%
  mutate(frequency = count / global_total)

#transform/ factor
data_all <- data_all %>%
  mutate(
    MajorCellType = factor(MajorCellType, levels = all_celltypes),
    InfectionStatus = factor(InfectionStatus, levels = all_conditions)
  )

# global scaling
global_max_freq <- max(data_all$frequency, na.rm = TRUE)


make_dotplot <- function(df, tissue_label) {
  df_tissue <- df %>% filter(Tissue == tissue_label)
  
  ggplot() +
    # plot only non-zero points, but provide the full y-axis range
    geom_point(
      data = df_tissue %>% filter(frequency > 0),
      aes(x = InfectionStatus, y = MajorCellType, size = frequency, color = MajorCellType)
    ) +
    scale_size_continuous(range = c(3, 10), limits = c(0, global_max_freq)) +
    scale_color_manual(values = colors, breaks = all_celltypes) +
    scale_y_discrete(drop = FALSE) +  # <- this ensures all MajorCellTypes are shown even if absent
    theme_minimal() +
    guides(alpha = "none", color = "none") +
    labs(title = tissue_label, size = "Frequency (of total cells per condition)") +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_text(size = 9)
    )
}
#Figure 2D and Figure S2G
dotplot_IEL <- make_dotplot(data_all, "Intestinal Epithelial Cells")
dotplot_LP  <- make_dotplot(data_all, "Lamina Propria")

dotplot_IEL + dotplot_LP + plot_layout(guides = "collect")
ggsave( "Ileum_bubbleplots_totalCondition_proportions_Fig2_Oct2025.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 11,  height = 10,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```


```{r}
colors_mln <- c("#CE864D", "#D07C3A","#D27227", "#D36813",  "#D55E00",   "#CC79A7", "#F0E442","#8DD4C0", "#71C9B1", "#55BEA1", "#38B492", "#1CA982", "#009E73", "#0072B2", "#014971")
MajorCellTypes_mln <- c("gd T Cells", "CD4 ab T Cells", "CD8 ab T Cells", "Innate Lymphoid Cells", "NKT Cells", "B Cells", "Plasma Cells", "Dendritic Cells", "Plasmacytoid DCs", "Macrophages", "Monocytes", "Neutrophils", "Basophils", "Mesenchymal",  "Endothelial")

```

```{r}
table(MLN_filt$MajorCellType)

data <- MLN_filt@meta.data %>%
  group_by(CoarseCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = (count / total)*100)
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))


stacked_bar <- ggplot(data, aes(x = InfectionStatus, y = count, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar

stacked_bar2 <- ggplot(data, aes(x = InfectionStatus, y = frequency, fill = CoarseCellType)) +
  geom_bar(aes(fill = CoarseCellType), stat = "identity") +
  theme_minimal() +
    scale_color_manual(values =c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Custom cluster colors
  scale_fill_manual(values = c("#CC79A7",  "#009E73", "#F0E442", "#0072B2", "#D55E00")) + # Fill color for hulls
  labs(title = "Cell Type Distribution by Condition", y = "Frequency") +
  theme( panel.background = element_rect(fill = "white", 
            colour = NA), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(), 
            strip.background = element_blank(), 
         axis.text.x = element_blank())
stacked_bar2


data <- MLN_filt@meta.data %>%
  group_by(MajorCellType, InfectionStatus) %>%
  summarise(count = n()) %>%
  ungroup()
# Calculate total number of cells per condition for normalization
total_counts <- data %>%
  group_by(InfectionStatus) %>%
  summarise(total = sum(count))
# Merge with original data to calculate frequency
data <- data %>%
  left_join(total_counts, by = "InfectionStatus") %>%
  mutate(frequency = count / total)

data$MajorCellType <- factor(data$MajorCellType , levels = MajorCellTypes_mln)
data$InfectionStatus <- factor(data$InfectionStatus, levels = c("Naive", "Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia"))
dotplot <- ggplot(data, aes(x = InfectionStatus, y = MajorCellType, size = frequency, color = MajorCellType)) +
  geom_point(aes(color = MajorCellType)) +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  scale_color_manual(values = colors_mln, breaks = MajorCellTypes_mln) + 
  #scale_fill_manual(values = colors, breaks = MajorCellTypes) +
   guides(alpha = "none", color = "none") +
  labs(size = "Frequency") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))
dotplot 


# Calculate total count by cell type
total_by_type <- data %>%
  group_by(MajorCellType) %>%
  summarise(total_count = sum(count))
total_by_type$MajorCellType <- factor(total_by_type$MajorCellType , levels = MajorCellTypes_mln)
# Horizontal bar chart # Will Need to be assigned new cell types and colors once we establish what we will show 
horizontal_bar <- ggplot(total_by_type, aes(x = total_count, y = MajorCellType)) +
  geom_bar(aes(fill = total_by_type$MajorCellType), stat = "identity") +
  theme_minimal() +scale_fill_manual(values = colors_mln, breaks = MajorCellTypes_mln) +
  labs(title = "Total Cell Type Counts", fill = "Major Cell Type", x = "Cell Number") +
  theme(axis.text.y = element_text(hjust = 1), legend.position = "none", 
        axis.title.y = element_blank())
horizontal_bar

#Figure 2C, Figure S2F
combined_plot <-(stacked_bar ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "MLN_filt_Total_fig2_Count.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

horizontal_bar
ggsave( "MLN_filt_fig2_Number.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 5,  height = 5,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

combined_plot <-(stacked_bar2 ) /( dotplot ) + plot_layout(heights = c(1, 2)) 
combined_plot 
ggsave( "MLN_filt_Total_fig2_Frequency.svg",  plot = last_plot() , device = NULL,  path = images,  scale = 1,  width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```


```{r}
# Finding Markers for Major cell types
library(presto)
Idents( Ileum_filt) <- "MajorCellType"
Ileum_markers <- FindAllMarkers(Ileum_filt, only.pos = TRUE, method = "wilcox")
Ileum_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 3) %>%
    ungroup() -> top3_ileum
DoHeatmap(subset(Ileum_filt, downsample = 100), features = top3_ileum$gene) + NoLegend()

Idents( MLN_filt) <- "MajorCellType"
MLN_markers <- FindAllMarkers(MLN_filt, only.pos = TRUE, method = "wilcox")
MLN_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 3) %>%
    ungroup() -> top3_MLN
DoHeatmap(subset(MLN_filt, downsample = 100), features = top3_MLN$gene) + NoLegend()

library(presto)
Idents( Ileum_filt) <- "MajorCellType"
prots <- unique(rownames(Ileum_filt@assays$ADT))
DoHeatmap(subset(Ileum_filt, downsample = 100), features =prots, assay = "ADT") + NoLegend()


#New Method to Find more Selective Markers for separation

avg_exp <- AverageExpression(
  Ileum_filt,
  group.by = "MajorCellType",
  assays = "RNA",
  slot = "data"
)$RNA  # matrix (genes x clusters)

Idents(Ileum_filt) <- "MajorCellType"

# Now compute pct_by_cluster
pct_by_cluster <- sapply(levels(Ileum_filt), function(clust) {
  cells_in_cluster <- WhichCells(Ileum_filt, idents = clust)
  expr_mat <- GetAssayData(Ileum_filt, assay = "RNA", slot = "data")[, cells_in_cluster]
  rowMeans(expr_mat > 0)  # fraction of cells expressing gene
})

pct_by_cluster <- as.matrix(pct_by_cluster)

#  Tau function to create a score for specificity
tau_score <- function(x) {
  if (max(x) == 0) return(NA)
  x <- x / max(x)
  tau <- sum(1 - x) / (length(x) - 1)
  return(tau)
}

tau_values <- apply(avg_exp, 1, tau_score)

tau_df <- data.frame(
  gene = rownames(avg_exp),
  tau = tau_values,
  max_cluster = apply(avg_exp, 1, function(x) colnames(avg_exp)[which.max(x)]),
  max_expr = apply(avg_exp, 1, max),
  stringsAsFactors = FALSE
)


tau_df$pct_in_max <- mapply(
  function(g, cl) pct_by_cluster[g, cl],
  tau_df$gene,
  tau_df$max_cluster
)

#  Filter to >50% expressed in max cluster
tau_df_filtered <- tau_df %>%
  filter(pct_in_max > 0.5)

# Top 10 per cluster
top_tau_genes <- tau_df_filtered %>%
  group_by(max_cluster) %>%
  slice_max(order_by = tau, n = 3) %>%
  arrange(max_cluster, desc(tau))

top_tau_genes




```


```{r}
#Advanced Heatmap Visualizations 
#Previously I used top3_ilem$gene and top3_ileum$cluster where now I am using top tau
library(ComplexHeatmap)
library(circlize)
avgexp_ileum = AggregateExpression(Ileum_filt, features = top_tau_genes$gene, group.by = "MajorCellType", assays = "RNA", return.seurat = T)
avgexp_ileum <- as.matrix(avgexp_ileum@assays$RNA$data )
desired_column_order <- c(
  "Enteric Nervous System", "Endothelial", "Mesenchymal", "Goblet Cells", 
  "Paneth Cells", "Enteroendocrine Cells", "Stem / Transit Amplifying Cells", 
  "Enterocytes", "Tuft Cells", "Innate Lymphoid Cells", "NKT Cells", "gd T Cells", 
  "CD4 ab T Cells", "CD8 ab T Cells", "B Cells", "Plasma Cells", "Plasmacytoid DCs", 
  "Mast Cells", "Neutrophils", "Dendritic Cells", "Macrophages", "Monocytes")

top_tau_genes$max_cluster <- factor(top_tau_genes$max_cluster, levels = desired_column_order)
top_tau_genes <- top_tau_genes[order(top_tau_genes$max_cluster), ]

desired_row_order <- top_tau_genes$gene

#Reorder the variables to your desired order
avgexp_ileum <- avgexp_ileum[, desired_column_order]
avgexp_ileum <- avgexp_ileum[desired_row_order, ]

cell_type_colors <- c(
  "gd T Cells" = "#CE864D", 
  "CD4 ab T Cells" = "#D07C3A",
  "CD8 ab T Cells" = "#D27227",
  "Innate Lymphoid Cells" = "#D36813",
  "NKT Cells" = "#D55E00",
  "B Cells" = "#CC79A7",
  "Plasma Cells" = "#F0E442",
  "Dendritic Cells" = "#8DD4C0",
  "Plasmacytoid DCs" = "#71C9B1",
  "Macrophages" = "#55BEA1",
  "Monocytes" = "#38B492",
  "Neutrophils" = "#1CA982",
  "Mast Cells" = "#009E73",
  "Enterocytes" = "#ABD7F1",
  "Goblet Cells" = "#9AD0EF",
  "Paneth Cells" = "#89C9EE",
  "Enteroendocrine Cells" = "#78C2EC",
  "Tuft Cells" = "#67BBEB",
  "Stem / Transit Amplifying Cells" = "#56B4E9",
  "Enteric Nervous System" = "#E60000",
  "Mesenchymal" = "#0072B2",
  "Endothelial" = "#014971"
)

annotation_col <- data.frame(cell_type = colnames(avgexp_ileum))
rownames(annotation_col) <- colnames(avgexp_ileum)

# Re-create the annotation object
ha <- HeatmapAnnotation(
  cell_type = annotation_col$cell_type,
  col = list(cell_type = cell_type_colors),
  show_annotation_name = FALSE,
  annotation_name_side = "left",
  show_legend = FALSE  #
)

ha_row <- rowAnnotation(
  cell_type = rownames(t(avgexp_ileum)),
  col = list(cell_type = cell_type_colors),
  show_annotation_name = FALSE, show_legend = FALSE
)
# Define the color scale for gene expression
color_range <- range(avgexp_ileum, na.rm = TRUE)
breaks <- seq(from = color_range[1], to = color_range[2], length.out = 9)
col_fun <- colorRamp2(
  breaks,
  brewer.pal(15, "Blues")
)

# Create the heatmap
ht <- Heatmap(
  matrix = t(avgexp_ileum),  
  name = "Expression",
  col = col_fun,
  left_annotation = ha_row,      
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  border = FALSE
)

draw(ht, heatmap_legend_side = "right")
#Figure S2D
# Save to SVG (properly structured)
svg(paste0(images, "/Ileum_MajorCellType_ComplexHeatmap_output_v2.svg"), width = 11, height = 5)
draw(ht, heatmap_legend_side = "right")
dev.off()





```

```{r}
#Adt  Ileum Heatmap

Ileum_filt[["ADT"]] <- JoinLayers(Ileum_filt[["ADT"]])
protMArkers <- FindAllMarkers(Ileum_filt, assay = "ADT", layer = "data", group.by = "MajorCellType", only.pos = T, 
                              min.pct = 0.25)
protMArkers$pct.diff <- protMArkers$pct.1-protMArkers$pct.2
protMArkers <-protMArkers[protMArkers$pct.diff >0.2,]

protMArkers <- protMArkers %>%
  group_by(cluster) %>%
  slice_max(order_by = pct.diff, n = 4) %>%
  arrange(cluster, desc(pct.diff))
avgexp_ileum = AverageExpression(Ileum_filt, features = protMArkers$gene, group.by = "MajorCellType", assays = "ADT", return.seurat = T)
avgexp_ileum <- as.matrix(avgexp_ileum@assays$ADT$data)
desired_column_order <- c(
  "Enteric Nervous System", "Endothelial", "Mesenchymal", "Goblet Cells", 
  "Paneth Cells", "Enteroendocrine Cells", "Stem / Transit Amplifying Cells", 
  "Enterocytes", "Tuft Cells", "Innate Lymphoid Cells", "NKT Cells", "gd T Cells", 
  "CD4 ab T Cells", "CD8 ab T Cells", "B Cells", "Plasma Cells", "Plasmacytoid DCs", 
  "Mast Cells", "Neutrophils", "Dendritic Cells", "Macrophages", "Monocytes")

annotation_col <- data.frame(cell_type = colnames(avgexp_ileum))
rownames(annotation_col) <- colnames(avgexp_ileum)

# Re-create the annotation object
ha <- HeatmapAnnotation(
  cell_type = annotation_col$cell_type,
  col = list(cell_type = cell_type_colors),
  show_annotation_name = FALSE,
  annotation_name_side = "left",
  show_legend = FALSE  #
)

ha_row <- rowAnnotation(
  cell_type = rownames(t(avgexp_ileum)),
  col = list(cell_type = cell_type_colors),
  show_annotation_name = FALSE, show_legend = FALSE
)
# Define the color scale for gene expression
color_range <- range(avgexp_ileum, na.rm = TRUE)
breaks <- seq(from = color_range[1], to = color_range[2], length.out = 9)
col_fun <- colorRamp2(
  breaks,
  brewer.pal(15, "Blues")
)

scaled_data <- scale(t(avgexp_ileum), center = TRUE, scale = TRUE)
library(ComplexHeatmap)
# Create the heatmap with scaled data
ht <- Heatmap(
  matrix = t(scaled_data),
  name = "Expression",
  col = colorRampPalette(brewer.pal(9, "Blues"))(100),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  border = FALSE
)

draw(ht, heatmap_legend_side = "right")

options(ComplexHeatmap.force_full_draw = TRUE)
options(ComplexHeatmap.raster_device = NULL)

svg(paste0(images, "/Ileum_MajorCellType_ComplexHeatmap_ADT_v2.svg"), width = 7, height = 8)
draw(ht, heatmap_legend_side = "right")
dev.off()
#Figure S2E
pdf(paste0(images, "/Ileum_MajorCellType_ComplexHeatmap_ADT_v2.pdf"), width = 7, height = 8, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right")
dev.off()


options(ComplexHeatmap.force_full_draw = FALSE)
```

```{r}
avgexp_MLN = AggregateExpression(MLN_filt, features = top3_MLN$gene, group.by = "MajorCellType", assays = "RNA", return.seurat = T)

avgexp_MLN <- as.matrix(avgexp_MLN@assays$RNA$data )
desired_column_order <- c(
   "Endothelial", "Mesenchymal",  "Innate Lymphoid Cells", "NKT Cells", "gd T Cells", 
  "CD4 ab T Cells", "CD8 ab T Cells", "B Cells", "Plasma Cells", "Plasmacytoid DCs", 
  "Basophils", "Neutrophils", "Dendritic Cells", "Macrophages", "Monocytes")

top3_MLN$cluster <- factor(top3_MLN$cluster, levels = desired_column_order)
top3_MLN <- top3_MLN[order(top3_MLN$cluster), ]

desired_row_order <- top3_MLN$gene


avgexp_MLN <- avgexp_MLN[, desired_column_order]
avgexp_MLN <- avgexp_MLN[desired_row_order, ]

cell_type_colors <- c(
  "gd T Cells" = "#CE864D", 
  "CD4 ab T Cells" = "#D07C3A",
  "CD8 ab T Cells" = "#D27227",
  "Innate Lymphoid Cells" = "#D36813",
  "NKT Cells" = "#D55E00",
  "B Cells" = "#CC79A7",
  "Plasma Cells" = "#F0E442",
  "Dendritic Cells" = "#8DD4C0",
  "Plasmacytoid DCs" = "#71C9B1",
  "Macrophages" = "#55BEA1",
  "Monocytes" = "#38B492",
  "Neutrophils" = "#1CA982",
  "Basophils" = "#009E73",
  "Mesenchymal" = "#0072B2",
  "Endothelial" = "#014971"
)

annotation_col <- data.frame(cell_type = colnames(avgexp_MLN))
rownames(annotation_col) <- colnames(avgexp_MLN)

# Re-create the annotation object
ha <- HeatmapAnnotation(
  cell_type = annotation_col$cell_type,
  col = list(cell_type = cell_type_colors),
  show_annotation_name = FALSE,
  annotation_name_side = "left",
  show_legend = FALSE )

# Define the color scale for gene expression
color_range <- range(avgexp_MLN, na.rm = TRUE)
breaks <- seq(from = color_range[1], to = color_range[2], length.out = 9)
col_fun <- colorRamp2(
  breaks,
  brewer.pal(9, "Blues")
)

# Create the heatmap
ht <- Heatmap(
  matrix = avgexp_MLN,
  name = "Expression",
  col = col_fun,
  top_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  border = FALSE
)
draw(ht, heatmap_legend_side = "right")
# Save to SVG 
svg(paste0(images, "/MLN_MajorCellType_ComplexHeatmap_output.svg"), width = 4, height = 8)
draw(ht, heatmap_legend_side = "right")
dev.off()
```

```{r}
# Using the improved method for the MLN 


#MLN 
avg_exp <- AverageExpression(
  MLN_filt,
  group.by = "MajorCellType",
  assays = "RNA",
  slot = "data"
)$RNA  # matrix (genes x clusters)

Idents(MLN_filt) <- "MajorCellType"

# Now compute pct_by_cluster
pct_by_cluster <- sapply(levels(MLN_filt), function(clust) {
  cells_in_cluster <- WhichCells(MLN_filt, idents = clust)
  expr_mat <- GetAssayData(MLN_filt, assay = "RNA", slot = "data")[, cells_in_cluster]
  rowMeans(expr_mat > 0)  # fraction of cells expressing gene
})

pct_by_cluster <- as.matrix(pct_by_cluster)

#  Tau function to create a score for specificity
tau_score <- function(x) {
  if (max(x) == 0) return(NA)
  x <- x / max(x)
  tau <- sum(1 - x) / (length(x) - 1)
  return(tau)
}

tau_values <- apply(avg_exp, 1, tau_score)

tau_df <- data.frame(
  gene = rownames(avg_exp),
  tau = tau_values,
  max_cluster = apply(avg_exp, 1, function(x) colnames(avg_exp)[which.max(x)]),
  max_expr = apply(avg_exp, 1, max),
  stringsAsFactors = FALSE
)


tau_df$pct_in_max <- mapply(
  function(g, cl) pct_by_cluster[g, cl],
  tau_df$gene,
  tau_df$max_cluster
)

#  Filter to >50% expressed in max cluster
tau_df_filtered <- tau_df %>%
  filter(pct_in_max > 0.5)

# Top 10 per cluster
top_tau_genes <- tau_df_filtered %>%
  group_by(max_cluster) %>%
  slice_max(order_by = tau, n = 3) %>%
  arrange(max_cluster, desc(tau))

top_tau_genes

avgexp_MLN = AggregateExpression(MLN_filt, features = top_tau_genes$gene, group.by = "MajorCellType", assays = "RNA", return.seurat = T)

avgexp_MLN <- as.matrix(avgexp_MLN@assays$RNA$data )
desired_column_order <- c(
   "Endothelial", "Mesenchymal",  "Innate Lymphoid Cells", "NKT Cells", "gd T Cells", 
  "CD4 ab T Cells", "CD8 ab T Cells", "B Cells", "Plasma Cells", "Plasmacytoid DCs", 
  "Basophils", "Neutrophils", "Dendritic Cells", "Macrophages", "Monocytes")

top_tau_genes$max_cluster <- factor(top_tau_genes$max_cluster, levels = desired_column_order)
top_tau_genes <- top_tau_genes[order(top_tau_genes$max_cluster), ]

desired_row_order <- top_tau_genes$gene


avgexp_MLN <- avgexp_MLN[, desired_column_order]
avgexp_MLN <- avgexp_MLN[desired_row_order, ]

cell_type_colors <- c(
  "gd T Cells" = "#CE864D", 
  "CD4 ab T Cells" = "#D07C3A",
  "CD8 ab T Cells" = "#D27227",
  "Innate Lymphoid Cells" = "#D36813",
  "NKT Cells" = "#D55E00",
  "B Cells" = "#CC79A7",
  "Plasma Cells" = "#F0E442",
  "Dendritic Cells" = "#8DD4C0",
  "Plasmacytoid DCs" = "#71C9B1",
  "Macrophages" = "#55BEA1",
  "Monocytes" = "#38B492",
  "Neutrophils" = "#1CA982",
  "Basophils" = "#009E73",
  "Mesenchymal" = "#0072B2",
  "Endothelial" = "#014971"
)

annotation_col <- data.frame(cell_type = colnames(avgexp_MLN))
rownames(annotation_col) <- colnames(avgexp_MLN)

# Re-create the annotation object
ha <- HeatmapAnnotation(
  cell_type = annotation_col$cell_type,
  col = list(cell_type = cell_type_colors),
  show_annotation_name = FALSE,
  annotation_name_side = "left",
  show_legend = FALSE )

# Define the color scale for gene expression
color_range <- range(avgexp_MLN, na.rm = TRUE)
breaks <- seq(from = color_range[1], to = color_range[2], length.out = 9)
col_fun <- colorRamp2(
  breaks,
  brewer.pal(9, "Blues")
)


ht <- Heatmap(
  matrix = avgexp_MLN,
  name = "Expression",
  col = col_fun,
  top_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  border = FALSE
)
draw(ht, heatmap_legend_side = "right")
# Save to SVG 
#Figure S2C
svg(paste0(images, "/MLN_MajorCellType_ComplexHeatmap_output_v2_oct2025.svg"), width = 4, height = 8)
draw(ht, heatmap_legend_side = "right")
dev.off()

```



#Transcripts per cell
```{r}

# Extract metadata into a data frame
library(forcats)
library(dplyr)

# Create combined dataset with proper sample names
combined_seurat <- as.data.frame(Ileum_filt@meta.data) %>%
  dplyr::select(Mouse, InfectionStatus, Tissue, nCount_RNA, nCount_ADT) %>%
  dplyr::rename(Sample = Mouse,
                Infection = InfectionStatus,
                nCount_RNA = nCount_RNA)

combined_seurat_mln <- as.data.frame(MLN_filt@meta.data) %>%
  dplyr::select(Mouse, InfectionStatus, Tissue, nCount_RNA, nCount_ADT) %>%
  dplyr::rename(Sample = Mouse,
                Infection = InfectionStatus,
                nCount_RNA = nCount_RNA)

combined_seurat <- rbind(combined_seurat, combined_seurat_mln)

# Factorize infection with custom order
infection_levels <- c("Naive", "Cryptosporidium", "Yersinia", "Candida",
                      "MNV", "SFB_YA", "Nippostrongylus")

combined_seurat$Infection <- factor(combined_seurat$Infection, levels = infection_levels)

# Calculate average nCount_RNA per sample
summary_data <- combined_seurat %>%
  group_by(Sample, Infection, Tissue) %>%
  summarise(mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
           .groups = 'drop') %>%
  dplyr::rename(Mean.Reads.per.Cell = mean_nCount_RNA)

# Tissue means for dashed red lines (using the summary data)
tissue_means <- summary_data %>%
  group_by(Tissue) %>%
  summarise(mean_reads = mean(Mean.Reads.per.Cell, na.rm = TRUE))

# Sample labels mapping
sample_labels <- setNames(as.character(summary_data$Infection), 
                        summary_data$Sample)
summary_data <- summary_data %>%
  arrange(Infection) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
# Plotting with averaged data
ReadNumber <- ggplot(summary_data,
                    aes(x = Sample, y = Mean.Reads.per.Cell, fill = Infection)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") +
  geom_hline(data = tissue_means, aes(yintercept = mean_reads),
             color = "red", linetype = "dashed") +
  scale_fill_manual(values = colors8, breaks = infection_levels) +
  facet_wrap(~Tissue, scales = "free_x") +
  scale_x_discrete(labels = sample_labels) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12),
    strip.background = element_rect(fill = "grey85", colour = "grey20"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ylab("Number of Reads")

ReadNumber
ggsave("nCount_RNA_BArchart_Supplement1.svg", plot = last_plot(), path = images, width = 11, height = 4)


# Calculate average nCount_RNA per sample
summary_data <- combined_seurat %>%
  group_by(Sample, Infection, Tissue) %>%
  summarise(mean_nCount_ADT = mean(nCount_ADT, na.rm = TRUE),
           .groups = 'drop') %>%
  dplyr::rename(Mean.Reads.per.Cell = mean_nCount_ADT)

# Tissue means for dashed red lines (using the summary data)
tissue_means <- summary_data %>%
  group_by(Tissue) %>%
  summarise(mean_reads = mean(Mean.Reads.per.Cell, na.rm = TRUE))

# Sample labels mapping
sample_labels <- setNames(as.character(summary_data$Infection), 
                        summary_data$Sample)
summary_data <- summary_data %>%
  arrange(Infection) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
# Plotting with averaged data
ReadNumber <- ggplot(summary_data,
                    aes(x = Sample, y = Mean.Reads.per.Cell, fill = Infection)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") +
  geom_hline(data = tissue_means, aes(yintercept = mean_reads),
             color = "red", linetype = "dashed") +
  scale_fill_manual(values = colors8, breaks = infection_levels) +
  facet_wrap(~Tissue, scales = "free_x") +
  scale_x_discrete(labels = sample_labels) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12),
    strip.background = element_rect(fill = "grey85", colour = "grey20"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ylab("Number of Reads")

ReadNumber
#Figure S1D
ggsave("nCount_ADT_BArchart_Supplement1.svg", plot = last_plot(), path = images, width = 11, height = 4)
```


#Figure Work - Graphing Subclustered Objects 
# Figure 2E and Figure S3A are made from the following sections  in which the subset objects were plotted
```{r}
#Load the Subclustered Objects 
Ileum_T_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Tcell_intermediate.qs2")
Ileum_Stromal_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Stromal_intermediate.qs2")
Ileum_Bcell_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Bcell_intermediate.qs2")
Ileum_Myeloid_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Myeloid_intermediate.qs2")
Ileum_Epithelial_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Epithelial_intermediate.qs2")
Ileum_ILCs_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_ILC_intermediate.qs2")
Ileum_CD4Tcell_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_CD4Tcell_intermediate.qs2")
Ileum_CD8Tcell_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_CD8Tcell_intermediate.qs2")
Ileum_Nervous_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Nervous_intermediate.qs2")
Ileum_Plasma_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Plasma_intermediate.qs2")
```


```{r}
cbf_18_colors <- c(
  "#1170AA",  
  "#FCB13F",  
  "#60BFC1",  
  "#EF6F6A",  
    "#937860", 
  "#D17A00",  
  "#B78CBA",  
  "#B3B3B3",  
    "#64B200", 
  "#D45E00",  
  "#7E8CCF",  
  "#E6A0C4",  
  "#568B3F",  
  "#C44E52",
   "#5FA2A3", 
  "#CCB974", 

  "#D0A6BE", 
  "#4E84C4"  
  
)

```

```{r}
#transfer Finest Cell Labels from ILC/CD4/CD8 groups to total T cells and examine
Ileum_T_filt$FinestCellType  <- "NA"
Ileum_T_filt <- AddMetaData(object = Ileum_T_filt, metadata = Ileum_CD4Tcell_filt$FinestCellType , col.name = "FinestCellType")
Ileum_T_filt <- AddMetaData(object = Ileum_T_filt, metadata = Ileum_CD8Tcell_filt$FinestCellType , col.name = "FinestCellType")
Ileum_T_filt <- AddMetaData(object = Ileum_T_filt, metadata = Ileum_ILCs_filt$FinestCellType , col.name = "FinestCellType")
DimPlot(Ileum_T_filt, reduction = "wnn.T.umap", group.by = "FinestCellType", label = T, repel = T, combine = T) + NoLegend()
```

```{r}
#Total T cell map
#Visualizing the UMAPS
umap_data <- Ileum_T_filt[["wnn.T.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_T_filt$FinestCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnntumap_1 = mean(wnntumap_1), wnntumap_2 = mean(wnntumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnntumap_1, y = wnntumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
  ylim(-15, 15) + xlim(-13, 10) +
  scale_color_viridis_d(option = "C") + 

  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_Tcell_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_Tcell_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_Tcell_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_T_filt[["wnn.T.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_T_filt$FinestCellType)
umap_data$condition <- Ileum_T_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_Tcell_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```


```{r}
#Cd4 t cells
umap_data <- Ileum_CD4Tcell_filt@reductions[["wnn.CD4Tcell.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_CD4Tcell_filt$FinestCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnncd4tcellumap_1 = mean(wnncd4tcellumap_1), wnncd4tcellumap_2 = mean(wnncd4tcellumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnncd4tcellumap_1, y = wnncd4tcellumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
    scale_color_manual(values = cbf_18_colors) + 
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_CD4Tcell_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_CD4Tcell_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_CD4Tcell_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_CD4Tcell_filt[["wnn.CD4Tcell.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_CD4Tcell_filt$FinestCellType)
umap_data$condition <- Ileum_CD4Tcell_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_CD4Tcell_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)


```


```{r}
#CD8 t cells
umap_data <- Ileum_CD8Tcell_filt@reductions[["wnn.CD8Tcell.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_CD8Tcell_filt$FinestCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnncd8tcellumap_1 = mean(wnncd8tcellumap_1), wnncd8tcellumap_2 = mean(wnncd8tcellumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnncd8tcellumap_1, y = wnncd8tcellumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
    scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_CD8Tcell_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_CD8Tcell_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_CD8Tcell_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_CD8Tcell_filt[["wnn.CD8Tcell.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_CD8Tcell_filt$FinestCellType)
umap_data$condition <- Ileum_CD8Tcell_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_CD8Tcell_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```


```{r}
umap_data <- Ileum_ILCs_filt[["wnn.ILC.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_ILCs_filt$FinestCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnilcumap_1 = mean(wnnilcumap_1), wnnilcumap_2 = mean(wnnilcumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnilcumap_1, y = wnnilcumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) + 
    scale_color_manual(values = cbf_18_colors) +  
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_ILCs_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_ILCs_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_ILCs_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_ILCs_filt[["wnn.ILC.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_ILCs_filt$FinestCellType)
umap_data$condition <- Ileum_ILCs_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_ILCs_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```


```{r}

umap_data <- Ileum_Bcell_filt[["wnn.Bcell.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_Bcell_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnbcellumap_1 = mean(wnnbcellumap_1), wnnbcellumap_2 = mean(wnnbcellumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnbcellumap_1, y = wnnbcellumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
    scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_Bcell_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_Bcell_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_Bcell_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_Bcell_filt[["wnn.Bcell.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_Bcell_filt$FineCellType)
umap_data$condition <- Ileum_Bcell_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_Bcell_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```

```{r}
umap_data <- Ileum_Myeloid_filt[["wnn.Myeloid.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_Myeloid_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnmyeloidumap_1 = mean(wnnmyeloidumap_1), wnnmyeloidumap_2 = mean(wnnmyeloidumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnmyeloidumap_1, y = wnnmyeloidumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
    scale_color_manual(values = cbf_18_colors) + 
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_Myeloid_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_Myeloid_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_Myeloid_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_Myeloid_filt[["wnn.Myeloid.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_Myeloid_filt$FineCellType)
umap_data$condition <- Ileum_Myeloid_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_Myeloid_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```


```{r}
umap_data <- Ileum_Epithelial_filt[["wnn.Epithelial.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_Epithelial_filt$FinestCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnepithelialumap_1 = mean(wnnepithelialumap_1), wnnepithelialumap_2 = mean(wnnepithelialumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnepithelialumap_1, y = wnnepithelialumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
   scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_Epithelial_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_Epithelial_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_Epithelial_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_Epithelial_filt[["wnn.Epithelial.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_Epithelial_filt$FinestCellType)
umap_data$condition <- Ileum_Epithelial_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_Epithelial_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```


```{r}
umap_data <- Ileum_Stromal_filt[["wnn.Stromal.umap"]]@cell.embeddings
cell_type_mapping <- Ileum_Stromal_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnstromalumap_1 = mean(wnnstromalumap_1), wnnstromalumap_2 = mean(wnnstromalumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnstromalumap_1, y = wnnstromalumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
    scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "Ileum_Stromal_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "Ileum_Stromal_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "Ileum_Stromal_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

umap_data <- as.data.frame(Ileum_Stromal_filt[["wnn.Stromal.umap"]]@cell.embeddings)
umap_data$cluster <- factor(Ileum_Stromal_filt$FineCellType)
umap_data$condition <- Ileum_Stromal_filt$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

cluster_condition_percent <- cluster_condition_percent[!cluster_condition_percent$cluster =="NA",]
 cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c(   "Yersinia", "SFB_YA", "Nippostrongylus",  "MNV", "Cryptosporidium", "Candida", "Naive" ))
# Create a stacked barplot of percent frequency of each cluster within each condition
ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Condition", y = "Percent Frequency", title = "Percent Frequency of Clusters Across Conditions") +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#117733", "#A27DF1", "#DDCC77",  "#AA4499",  "#882255")) +
  theme_minimal(base_size = 15) +  
  theme( axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    legend.position = "right"  ) + coord_flip()
ggsave( "Ileum_Stromal_Frequency_barchart.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```


#Figure Work Graphing Subclustered MLN
#Figure S3A
```{r}
#Load the Subclustered Objects 
MLN_Bcell_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/MLN_Bcell_intermediate.qs2")
MLN_Stromal_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/MLN_Stromal_intermediate.qs2")
MLN_Myeloid_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/MLN_Myeloid_intermediate.qs2")
MLN_CD4Tcell_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/MLN_CD4Tcell_intermediate.qs2")
MLN_CD8Tcell_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/MLN_CD8Tcell_intermediate.qs2")
MLN_T_filt <- qs_read("/path/to/analysis/directory/Seurat_Files/MLN_Tcell_intermediate.qs2")
```


```{r}
cbf_18_colors <- c(
  "#1170AA",  
  "#FCB13F",  
  "#60BFC1",  
  "#EF6F6A",  
    "#937860", 
  "#D17A00",  
  "#B78CBA",  
  "#B3B3B3",  
    "#64B200", 
  "#D45E00",  
  "#7E8CCF",  
  "#E6A0C4",  
  "#568B3F",  
  "#C44E52",
   "#5FA2A3", 
  "#CCB974", 

  "#D0A6BE", 
  "#4E84C4"  
  
)
ncol(MLN_Bcell_filt)
unique(MLN_CD4Tcell_filt$FineCellType)
unique(MLN_CD8Tcell_filt$FineCellType)
unique(MLN_Myeloid_filt$FineCellType)
unique(MLN_Stromal_filt$FineCellType)
unique(MLN_Bcell_filt$FineCellType)
```


```{r}
#Total T cell map
#Visualizing the UMAPS
umap_data <- MLN_CD4Tcell_filt[["wnn.CD4Tcell.umap"]]@cell.embeddings
cell_type_mapping <- MLN_CD4Tcell_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnncd4tcellumap_1 = mean(wnncd4tcellumap_1), wnncd4tcellumap_2 = mean(wnncd4tcellumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnncd4tcellumap_1, y = wnncd4tcellumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) + 
  ylim(-15, 15) + xlim(-13, 10) +
  scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "MLN_CD4Tcell_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "MLN_CD4Tcell_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "MLN_CD4Tcell_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```
```{r}
#CD8 T cells MLN
umap_data <- MLN_CD8Tcell_filt[["wnn.CD8Tcell.umap"]]@cell.embeddings
cell_type_mapping <- MLN_CD8Tcell_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnncd8tcellumap_1 = mean(wnncd8tcellumap_1), wnncd8tcellumap_2 = mean(wnncd8tcellumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnncd8tcellumap_1, y = wnncd8tcellumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) + 
  ylim(-15, 15) + xlim(-13, 10) +
   scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "MLN_CD8Tcell_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "MLN_CD8Tcell_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "MLN_CD8Tcell_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```
```{r}
# MLN Stromal
umap_data <- MLN_Stromal_filt[["wnn.Stromal.umap"]]@cell.embeddings
cell_type_mapping <- MLN_Stromal_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnstromalumap_1 = mean(wnnstromalumap_1), wnnstromalumap_2 = mean(wnnstromalumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnstromalumap_1, y = wnnstromalumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) + 
  ylim(-15, 15) + xlim(-13, 10) +
   scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "MLN_Stromal_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "MLN_Stromal_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "MLN_Stromal_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```
```{r}
#MLn B cell 
umap_data <- MLN_Bcell_filt[["wnn.Bcell.umap"]]@cell.embeddings
cell_type_mapping <- MLN_Bcell_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnbcellumap_1 = mean(wnnbcellumap_1), wnnbcellumap_2 = mean(wnnbcellumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnbcellumap_1, y = wnnbcellumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
  ylim(-15, 15) + xlim(-13, 10) +
   scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "MLN_Bcell_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "MLN_Bcell_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "MLN_Bcell_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```
```{r}
#MLn B cell 
umap_data <- MLN_Myeloid_filt[["wnn.Myeloid.umap"]]@cell.embeddings
cell_type_mapping <- MLN_Myeloid_filt$FineCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnmyeloidumap_1 = mean(wnnmyeloidumap_1), wnnmyeloidumap_2 = mean(wnnmyeloidumap_2))

library(ggrepel)
plot <- ggplot(umap_data, aes(x = wnnmyeloidumap_1, y = wnnmyeloidumap_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  
  ylim(-15, 15) + xlim(-13, 10) +
   scale_color_manual(values = cbf_18_colors) +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 
plot 
ggsave( "MLN_Myeloid_subclusters_UMAP_NoLabels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0)
ggsave( "MLN_Myeloid_subclusters_UMAP_Labels.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 9,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
plot +  geom_point(alpha = 1) + geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
                   min.segment.length = 0) + theme(legend.position = "right")
ggsave( "MLN_Myeloid_subclusters_legend.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 11,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```

# Milo based neighborhood calculation and differential abundance 
```{r}
#additional Libraries
library(miloR)
library(scater)
library(patchwork)
library(igraph)
library(ggbeeswarm)
```

```{r}
#This component was performed in a separate script
#Convert Each Object to sce and perform the Milo workflow - make a nn graph
Ileum_CD4Tcell_filt[["RNA"]] <- JoinLayers(Ileum_CD4Tcell_filt[["RNA"]])
umap_embed <- Embeddings(Ileum_CD4Tcell_filt, reduction = "wnn.CD4Tcell.umap")
rna_embed <- Embeddings(Ileum_CD4Tcell_filt, "integrated.harmony.CD4Tcell")[, 1:25]
adt_embed <- Embeddings(Ileum_CD4Tcell_filt, "integrated.ADT.harmony.CD4Tcell")[, 1:20]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(Ileum_CD4Tcell_filt, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 45, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 45, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/CD4Tcell_milo.rds")




Ileum_ILC_filt[["RNA"]] <- JoinLayers(Ileum_ILC_filt[["RNA"]])
umap_embed <- Embeddings(Ileum_ILC_filt, reduction = "wnn.ILC.umap")
rna_embed <- Embeddings(Ileum_ILC_filt, "integrated.harmony.ILC")[, 1:25]
adt_embed <- Embeddings(Ileum_ILC_filt, "integrated.ADT.harmony.ILC")[, 1:20]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(Ileum_ILC_filt, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 45, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 45, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/ILC_milo.rds")


Ileum_Bcell_filt[["RNA"]] <- JoinLayers(Ileum_Bcell_filt[["RNA"]])
umap_embed <- Embeddings(Ileum_Bcell_filt, reduction = "wnn.Bcell.umap")
rna_embed <- Embeddings(Ileum_Bcell_filt, "integrated.harmony.Bcell")[, 1:25]
adt_embed <- Embeddings(Ileum_Bcell_filt, "integrated.ADT.harmony.Bcell")[, 1:20]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(Ileum_Bcell_filt, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 45, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 45, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/Bcell_milo.rds")


Ileum_Stromal_filt[["RNA"]] <- JoinLayers(Ileum_Stromal_filt[["RNA"]])
umap_embed <- Embeddings(Ileum_Stromal_filt, reduction = "wnn.Stromal.umap")
rna_embed <- Embeddings(Ileum_Stromal_filt, "integrated.harmony.Stromal")[, 1:25]
adt_embed <- Embeddings(Ileum_Stromal_filt, "integrated.ADT.harmony.Stromal")[, 1:20]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(Ileum_Stromal_filt, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 45, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 45, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/Stromal_milo.rds")



Ileum_Epithelial_filt[["RNA"]] <- JoinLayers(Ileum_Epithelial_filt[["RNA"]])
umap_embed <- Embeddings(Ileum_Epithelial_filt, reduction = "wnn.Epithelial.umap")
rna_embed <- Embeddings(Ileum_Epithelial_filt, "integrated.harmony.Epithelial")[, 1:25]
adt_embed <- Embeddings(Ileum_Epithelial_filt, "integrated.ADT.harmony.Epithelial")[, 1:20]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(Ileum_Epithelial_filt, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 45, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 45, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/Epithelial_milo.rds")



Ileum_Myeloid_filt[["RNA"]] <- JoinLayers(Ileum_Myeloid_filt[["RNA"]])
umap_embed <- Embeddings(Ileum_Myeloid_filt, reduction = "wnn.Myeloid.umap")
rna_embed <- Embeddings(Ileum_Myeloid_filt, "integrated.harmony.Myeloid")[, 1:25]
adt_embed <- Embeddings(Ileum_Myeloid_filt, "integrated.ADT.harmony.Myeloid")[, 1:20]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(Ileum_Myeloid_filt, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 45, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 45, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/Myeloid_milo.rds")


Ileum_CD8Tcell_filt[["RNA"]] <- JoinLayers(Ileum_CD8Tcell_filt[["RNA"]])
umap_embed <- Embeddings(Ileum_CD8Tcell_filt, reduction = "wnn.CD8Tcell.umap")
rna_embed <- Embeddings(Ileum_CD8Tcell_filt, "integrated.harmony.CD8Tcell")[, 1:25]
adt_embed <- Embeddings(Ileum_CD8Tcell_filt, "integrated.ADT.harmony.CD8Tcell")[, 1:15]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(Ileum_CD8Tcell_filt, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 40, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 40, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/CD8Tcell_milo.rds")



#Plasma cell 

plasma_ileum <- qs2::qs_read("/path/to/analysis/directory/Seurat_Files/Ileum_Plasma_intermediate.qs2")
plasma_ileum[["RNA"]] <- JoinLayers(plasma_ileum[["RNA"]])
umap_embed <- Embeddings(plasma_ileum, reduction = "wnn.Plasma.umap")
rna_embed <- Embeddings(plasma_ileum, "integrated.harmony.Plasma")[, 1:25]
adt_embed <- Embeddings(plasma_ileum, "integrated.ADT.harmony.Plasma")[, 1:15]
multi_modal_embed <- cbind(rna_embed, adt_embed)
sce <- as.SingleCellExperiment(plasma_ileum, assay = "RNA")
reducedDims(sce)$UMAP <- umap_embed
reducedDims(sce)$multi_modal <- multi_modal_embed
milo <- Milo(sce)
print(Sys.time())
milo <- buildGraph(milo, k = 30, d = 40, reduced.dim = "multi_modal")
print(Sys.time())
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 40, reduced_dims = "multi_modal",  refined = TRUE)
milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
saveRDS(milo, file = "/path/to/analysis/directory/Seurat_Files/Plasma_milo.rds")
DimPlot(plasma_ileum, reduction = "wnn.Plasma.umap", group.by = "FineCellType")
FeaturePlot(plasma_ileum, "Il10", reduction = "wnn.Plasma.umap")
```

```{r}
#Loop through all the milo results and plot the initial plots for visualization
milo_files <- list(
  "CD4Tcell" = "/path/to/analysis/directory/Seurat_Files/CD4Tcell_milo.rds",
  "CD8Tcell" = "/path/to/analysis/directory/Seurat_Files/CD8Tcell_milo.rds",
  "ILC" = "/path/to/analysis/directory/Seurat_Files/ILC_milo.rds",
  "Epithelial" = "/path/to/analysis/directory/Seurat_Files/Epithelial_milo.rds",
  "Myeloid" = "/path/to/analysis/directory/Seurat_Files/Myeloid_milo.rds",
  "Stromal" = "/path/to/analysis/directory/Seurat_Files/Stromal_milo.rds",
  "Bcell" = "/path/to/analysis/directory/Seurat_Files/Bcell_milo.rds"
)

infections <- c("Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia")

sample_metadata <- Ileum_filt@meta.data


run_milo_DA <- function(file, label) {
  message("Processing: ", label)
  
  milo <- readRDS(file)
  milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")

  samples_in_milo <- colnames(nhoodCounts(milo))
  
  design.df <- sample_metadata %>%
    filter(Mouse %in% samples_in_milo) %>%
    distinct(Mouse, .keep_all = TRUE) %>%
    arrange(factor(Mouse, levels = samples_in_milo))
  rownames(design.df) <- design.df$Mouse
  
  design.df$InfectionStatus <- factor(
    design.df$InfectionStatus,
    levels = c("Naive", infections)
  )

  model <- model.matrix(~ 0 + InfectionStatus, data = design.df)
  colnames(model) <- gsub("InfectionStatus", "", colnames(model))
  contrasts <- paste0(infections, " - Naive")
  contrast.matrix <- makeContrasts(contrasts = contrasts, levels = model)

  contrast_results <- testNhoods(
    milo,
    design = model,
    design.df = design.df,
    model.contrasts = contrasts,
    reduced.dim = "multi_modal",
    fdr.weighting = "graph-overlap",
    norm.method = "TMM"
  )

  # Annotate with cell types
  if ("FinestCellType" %in% colnames(colData(milo))) {
    contrast_results <- annotateNhoods(milo, contrast_results, coldata_col = "FinestCellType")
    group_col <- "FinestCellType"
  } else {
    contrast_results <- annotateNhoods(milo, contrast_results, coldata_col = "FineCellType")
    group_col <- "FineCellType"
  }

 
  outdir <- file.path("/path/to/analysis/directory/Images/Milo_Plots", label)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  for (inf in infections) {
    logFC_col <- paste0("logFC.", inf, "...Naive")
    if (!logFC_col %in% colnames(contrast_results)) next
    
    res <- contrast_results[, c("Nhood", group_col, "SpatialFDR", logFC_col)]
    colnames(res)[colnames(res) == logFC_col] <- "logFC"
    
    p <- plotDAbeeswarm(res, group.by = group_col, alpha = 0.1) +
      theme_minimal(base_size = 13) +
      labs(title = paste(label, "- Neighborhood DA:", inf, "vs Naive")) +
      theme(
        axis.title.y = element_blank(),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "bottom"
      )
    
    ggsave(
      filename = file.path(outdir, paste0(label, "_", inf, "_DAbeeswarm.png")),
      plot = p, width = 8, height = 6, dpi = 300
    )
  }
  
  return(contrast_results)
}



# List of infections to loop over
infections <- c("Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia")

# Loop through each Milo result object in results_list
for (label in names(results_list)) {
  contrast_results <- results_list[[label]]
  
  if (is.null(contrast_results)) next  # skip if analysis failed
  
  # Pick the right grouping column
  group_col <- if ("FinestCellType" %in% colnames(contrast_results)) "FinestCellType" else "FineCellType"
  
  # Generate one beeswarm per infection
  plots <- lapply(infections, function(inf) {
    logFC_col <- paste0("logFC.", inf, "...Naive")
    if (!logFC_col %in% colnames(contrast_results)) return(NULL)
    
    res <- contrast_results[, c("Nhood", group_col, "SpatialFDR", logFC_col)]
    colnames(res)[colnames(res) == logFC_col] <- "logFC"
    
    p <- plotDAbeeswarm(res, group.by = group_col, alpha = 0.1) +
      theme_minimal(base_size = 12) +
      labs(title = paste(inf, "vs Naive")) +
      theme(
        axis.title.y = element_blank(),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "bottom"
      )
    return(p)
  })
  
  # Remove any NULL plots (missing infections)
  plots <- Filter(Negate(is.null), plots)
  
  # Combine all infection plots for this cell type
  combined_plot <- plot_grid(plotlist = plots, ncol = 3, align = "hv")
  
  # Save combined image
  outdir <- file.path("/path/to/analysis/directory/Images/Milo_Plots", label)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = file.path(outdir, paste0(label, "_Combined_Beeswarm.png")),
    plot = combined_plot, width = 15, height = 8, dpi = 300
  )
  
  message("✅ Saved combined plot for ", label)
}


# Remove NULLs and arrange all plots together
plots <- Filter(Negate(is.null), plots)
combined_plot <- plot_grid(plotlist = plots, ncol = 3, align = "hv")

print(combined_plot)
ggsave("/path/to/analysis/directory/Images/Milo_Plots/Combined_Beeswarm_AllInfections.png",
       combined_plot, width = 15, height = 8, dpi = 300)

```

```{r}
#Figure 3A and Figure 3B can be generated below. For 3A, the Epithelial milo object must be used. For 3B the Myeloid Milo object must be loaded
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/CD4Tcell_milo.rds")
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/ILC_milo.rds")
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/Epithelial_milo.rds")
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/Myeloid_milo.rds")
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/Stromal_milo.rds")
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/Bcell_milo.rds")
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/CD8Tcell_milo.rds")
plotNhoodSizeHist(milo)

```

```{r}

milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
#Get correct sample names from milo
samples_in_milo <- colnames(nhoodCounts(milo))
unique(samples_in_milo)
#  Filter and order sample metadata
sample_metadata<- Ileum_filt@meta.data

design.df <- sample_metadata %>%
  filter(Mouse %in% samples_in_milo) %>%
  distinct(Mouse, .keep_all = TRUE) %>%
  arrange(factor(Mouse, levels = samples_in_milo))

#  Set rownames
rownames(design.df) <- design.df$Mouse
design.df <- design.df[samples_in_milo, , drop = FALSE]  # enforce order

# Ensure InfectionStatus is a factor with all levels
design.df$InfectionStatus <- factor(
  design.df$InfectionStatus,
  levels = c("Naive", "Candida", "Cryptosporidium", "MNV", 
             "Nippostrongylus", "SFB_YA", "Yersinia")
)

# Ensure Mouse is a factor and matches the order of sample names from Milo
design.df <- design.df[match(samples_in_milo, design.df$Mouse), , drop = FALSE]
rownames(design.df) <- design.df$Mouse

# Build model matrix (now should be non-zero)
model <- model.matrix(~ 0 + InfectionStatus, data = design.df)
colnames(model) <- gsub("InfectionStatus", "", colnames(model))

# Define contrasts
library(limma)
infections <- c("Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia")
contrasts <- paste0(infections, " - Naive")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = model)

#  Run testNhoods
contrast_results <- testNhoods(
  milo,
  design = model,
  design.df = design.df,
  model.contrasts = contrasts,
  reduced.dim = "multi_modal",
  fdr.weighting = "graph-overlap",
  norm.method = "TMM"
)

contrast_results <- annotateNhoods(milo, contrast_results, coldata_col = "FineCellType")
#contrast_results <- annotateNhoods(milo, contrast_results, coldata_col = "FinestCellType")
```


```{r}
# Rename for compatibility with plotDAbeeswarm
candida_res <- contrast_results[, c("Nhood", "FinestCellType", "SpatialFDR", "logFC.Candida...Naive")]
colnames(candida_res)[colnames(candida_res) == "logFC.Candida...Naive"] <- "logFC"
plotDAbeeswarm(candida_res, group.by = "FinestCellType")

Yersinia_res <- contrast_results[, c("Nhood", "FineCellType", "SpatialFDR", "logFC.Yersinia...Naive")]
colnames(Yersinia_res)[colnames(Yersinia_res) == "logFC.Yersinia...Naive"] <- "logFC"
plotDAbeeswarm(Yersinia_res, group.by = "FineCellType")

MNV_res <- contrast_results[, c("Nhood", "FinestCellType", "SpatialFDR", "logFC.MNV...Naive")]
colnames(MNV_res)[colnames(MNV_res) == "logFC.MNV...Naive"] <- "logFC"
plotDAbeeswarm(MNV_res, group.by = "FinestCellType")

Cryptosporidium_res <- contrast_results[, c("Nhood", "FinestCellType", "SpatialFDR", "logFC.Cryptosporidium...Naive")]
colnames(Cryptosporidium_res)[colnames(Cryptosporidium_res) == "logFC.Cryptosporidium...Naive"] <- "logFC"
plotDAbeeswarm(Cryptosporidium_res, group.by = "FinestCellType")

Cryptosporidium_res <- contrast_results[, c("Nhood", "FineCellType", "SpatialFDR", "logFC.Cryptosporidium...Naive")]
colnames(Cryptosporidium_res)[colnames(Cryptosporidium_res) == "logFC.Cryptosporidium...Naive"] <- "logFC"
plotDAbeeswarm(Cryptosporidium_res, group.by = "FineCellType")

SFB_YA_res <- contrast_results[, c("Nhood", "FinestCellType", "SpatialFDR", "logFC.SFB_YA...Naive")]
colnames(SFB_YA_res)[colnames(SFB_YA_res) == "logFC.SFB_YA...Naive"] <- "logFC"
plotDAbeeswarm(SFB_YA_res, group.by = "FinestCellType")

Nippostrongylus_res <- contrast_results[, c("Nhood", "FinestCellType", "SpatialFDR", "logFC.Nippostrongylus...Naive")]
colnames(Nippostrongylus_res)[colnames(Nippostrongylus_res) == "logFC.Nippostrongylus...Naive"] <- "logFC"
plotDAbeeswarm(Nippostrongylus_res, group.by = "FinestCellType", alpha = 0.1) + theme_minimal(base_size = 13) + 
    labs(title = "Neighborhood DA: Nippostrongylus vs Naive") +theme() + 
  theme(axis.title.y = element_blank(), 
        text = element_text(color = "black"), 
        axis.text = element_text(color = "black"), 
        legend.position = "bottom")

```



```{r}
library(ggbeeswarm)
#Run for Epithelial

Nippostrongylus_res$FDR_signFC <- ifelse(Nippostrongylus_res$logFC < 0, -(1-Nippostrongylus_res$SpatialFDR), 1-(Nippostrongylus_res$SpatialFDR))
Nippostrongylus_res$FDR_FC <- ifelse(Nippostrongylus_res$logFC < 0, (1-Nippostrongylus_res$SpatialFDR), 1-(Nippostrongylus_res$SpatialFDR))
custom_colors <- c(
  "#1C75BC",  # dark blue at -1
  "gray80",   # gray at -0.9 (start of middle flat zone)
  "gray80",   # gray at 0.9 (end of middle flat zone)
  "#F7941D"   # dark red at 1
)


color_breaks <- scales::rescale(c(-1, -0.5,0.5, 1)) #only give color if FDR < 0.2 (in transformed 1-FDR)

ggplot(Nippostrongylus_res, aes(x = logFC, y = FinestCellType)) +
  geom_quasirandom(aes(color = FDR_signFC, alpha = FDR_FC), groupOnX = TRUE, size = 1.5) +
  scale_color_gradientn(
    colors = custom_colors,
    values = color_breaks,
    limits = c(-1, 1),
    oob = scales::squish,
    guide = "colorbar",
    name = "1 - FDR"
  ) +
 # scale_alpha_continuous(name = "1 - SpatialFDR", range = c(0.2, 1)) +
  theme_minimal(base_size = 13) +
  labs( title = "Neighborhood DA: Nippo", x = "logFC") +
  theme(
    axis.title.y = element_blank(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )

ggsave( "neighborhoodDA_Nippo_Epithelial.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)

```

```{r}
#Figure 3 demonstrating that Nippo has expansion fo Goblet cells and Tuft cells
library(cowplot)
library(scales)

# Filter to only Tuft and Goblet cells
unique(Ileum_filt$FineCellType)
cell_subset <- c("Tuft Cells", "Goblet Cells")

infections <- c("Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia")
custom_colors <- c("#1C75BC", "gray80", "gray80", "#F7941D")
color_breaks <- rescale(c(-1, -0.5, 0.5, 1))  # symmetric gradient around 0

# Determine global x-axis limits across all logFC cols
logfc_cols <- grep("^logFC\\.", colnames(contrast_results), value = TRUE)
global_range <- range(contrast_results[, logfc_cols], na.rm = TRUE)

plots <- lapply(seq_along(infections), function(i) {
  inf <- infections[i]
  logFC_col <- paste0("logFC.", inf, "...Naive")
  if (!logFC_col %in% colnames(contrast_results)) return(NULL)
  
  res <- contrast_results[, c("Nhood", "FineCellType", "SpatialFDR", logFC_col)]
  colnames(res)[colnames(res) == logFC_col] <- "logFC"

  res$FDR_signFC <- ifelse(res$logFC < 0, -(1 - res$SpatialFDR), 1 - res$SpatialFDR)
  res$FDR_FC <- 1 - res$SpatialFDR
  
  res <- res[res$FineCellType %in% cell_subset, ]
  if (nrow(res) == 0) return(NULL)

  color_guide <- if (i == length(infections)) "colorbar" else "none"
  
  ggplot(res, aes(x = logFC, y = FineCellType)) +
    geom_quasirandom(aes(color = FDR_signFC, alpha = 0.5),
                     groupOnX = TRUE, size = 1.25) +
    scale_color_gradientn(
      colors = custom_colors,
      values = color_breaks,
      limits = c(-1, 1),
      oob = scales::squish,
      guide = color_guide,
      name = "1 - FDR"
    ) +
    scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
    scale_x_continuous(limits = global_range) +     # ★★ FIXED X-AXIS RANGE ★★
    theme_minimal(base_size = 13) +
    labs(title = paste0(inf, " vs Naive"), x = "logFC") +
    theme(
      axis.title.y = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (i == length(infections)) "right" else "none",
      axis.text.y = if (i == 1) element_text() else element_blank()
    )
})


# Combine all infection plots horizontally
combined_plot <- plot_grid(plotlist = plots, ncol = 6, align = "hv")
combined_plot
# Save the combined plot
#Figure 3A
ggsave(
  filename = file.path(images, "Tuft_Goblet_AllInfections_Beeswarm.svg"),
  plot = combined_plot,
  width = 22, height = 5, dpi = 600, bg = "white"
)

```

```{r}
library(ggbeeswarm)
#Run for Epithelial
Cryptosporidium_res$FDR_signFC <- ifelse(Cryptosporidium_res$logFC < 0, -(1-Cryptosporidium_res$SpatialFDR), 1-(Cryptosporidium_res$SpatialFDR))
Cryptosporidium_res$FDR_FC <- ifelse(Cryptosporidium_res$logFC < 0, (1-Cryptosporidium_res$SpatialFDR), 1-(Cryptosporidium_res$SpatialFDR))
custom_colors <- c(
  "#1C75BC",  # dark blue at -1
  "gray80",   # gray at -0.9 (start of middle flat zone)
  "gray80",   # gray at 0.9 (end of middle flat zone)
  "#F7941D"   # dark red at 1
)

color_breaks <- scales::rescale(c(-1, -0.5,0.5, 1)) #only give color if FDR < 0.2 (in transformed 1-FDR)

ggplot(Cryptosporidium_res, aes(x = logFC, y = FineCellType)) +
  geom_quasirandom(aes(color = FDR_signFC, alpha = FDR_FC), groupOnX = TRUE, size = 1.5) +
  scale_color_gradientn(
    colors = custom_colors,
    values = color_breaks,
    limits = c(-1, 1),
    oob = scales::squish,
    guide = "colorbar",
    name = "1 - FDR"
  ) +
 # scale_alpha_continuous(name = "1 - SpatialFDR", range = c(0.2, 1)) +
  theme_minimal(base_size = 13) +
  labs( title = "Neighborhood DA: Nippo", x = "logFC") +
  theme(
    axis.title.y = element_blank(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )

ggsave( "neighborhoodDA_Crypto_Epithelial.svg",  plot = last_plot() , device = NULL, 
        path = images, width = 6,  height = 9,  units = c("in"),  dpi = 600,  limitsize = TRUE,  bg = NULL)
```

```{r}
library(scales)
milo <- readRDS("/path/to/analysis/directory/Seurat_Files/Myeloid_milo.rds")
# transfer labels
new_meta <- Ileum_filt@meta.data
new_meta$cellID <- rownames(new_meta)
cells_in_milo <- colnames(milo)
any(cells_in_milo %in% rownames(new_meta))
new_meta <- new_meta[new_meta$cellID%in% cells_in_milo, ]  

new_meta <- new_meta %>% dplyr::select(secondlevel)
colData(milo) <- cbind(colData(milo), new_meta)


milo <- countCells(milo, meta.data = colData(milo), sample = "Mouse")
samples_in_milo <- colnames(nhoodCounts(milo))
sample_metadata<- Ileum_filt@meta.data

design.df <- sample_metadata %>%
  filter(Mouse %in% samples_in_milo) %>%
  distinct(Mouse, .keep_all = TRUE) %>%
  arrange(factor(Mouse, levels = samples_in_milo))

rownames(design.df) <- design.df$Mouse
design.df <- design.df[samples_in_milo, , drop = FALSE]  # enforce order

# InfectionStatus is a factor with all levels
design.df$InfectionStatus <- factor(
  design.df$InfectionStatus,
  levels = c("Naive", "Candida", "Cryptosporidium", "MNV", 
             "Nippostrongylus", "SFB_YA", "Yersinia")
)

# Ensure Mouse is a factor and matches the order of sample names from Milo
design.df <- design.df[match(samples_in_milo, design.df$Mouse), , drop = FALSE]
rownames(design.df) <- design.df$Mouse

# Model matric
model <- model.matrix(~ 0 + InfectionStatus, data = design.df)
colnames(model) <- gsub("InfectionStatus", "", colnames(model))

#Define contrasts
library(limma)
infections <- c("Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia")
contrasts <- paste0(infections, " - Naive")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = model)

# Run testNhoods
contrast_results <- testNhoods(
  milo,
  design = model,
  design.df = design.df,
  model.contrasts = contrasts,
  reduced.dim = "multi_modal",
  fdr.weighting = "graph-overlap",
  norm.method = "TMM"
)

contrast_results <- annotateNhoods(milo, contrast_results, coldata_col = "secondlevel")



unique(Ileum_filt$FineCellType)
cell_subset <- c("Monocytes", "Neutrophils")

infections <- c("Candida", "Cryptosporidium", "MNV", "Nippostrongylus", "SFB_YA", "Yersinia")
custom_colors <- c("#1C75BC", "gray80", "gray80", "#F7941D")
color_breaks <- rescale(c(-1, -0.5, 0.5, 1))  # symmetric gradient around 0


# Determine global logFC limits across all contrast columns
logfc_cols <- grep("^logFC\\.", colnames(contrast_results), value = TRUE)
global_range <- range(contrast_results[, logfc_cols], na.rm = TRUE)

plots <- lapply(seq_along(infections), function(i) {
  inf <- infections[i]
  logFC_col <- paste0("logFC.", inf, "...Naive")
  if (!logFC_col %in% colnames(contrast_results)) return(NULL)
  
  res <- contrast_results[, c("Nhood", "secondlevel", "SpatialFDR", logFC_col)]
  colnames(res)[colnames(res) == logFC_col] <- "logFC"
  
  # Compute transformed values for color & alpha
  res$FDR_signFC <- ifelse(res$logFC < 0, -(1 - res$SpatialFDR), 1 - res$SpatialFDR)
  res$FDR_FC     <- 1 - res$SpatialFDR
  
  # Filter to Monocytes / Neutrophils or other selected cells
  res <- res[res$secondlevel %in% cell_subset, ]
  if (nrow(res) == 0) return(NULL)

  # Only last plot gets the colorbar
  color_guide <- if (i == length(infections)) "colorbar" else "none"
  
  ggplot(res, aes(x = logFC, y = secondlevel)) +
    geom_quasirandom(
      aes(color = FDR_signFC, alpha = 0.5),
      groupOnX = TRUE, size = 1.25
    ) +
    scale_color_gradientn(
      colors = custom_colors,
      values = color_breaks,
      limits = c(-1, 1),
      oob = scales::squish,
      guide = color_guide,
      name = "1 - FDR"
    ) +
    scale_x_continuous(limits = global_range) +   # ★ shared x-axis ★
    scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
    theme_minimal(base_size = 13) +
    labs(title = paste0(inf, " vs Naive"), x = "logFC") +
    theme(
      axis.title.y = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (i == length(infections)) "right" else "none",
      axis.text.y = if (i == 1) element_text() else element_blank()
    )
})


# Combine all infection plots horizontally
combined_plot <- plot_grid(plotlist = plots, ncol = 6, align = "hv")
combined_plot
# Save the combined plot
# Figure 3B
ggsave(
  filename = file.path(images, "Monocyte_Neutrophils_Neighborhood_milo_Beeswarm.svg"),
  plot = combined_plot,
  width = 22, height = 5, dpi = 600, bg = "white"
)

```

```{r}

qs_save(Ileum_filt, "/path/to/analysis/directory/Seurat_Files/Analysis2.5_Ileum.qs2")

```
