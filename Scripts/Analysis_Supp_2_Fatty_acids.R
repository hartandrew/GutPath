# Description. This analysis will load the weight normalized abundance of Fatty acids measured in Nippostrongylus infected mouse samples and naive mouse samples. Samples include fecal/luminal contents of the colon, fecal/luminal contents of the distal small intestine, and scrapings of the distal small intestinal (likely capturing epithelial cells, mucus, etc)
# Figure 5G and Figure S6C data are generated in these scripts

# Load the Libraries----
library(tidyverse)
library(Matrix)
library(scales)
library(rjson)
library(DT)
library(readxl)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(ggbreak)
library(rstatix)
library(ggpubr)
library(vegan)
library(pheatmap)
library(RColorBrewer)


# Read data -----
data <- read_excel("/path/to/directory/FA_results_R.xlsx")
colnames(data)

out <- "/path/to/directory/Output"
seurat <- "/path/to/directory/Seurat_Files"
images <- "/path/to/directory/Images"
CSV <- "/path/to/directory/CSV"

# transform to long data format with analytes in one column and their values in another 
data_long <- data %>%
  pivot_longer(
    cols = Palmitic_acid:FA_24_0, 
    names_to = "Analyte", 
    values_to = "Value"
  )

# Add Common Names. Notably, the common names are only best guesses for the fatty acids with multiple unsaturations. The saturated fatty acids are more definite as there aren't often alternate conformations 
data_long$Common_Analyte <- data_long$Analyte
data_long$Common_Analyte[data_long$Common_Analyte == "FA_18_2"] <- "Linoleic_acid"
data_long$Common_Analyte[data_long$Common_Analyte == "FA_20_0"] <- "Arachidic_acid"
data_long$Common_Analyte[data_long$Common_Analyte == "FA_20_3"] <- "DGLA"
data_long$Common_Analyte[data_long$Common_Analyte == "FA_20_4"] <- "Arachidonic_acid"
data_long$Common_Analyte[data_long$Common_Analyte == "FA_22_5"] <- "Docosapentaenoic_acid"
data_long$Common_Analyte[data_long$Common_Analyte == "FA_24_0"] <- "Lignoceric_acid"



# Function for graphing all the data, separating by type - which type of sample
make_fa_plots <- function(current_type) {
  
  # Filter for the sample type
  type_data <- data_long %>% filter(Type == current_type)
  
  ggplot(type_data, aes(x = Condition, y = Value, fill = Condition)) +
    stat_summary(fun = "mean", geom = "bar", alpha = 0.7, color = "black") +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
    scale_fill_manual(values = c("#DBB075", "#D3D3D3"), breaks = c("Nippo", "Naive"))+
    # jitter the overlayed dots
    geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.8) +
    # wrap by fatty acid
    facet_wrap(~Analyte, scales = "free_y") + 
    theme_classic() +
    labs(
      title = paste("Fatty Acid Levels -", current_type),
      y = "Concentration (nmol/mg) ",
      x = "Condition"
    ) +
    theme(legend.position = "none")
}

# Plot each
unique_types <- unique(data$Type)

for (t in unique_types) {
  print(make_fa_plots(t))
}

unique(data_long$Analyte)

# Make a long plot showing all the analytes for each sample type

make_combined_plot <- function(current_type) {
  
  type_data <- data_long %>% filter(Type == current_type)
  
  ggplot(type_data, aes(x = Analyte, y = Value, fill = Condition)) +
    stat_summary(fun = "mean", geom = "bar", 
                 position = position_dodge(width = 0.8), 
                 alpha = 0.7, color = "black", width = 0.7) +
    # Error bars to dodged bars
    stat_summary(fun.data = "mean_se", geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
               size = 1.5, color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("#DBB075", "#D3D3D3"), breaks = c("Nippo", "Naive"))+
    theme_classic() +
    labs(
      title = paste("Fatty Acid Profile -", current_type),
      y = "Value",
      x = "Fatty Acid Analyte"
    ) +
    # Rotate X-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


unique_types <- unique(data$Type)

for (t in unique_types) {
  print(make_combined_plot(t))
}




# Combine the plots into one 
plot_list <- list()
unique_types <- unique(data$Type)

# Loop through
for (i in seq_along(unique_types)) {
  current_type <- unique_types[i]
  
  p <- data_long %>% 
    filter(Type == current_type) %>%
    ggplot(aes(x = Analyte, y = Value, fill = Condition)) +
    scale_fill_manual(values = c("#DBB075", "#D3D3D3"), breaks = c("Nippo", "Naive"))+
    stat_summary(fun = "mean", geom = "bar", 
                 position = position_dodge(width = 0.8), 
                 alpha = 0.7, color = "black", width = 0.7) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
               size = 1.2, color = "black", alpha = 0.8) +
    theme_classic() +
    labs(subtitle = paste("Type:", current_type), y = "Value", x = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Only show the legend on the top plot
  if(i > 1) p <- p + theme(legend.position = "none")
  
  plot_list[[i]] <- p
}
combined_stack <- wrap_plots(plot_list, ncol = 1) + 
  plot_annotation(title = "Fatty Acid Profiles Across Types")


print(combined_stack)

# Separate them out and customize

a <- plot_list[[1]] + theme(axis.text.x = element_blank())
b <- plot_list[[2]] + theme(axis.text.x = element_blank())
c <- plot_list[[3]]

a / b/ c
ggsave("BArchart_Fatty_Acids_All_Samples.svg", plot = last_plot(), path = images, width = 7, height = 10)

# Perform Stats and correct within each sample type
stats_results <- data_long %>%
  group_by(Type, Analyte) %>%  
  t_test(Value ~ Condition) %>%
  group_by(Type) %>%           
  adjust_pvalue(method = "BH") %>%
  add_significance()




# Many of the individual comparisons seem significant but do not survive multiple comparisons corrections. However, there is a shared trend across all FAs and we can test whether this is significant 



tissue_data_clean <- data[data$Type == "Ileum_S",] %>%
  drop_na(Palmitic_acid:FA_24_0)

fatty_acids_matrix <- tissue_data_clean %>% select(Palmitic_acid:FA_24_0)

permanova_result <- adonis2(fatty_acids_matrix ~ Condition, 
                            data = tissue_data_clean, 
                            method = "euclidean")


print(permanova_result)

# This PERMANOVA looks at the difference in composition across the entire input matrix - in this case asking a broader question - are FAs different in Ileum scrapings of Nippo infected mice compared to Naive
# The results suggest that r2 = 0.47 which is a large association of nearly 50% of variation explained by the conditions. the significant p value also supports the conclusion that Nippo affects the fatty acids in the mice. 



# Graph the relationships of Fatty acids to Worm counts 

data_nippo  <- data[data$Condition == "Nippo",]
data_nippo$Worm_Count <- c(41, 41, 41, 32, 32, 32, 59, 59, 59, 68, 68, 68) # TWorm counts were alculated on the day of tissue sampling by APH prior to reporting of fatty acid results


ggplot(data_nippo[data_nippo$Type == "Ileum_S",], aes(x =Oleic_acid, y = Worm_Count)) +
  geom_point(aes(color = Type), size = 5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  # Add correlation coefficient
  #stat_cor(method = "pearson") + 
  theme_classic() +
  scale_color_manual(values = "#DBB075") +
  labs(
    title = "Oleic Acid vs. Parasite Burden",
    x = "Oleic Acid Concentration",
    y = "Worm Count",
    color = "Tissue Type"
  )
ggsave("Worm_vs_Oleic_acid_line.svg", plot = last_plot(), path = images, width = 7, height = 7)

ggplot(data_nippo[data_nippo$Type == "Ileum_S",], aes(x =Palmitic_acid, y = Worm_Count)) +
  geom_point(aes(color = Type), size = 5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = "#DBB075") +
  labs(
    title = "Palmitic Acid vs. Parasite Burden",
    x = "Palmitic Acid Concentration",
    y = "Worm Count",
    color = "Tissue Type"
  )
ggsave("Worm_vs_Palmitic_acid_line.svg", plot = last_plot(), path = images, width = 7, height = 7)


