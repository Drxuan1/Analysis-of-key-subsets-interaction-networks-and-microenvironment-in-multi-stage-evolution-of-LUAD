library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(tidyverse)

setwd("")

pbmc <- readRDS("low_con_scobjendann2.rds")

pbmc$sample_id <- pbmc@meta.data$orig.ident

pbmc$group <- ifelse(grepl("IAC", pbmc$sample_id), "IAC", "AIS")

pbmc$celltype <- pbmc@active.ident

cell_type_mapping <- list(
  immune = c("T cells", "NK cells", "B cells", "Macrophage cells", "Neutrophil cells"),
  non_immune = c("Epithelial cells", "Endothelial cells", "Fibroblast cells")
)

cell_proportions <- pbmc@meta.data %>%
  group_by(group, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  rename(Celltype = celltype) %>%
  mutate(
    Cell_class = case_when(
      Celltype %in% cell_type_mapping$immune ~ "Immune cells",
      Celltype %in% cell_type_mapping$non_immune ~ "Non-immune cells",
      TRUE ~ "Other"
    )
  ) %>%

  filter(Cell_class != "Other") %>%

  group_by(group, Cell_class) %>%
  mutate(
    total_in_class = sum(count),  
    proportion_percent = (count / total_in_class) * 100  
  ) %>%
  ungroup() %>%

  mutate(
    Celltype = factor(
      Celltype,
      levels = c(cell_type_mapping$immune, cell_type_mapping$non_immune)
    )
  )

ais_samples <- str_sort(unique(filter(pbmc@meta.data, group == "AIS")$sample_id), numeric = TRUE)

iac_samples <- str_sort(unique(filter(pbmc@meta.data, group == "IAC")$sample_id), numeric = TRUE)

all_samples <- c(ais_samples, iac_samples)

cell_proportions_sample <- pbmc@meta.data %>%

  group_by(celltype, group, sample_id) %>%
  
  summarise(sample_count = n(), .groups = "drop") %>%
  
  mutate(
    Cell_class = case_when(
      celltype %in% cell_type_mapping$immune ~ "Immune cells",
      celltype %in% cell_type_mapping$non_immune ~ "Non-immune cells",
      TRUE ~ "Other"
    )
  ) %>%  

  group_by(group) %>%
  mutate(
    group_total = sum(sample_count),  
    sample_proportion = (sample_count / group_total) * 100 
  )%>%
  ungroup()%>% 
  mutate(
    celltype = factor(
      celltype,
      levels = c(cell_type_mapping$immune, cell_type_mapping$non_immune)
    )
  )

cell_proportions_sample$sample_id <- as.factor(as.character(cell_proportions_sample$sample_id))

max_immune <- max(filter(cell_proportions, Cell_class == "Immune cells")$proportion_percent)

max_non_immune <- max(filter(cell_proportions, Cell_class == "Non-immune cells")$proportion_percent)

y_limit_immune <- max_immune * 1.15

y_limit_non_immune <- max_non_immune * 1.15

group_colors <- c("AIS" = "#3B9AB2", "IAC" = "#F21A00")

ais_pal <- colorRampPalette(c(group_colors["AIS"], "#A6D8F0"))(length(ais_samples))#"#A6D8F0"

iac_pal <- colorRampPalette(c(group_colors["IAC"], "#F9B2A8"))(length(iac_samples))#"#F9B2A8"

sample_colors <- c(setNames(ais_pal, ais_samples), setNames(iac_pal, iac_samples))

cell_proportions_group <- cell_proportions_sample %>%
  group_by(celltype, group) %>%
  summarise(proportion_percent = sum(sample_proportion), .groups = "drop") %>%
  mutate(cell_group = interaction(celltype, group, sep = "\n"))%>%
  mutate(
    Cell_class = case_when(
      celltype %in% cell_type_mapping$immune ~ "Immune cells",
      celltype %in% cell_type_mapping$non_immune ~ "Non-immune cells",
      TRUE ~ "Other"
    ))

p_immune <- ggplot(filter(cell_proportions_sample, Cell_class == "Immune cells"),
                   aes(x = celltype, y = sample_proportion, fill = sample_id))+

  geom_bar(stat = "identity", position = "stack", width = 0.8)+

  facet_wrap(~group, ncol = 2)+

  scale_fill_brewer(palette = "Set1")+ 
  
  labs(title = "Cell Type Proportions by Group",
       x = "Cell Type", y = "Proportion", fill = "Sample ID")+
  
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#f0f0f0"))+
  
  geom_label( 
    data = filter(cell_proportions_group[which(cell_proportions_group$Cell_class == "Immune cells"), ]), 
    aes(x = celltype, y = proportion_percent, 
        label = paste0(round(proportion_percent, 1), "%")), 
    position = position_dodge(width = 0.8), 
    vjust = -0.001, 
    hjust = 0.5, 
    size = 3.5, 
    color = "#333333", 
    fill = "white", 
    fontface = "bold",
    label.size = NA,
    show.legend = FALSE 
  )+

  labs(
    x = " ", y = "Proportion within Group (%)",
    title = "Immune cell composition comparison",
    fill = "Sample ID"
  )+

  scale_fill_manual(values = sample_colors)+
  
  geom_vline(xintercept = seq(1.5, length(cell_type_mapping$immune)-0.5, 1), 
             color = "#CCCCCC", linetype = "solid", size = 1)+

  theme_minimal()+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#CCCCCC", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "#333333"),
    axis.text.y = element_text(size = 12, color = "#333333"),
    axis.title = element_text(size = 14, face = "bold", color = "#333333"),
    axis.ticks = element_line(color = "#CCCCCC"),
    legend.position = "top",  
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 40, 10),  
    strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC"),
    strip.text = element_text(size = 14, face = "bold", color = "#333333")
  )
#
p_non_immune <- ggplot(filter(cell_proportions_sample, Cell_class == "Non-immune cells"),
                       aes(x = celltype, y = sample_proportion, fill = sample_id)) +
  
  geom_bar(stat = "identity", position = "stack", width = 0.8)+

  facet_wrap(~group, ncol = 2)+

  scale_fill_brewer(palette = "Set1")+ 
  
  labs(title = "Cell Type Proportions by Group",
       x = "Cell Type", y = "Proportion", fill = "Sample ID") +
  
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#f0f0f0"))+
  
  geom_label( 
    data = filter(cell_proportions_group[which(cell_proportions_group$Cell_class == "Non-immune cells"), ]), 
    aes(x = celltype, y = proportion_percent, 
        label = paste0(round(proportion_percent, 1), "%")), 
    position = position_dodge(width = 0.8), 
    vjust = -0.001, 
    hjust = 0.5, 
    size = 3.5, 
    color = "#333333", 
    fill = "white",
    fontface = "bold",
    label.size = NA,
    show.legend = FALSE 
  )+
  
  labs(
    x = " ", y = "Proportion within Group (%)",
    title = "Non-immune cell composition comparison",
    fill = "Sample ID"
  )+

  scale_fill_manual(values = sample_colors)+
  
  geom_vline(xintercept = seq(1.5, length(cell_type_mapping$non_immune)-0.5, 1), 
             color = "#CCCCCC", linetype = "solid", size = 1) +

  theme_minimal()+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#CCCCCC", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "#333333"),
    axis.text.y = element_text(size = 12, color = "#333333"),
    axis.title = element_text(size = 14, face = "bold", color = "#333333"),
    axis.ticks = element_line(color = "#CCCCCC"),
    legend.position = "none",  
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 40, 10),  
    strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC"),
    strip.text = element_text(size = 14, face = "bold", color = "#333333")
  )

combined_plot <- p_immune / p_non_immune

print(combined_plot)

