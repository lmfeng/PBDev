library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(ggrepel)
library(dplyr)

seurat_obj <-  readRDS("scrds.rds")
geneset_file <- "geneset_pig_Final.gs"
genesets <- read.delim(geneset_file, header = TRUE, stringsAsFactors = FALSE)

time_order <- c('E45', 'E55', 'E66', 'E76', 'E85', 'E94', 'E104', 'E109', 'P0', 'P3')
regions <- c('PFC', 'STR', 'THA')
disease_names <- c("MDD_all", "SCZ_all", "NEUROTICISM_all")

seurat_obj@meta.data$time <- factor(seurat_obj@meta.data$time, levels = time_order)
seurat_obj@meta.data$region <- factor(seurat_obj@meta.data$region, levels = regions)

disease_genesets <- genesets %>%
  dplyr::filter(TRAIT %in% disease_names) %>%
  dplyr::select(TRAIT, GENESET) %>%
  tidyr::separate_rows(GENESET, sep = ",") %>%
  tidyr::separate(GENESET, into = c("gene", "weight"), sep = ":", convert = TRUE) %>%
  dplyr::rename(disease = TRAIT) %>%
  dplyr::mutate(weight = 1) %>%
  dplyr::filter(!is.na(gene), gene != "")


risk_genes <- unique(disease_genesets$gene)
present_genes <- intersect(risk_genes, rownames(seurat_obj))


avg_expr <- AverageExpression(seurat_obj,
                              features = present_genes,
                              group.by = c("time", "region"),
                              slot = "data")$RNA


expr_df <- as.data.frame(avg_expr) %>%
  rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols = -gene, names_to = "group", values_to = "expression") %>%
  tidyr::separate(group, into = c("time", "region"), sep = "_") %>%
  dplyr::mutate(
    time = factor(time, levels = time_order),
    region = factor(region, levels = regions)
  ) %>%
  dplyr::inner_join(disease_genesets, by = "gene", relationship = "many-to-many")

gene_specificity <- expr_df %>%
  group_by(gene, disease) %>%
  summarize(

    THA_E76_E104_mean = mean(expression[region == "THA" & time %in% c("E76", "E85", "E94", "E104")]),

    THA_other_mean = mean(expression[region == "THA" & !time %in% c("E76", "E85", "E94", "E104")]),

    other_regions_E76_E104_mean = mean(expression[region != "THA" & time %in% c("E76", "E85", "E94", "E104")]),
    
    other_regions_other_mean = mean(expression[region != "THA" & !time %in% c("E76", "E85", "E94", "E104")]),

    THA_specificity = THA_E76_E104_mean / (THA_other_mean + 0.01),  
    region_specificity = THA_E76_E104_mean / (other_regions_E76_E104_mean + 0.01),
    overall_specificity = THA_specificity * region_specificity,
    
    .groups = "drop"
  )

specific_genes <- gene_specificity %>%
  filter(
    THA_E76_E104_mean > 1.5,  
    THA_specificity > 1  ,    
    region_specificity > 1  
  ) %>%
  arrange(desc(overall_specificity))

print(specific_genes)

specific_genes_expr <- expr_df %>%
  filter(gene %in% specific_genes$gene)


specific_genes_expr <- specific_genes_expr %>%
  left_join(
    specific_genes %>% select(gene, disease, overall_specificity),
    by = c("gene", "disease")
  )


heatmap_data <- specific_genes_expr %>%
  group_by(gene, time, region, disease, overall_specificity) %>%
  summarize(mean_expr = mean(expression), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(z_score = scale(mean_expr)) %>%
  ungroup()


heatmap_data <- heatmap_data %>%
  group_by(disease) %>%
  mutate(
    gene = factor(gene, levels = unique(gene[order(-overall_specificity)]))
  ) %>%
  ungroup()


p_heatmap <- ggplot(heatmap_data, aes(x = interaction(time, region), y = gene, fill = z_score)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "magma", name = "Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(face = "italic", size = 8),
    strip.text = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  facet_wrap(~ disease, scales = "free_y", ncol = 1)

print(p_heatmap)


bubble_data_raw <- specific_genes_expr %>%
  group_by(gene, disease, region, time) %>%
  summarize(
    mean_expr = mean(expression),
    .groups = "drop"
  ) %>%
  inner_join(
    specific_genes %>% select(gene, disease, overall_specificity, region_specificity),
    by = c("gene", "disease")
  ) %>%
  group_by(time, region) %>%
  mutate(
    is_empty = all(mean_expr <= 0)  
  ) %>%
  ungroup()

bubble_data <- bubble_data_raw %>%
  filter(is_empty == FALSE) %>%  
  
  mutate(
    time = factor(time, levels = time_order[time_order %in% unique(time)])
  ) %>%
  
  mutate(
    region = factor(region, levels = regions[regions %in% unique(region)])
  ) %>%
  
  group_by(disease, region) %>%
  arrange(desc(region_specificity), .by_group = TRUE) %>%
  
  mutate(gene = factor(gene, levels = unique(gene))) %>%
  ungroup() %>%
  
  group_by(gene) %>%
  mutate(z_score = as.vector(scale(mean_expr))) %>%
  ungroup()

p_bubble <- ggplot(bubble_data, 
                   aes(x = time, y = gene, 
                       size = mean_expr,
                       color = z_score)) +
  geom_point(alpha = 0.8) +
  
  facet_grid(disease ~ region, scales = "free_y", space = "free_y") +
  
  scale_color_gradientn(
    colors = c("#3B4992", "#FFFFFF", "#E41A1C"),
    name = "Z-score",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  
  scale_size_continuous(
    name = "Average expression",
    range = c(1, 8),
    guide = guide_legend(title.position = "top", title.hjust = 0.5),
    breaks = scales::breaks_pretty(n = 4)
  ) +

  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 10, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, 
                              margin = margin(b = 15)),
    axis.title.x = element_text(size = 10, margin = margin(t = 10)),
    axis.title.y = element_text(size = 10, margin = margin(r = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, face = "italic", color = "black"),
    strip.text.x = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    strip.background = element_rect(fill = "gray90", color = NA),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = -10),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.8, "cm"),
    panel.spacing.y = unit(0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  
  labs(
    title = "Risk genes with specific high expression during the E76-E104 period in the THA region",
    x = "Time",
    y = "Gene"
  )

print(p_bubble)

ggsave(
  filename = "THA_specific_genes_bubbleplot.pdf",
  plot = p_bubble,
  width = 14 + length(unique(bubble_data$region)) * 0.5,
  height = 5 + nrow(specific_genes) * 0.15,
  dpi = 300,
  bg = "white",
  device = cairo_pdf
)
