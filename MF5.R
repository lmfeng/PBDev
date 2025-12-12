library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(Seurat)
library(tidyverse)

Time_order <- c("E45", "E55", "E66", "E76", "E85", "E94", "E104", "E109", "P0/P3")
region_order <- c("PFC", "STR", "THA")
subclass_order <- c("OPC cycle", "OPC", "COP", "Oligo")

cols_celltype <- c(
  COP = "#1f77b4", 
  OPC = "#ff7f0e", 
  "OPC cycle" = "#2ca02c", 
  Oligo = "#d62728"
)

cols_time <- c(
  E45 = "#edc815", 
  E55 = "#35459d", 
  E66 = "#b12524", 
  E76 = "#ef8263", 
  E85 = "#96cecd", 
  E94 = "#9682bd", 
  E104 = "#0e8c45", 
  E109 = "#e18a3b", 
  "P0/P3" = "#dd7694"
)

cols_region <- c(
  PFC = '#1F77D4',
  STR = '#FF7F1E',
  THA = '#2CA06C'
)

pig_opc <- readRDS("pig_opc.rds")

cell_metadata <- final_pig_opc_final@meta.data

cell_metadata$Time <- ifelse(cell_metadata$time %in% c("P0", "P3"), "P0/P3", as.character(cell_metadata$time))
cell_metadata$Time <- factor(cell_metadata$Time, levels = Time_order)


cell_counts <- table(
  cell_metadata$subclass,
  cell_metadata$region, 
  cell_metadata$Time
)

cell_counts_long <- as.data.frame.table(cell_counts) %>%
  dplyr::rename(Subclass = Var1, Region = Var2, Time = Var3, CellCount = Freq) %>%
  mutate(
    Time = factor(Time, levels = Time_order),
    Region = factor(Region, levels = region_order),
    Subclass = factor(Subclass, levels = subclass_order)
  )

cell_proportions <- cell_counts_long %>%
  group_by(Region, Time) %>%
  summarise(TotalCells = sum(CellCount), .groups = "drop") %>%
  mutate(Proportion = TotalCells / sum(TotalCells))

rownames(regulonActivity_byGroup_Scaled) <- gsub("\\(\\+\\)", "", rownames(regulonActivity_byGroup_Scaled))
tf_activity <- regulonActivity_byGroup_Scaled

subclass_mapping <- cell_metadata %>%
  select(subclass, region, Time) %>%
  distinct() %>%
  mutate(Group = paste(region, Time, sep = "_"))

all_region_time_combinations <- expand.grid(Region = region_order, Time = Time_order) %>%
  mutate(Group = paste(Region, Time, sep = "_")) %>%
  arrange(factor(Region, levels = region_order), factor(Time, levels = Time_order))

tf_activity_region_time <- matrix(NA, 
                                  nrow = nrow(tf_activity), 
                                  ncol = nrow(all_region_time_combinations),
                                  dimnames = list(rownames(tf_activity), all_region_time_combinations$Group))

for (i in 1:nrow(all_region_time_combinations)) {
  group <- all_region_time_combinations$Group[i]
  region <- all_region_time_combinations$Region[i]
  time <- all_region_time_combinations$Time[i]
  
  subclasses_in_group <- subclass_mapping %>%
    filter(region == !!region, Time == !!time) %>%
    pull(subclass)
  
  if (length(subclasses_in_group) > 0) {
    group_activity <- rep(0, nrow(tf_activity))
    total_weight <- 0
    
    for (subclass in subclasses_in_group) {
      if (subclass %in% colnames(tf_activity)) {
        subclass_count <- cell_counts_long %>%
          filter(Subclass == subclass, Region == region, Time == time) %>%
          pull(CellCount)
        
        if (length(subclass_count) > 0 && subclass_count > 0) {
          group_activity <- group_activity + tf_activity[, subclass] * subclass_count
          total_weight <- total_weight + subclass_count
        }
      }
    }
    
    if (total_weight > 0) {
      tf_activity_region_time[, group] <- group_activity / total_weight
    }
  }
}

tf_activity_original <- tf_activity_region_time

tf_activity_zscore <- t(apply(tf_activity_original, 1, function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  if (sd_val == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - mean_val) / sd_val)
  }
}))

tf_subclass_zscore <- t(apply(regulonActivity_byGroup_Scaled, 1, function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  if (sd_val == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - mean_val) / sd_val)
  }
}))

activity_threshold <- 0.5 

tf_subclass_comb <- apply(tf_subclass_zscore, 1, function(z) {
  high_subclasses <- colnames(tf_subclass_zscore)[z > activity_threshold]
  if (length(high_subclasses) == 0) {
    max_subclass <- colnames(tf_subclass_zscore)[which.max(z)]
    return(max_subclass)
  } else {
    high_subclasses_sorted <- high_subclasses[order(match(high_subclasses, subclass_order))]
    return(paste(high_subclasses_sorted, collapse = "+"))
  }
})


tf_comb_mapping <- data.frame(
  TF = names(tf_subclass_comb),
  Subclass_Combination = unlist(tf_subclass_comb),
  stringsAsFactors = FALSE
)


all_combinations <- unique(tf_comb_mapping$Subclass_Combination)

sorted_combinations <- all_combinations[order(
  sapply(all_combinations, function(comb) {
    subclasses <- strsplit(comb, "\\+")[[1]]
    length(subclasses)
  }),
  sapply(all_combinations, function(comb) {
    subclasses <- strsplit(comb, "\\+")[[1]]
    paste(match(subclasses, subclass_order), collapse = ",")
  })
)]

tf_comb_factor <- factor(
  tf_comb_mapping$Subclass_Combination,
  levels = sorted_combinations
)
names(tf_comb_factor) <- tf_comb_mapping$TF


valid_tfs <- intersect(rownames(tf_activity_zscore), names(tf_comb_factor))
tf_activity_filtered <- tf_activity_zscore[valid_tfs, , drop = FALSE]
tf_comb_filtered <- tf_comb_factor[valid_tfs]


cell_prop_anno <- cell_proportions %>%
  mutate(Group = paste(Region, Time, sep = "_")) %>%
  select(Group, Proportion) %>%
  right_join(all_region_time_combinations, by = "Group") %>%
  mutate(Proportion = ifelse(is.na(Proportion), 0, Proportion)) %>%
  column_to_rownames("Group")

cell_prop_anno <- cell_prop_anno[colnames(tf_activity_filtered), , drop = FALSE]

col_anno <- data.frame(
  Region = sapply(strsplit(colnames(tf_activity_filtered), "_"), `[`, 1),
  Time = sapply(strsplit(colnames(tf_activity_filtered), "_"), `[`, 2),
  stringsAsFactors = FALSE
)
rownames(col_anno) <- colnames(tf_activity_filtered)
col_anno$Time <- factor(col_anno$Time, levels = Time_order)
col_anno$Region <- factor(col_anno$Region, levels = region_order)


zscore_range <- max(abs(tf_activity_filtered), na.rm = TRUE)
zscore_range <- max(zscore_range, 2)  

activity_colors <- colorRamp2(
  seq(-zscore_range, zscore_range, length = 11), 
  rev(brewer.pal(11, "RdBu"))
)

comb_colors <- c()
for (comb in sorted_combinations) {
  subclasses <- strsplit(comb, "\\+")[[1]]
  if (length(subclasses) == 1) {
    comb_colors[comb] <- cols_celltype[subclasses]
  } else {
    if (length(subclasses) == 2) {
      color1 <- cols_celltype[subclasses[1]]
      color2 <- cols_celltype[subclasses[2]]
      comb_colors[comb] <- colorRampPalette(c(color1, color2))(3)[2]  
    } else {
      comb_colors[comb] <- "#808080"
    }
  }
}

row_anno <- rowAnnotation(
  SubclassCombination = anno_text(
    x = as.character(tf_comb_filtered),
    which = "row",
    gp = gpar(fontsize = 8, col = "black"),
    location = 0.5,
    just = "center"
  ),
  ColorBlock = anno_simple(
    as.character(tf_comb_filtered),
    col = comb_colors,
    width = unit(0.5, "cm")  
  ),
  annotation_name_rot = 90,
  annotation_name_side = "top",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  gap = unit(2, "mm")
)


top_anno <- HeatmapAnnotation(
  CellProportion = anno_barplot(
    cell_prop_anno$Proportion,
    bar_width = 0.7,
    gp = gpar(fill = "darkgreen", col = "darkgreen"),
    axis_param = list(
      gp = gpar(fontsize = 8),
      at = c(0, max(cell_prop_anno$Proportion, na.rm = TRUE)),
      labels = c("0", "Max")
    ),
    height = unit(2, "cm")
  ),
  Region = col_anno$Region,
  Time = col_anno$Time,
  col = list(
    Region = cols_region,
    Time = cols_time
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Region = list(nrow = 1),
    Time = list(nrow = 2)
  )
)


mean_activity_by_tf <- rowMeans(tf_activity_original[valid_tfs, ], na.rm = TRUE)

top_tf_labels <- lapply(levels(tf_comb_filtered), function(comb) {
  tfs_in_comb <- names(tf_comb_filtered)[tf_comb_filtered == comb]
  if (length(tfs_in_comb) == 0) return(character(0))
  
  comb_activity <- mean_activity_by_tf[tfs_in_comb]
  threshold <- quantile(comb_activity, probs = 0.0, na.rm = TRUE)
  top_tfs <- names(comb_activity)[comb_activity >= threshold]
  return(top_tfs)
})

top_tf_labels <- unlist(top_tf_labels)

row_labels <- ifelse(
  rownames(tf_activity_filtered) %in% top_tf_labels,
  rownames(tf_activity_filtered),
  ""
)

main_heatmap <- Heatmap(
  tf_activity_filtered,
  name = "TF Activity\n(z-score)",
  col = activity_colors,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = col_anno$Region,
  column_order = order(col_anno$Region, col_anno$Time),
  row_split = tf_comb_filtered,
  row_labels = row_labels,
  row_title = "Subclass Combination",
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 12),
  column_title = "Brain Region",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 8),
  top_annotation = top_anno,
  right_annotation = row_anno, 
  heatmap_legend_param = list(
    title_position = "leftcenter-rot",
    legend_height = unit(4, "cm")
  )
)

bottom_anno <- HeatmapAnnotation(
  Time = anno_text(
    col_anno$Time,
    location = 0.5,
    rot = 45,
    just = "center",
    gp = gpar(fontsize = 8, col = cols_time[as.character(col_anno$Time)])
  ),
  show_annotation_name = FALSE
)

combined_heatmap <- main_heatmap %v% bottom_anno

subclass_legend <- Legend(
  labels = names(comb_colors),
  legend_gp = gpar(fill = comb_colors),
  title = "Subclass Combination",
  ncol = 1,
  title_position = "topleft"
)

pdf("Final_optimized_tf_activity_heatmap_with_all_subclass_combinations.pdf", width = 16, height = 14)
draw(
  combined_heatmap,
  column_title = "TF Activity Across Brain Regions and Developmental Stages\n(Grouped by Subclass Combinations, Z-score Normalized by TF)",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(subclass_legend)
)
dev.off()
