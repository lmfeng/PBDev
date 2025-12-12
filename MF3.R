source("~/R/basic code/library.R")
load("/Project_PBD/Rds/final/V2/ExN_L/Monocle3/cds.ExN_2nd.RData")
cds.ExN_2nd <- order_cells(cds.ExN_2nd, root_pr_nodes=get_earliest_principal_node2(cds.ExN_2nd,age_bin = "E45",celltype = "enIPC cycle"))

pseudotime_df <- data.frame(cell = names(cds.ExN_2nd@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]),
                            pseudotime = cds.ExN_2nd@principal_graph_aux@listData[["UMAP"]][["pseudotime"]],
                            stringsAsFactors = FALSE)

rownames(pseudotime_df) <- pseudotime_df$cell
ExN_2nd@meta.data$m3_pseudotime <- pseudotime_df[rownames(ExN_2nd@meta.data), "pseudotime"]

ExN_2nd@meta.data$Pbin <- cut(ExN_2nd@meta.data$m3_pseudotime, 
                              breaks = 1000, 
                              labels = FALSE,
                              include.lowest = TRUE)

meta_df <- ExN_2nd@meta.data[, c("Pbin", "subclass_ord")]
count_df <- meta_df %>%
  group_by(Pbin, subclass_ord) %>%
  summarise(cell_count = n(), .groups = 'drop')

bin_totals <- count_df %>%
  group_by(Pbin) %>%
  summarise(total_cells = sum(cell_count), .groups = 'drop')

percentage_df <- count_df %>%
  left_join(bin_totals, by = "Pbin") %>%
  mutate(percentage = (cell_count / total_cells) * 100)

###
threshold_percentage <- 1  
selected_subclasses <- c("enIPC", "mExN1", "mExN2","UpL ExN","DpL ExN","L6 CT","L6 B")
smooth_method <- "loess"  #"gam"/"loess"

process_percentage_data <- function(df, threshold) {
  
  complete_grid <- expand.grid(
    Pbin = unique(df$Pbin),
    subclass_ord = selected_subclasses,
    stringsAsFactors = FALSE
  )
  
  df %>%
    group_by(Pbin) %>%
    mutate(
      percentage_filtered = ifelse(percentage < threshold, 0, percentage),
      new_total = sum(percentage_filtered),
      percentage_normalized = ifelse(
        new_total == 0,
        0,
        percentage_filtered * (100 / new_total)
      )
    ) %>%
    ungroup() %>%
    select(-new_total) %>%
    right_join(complete_grid, by = c("Pbin", "subclass_ord")) %>%
    mutate(
      across(c(cell_count, percentage, percentage_filtered, percentage_normalized),
             ~ ifelse(is.na(.), 0, .))
    ) %>%
    mutate(subclass_ord = factor(subclass_ord, levels = selected_subclasses)) %>%
    
    arrange(Pbin, subclass_ord)
}

filtered_data <- process_percentage_data(percentage_df, threshold_percentage) %>%
  filter(subclass_ord %in% selected_subclasses)

plot_data <- filtered_data
plot_data$percentage_filtered[plot_data$subclass_ord=="DpL ExN"]

filter_consecutive_values <- function(x, min_consecutive = 5) {
  n <- length(x)
  keep <- rep(FALSE, n)
  
  non_zero <- x != 0
  
  rle_result <- rle(non_zero)
  starts <- cumsum(c(1, rle_result$lengths[-length(rle_result$lengths)]))
  ends <- cumsum(rle_result$lengths)
  
  for (i in 1:length(rle_result$values)) {
    if (rle_result$values[i] && rle_result$lengths[i] >= min_consecutive) {
      keep[starts[i]:ends[i]] <- TRUE
    }
  }
  
  return(ifelse(keep, x, 0))
}

plot_data <- plot_data %>%
  group_by(subclass_ord) %>%
  mutate(
    percentage_filtered = filter_consecutive_values(percentage_filtered, min_consecutive = 5)
  ) %>%
  ungroup()

plot_data$percentage_filtered[plot_data$subclass_ord=="DpL ExN"]

create_plot <- function(data, use_filtered = TRUE) {
  filtered_cols <- subclass_cols[names(subclass_cols) %in% selected_subclasses]
  
  if(smooth_method == "gam") {
    smooth_layer <- geom_smooth(
      aes(y = .data[["percentage_filtered"]], color = subclass_ord),
      method = "gam",
      formula = y ~ s(x, bs = "tp", k = 5),
      se = F,
      linewidth = 2.5
    )
  } else if(smooth_method == "loess") {
    smooth_layer <- geom_smooth(
      aes(y = .data[["percentage_filtered"]], color = subclass_ord),
      method = "loess",
      span = 0.8,
      se = T,
      linewidth = 3
    )
  }
  
  ggplot(data, aes(x = Pbin)) +
    geom_point(
      aes(y = .data[["percentage_filtered"]], color = subclass_ord,fill = subclass_ord), 
      alpha = 1, 
      size = 0.9
    ) +
    smooth_layer +
    scale_color_manual(values = filtered_cols) +
    scale_fill_manual(values = filtered_cols) +
    labs(
      x = "Pseudotime",
      y = "Percentage",
      color = "subclass",
      title = paste()
    ) +
    theme_bw(base_size = 16) +
    scale_x_continuous(breaks = seq(0, 1000, by = 200), limits = c(1, 1000)) +
    scale_y_continuous(limits = c(0, NA))
}

final_plot <- create_plot(plot_data, use_filtered = TRUE)

EXN_m <- readRDS("/Project_PBD/Rds/final/V2/ExN_L/WGCNA/EXN_m.rds")
DefaultAssay(EXN_m) <- "RNA"
EXN_m@assays$RNA@data <- EXN_m@assays$RNA@counts
EXN_m <- NormalizeData(EXN_m)
EXN_m <- ScaleData(EXN_m, features=VariableFeatures(EXN_m))

PlotDendrogram(EXN_m, main='Dendrogram')
load("/Project_PBD/Rds/final/V2/ExN_L/Monocle3/cds.ExN_per.RData")
cds.ExN_per <- order_cells(cds.ExN_per, root_pr_nodes=get_earliest_principal_node2(cds.ExN_per,age_bin = "E76",celltype = "mExN1"))

pseudotime_df <- data.frame(cell = names(cds.ExN_per@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]),
                            pseudotime = cds.ExN_per@principal_graph_aux@listData[["UMAP"]][["pseudotime"]],
                            stringsAsFactors = FALSE)
rownames(pseudotime_df) <- pseudotime_df$cell
ExN_per@meta.data$m3_pseudotime <- pseudotime_df[rownames(ExN_per@meta.data), "pseudotime"]

MEs <- GetMEs(EXN_m, harmonized = FALSE)
MEs$grey <- NULL

extra_info <- data.frame(
  pseudotime = numeric(nrow(MEs)),
  Age = character(nrow(MEs)),     
  subclass = character(nrow(MEs)),  
  row.names = rownames(MEs)
)

for (cell_barcode in rownames(MEs)) {
  cell_index <- which(rownames(ExN_per@meta.data) == cell_barcode)
  if (length(cell_index) > 0) {
    extra_info[cell_barcode, "pseudotime"] <- as.numeric(ExN_per@meta.data[cell_index, "m3_pseudotime"])
    extra_info[cell_barcode, "Age"] <- as.character(ExN_per@meta.data[cell_index, "time"])
    extra_info[cell_barcode, "subclass"] <- as.character(ExN_per@meta.data[cell_index, "subclass"]) 
  } else {
    warning(paste("Cell", cell_barcode, "not found in ExN_per meta.data"))
  }
}

MEs_info <- cbind(extra_info, MEs)
MEs_long <- reshape2::melt(
  MEs_info,
  id.vars = c("pseudotime", "Age", "subclass"),
  variable.name = "module",
  value.name = "Eigengene"
)
MEs_long$Days[MEs_long$Age == "E76"] = 76
MEs_long$Days[MEs_long$Age == "E85"] = 85
MEs_long$Days[MEs_long$Age == "E94"] = 94
MEs_long$Days[MEs_long$Age == "E104"] = 104
MEs_long$Days[MEs_long$Age == "E109"] = 109
MEs_long$Days[MEs_long$Age == "P3"] = 118

MEs_long$module <- factor(MEs_long$module, levels = c("M1", "M2" ,"M3","M4","M5"))
MEs_long$Age <- factor(MEs_long$Age, levels = age_ref)
MEs_long$subclass <- factor(MEs_long$subclass, levels = subclass_ref)
new_colors <- c("#f47983","#1685a9","#21a675","#725e82","#b0a4e3")

ggplot(MEs_long, aes(x = pseudotime, y = Eigengene, color = module, fill = module)) +
  geom_smooth(se = TRUE, alpha = 0.7 ,level = 0.95, method = "gam",formula = y ~ s(x, k = 5), linewidth = 1.5) +
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    linewidth = 1.5,
    color = "gray40"
  ) +
  labs(
    title = "",
    x = "Pseudotime",
    y = "Eigengene"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  )+
  scale_color_manual(values = new_colors) +
  scale_fill_manual(values = new_colors)

ggplot(MEs_long, aes(x = Days, y = Eigengene, color = module, fill = module)) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = T,
    alpha = 0.6,
    linewidth = 2
  ) +
  geom_vline(
    xintercept = 115,
    linetype = "solid",
    linewidth = 1,
    color = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 1.5,
    color = "gray40"
  ) +
  labs(
    title = "",
    x = "Age",
    y = "Eigengene"
  ) +
  scale_x_continuous(
    breaks = unique(MEs_long$Days),
    labels = unique(MEs_long$Age)
  ) +
  scale_color_manual(values = new_colors) +
  scale_fill_manual(values = new_colors) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  )

species_cols <- c(
  "Human" = "#007947",
  "Macaque" = "#f173ac",
  "Marmoset" = "#009ad6",
  "Pig" = "#ed1941",
  "Mouse" = "#f58220"
)

newtype_cols <- c(
  "mExN1" = "#FFB6C1", 
  "L2_3 IT" = "#FF4757", 
  "L3_5 IT" = "#E84148", 
  "L5_6 IT" = "#C23616",     
  
  "L5 ET" = "#6A5ACD",       
  "L5_6 NP" = "#5D478B",    
  "L6 CT" = "#4B0082",       
  "L6 B" = "#2E0854"         
)

species_ord <- c("Human","Macaque","Marmoset", "Pig", "Mouse")
# age_ord <- c(
#   "Hum-GW27","Hum-GW30","Hum-GW34","Hum-GW36","Hum-GW37","Hum-GW37+5d","Mac-E127","Mac-E147","Mac-E155","Mar-E135","Mar-P0","Pig-E76","Pig-E85",
#   "Pig-E94","Pig-E104","Pig-E109","Pig-P3","Mou-P4","Mou-P7"
# )
newtype_ord <- c("mExN1","L2_3 IT","L3_5 IT","L5_6 IT","L5 ET","L5_6 NP","L6 CT","L6 B")

gexpr_per <- readRDS("/Project_PBD/SpeciesIntegration/Data/Integration2/gexpr_per.ori.rds")
DefaultAssay(gexpr_per) <- "RNA"
gexpr_per@assays$RNA@data <- gexpr_per@assays$RNA@counts
gexpr_per <- NormalizeData(gexpr_per)
gexpr_per <- ScaleData(gexpr_per, features=VariableFeatures(gexpr_per))

gexpr_per@meta.data$species_ord <- gexpr_per$species
gexpr_per$species_ord <- factor(gexpr_per$species,levels = c("Human","Macaque","Marmoset","Pig","Mouse"))
gexpr_per@meta.data$newtype_ord <- gexpr_per$newtype
gexpr_per$newtype_ord <- factor(gexpr_per$newtype_ord,levels = newtype_ord)

DefaultAssay(gexpr_per) <- "RNA"

umap_coords <- read.csv("/Project_PBD/SpeciesIntegration/File/V2_per_scvi_umap_7_70.csv", 
                        row.names = 1)
gexpr_per[["scvi"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_coords),
  key = "scvi_",
  assay = "RNA"
)


DimPlot(gexpr_per, reduction = "scvi",
        group.by = "newtype",
        split.by = "species_ord",
        label =T , label.size = 5,raster = F,pt.size = 0.2,shuffle = T,
        cols = newtype_cols
) +
  theme(text = element_text(size = 12 , face = "bold")) +
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = arrow(length = unit(0.2,'inches'), type = "closed"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))+
  labs(title = '')

umap_coords <- read.csv("/Project_PBD/SpeciesIntegration/File/V2_per_scvi_umap_7_70.csv", 
                        row.names = 1)
gexpr_per[["scvi"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_coords),
  key = "scvi_",
  assay = "RNA"
)

DimPlot(gexpr_per, reduction = "scvi",
        group.by = "newtype",
        split.by = "species_ord",
        label =T , label.size = 5,raster = F,pt.size = 0.4,shuffle = T,
        cols = newtype_cols
) +
  theme(text = element_text(size = 12 , face = "bold")) +
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = arrow(length = unit(0.2,'inches'), type = "closed"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))+
  labs(title = '')


DefaultAssay(gexpr_per) <- "RNA"
cds <- run_monocle3(gexpr_per,use_partition = T,umap = "scvi")
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node3(cds,genes = c("UNC5D"),newtype = "mExN1"))
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 2.5,
           show_trajectory_graph = F,
           # trajectory_graph_color = "#ffc20e",
           label_roots = T) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

cds_sub = cds
cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

pseudotime_values <- cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
cell_groups <- colData(cds_sub)$newtype
non_target_cells <- !cell_groups %in% c("mExN1")
pseudotime_values[non_target_cells] <- Inf
cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] <- pseudotime_values
head(pseudotime_values)

plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 2.5,
           trajectory_graph_color = "#ffc20e",
           label_roots = T) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

pseudotime_values <- cds_sub@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
cell_groups <- colData(cds_sub)$newtype
target_cells <- cell_groups %in% c("mExN1")
pseudotime_values <- pseudotime_values[target_cells]

pt <- as.data.frame(pseudotime_values)
pt <- cbind(pt, 
            colData(cds_sub)[target_cells,]$newtype,
            colData(cds_sub)[target_cells,]$species)
colnames(pt) <- c("pt_value", "newtype", "species")

pt_filtered <- pt

species_order <- c("Human", "Macaque", "Marmoset", "Pig", "Mouse")
species_cols <- c(
  "Human" = "#007947",
  "Macaque" = "#f173ac",
  "Marmoset" = "#009ad6",
  "Pig" = "#ed1941",
  "Mouse" = "#f58220"
)

pt_filtered$pseudotime_bin <- cut(pt_filtered$pt_value, 
                                  breaks = 50, 
                                  labels = FALSE,
                                  include.lowest = TRUE)

bin_counts <- pt_filtered %>%
  group_by(species, pseudotime_bin) %>%
  summarise(cell_count = n(), .groups = "drop")

bin_counts <- bin_counts %>%
  group_by(species) %>%
  mutate(normalized_count = cell_count / max(cell_count)) %>%
  ungroup()

all_combinations <- expand.grid(
  species = factor(species_order, levels = species_order),
  pseudotime_bin = 1:50
)

bin_counts_complete <- all_combinations %>%
  left_join(bin_counts, by = c("species", "pseudotime_bin")) %>%
  mutate(
    cell_count = ifelse(is.na(cell_count), 0, cell_count),
    normalized_count = ifelse(is.na(normalized_count), 0, normalized_count)
  )

get_species_gradient <- function(species, value) {
  base_color <- species_cols[species]
  adjusted_value <- value^0.7
  
  rgb_color <- col2rgb(base_color)
  colorRampPalette(c("white", base_color))(101)[round(adjusted_value * 100) + 1]
}

bin_counts_complete <- bin_counts_complete %>%
  rowwise() %>%
  mutate(
    fill_color = get_species_gradient(as.character(species), normalized_count)
  ) %>%
  ungroup()

ggplot(bin_counts_complete, aes(x = pseudotime_bin, y = species, fill = fill_color)) +
  geom_tile(color = NA) + 
  scale_fill_identity() +
  scale_y_discrete(limits = rev(species_order)) +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
  ) +
  guides(fill = guide_colorbar(title = "")) +
  scale_x_continuous(expand = c(0, 0))

xx <- load("/Project_PBD/SpeciesIntegration/File/V2/species_specific_newtype_DEGs.RData")
xx
filter_deg <- function(deg_df,n = 1, m = 0) {
  deg_df %>%
    mutate(diff_pct = pct.1 - pct.2) %>%
    filter(p_val_adj < 0.05,
           avg_log2FC > log2(n)
           # ,diff_pct > m
    )
}
hum_deg <- filter_deg(all_species_degs[["Human"]])
# hum_deg <- all_species_degs[["Human"]]
hum_deg$cluster <- paste0("Human-",hum_deg$cluster)
hum_deg <- hum_deg[hum_deg$cluster == "Human-mExN1",]
table(hum_deg$cluster)

# pig_deg <- filter_deg(all_species_degs[["Pig"]])
pig_deg <- all_species_degs[["Pig"]]
pig_deg$cluster <- paste0("Pig-",pig_deg$cluster)
pig_deg <- pig_deg[pig_deg$cluster == "Pig-mExN1",]
table(pig_deg$cluster)

# mar_deg <- filter_deg(all_species_degs[["Marmoset"]])
mar_deg <- all_species_degs[["Marmoset"]]
mar_deg$cluster <- paste0("Marmoset-",mar_deg$cluster)
mar_deg <- mar_deg[mar_deg$cluster == "Marmoset-mExN1",]
table(mar_deg$cluster)

# mou_deg <- filter_deg(all_species_degs[["Mouse"]])
mou_deg <- all_species_degs[["Mouse"]]
mou_deg$cluster <- paste0("Mouse-",mou_deg$cluster)
mou_deg <- mou_deg[mou_deg$cluster == "Mouse-mExN1",]
table(mou_deg$cluster)

# mac_deg <- filter_deg(all_species_degs[["Macaque"]])
mac_deg <- all_species_degs[["Macaque"]]
mac_deg$cluster <- paste0("Macaque-",mac_deg$cluster)
mac_deg <- mac_deg[mac_deg$cluster == "Macaque-mExN1",]
table(mac_deg$cluster)

species_conserved <- Reduce(intersect,list(hum_deg$gene,mar_deg$gene,mac_deg$gene,pig_deg$gene,mou_deg$gene))
print(species_conserved)

slc.g <- species_conserved

subset_cells <- subset(gexpr_per, subset = newtype %in% c("mExN1"))
avg_expr_list <- list()
for (sp in unique(subset_cells$species)) {
  sp_cells <- subset(subset_cells, subset = species == sp)
  avg_expr_list[[sp]] <- rowMeans(GetAssayData(sp_cells, slot = "data")[slc.g,])
}
avg_expr_matrix <- do.call(cbind, avg_expr_list)
cor_matrix <- cor(avg_expr_matrix, method = "pearson")

col_fun <- colorRamp2(
  breaks = c(min(cor_matrix), mean(cor_matrix), max(cor_matrix)),
  colors = rev(brewer.pal(11, "RdBu"))[c(3,6,9)]
)

Heatmap(
  cor_matrix, col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE,
  name = "Similarity", 
  heatmap_height = unit(10, "cm"), 
  heatmap_width = unit(10, "cm"),
  show_row_names = T, 
  show_column_names = T,
  border = T,
  border_gp = gpar(col = "black")
)

cell_types <- c("IT_ExN", "L5_6 NP", "L6 CT", "L6 B")
for (ct in cell_types) {
  cat("Processing cell type:", ct, "\n")
  
  if (ct == "IT_ExN") {
    subset_cells <- subset(gexpr_per, subset = newtype %in% c("L2_3 IT", "L3_5 IT", "L5_6 IT"))
    
    conserved_genes_list <- list()
    for (sp in c("Human", "Pig", "Marmoset", "Mouse", "Macaque")) {
      it_genes <- c()
      for (it_ct in c("L2_3 IT", "L3_5 IT", "L5_6 IT")) {
        deg_df <- all_species_degs[[sp]]
        deg_df <- deg_df[deg_df$cluster == it_ct, ]
        it_genes <- union(it_genes, deg_df$gene)
      }
      conserved_genes_list[[sp]] <- it_genes
    }
    species_conserved <- Reduce(intersect, conserved_genes_list)
    
  } else {
    subset_cells <- subset(gexpr_per, subset = newtype == ct)
    
    conserved_genes <- list()
    for (sp in c("Human", "Pig", "Marmoset", "Mouse", "Macaque")) {
      deg_df <- all_species_degs[[sp]]
      deg_df <- deg_df[deg_df$cluster == ct, ]
      conserved_genes[[sp]] <- deg_df$gene
    }
    species_conserved <- Reduce(intersect, conserved_genes)
  }
  
  avg_expr_list <- list()
  for (sp in unique(subset_cells$species)) {
    sp_cells <- subset(subset_cells, subset = species == sp)
    avg_expr_list[[sp]] <- rowMeans(GetAssayData(sp_cells, slot = "data")[species_conserved, ])
  }
  avg_expr_matrix <- do.call(cbind, avg_expr_list)
  cor_matrix <- cor(avg_expr_matrix, method = "pearson")
  
  col_fun <- colorRamp2(
    breaks = c(min(cor_matrix), mean(cor_matrix), max(cor_matrix)),
    colors = rev(brewer.pal(11, "RdBu"))[c(3,6,9)]
  )
  
  draw(Heatmap(
    cor_matrix, col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE,
    name = ct, 
    heatmap_height = unit(10, "cm"), 
    heatmap_width = unit(10, "cm"),
    show_row_names = T, 
    show_column_names = T,
    border = T,
    border_gp = gpar(col = "black")
  ))
}
