```{r}

# load libraries
library(tidyverse)
library(data.table)
library(Matrix)
library(Rfast)
library(matrixStats)
library(reticulate)
library(anndata)
library(future.apply)
library(gtools)
library(ggridges)
library(scales)
library(ComplexHeatmap)
library(forcats)
library(igraph)
library(mclust)
library(UpSetR)
library(patchwork)
library(circlize)
library(imager)

# source functions
source("/inkwell05/ameer/functions/0_source_functions.R")

```


```{r}

# speed and RAM benchmarking
benchmark = as.data.frame(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/speed_memory_benchmark.csv"))
benchmark = benchmark[!(benchmark$ram_MB == -1 | benchmark$time_sec == -1),]

# plot together
result = benchmark
result$ram_MB   = result$ram_MB/1024
result$time_sec = result$time_sec/60
colnames(result)[3:4] = c("Memory","Time")
result = pivot_longer(data = result, cols = 3:4, names_to = "type", values_to = "vals") |> as.data.frame()

# time
X = result[result$type == "Time",]
ggplot(X) + aes(x = cells, y = vals, color = factor(genes)) + geom_point() + geom_line() + theme_classic() +
  scale_x_continuous(breaks = benchmark$cells, labels = label_number(scale = 1e-3, suffix = "k"), transform = "log10") +
  labs(x = "Number of cells", y = "Time (minutes)", color = "Gene panel", linetype = "") +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 15), legend.text = element_text(size = 12),
        plot.subtitle = element_text(size = 15), legend.title = element_text(size = 15))

# memory
X = result[result$type == "Memory",]
ggplot(X) + aes(x = cells, y = vals, color = factor(genes)) + geom_point() + geom_line() + theme_classic() +
  scale_x_continuous(breaks = benchmark$cells, labels = label_number(scale = 1e-3, suffix = "k"), transform = "log10") +
  labs(x = "Number of cells", y = "Peak memory usage (GB)", color = "Gene panel", linetype = "") +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 15), legend.text = element_text(size = 12),
        plot.subtitle = element_text(size = 15), legend.title = element_text(size = 15))

```


```{r}

# total number of cells and genes in each dataset
data = as.data.frame(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_total_cells_genes.csv"))

# change names
#data$data[data$data == "Zhuang_2023_MERFISH"] = "Zhuang_2023a_MERFISH"
#data$data[data$data == "Zhuang_2024_MERFISH"] = "Zhuang_2023b_MERFISH"

# cells vs datasets
ggplot(data, aes(x = reorder(data, cells), y = cells)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  scale_y_log10(
    breaks = c(10, 100, 1e3, 1e4, 1e5, 1e6, 1e7),
    labels = label_number(scale_cut = cut_short_scale())) +
  coord_flip() +
  geom_hline(yintercept = median(data$cells), color = "brown3", linetype = "dashed") +
  labs(x = "", y = "Number of cells (post-QC)",
       subtitle = str_c("Median (dashed line) = ", comma(median(data$cells)), " cells")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), plot.subtitle = element_text(size = 15))

ggsave("/inkwell05/images/cell_counts.svg", device = "svg", dpi = 600, units = "in", height = 8, width = 8)

# unique genes as a function of number of samples included (100 iterations)
# represent this as one figure (with axes as proportions)
output = as.data.frame(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_predictions_vs_slices_cross_validation.csv"))
output = output[!str_detect(output$data, "Zhuang_2024_Mouse_MERFISH"),] # remove dataset due to age vs young samples
output$data = str_remove(output$data, "_Mouse")

for (i in 1:length(unique(output$data))){output$slices_included[output$data %in% unique(output$data)[i]] = scale_01(output$slices_included[output$data %in% unique(output$data)[i]])} # normalize slices included b/w 0-1

# remove datasets with fewer than N slices
to_remove = table(output$data) < 15
output = output[!(output$data %in% names(to_remove)[to_remove]),]
output$data = str_remove(output$data, "Mouse_")

# assign colors by dataset source, not by brain
source_data = vector(mode = "character", length = nrow(output))
source_data[output$data %in% "Gillis_Unpublished_BARseq"] = "Gillis_Unpublished_BARseq"
source_data[output$data %in% "Wang_2023_STARmap"] = "Wang_2023_STARmap"
source_data[str_detect(output$data, "Zador_2024_BARseq")] = "Zador_2024_BARseq"
source_data[output$data %in% "Zeng_2023_MERSCOPE"] = "Zeng_2023_MERSCOPE"
source_data[output$data %in% "Zeng_2023_MERSCOPE"] = "Zeng_2023_MERSCOPE"
source_data[str_detect(output$data, "Zhuang_2023_MERFISH")] = "Zhuang_2023_MERFISH"
output = data.frame(source_data, output)

# change names
output$source_data = "Zhuang_2023a_MERFISH"
output$data = str_replace_all(output$data, "Zhuang_2023", "Zhuang_2023a")

# make plot
ggplot(output) + aes(x = slices_included, y = auroc, color = data) + geom_point() +
  theme_minimal() + geom_smooth(method = "loess") +
  guides(color = guide_legend(override.aes = list(size = 3)), shape = "none") +
  theme(legend.position = "right", legend.text = element_text(size = 10), axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  labs(x = "Proportion of samples (slices) used for prediction\n(minimum samples required = 15)", y = "AUROC\n(prediction of held-out sample's cross-expression)", color = "")
ggsave("/inkwell05/images/pvalues_correlation_predictions.png", device = "png", width = 6, height = 4, dpi = 300)

# cross-expression similarity vs distance between slices across datasets
target = as.data.frame(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_similarity_vs_distance.csv"))
target = target[!is.na(target$distance),]
target$brain = str_remove(target$brain, "_Mouse")
target = target[!str_detect(target$brain, "sagittal_3"),] # remove sagittal 3 brain due to lack of samples

for (i in 1:length(unique(target$brain))){
  #target$distance[target$brain %in% unique(target$brain)[i]]   = scale_01(target$distance[target$brain %in% unique(target$brain)[i]])
  #target$same_genes[target$brain %in% unique(target$brain)[i]] = scale_01(target$same_genes[target$brain %in% unique(target$brain)[i]])
} # normalize distance and same genes b/w brains to 0 and 1

# correlation
temp_target = target[, c(1:2,5)]
temp_target = temp_target[!(is.na(temp_target[,3]) | is.nan(temp_target[,3])),]

x = temp_target |> group_by(brain) |> summarise(corr = cor(distance, correlation, method = "spearman")) |> mutate(label = paste0(brain, "\nSpearman's rho = ", round(corr, 2))) |> as.data.frame()
temp_target = left_join(temp_target, x, by = "brain")

ggplot(temp_target, aes(x = distance, y = correlation, color = brain)) +
  geom_point(size = 0) + facet_wrap(~ label, scales = "free") + theme_minimal() +
  labs(x = "Distance between slices (separation by slice order)", y = "Correlation between slices' cross-expression profiles") +
  theme(legend.position = "none")

ggplot(temp_target, aes(x = distance, y = correlation, color = brain)) +
  geom_point(size = 0) + theme_minimal() +
  labs(x = "Distance between slices (separation by slice order)", y = "Correlation between slices' cross-expression profiles", color = "",
       subtitle = str_c("Spearman's rho = ", signif(cor(temp_target$distance, temp_target$correlation, method = "spearman"), digits = 2))) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "right")

# same genes
temp_target = target[, c(1:2,4)]
temp_target = temp_target[!(is.na(temp_target[,3]) | is.nan(temp_target[,3])),]

ggplot(temp_target, aes(x = distance, y = same_genes, color = brain)) +
  geom_point(size = 0) + facet_wrap(~brain, scales = "free") + theme_minimal() +
  labs(x = "Distance between slices (in terms of slice order)", y = "Number of cross-expressing genes between slices") +
  theme(legend.position = "none")

```


```{r}

# overlap in gene panels
overlap = as.matrix(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_gene_panel_overlap.csv"))
rownames(overlap) = colnames(overlap)
overlap[lower.tri(overlap)] = NA  # Keep only upper triangle

# store diagonal values separately and temporarily make them the minimum value
diagonal_values = diag(overlap)
min_off_diag = min(overlap, na.rm = TRUE)
diag(overlap) = min_off_diag

# define color breaks for off-diagonal elements only
breaks = seq(min_off_diag, max(overlap, na.rm = TRUE), length.out = 50)
color_palette = colorRampPalette(c("yellow", "red"))(50)

# manually insert diagonal values back into display_numbers
display_matrix = formatC(overlap, format = "f", digits = 0)
diag(display_matrix) = formatC(diagonal_values, format = "f", digits = 0)  # Restore correct diagonal numbers

# plot heatmap
pheatmap(overlap, 
         display_numbers = display_matrix,  # manually insert correct diagonal values
         fontsize_number = 10,
         fontsize_row = 12,
         fontsize_col = 12,
         cluster_rows = FALSE, cluster_cols = FALSE,
         angle_col = "315", legend = FALSE,
         number_color = "blue",
         color = color_palette,
         breaks = breaks,  # off-diagonal scaling only
         main = "Number of overlapping genes between panels", fontsize = 9)

# process data to create an upset plot

# list of background genes

dirr = "/inkwell05/ameer/databases/spatial/"
dirr = str_c(dirr, list.files(dirr))

background_genes = c()

for (i in 1:length(dirr)){
  
  temp_dirr = str_c(dirr[i], "/", list.files(dirr[i]))
  temp_dirr = temp_dirr[str_detect(temp_dirr, "expression")]
  temp_dirr = str_c(temp_dirr, "/", list.files(temp_dirr))
  exp_dirr  = temp_dirr[1]
  
  temp_dirr = str_c(dirr[i], "/", list.files(dirr[i]))
  temp_dirr = temp_dirr[str_detect(temp_dirr, "metadata")]
  temp_dirr = str_c(temp_dirr, "/", list.files(temp_dirr))
  md_dirr   = temp_dirr[1]
  
  temp_dirr = spatial_QC(path_to_expression = exp_dirr, path_to_metadata = md_dirr)
  temp_dirr = colnames(temp_dirr$data)
  
  background_genes = c(background_genes, temp_dirr)
}

# find gene symbols
master_list = as.data.frame(fread("/inkwell05/ameer/databases/scRNAseq/Zeng_2023_Mouse_10x/metadata/gene_master_list.csv"))

x = master_list[match(background_genes, master_list$gene_identifier), "gene_symbol"]
y = master_list[match(background_genes, master_list$gene_symbol), "gene_symbol"]

x[is.na(x)] = y[is.na(x)]
y[is.na(y)] = x[is.na(y)]

background_genes = unique(c(x,y))

# assess whether dataset-specific panels include background genes
hits_background = matrix(data = 0, nrow = length(background_genes), ncol = length(dirr))
rownames(hits_background) = background_genes
colnames(hits_background) = str_remove(dirr, "/inkwell05/ameer/databases/spatial/")

for (i in 1:ncol(hits_background)){
  
  temp_dirr = str_c(dirr[i], "/", list.files(dirr[i]))
  temp_dirr = temp_dirr[str_detect(temp_dirr, "expression")]
  temp_dirr = str_c(temp_dirr, "/", list.files(temp_dirr))
  exp_dirr  = temp_dirr[1]
  
  temp_dirr = str_c(dirr[i], "/", list.files(dirr[i]))
  temp_dirr = temp_dirr[str_detect(temp_dirr, "metadata")]
  temp_dirr = str_c(temp_dirr, "/", list.files(temp_dirr))
  md_dirr   = temp_dirr[1]
  
  temp_dirr = spatial_QC(path_to_expression = exp_dirr, path_to_metadata = md_dirr)
  temp_dirr = colnames(temp_dirr$data)
  
  x = master_list[match(temp_dirr, master_list$gene_identifier), "gene_symbol"]
  y = master_list[match(temp_dirr, master_list$gene_symbol), "gene_symbol"]

  x[is.na(x)] = y[is.na(x)]
  y[is.na(y)] = x[is.na(y)]
  
  temp_dirr = unique(c(x,y))
  
  hits_background[,i] = as.numeric(background_genes %in% temp_dirr)
}

# create upset plot with intersection
hits_background = as.data.frame(hits_background)
colnames(hits_background) = str_remove(colnames(hits_background), "Mouse_")

upset(hits_background, nintersects = 25,
      sets = colnames(hits_background),
      order.by = "freq",
      text.scale = c(0, 0, 1.2, 1, 0.8, 1),
      sets.x.label = "Gene panel size",
      set_size.show = FALSE)

# create upset plot without self-intersection
hits_background = as.data.frame(hits_background)
colnames(hits_background) = str_remove(colnames(hits_background), "Mouse_")

wo_intersection = hits_background[rowSums(hits_background) > 1,]

upset(wo_intersection, nintersects = 25,
      sets = colnames(hits_background),
      order.by = "freq",
      text.scale = c(0, 0, 1.2, 1, 0.8, 1),
      sets.x.label = "Gene panel size",
      set_size.show = FALSE)

# create gene panel size (w/ same order as upset plot)
panel_size = data.frame(data = names(diag(overlap)), size = as.numeric(diag(display_matrix)))
correct_order = c("Resolve_2021_Molecular_Cartography", "Zador_2024_BARseq", "Gillis_Unpublished_BARseq", "10x_replicates4_2023_Xenium", "10x_AD-controls_2023_Xenium", "Zhuang_2024_MERFISH", 
                  "Linnarsson_2022_EEL_FISH", "Vizgen_2022_MERSCOPE", "Zeng_2023_MERSCOPE", "NanoString_CosMx", "Wang_2023_STARmap", "Zhuang_2023_MERFISH", "10x_5K_panel_2024_Xenium")

panel_size$data = factor(panel_size$data, levels = correct_order)
panel_size = panel_size[order(panel_size$data), ]

ggplot(panel_size, aes(x = size, y = data)) +
  geom_bar(stat = "identity", fill = "#333333", width = 0.6) +
  labs(x = "Gene panel size", y = NULL) +
  theme_classic() +
  scale_x_reverse() +
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(hjust = 1),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.y.right = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15))

# upset plots with (non-unique) intersections
colnames(hits_background)[colnames(hits_background) == "Zhuang_2023_MERFISH"] = "Zhuang_2023a_MERFISH"
colnames(hits_background)[colnames(hits_background) == "Zhuang_2024_MERFISH"] = "Zhuang_2023b_MERFISH"

comb_mat = make_comb_mat(hits_background, mode = "intersect")
comb_mat = comb_mat[, str_count(comb_name(comb_mat), "1") != 1 & str_detect(comb_name(comb_mat), "^[01]+$")] # remove self-intersections

top_n = 50  # number of intersections to show
comb_mat = comb_mat[, order(comb_size(comb_mat), decreasing = TRUE)[1:top_n]]

anno_left  = rowAnnotation(set_name = anno_text(set_name(comb_mat), location = 0.5, just = "center", gp = gpar(fontsize = 8)))
anno_right = rowAnnotation("Genes" = anno_barplot(set_size(comb_mat), border = FALSE, gp = gpar(fill = "black")))
anno_top   = HeatmapAnnotation("Gene overap\n(non-unique)" = anno_barplot(comb_size(comb_mat), ylim = c(0, max(comb_size(comb_mat))*1.1), border = FALSE, gp = gpar(fill = "black"), height = unit(4, "cm")),
                               annotation_name_side = "left", annotation_name_rot = 90)

UpSet(comb_mat, 
      comb_order = order(comb_size(comb_mat), decreasing = TRUE),
      top_annotation = anno_top,
      right_annotation = anno_right,
      left_annotation = anno_left,
      show_row_names = FALSE)

# save upset plot
svg("/inkwell05/images/gene_upset.svg", width = 9, height = 5)
UpSet(comb_mat, 
      comb_order = order(comb_size(comb_mat), decreasing = TRUE),
      top_annotation = anno_top,
      right_annotation = anno_right,
      left_annotation = anno_left,
      show_row_names = FALSE)
dev.off()

# median number of intersecting genes
print(str_c("Medium number of intersecting genes between a pair of datasets = ", median(upper_tri(overlap))))

```


```{r}

# replicability (correlation used to predict significant genes -- AUROC) within and between datasets
output = as.data.frame(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_replicability_across_samples_and_datasets.csv"))

result_corr = matrix(data = 0, nrow = length(unique(output$data1)), ncol = length(unique(output$data2)),
                     dimnames = list(str_remove(unique(output$data1), "Mouse_"), str_remove(unique(output$data1), "Mouse_")))

for (i in 1:nrow(result_corr)){
  for (j in 1:ncol(result_corr)){
    x = output[(output$data1 %in% unique(output$data1)[i]) & (output$data2 %in% unique(output$data2)[j]),]
    x = x$auroc_corr
    x = x[!is.na(x)]
    result_corr[i,j] = median(x)
  }
}

# change names
rownames(result_corr)[rownames(result_corr) == "Zhuang_2023_MERFISH"] = "Zhuang_2023a_MERFISH"
colnames(result_corr) = rownames(result_corr)

rownames(result_corr)[rownames(result_corr) == "Zhuang_2024_MERFISH"] = "Zhuang_2023b_MERFISH"
colnames(result_corr) = rownames(result_corr)

# heatmap with max. AUROC
output = result_corr
#output = (t(output) + output)/2
#output[lower.tri(output)] = NA
#output[is.nan(output)] = NA
#dimnames(output) = dimnames(result_corr)

# plot heatmap
pheatmap(output, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_color = "blue", color = colorRampPalette(c("yellow", "red"))(50), legend = FALSE,
         main = "Median AUROC", fontsize = 9, angle_col = "315", fontsize_number = 10, fontsize_row = 10, fontsize_col = 10)

# save heatmap
svg("/inkwell05/images/AUROC_heatmap.svg", width = 10, height = 6)
pheatmap(output, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_color = "blue", color = colorRampPalette(c("yellow", "red"))(50), legend = FALSE,
         main = "", fontsize = 9, angle_col = "315", fontsize_number = 12, fontsize_row = 10, fontsize_col = 10)
dev.off()

```


```{r}

# union meta-analytic network: genes significant b/w all pairs of studies, including those appearing in a single study
# intersection meta-analytic network: genes significant b/w pairs of studies with good AUROC (above median), excluding those appearing in a single study
network = "union"
network = "intersection"

# load significant genes
sig_genes = as.data.frame(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_significant_genes_across_samples_and_datasets.csv"))
sig_genes$data1 = str_remove(sig_genes$data1, "Mouse_")
sig_genes$data2 = str_remove(sig_genes$data2, "Mouse_")

if (network == "intersection"){
  
  # remove genes appearing in a single study
  sig_genes = sig_genes[sig_genes$data1 != sig_genes$data2,]
  
  # keep pairs of studies with good AUROC (median = approx. 0.7)
  study_pairs = output; study_pairs[is.na(study_pairs)] = 0
  study_pairs = study_pairs >= median(as.vector(output), na.rm = TRUE)
  threshold   = vector(mode = "numeric", length = nrow(sig_genes))
  
  for (i in 1:nrow(study_pairs)){for (j in 1:ncol(study_pairs)){threshold[(sig_genes$data1 %in% rownames(study_pairs)[i] | sig_genes$data2 %in% rownames(study_pairs)[i]) & (sig_genes$data1 %in% colnames(study_pairs)[j] | sig_genes$data2 %in% colnames(study_pairs)[j])] = as.numeric(study_pairs[i,j])}}
  
  sig_genes = sig_genes[as.logical(threshold),]
}

# consider both directions and remove duplicated genes
pair1 = str_c(sig_genes$gene1_ensembl, "_", sig_genes$gene2_ensembl)
pair2 = str_c(sig_genes$gene2_ensembl, "_", sig_genes$gene1_ensembl)
sig_genes$pair = pmin(pair1, pair2)
sig_genes = sig_genes[!(sig_genes$gene1_symbol == "" | sig_genes$gene2_symbol == ""),]

sig_genes = sig_genes[!duplicated(sig_genes$pair), ]
sig_genes$pair = NULL

# remove dataset information
sig_genes$data1 = NULL; sig_genes$data2 = NULL

# node degree
node = table(c(sig_genes$gene1_symbol, sig_genes$gene2_symbol))
node = data.frame(node = names(node), degree = as.numeric(node))

# initialize node community membership
graph = sig_genes[,c(3,4)]
graph = graph_from_data_frame(graph, directed = FALSE)
graph = igraph::simplify(graph, remove.loops = TRUE)

# Leiden clustering nodes multiple times
plan(multisession, workers = parallel::detectCores() - 1)
iter = 100

cluster_list = future_lapply(1:iter, \(i)
  as.numeric(membership(cluster_leiden(graph, objective_function = "modularity", n_iterations = 100))), future.seed = TRUE)

clusters = do.call(cbind, cluster_list)
rownames(clusters) = names(membership(cluster_leiden(graph, objective_function = "modularity", n_iterations = 1)))
colnames(clusters) = 1:iter

# compute adjusted rand indices b/w various partitions
ARI = matrix(data = 0, nrow = iter, ncol = iter)
for (i in 1:iter){for (j in 1:iter){ARI[i,j] = adjustedRandIndex(clusters[,i], clusters[,j])}}
clusters = clusters[, sample(which(colMeans(ARI) == max(colMeans(ARI))), size = 1)]

node$community = as.numeric(clusters[match(node$node, names(clusters))])
node = node[order(node$community),]

```


```{r}

# background genes (overlap of gene panels across all datasets) for GO enrichment analysis
dirr = "/inkwell05/ameer/databases/spatial/"
dirr = str_c(dirr, list.files(dirr))

background_genes = c()

for (i in 1:length(dirr)){
  
  temp_dirr = str_c(dirr[i], "/", list.files(dirr[i]))
  temp_dirr = temp_dirr[str_detect(temp_dirr, "expression")]
  temp_dirr = str_c(temp_dirr, "/", list.files(temp_dirr))
  exp_dirr  = temp_dirr[1]
  
  temp_dirr = str_c(dirr[i], "/", list.files(dirr[i]))
  temp_dirr = temp_dirr[str_detect(temp_dirr, "metadata")]
  temp_dirr = str_c(temp_dirr, "/", list.files(temp_dirr))
  md_dirr   = temp_dirr[1]
  
  temp_dirr = spatial_QC(path_to_expression = exp_dirr, path_to_metadata = md_dirr)
  temp_dirr = colnames(temp_dirr$data)
  
  background_genes = c(background_genes, temp_dirr)
}

# find gene symbols
master_list = as.data.frame(fread("/inkwell05/ameer/databases/scRNAseq/Zeng_2023_Mouse_10x/metadata/gene_master_list.csv"))

x = master_list[match(background_genes, master_list$gene_identifier), "gene_symbol"]
y = master_list[match(background_genes, master_list$gene_symbol), "gene_symbol"]

x[is.na(x)] = y[is.na(x)]
y[is.na(y)] = x[is.na(y)]

background_genes = unique(c(x,y))

```


```{r}

# GO enrichment label of each node's community

# GO enrichment using union of gene panels as background and genes in each community as target
plan(multisession, workers = parallel::detectCores() - 1)

run_GO_enrichment = function(i) {
  community_id = unique(node$community)[i]
  test_genes = node$node[node$community == community_id]
  
  GO = GO_enrichment(test_genes = test_genes, background_genes = background_genes)
  
  GO_pvalues = GO$p_values[as.logical(GO$p_values$sig_GO),]
  GO_pvalues = GO_pvalues[order(GO_pvalues$pvals_FDR),]
  
  if (nrow(GO_pvalues) == 0) return(NULL)  # skip empty results
  
  return(data.frame(GO_pvalues, community = community_id))
}

results_list = future_lapply(1:length(unique(node$community)), run_GO_enrichment)
results_GO = do.call(rbind, results_list)
plan(sequential)

# load semantic similarity b/w GO terms of each community from REVIGO's outputs
annotations = as.data.frame(fread("/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_meta_analytic_network_GO_annotations.csv"))

# choose the best annotations for each aspect (biological process, cellular component, and molecular function) per community
annot = data.frame(community = c(1, 2, 3, 4, 5, 6, 7, 11),
                   BP = c("Cell adhesion (junction)/ Immune signaling",
                          "Extracellular (proteoglycan) and vascular structure development",
                          "G protein coupled receptor signaling (small molecules: nitric oxide, amines, hormones, peptides)",
                          "Synapse organization",
                          "Neurotransmitter/ Neuropeptide regulation and signaling",
                          "Homeostasis: circulatory, renal, endocrine, respiration, metabolism, muscle contraction",
                          "CNS and mesenchyme development and differentiation",
                          "Amyloid precursor protein metabolic process"),
                   CC = c("Cell adhesion/ Integrin complex",
                          "External encapsulating structure",
                          "Receptor complex/ Synaptic membrane",
                          "Postsynaptic membrane",
                          "Presynapse and transport vesicle",
                          "Receptor complex/ Membrane microdomain",
                          "Presynapse",
                          "External encapsulating structure/ Extracellular matrix"),
                   MF = c("Cytokine binding/ Cargo receptor",
                          "Integrin binding/ Glycan binding",
                          "G protein coupled receptor activity",
                          NA,
                          "Neuropeptide receptor binding and activity",
                          "G protein coupled receptor binding and activity",
                          "Transcription factor activity",
                          "Amyloid-beta binding"))

# append GO information of each community
node_annot = merge(node, annot, by = "community", all.x = TRUE)

# save nodes and edges
fwrite(node_annot, str_c("/inkwell05/images/network_", network, "_nodes.csv"))
fwrite(sig_genes[,3:4], str_c("/inkwell05/images/network_", network, "_edges.csv"))

```


```{r}

# cross-expression example 1: Cst3 and Itm2b genes in Alzheimer's disease (from 10x_AD-controls_2023_Mouse_Xenium)

# load expression matrices and metadata
dirr = "/inkwell05/ameer/databases/spatial/10x_AD-controls_2023_Mouse_Xenium/"
dirr = str_c(dirr, list.files(dirr))

metadata = dirr[str_detect(dirr, "metadata")]
metadata = str_c(metadata, "/", list.files(metadata))
metadata = mixedsort(metadata)

data = dirr[str_detect(dirr, "expression")]
data = str_c(data, "/", list.files(data))
data = mixedsort(data)

nam  = str_remove(data, "/inkwell05/ameer/databases/spatial/10x_AD-controls_2023_Mouse_Xenium/expression/10x_AD-controls_2023_Mouse_Xenium_brain_")
nam  = str_remove(nam, "_slice_1_expression.h5ad")

df = vector(mode = "list", length = length(data))
md = vector(mode = "list", length = length(data))
names(df) = nam; names(md) = nam

master_list = as.data.frame(fread("/inkwell05/ameer/databases/scRNAseq/Zeng_2023_Mouse_10x/metadata/gene_master_list.csv"))

for (i in 1:length(df)){
  
  x = spatial_QC(path_to_expression = data[i], path_to_metadata = metadata[i])
  df[[i]] = x$data
  md[[i]] = x$metadata
  
  x = master_list[match(colnames(df[[i]]), master_list$gene_identifier), "gene_symbol"]
  y = master_list[match(colnames(df[[i]]), master_list$gene_symbol), "gene_symbol"]

  x[is.na(x)] = y[is.na(x)]
  y[is.na(y)] = x[is.na(y)]

  colnames(df[[i]]) = unique(c(x,y))
}

data = df; metadata = md; remove(df, md)

# we looked at the genes in community 11 (Alzheimer's) and found Cst3 and Itm2b interesting
# investigate genes of interest
genes  = c("Cst3","Itm2b")
counts = data.frame()

thm = theme(legend.position = "none", axis.text = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())

for (i in 1:3){
  
  # subset brains
  df1 = data[[i]][,genes]          # AD data
  md1 = metadata[[i]][,c("x","y")]
  
  df2 = data[[i + 3]][,genes]      # WT data
  md2 = metadata[[i + 3]][,c("x","y")]
  
  # cross-expression plots
  pval = signif(cross_expression(data = df1, locations = md1)$cross_pvalue, digits = 2)
  p1 = tissue_expression_plot(data = df1, locations = md1, gene1 = genes[1], gene2 = genes[2], cross_expression = TRUE) + thm + labs(subtitle = str_c(names(data)[i], ", p-value = ", pval))

  pval = signif(cross_expression(data = df2, locations = md2)$cross_pvalue, digits = 2)
  p2 = tissue_expression_plot(data = df2, locations = md2, gene1 = genes[1], gene2 = genes[2], cross_expression = TRUE) + thm + labs(subtitle = str_c(names(data)[i + 3], ", p-value = ", pval))
  
  # view both plots
  print(p1 + p2)
  
  # DE between genes: test that WT expression is larger than AD expression
  if (i == 1){age = "Young, 2.5 months"}; if (i == 2){age = "Middle, 5.7 months"}; if (i == 3){age = "Old, 13.4+ months"}
  df1 = data.frame(as.matrix(df1), age, status = "AD")
  df2 = data.frame(as.matrix(df2), age, status = "WT")
  dff = rbind(df1, df2)
  dff = as.data.frame(pivot_longer(dff, cols = 1:2, values_to = "counts", names_to = "genes"))
  
  counts = rbind(counts, dff)

  print(i)
}

# make DE plot
counts$age   = factor(counts$age,   levels = c("Young, 2.5 months", "Middle, 5.7 months", "Old, 13.4+ months"))
counts$genes = factor(counts$genes, levels = c("Itm2b", "Cst3"))

ggplot(counts, aes(x = counts, color = status, y = after_stat(scaled))) +
  geom_density(alpha = 0.5) + facet_grid(genes ~ age, scales = "free", switch = "y") +
  labs(x = "Gene counts per cell", y = "", color = "") + theme_minimal() +
  theme(legend.position = "top", axis.text.y = element_blank(), legend.text = element_text(size = 10), strip.text = element_text(size = 10))

# conclusion: Cross-expression ceases to exist in Old AD mice but is otherwise present
# this likely occurs because Cst3 inhibits proteases that clear Apoe tangles and Itm2b prevents Apoe tangle formation
# so, we get higher Cst3 (more tangles) and Itm2b not in adjacent cells (less prevention) so both cause more AD

# NOTE: only one other dataset (NanoString_Mouse_CosMx), which has just one slice, contains these two genes and cross-expression is non-significant
# so, overall, we do not have enough datasets here to robustly evaluate this claim...

```


```{r}

# we consider Parkinson's genes Drd1 and Gpr6 
# we compute their cross-expression across numerous slices ordered in anterior-to-posterior or medial-to-lateral directions
dirr = "/inkwell05/ameer/databases/spatial/"
dirr = str_c(dirr, list.files(dirr))

master_list = as.data.frame(fread("/inkwell05/ameer/databases/scRNAseq/Zeng_2023_Mouse_10x/metadata/gene_master_list.csv"))
genes = c("Drd1","Gpr6")
data  = list()
metadata = list()

# datasets that contain these gene pairs
dirr = dirr[dirr %in% c("/inkwell05/ameer/databases/spatial/Zhuang_2023_Mouse_MERFISH", "/inkwell05/ameer/databases/spatial/Vizgen_2022_Mouse_MERSCOPE")]

idx = 1 # index for storing results

for (i in 1:length(dirr)){
  
  # load specific dataset
  temp_dirr = str_c(dirr[i], "/")
  temp_dirr = str_c(temp_dirr, list.files(temp_dirr))
  
  # expression matrices
  exp_dirr = temp_dirr[str_detect(temp_dirr, "expression")]
  exp_dirr = str_c(exp_dirr, "/")
  exp_dirr = str_c(exp_dirr, list.files(exp_dirr))
  exp_dirr = mixedsort(exp_dirr)
  
  # metadata
  md_dirr = temp_dirr[str_detect(temp_dirr, "metadata")]
  md_dirr = str_c(md_dirr, "/")
  md_dirr = str_c(md_dirr, list.files(md_dirr))
  md_dirr = mixedsort(md_dirr)
  
  # check if genes of interest are present in the dataset
  z = spatial_QC(path_to_expression = exp_dirr[1], path_to_metadata = md_dirr[1])
  z = z$data
  x = master_list[match(colnames(z), master_list$gene_identifier), "gene_symbol"]
  y = master_list[match(colnames(z), master_list$gene_symbol), "gene_symbol"]
  
  x[is.na(x)] = y[is.na(x)]
  y[is.na(y)] = x[is.na(y)]
  
  colnames(z) = dplyr::coalesce(x, y)
  
  if (sum(colnames(z) %in% genes) != 2){next}
  
  # load datasets
  for (j in 1:length(exp_dirr)){
    
    # load expression matrix and metadata
    x = spatial_QC(path_to_expression = exp_dirr[j], path_to_metadata = md_dirr[j])
    
    # columns to extract
    to_extract = c("x","y","CCFx","CCFy","CCFz")
    
    # store results
    data[[idx]]     = x$data
    metadata[[idx]] = x$metadata[,colnames(x$metadata) %in% to_extract]
    
    # ensure appropriate gene symbols
    x = master_list[match(colnames(data[[idx]]), master_list$gene_identifier), "gene_symbol"]
    y = master_list[match(colnames(data[[idx]]), master_list$gene_symbol), "gene_symbol"]
    
    x[is.na(x)] = y[is.na(x)]
    y[is.na(y)] = x[is.na(y)]
    
    colnames(data[[idx]]) = dplyr::coalesce(x, y)
    
    # subset to genes of interest
    data[[idx]] = data[[idx]][, colnames(data[[idx]]) %in% genes]
    
    # name datasets
    nam = str_remove(exp_dirr, "/inkwell05/ameer/databases/spatial/")
    nam = str_remove(nam, "_expression.h5ad")
    nam = str_remove(nam, "^.*?/expression/")
    names(data)[idx] = nam[j]
    names(metadata)[idx] = nam[j]
    
    # update index
    idx = idx + 1
  }
  print(i/length(dirr))
}

# re-order Vizgen slices across brains to consider it a single coronal dataset
order_idx = c(9, 3, 6, 2, 5, 8, 7, 4, 1)

data[1:9]     = data[order_idx]
metadata[1:9] = metadata[order_idx]

names(data)[1:9]     = names(data)[order_idx]
names(metadata)[1:9] = names(metadata)[order_idx]

# re-order two sagittal brains (sagittal 3 and sagittal 23) into a single dataset using their CCF coordinates (z-direction, medial to lateral)
x  = which(str_detect(names(data), "sagittal"))
dd = data.frame()

for (i in 1:length(x)){
  if (str_detect(names(data)[x], "sagittal_23")[i]){brain = "sagittal_23"} else {brain = "sagittal_3"}
  temp_dd = data.frame(brain, slice = names(data)[x[i]], z = median(metadata[[x[i]]][,"CCFz"], na.rm = TRUE), idx = x[i])
  dd = rbind(dd, temp_dd)
}

dd$z[nrow(dd)] = 1 # last slice is not registered to assign it a dummy value
dd = dd[order(dd$z, decreasing = TRUE),]

data[x]     = data[dd$idx]
metadata[x] = metadata[dd$idx]

names(data)[x]     = names(data)[dd$idx]
names(metadata)[x] = names(metadata)[dd$idx]

# re-order two coronal brains (coronal 66 and coronal 147) into a single dataset using their CCF coordinates (x-direction, anterior to posterior)
x  = which(str_detect(names(data), "coronal"))
dd = data.frame()

for (i in 1:length(x)){
  if (str_detect(names(data)[x], "coronal_147")[i]){brain = "coronal_147"} else {brain = "coronal_66"}
  temp_dd = data.frame(brain, slice = names(data)[x[i]], x = median(metadata[[x[i]]][,"CCFx"], na.rm = TRUE), idx = x[i])
  dd = rbind(dd, temp_dd)
}

# assign indices/ positions to slices
dd$x[1:3]   = -12:-10 # olfactory bulb slices of coronal 66 brain
dd$x[67:75] = -9:-1   # olfactory bulb slices of coronal 147 brain

dd$x[58:66]   = 13:21 # cerebellum of coronal 66 brains
dd$x[205:213] = 22:30 # cerebellum of coronal 147 brains

# notes on ordering: olfactory bulb
# olfactory bulb slices of coronal 66 are worse than those of coronal 147
# so, coronal 66 olfactory slices are more anterior and then it's the coronal 147 slices
# rest of the CCF registered slices are posterior to these slices

# notes on ordering: cerebellum
# cerebellum slices of coronal 66 brains are better than those of coronal 147
# so, coronal 66 cerebellum slices are more anterior and then it's the coronal 147 slices
# rest of the CCF registered slices are more anterior to these slices

# re-order dataset
dd = dd[order(dd$x),]

data[x]     = data[dd$idx]
metadata[x] = metadata[dd$idx]

names(data)[x]     = names(data)[dd$idx]
names(metadata)[x] = names(metadata)[dd$idx]

# arrange the Vizgen slices within the coronal sections
from = 1:9                                         # first 9 datasets needing new positions/ slots
to   = c(84, 85, 93, 119, 123, 124, 149, 150, 151) # positions where these datasets are slotted

# Function to insert an element into a list
insert_at = function(lst, value, name, position) {
  append(lst[1:(position - 1)], setNames(list(value), name), after = position - 1) |>
    append(lst[position:length(lst)])
}

# apply mapping in reverse order to preserve correct positions during shifting
for (i in length(from):1) {
  idx_from = from[i]
  idx_to   = to[i]

  # insert into `data`
  data   = insert_at(data[-idx_from], data[[idx_from]], names(data)[idx_from], idx_to)

  # insert into `metadata`
  metadata = insert_at(metadata[-idx_from], metadata[[idx_from]], names(metadata)[idx_from], idx_to)
}

# compute p-value for each slice
pvals = c()
for (i in 1:length(data)){pvals = c(pvals, cross_expression(data = data[[i]], locations = metadata[[i]][,c("x","y")])$cross_pvalue); print(i/length(data))}
pvals = p.adjust(pvals, method = "BH")

# plot all slices' cross-expression
if (FALSE){
  for (i in 1:length(data)){
    p = tissue_expression_plot(data = data[[i]], locations = metadata[[i]][,c("x","y")], gene1 = "Gpr6", gene2 = "Drd1", cross_expression = TRUE) + labs(subtitle = str_c(names(data)[i], "\np-value = ", signif(pvals[i], 2), "; id = ", i)); print(p); print(i/length(data))
  }
}

# plot significance against ordered slices
slice_type = vector(mode = "character", length   = length(data))
slice_type[str_detect(names(data), "sagittal")]  = "Sagittal slices, medial-to-lateral order"
slice_type[!str_detect(names(data), "sagittal")] = "Coronal slices, anterior-to-posterior order"

slice_order = vector(mode = "numeric", length = length(data))
slice_order[which(str_detect(names(data), "sagittal"))]  = 1:sum(str_detect(names(data), "sagittal"))
slice_order[which(!str_detect(names(data), "sagittal"))] = 1:sum(!str_detect(names(data), "sagittal"))

threshold = 5*10^-8; threshold = 0.05
temp_pvals = slice_type
temp_pvals[pvals <= threshold] = "Cross-expression"
temp_pvals[pvals >  threshold] = "No cross-expression"

dd = data.frame(slice_order, slice_type, pvals = -log10(pvals), significant = temp_pvals)

ggplot(dd) + aes(x = slice_order, y = pvals, color = significant) + geom_point(size = 1) +
  facet_wrap(~slice_type, scales = "free") + theme_minimal() +
  labs(x = "Slice order", y = "−log10 p-value", color = "") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = c("Cross-expression" = "brown3", "No cross-expression" = "#bbbbbb")) +
  theme(legend.position = "top", axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        legend.text = element_text(size = 10), strip.text = element_text(size = 10))

# plot each facet in turn and save it separately
for (i in 1:length(unique(dd$slice_type))){
  
  temp_dd = dd[dd$slice_type %in% unique(dd$slice_type)[i],]
  
  p = ggplot(dd[dd$slice_type %in% unique(dd$slice_type)[i],]) + aes(x = slice_order, y = pvals, color = significant) + geom_point(size = 1.5) +
    theme_minimal() + labs(x = "Slice order", y = "−log10 p-value", color = "", subtitle = unique(dd$slice_type)[i]) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_color_manual(values = c("Cross-expression" = "brown3", "No cross-expression" = "#bbbbbb")) +
    theme(legend.position = "right", axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          legend.text = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5, size = 12))
  print(p)

  ggsave(str_c("/inkwell05/images/PD_GWAS_", i, ".png"), device = "png", plot = p, width = 6, height = 4, dpi = 300)
}

# make joint plots
coronal  = dd$pvals[str_detect(dd$slice_type, "Coronal")]
sagittal = dd$pvals[str_detect(dd$slice_type, "Sagittal")]

Z = outer(coronal, sagittal, FUN = "*")

# convolve the outer product
Z = as.cimg(Z)
Z = isoblur(Z, sigma = 3)
Z = as.matrix(Z)

# heatmap annotations
row_origin = names(data)[str_detect(names(data), "sagittal")]
row_origin = str_extract(row_origin, "sagittal_\\d+")
row_origin = tools::toTitleCase(row_origin)

col_origin = names(data)[!str_detect(names(data), "sagittal")]
col_origin[str_detect(col_origin, "brain_1")] = "Replicate_1"
col_origin[str_detect(col_origin, "brain_2")] = "Replicate_2"
col_origin[str_detect(col_origin, "brain_3")] = "Replicate_3"
col_origin[str_detect(col_origin, "coronal")] = str_extract(col_origin[str_detect(col_origin, "coronal")], "coronal_\\d+")
col_origin = tools::toTitleCase(col_origin)

# define a common color mapping for shared annotations
brain_levels = unique(c(row_origin, col_origin))
brain_colors = setNames(RColorBrewer::brewer.pal(length(brain_levels), "Set1"), brain_levels)
brain_col_fun = list(Brain = brain_colors)

row_anno = HeatmapAnnotation(Brain = factor(row_origin, levels = brain_levels), col = brain_col_fun, which = "row", show_legend = FALSE, show_annotation_name = FALSE)
col_anno = HeatmapAnnotation(Brain = factor(col_origin, levels = brain_levels), col = brain_col_fun, which = "column", show_legend = TRUE,  show_annotation_name = FALSE, annotation_legend_param = list(title = "Brain"))

# draw heatmap
col_fun = colorRamp2(c(0, max(Z, na.rm = TRUE)), c("white", "red"))
ht = Heatmap(t(Z), cluster_rows = FALSE, cluster_columns = FALSE, heatmap_legend_param = list(title = "Cross-expression\nscore"), bottom_annotation = col_anno, left_annotation = row_anno, col = col_fun,
             row_title = str_c("Lateral-to-medial (", ncol(Z), " slices)"), column_title = str_c("Anterior-to-posterior (", nrow(Z), " slices)"), column_title_side = "top", row_title_side = "left")

draw(ht, merge_legend = TRUE, annotation_legend_side = "right", heatmap_legend_side = "right")

svg("/inkwell05/images/PD_cross_expression_heatmap.svg", width = 10, height = 6)
draw(ht, merge_legend = TRUE, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

```

```{r}




```


```{r}

# compare bivariate global Moran's I with our correlation-based approach to assess if they are similar (one example slice)
example_data = spatial_QC(path_to_expression = "/inkwell05/ameer/databases/spatial/Vizgen_2022_Mouse_MERSCOPE/expression/Vizgen_2022_Mouse_MERSCOPE_brain_2_slice_2_expression.h5ad",
                          path_to_metadata   = "/inkwell05/ameer/databases/spatial/Vizgen_2022_Mouse_MERSCOPE/metadata/Vizgen_2022_Mouse_MERSCOPE_brain_2_slice_2_metadata.csv")
example_meta = example_data$metadata
example_data = example_data$data

neighbors = RANN::nn2(data = example_meta[,c("x","y")], query = example_meta[,c("x","y")], k = 2)$nn.idx[,2] # neighbor indices

# correlation-based approach
mask_data = example_data
mask_data[mask_data > 0] = 1

mask_data_temp = example_data[neighbors,]
mask_data_temp[mask_data_temp > 0] = 1

X = mask_data * (1 - mask_data_temp)
Y = mask_data_temp * (1 - mask_data)

X = X * example_data
Y = Y * example_data[neighbors,]

corr = correlation(X, Y)

# bivariate global Moran's I approach
varX = example_data
varY = example_data

mat_dist = sparseMatrix(i = seq_len(length(neighbors)), j = neighbors, x = 1, dims = c(length(neighbors), length(neighbors)))
#mat_dist = mat_dist + t(mat_dist)   # commenting this off yields our correlation-based approach, where we use asymmetric relations
#mat_dist@x[] = 1

result   = mat_dist %*% varY           # same as re-ordering the expression matrix by neighbor indices

mask_data = varX
mask_data[mask_data > 0] = 1

mask_data_temp = result
mask_data_temp[mask_data_temp > 0] = 1

X = mask_data * (1 - mask_data_temp)
Y = mask_data_temp * (1 - mask_data)

X = X * varX
Y = Y * result

varX = scale(X)
varY = scale(Y)

bv_moran = (t(varX) %*% varY) / (nrow(varX) - 1)

```
