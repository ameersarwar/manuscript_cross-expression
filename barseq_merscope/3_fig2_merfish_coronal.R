```{r}

# MERFISH brain receptor map
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))
registerDoParallel(detectCores() - 1)
plan("multisession", workers = 8)

# choose slice and replicate (mouse)
slice = 2; replicate = 2; cluster = FALSE

# load dataset
data     <- fread(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice",slice,"_Replicate",replicate,"_cell_by_gene.csv"))
metadata <- fread(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice",slice,"_Replicate",replicate,"_cell_metadata.csv"))
CCFnames <- fread(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice",slice,"_Replicate",replicate,"_CCFnames.csv"))
colnames(data)[1] <- "sample_id"; colnames(metadata)[1] <- "sample_id"; data <- as.data.frame(data); colnames(metadata)[4:5] <- c("pos_x","pos_y")
data <- data[,which(!str_detect(colnames(data), "^Blank-"))]; metadata <- as.data.frame(metadata); CCFnames <- as.data.frame(CCFnames)

# perform cell type clustering and plot on brain slice
if (cluster){
  data_matrix   <- t(Matrix(as.matrix(data[,2:ncol(data)]), sparse = TRUE)); colnames(data_matrix) <- 1:nrow(data)
  seurat_object <- CreateSeuratObject(counts = data_matrix)
  seurat_object <- SCTransform(seurat_object, clip.range = c(-10,10), )
  seurat_object <- RunPCA(seurat_object, npcs = 30, features = colnames(data))
  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:30)
  seurat_object <- FindClusters(seurat_object, resolution = 0.3)
  clusters      <- seurat_object@meta.data[["seurat_clusters"]]
  metadata$clusters <- clusters
}

# process the data
metadata <- metadata[,which(!colnames(metadata) %in% c("fov","volume","min_x","max_x","min_y","max_y"))]
metadata <- data.frame(metadata, CCFnames)
metadata <- metadata[which(!is.na(metadata$CCFparentname)),]; data <- data[which(data$sample_id %in% metadata$sample_id),]
data <- data[which(rowSums(data[,2:ncol(data)]) >= 50),]; metadata <- metadata[which(metadata$sample_id %in% data$sample_id),]

# rotate slice (counter-clockwise) and center to origin
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
degs <- 40
rot  <- rotate_coordinates(x=metadata$pos_x, y=metadata$pos_y, n_degrees = degs)
metadata$pos_x <- rot$pos_x; metadata$pos_y <- rot$pos_y

```


```{r}

# find gene pairs with regional cross-expression
# choose great looking gene pairs
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
cross = cross_expression(data = data[,2:ncol(data)], locations = metadata[,c("pos_x","pos_y")], alpha_co = 0.01, alpha_cross = 0.01)
cross = cross[cross$gene1 %in% "Lgr6" | cross$gene2 %in% "Lgr6",]
cross = cross[cross$cross_sig == 1, ]

for (i in 1:nrow(cross)){
  
  # genes' cell-neighbor pairs
  df = data[,c("sample_id", cross$gene1[i], cross$gene2[i])]
  locations = metadata[,c("pos_x","pos_y")]
  distances = RANN::nn2(locations, locations, k = 2, searchtype = "priority")
  distances = distances$nn.idx[,2]
  df_neig   = df[distances,]
  
  cells     = df$sample_id
  neigs     = df$sample_id[distances]
  
  df      = df[,2:ncol(df)]
  df_neig = df_neig[,2:ncol(df_neig)]
  
  # cell-neighbor pairs to light up
  pair1 = df[,1] > 0 & df[,2] == 0 & df_neig[,2] > 0 & df_neig[,1] == 0
  pair2 = df[,2] > 0 & df[,1] == 0 & df_neig[,1] > 0 & df_neig[,2] == 0
  
  # cells lit up with genes 1 or 2
  gene1 = unique(c(cells[pair1], neigs[pair2]))
  gene2 = unique(c(cells[pair2], neigs[pair1]))
  
  # vector for lighting up cells with genes 1 or 2
  light = vector(mode = "character", length = nrow(metadata))
  light[metadata$sample_id %in% gene1] = cross$gene1[i]
  light[metadata$sample_id %in% gene2] = cross$gene2[i]
  light[light == ""] = "Neither"
  
  # process for plotting
  meta = data.frame(metadata, light)
  meta$light = factor(meta$light)
  
  meta_neither = meta[meta$light == "Neither", ]
  meta_others  = meta[meta$light != "Neither", ]
  meta_reordered = rbind(meta_neither, meta_others)
  
  # set color scheme
  colors = c("Neither" = "gray88")
  levels_without_neither = setdiff(levels(meta$light), "Neither")
  default_colors = c("brown3", "deepskyblue4")
  names(default_colors) = levels_without_neither
  final_colors = c(colors, default_colors)
  
  # plot expression on tissue
  p = ggplot(meta_reordered) + aes(x = pos_x, y = pos_y, color = light) +
    geom_point(size = 0.1) +
    scale_color_manual(values = final_colors) +
    guides(color = guide_legend(override.aes = list(size = 3), reverse = TRUE)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "none",
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_blank())
  print(p)
  ggsave("/Users/AmeerSarwar/Downloads/mask.png", device = "png", width = 6, height = 4, dpi = 300)
  
  # plot cells colored by expression
  light = vector(mode = "character", length = nrow(locations))
  light[df[,1] > 0 & df[,2] == 0] = colnames(df)[1]
  light[df[,1] == 0 & df[,2] > 0] = colnames(df)[2]
  light[df[,1] > 0 & df[,2] > 0]  = "Both"
  light[light == ""] = "Neither"
  
  # process for plotting
  meta = data.frame(locations, light)
  meta$light = factor(meta$light)
  
  meta_neither = meta[meta$light == "Neither", ]
  meta_others  = meta[meta$light != "Neither", ]
  meta_reordered = rbind(meta_neither, meta_others)
  
  # set color scheme
  colors = c("Neither" = "gray88", "Both" = "chartreuse3")
  levels_without_neither = setdiff(levels(meta$light), c("Neither","Both"))
  default_colors = c("brown3", "deepskyblue4")
  names(default_colors) = levels_without_neither
  final_colors = c(colors, default_colors)
  
  # plot tissue expression
  p = ggplot(meta_reordered) + aes(x = pos_x, y = pos_y, color = light) +
    geom_point(size = 0.1) +
    scale_color_manual(values = final_colors) +
    guides(color = guide_legend(override.aes = list(size = 3), reverse = TRUE)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "none",
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_blank())
  print(p)
  ggsave("/Users/AmeerSarwar/Downloads/mask.png", device = "png", width = 6, height = 4, dpi = 300)
  
  # plot gene vs gene
  colnames(df) = c("gene1","gene2")
  p = ggplot(df) + aes(x = gene1, y = gene2) +
    geom_jitter(size = 0.1) + labs(x = cross$gene1[i], y = cross$gene2[i],
                         subtitle = str_c("Spearman's rho = ",
                                          signif(cor(df[,1], df[,2], method = "spearman"), digits = 2))) +
    theme_bw()
  print(p)
  
  print(i/nrow(cross))
}

```


```{r}

# we see genes cross-expressed with Lgr6 in the thalamus
# check if they are co-expressed with each other in the thalamus
# co-expression of other genes
dat = data[,unique(c(cross$gene1, cross$gene2))]
dat = dat[,!colnames(dat) %in% "Lgr6"]

# significant co-expression
dat_cross = cross_expression(data = dat, locations = metadata[,c("pos_x","pos_y")])

for (i in 1:(ncol(dat)-1)) {
  
  for (j in (i+1):ncol(dat)) {
    
    dd = dat[,c(i,j)]
    locations = metadata[,c("pos_x","pos_y")]
    
    # gene-gene co-expression
    colnames(dd) = c("gene1","gene2")
    p = ggplot(dd) + aes(x = gene1, y = gene2) +
      geom_jitter(size = 0.1) + labs(x = colnames(dat)[i], y = colnames(dat)[j],
                         subtitle = str_c("Spearman's rho = ",
                                          signif(cor(dd[,1], dd[,2], method = "spearman"), digits = 2))) +
      theme_bw()
    print(p)
    
    # plot cells colored by expression
    dd = dat[,c(i,j)]
    light = vector(mode = "character", length = nrow(locations))
    light[dd[,1] > 0 & dd[,2] == 0] = colnames(dd)[1]
    light[dd[,1] == 0 & dd[,2] > 0] = colnames(dd)[2]
    light[dd[,1] > 0 & dd[,2] > 0]  = "Both"
    light[light == ""] = "Neither"
  
    # process for plotting
    meta = data.frame(locations, light)
    meta$light = factor(meta$light)
  
    meta_neither = meta[meta$light == "Neither", ]
    meta_others  = meta[meta$light != "Neither", ]
    meta_reordered = rbind(meta_neither, meta_others)
  
    # set color scheme
    colors = c("Neither" = "gray88", "Both" = "chartreuse3")
    levels_without_neither = setdiff(levels(meta$light), c("Neither","Both"))
    default_colors = c("brown3", "deepskyblue4")
    names(default_colors) = levels_without_neither
    final_colors = c(colors, default_colors)
  
    # plot tissue expression
    p = ggplot(meta_reordered) + aes(x = pos_x, y = pos_y, color = light) +
      geom_point(size = 0.1) +
      scale_color_manual(values = final_colors) +
      guides(color = guide_legend(override.aes = list(size = 3), reverse = TRUE)) +
      theme_classic() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            #legend.position = "none",
            legend.background = element_blank(),
            legend.key = element_blank(),
            panel.border = element_blank())
    print(p)
    ggsave("/Users/AmeerSarwar/Downloads/mask.png", device = "png", width = 6, height = 4, dpi = 300)
    
  }
  print(i/(ncol(dat)-1))
}

```
```{r}

# spatial enrichment of cross-expression and co-expression
sp = data.frame(gene1 = c("Lgr6","Lgr6"), gene2 = c("Adra2b","Ret"))

for (i in 1:nrow(sp)){
  
  sp_enrich = spatial_enrichment(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = sp$gene1[i], gene2 = sp$gene2[i])
  plot_enrich = data.frame(vals = c(sp_enrich$target, sp_enrich$null), type = rep(c("Cross-expressing pairs","Random pairs"), times = c(length(sp_enrich$target), length(sp_enrich$null))))
  plot_enrich$type <- factor(plot_enrich$type, levels = c("Random pairs", "Cross-expressing pairs"))
  
  p = ggplot(plot_enrich) + aes(x = vals, y = after_stat(scaled), fill = type) +
    geom_density(aes(color = type), alpha = 0.5) +
    labs(subtitle = str_c("p-value = ", signif(sp_enrich$pvalue, digits = 2)),
         x = "Eucliean distance (µm)", y = "Density", fill ="") +
    scale_fill_manual(values = c("Cross-expressing pairs" = "steelblue", "Random pairs" = "gray88")) +
    scale_color_manual(values = c("Cross-expressing pairs" = "darkblue", "Random pairs" = "gray")) +
    guides(color = "none") +
    theme_classic() +
    theme(legend.position = "top",
          legend.text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.subtitle = element_text(size = 12))
  print(p)
}

```


```{r}

# some genes cross-express with Lgr6
# they are co-expressed with each other
# show higher co-expression of these genes in neighbors of Lgr6 cells vs other cells

# expression matrices
dat  = data[,unique(c(cross$gene1, cross$gene2))]
lgr6 = data.frame(Lgr6 = dat$Lgr6)
dat  = dat[,!colnames(dat) %in% "Lgr6"]
locations = metadata[,c("pos_x","pos_y")]

distances = RANN::nn2(locations, locations, k = 2, searchtype = "priority")
distances = distances$nn.idx[,2]
dat_neig  = dat[distances,]

# Lgr6 neighbors vs other cells co-expression
neigs = dat_neig[which(lgr6 > 0),]
down  = nrow(neigs)
cells = dat_neig[which(lgr6 == 0), ]

neigs = cor(neigs, method = "spearman")
neigs = upper_tri(neigs)
neigs = data.frame(id = 1:length(neigs), coexp = neigs, class = "Neighbors")

# downsample cells and take average co-expression
corr_cell = 0

for (i in 1:100){
  cell = cells[sample(1:nrow(cells), size = down, replace = FALSE), ]
  cell = cor(cell, method = "spearman")
  cell = upper_tri(cell)
  corr_cell = corr_cell + cell
}

corr_cell = corr_cell / 100
cells = data.frame(id = 1:length(corr_cell), coexp = corr_cell, class = "Non-neighbors")

# process co-expression for plotting
coexp = rbind(cells, neigs)
coexp$class = factor(coexp$class, levels = c("Neighbors","Non-neighbors"))
neigs = neigs$coexp
cells = cells$coexp
pval  = wilcox.test(neigs, cells)

# show boxplots
ggplot(coexp) + aes(y = class, x = coexp, group = id) +
  geom_violin(aes(group = class), fill = "grey", alpha = 0.5) +
  geom_point(aes(color = as.factor(id)), size = 2) +
  geom_line(aes(color = as.factor(id)), linewidth = 0.2) +
  labs(title = "Whole-slice tissue",
       x = "Co-expression among genes cross-expressed with Lgr6",
       y = "Cell location in relation to Lgr6-expressing cells",
       subtitle = str_c("P-value = ", signif(pval$p.value, digits = 2))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 90, margin = margin(l = 5, t = 10), size = 10),
        axis.text.x = element_text(margin = margin(b = 2)))
ggsave("/Users/AmeerSarwar/Downloads/violin.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# co-expression differences above might be due to sampling from thalamic vs non-thalamic tissue
# restrict to the thalamus and see if neighbors of Lgr6-expressing cells have higher co-expression than non-neighbors

# thalamic regions
regions = c("ATN","ILM","LAT","MED","MTN","VENT","VP") # missing for mask DORpm
regions = metadata[metadata$CCFparentname %in% regions, ]

# expression matrices
dat  = data[data$sample_id %in% regions$sample_id, unique(c(cross$gene1, cross$gene2))]
lgr6 = data.frame(Lgr6 = dat$Lgr6)
dat  = dat[,!colnames(dat) %in% "Lgr6"]
locations = regions[,c("pos_x","pos_y")]

distances = RANN::nn2(locations, locations, k = 2, searchtype = "priority")
distances = distances$nn.idx[,2]
dat_neig  = dat[distances,]

# Lgr6 neighbors vs other cells co-expression
neigs = dat_neig[which(lgr6 > 0),]
down  = nrow(neigs)
cells = dat_neig[which(lgr6 == 0), ]

neigs = cor(neigs, method = "spearman")
neigs = upper_tri(neigs)
neigs = data.frame(id = 1:length(neigs), coexp = neigs, class = "Neighbors")

# downsample cells and take average co-expression
corr_cell = 0

for (i in 1:100){
  cell = cells[sample(1:nrow(cells), size = down, replace = FALSE), ]
  cell = cor(cell, method = "spearman")
  cell = upper_tri(cell)
  corr_cell = corr_cell + cell
}

corr_cell = corr_cell / 100
cells = data.frame(id = 1:length(corr_cell), coexp = corr_cell, class = "Non-neighbors")

# process co-expression for plotting
coexp = rbind(cells, neigs)
coexp$class = factor(coexp$class, levels = c("Neighbors","Non-neighbors"))
neigs = neigs$coexp
cells = cells$coexp
pval  = wilcox.test(neigs, cells)

# show boxplots
ggplot(coexp) + aes(y = class, x = coexp, group = id) +
  geom_violin(aes(group = class), fill = "grey", alpha = 0.5) +
  geom_point(aes(color = as.factor(id)), size = 2) +
  geom_line(aes(color = as.factor(id)), linewidth = 0.2) +
  labs(title = "Thalamic tissue",
       x = "Co-expression among genes cross-expressed with Lgr6",
       y = "Cell location in relation to Lgr6-expressing cells",
       subtitle = str_c("P-value = ", signif(pval$p.value, digits = 2))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 90, margin = margin(l = 5, t = 10), size = 10),
        axis.text.x = element_text(margin = margin(b = 2)))
ggsave("/Users/AmeerSarwar/Downloads/violin.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# co-expression between neighbors and non-neighbors of Lgr6-expressing cells in the thalamus is similar
# but Lgr6 is expressed in non-thalamic regions, too
# here, we compare neighbors vs non-neighbors of Lgr6-expressing cells in non-thalamic regions

# thalamic regions
regions = c("ATN","ILM","LAT","MED","MTN","VENT","VP") # missing for mask DORpm
regions = metadata[!metadata$CCFparentname %in% regions, ]

# expression matrices
dat  = data[data$sample_id %in% regions$sample_id, unique(c(cross$gene1, cross$gene2))]
lgr6 = data.frame(Lgr6 = dat$Lgr6)
dat  = dat[,!colnames(dat) %in% "Lgr6"]
locations = regions[,c("pos_x","pos_y")]

distances = RANN::nn2(locations, locations, k = 2, searchtype = "priority")
distances = distances$nn.idx[,2]
dat_neig  = dat[distances,]

# Lgr6 neighbors vs other cells co-expression
neigs = dat_neig[which(lgr6 > 0),]
down  = nrow(neigs)
cells = dat_neig[which(lgr6 == 0), ]

neigs = cor(neigs, method = "spearman")
neigs = upper_tri(neigs)
neigs = data.frame(id = 1:length(neigs), coexp = neigs, class = "Neighbors")

# downsample cells and take average co-expression
corr_cell = 0

for (i in 1:100){
  cell = cells[sample(1:nrow(cells), size = down, replace = FALSE), ]
  cell = cor(cell, method = "spearman")
  cell = upper_tri(cell)
  corr_cell = corr_cell + cell
}

corr_cell = corr_cell / 100
cells = data.frame(id = 1:length(corr_cell), coexp = corr_cell, class = "Non-neighbors")

# process co-expression for plotting
coexp = rbind(cells, neigs)
coexp$class = factor(coexp$class, levels = c("Neighbors","Non-neighbors"))
neigs = neigs$coexp
cells = cells$coexp
pval  = wilcox.test(neigs, cells)

# show boxplots
ggplot(coexp) + aes(y = class, x = coexp, group = id) +
  geom_violin(aes(group = class), fill = "grey", alpha = 0.5) +
  geom_point(aes(color = as.factor(id)), size = 2) +
  geom_line(aes(color = as.factor(id)), linewidth = 0.2) +
  labs(title = "Non-thalamic tissue",
       x = "Co-expression among genes cross-expressed with Lgr6",
       y = "Cell location in relation to Lgr6-expressing cells",
       subtitle = str_c("P-value = ", signif(pval$p.value, digits = 2))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 90, margin = margin(l = 5, t = 10), size = 10),
        axis.text.x = element_text(margin = margin(b = 2)))
ggsave("/Users/AmeerSarwar/Downloads/violin.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# we find Lgr6-Adra2b and Lgr6-Ret as showing cross-expression in the thalamus
# even though they are expressed across the brain
# compute % of cells cross-expressing in the thalamus vs non-thalamus
# compute % of cells expressing each gene in thalamus vs non-thalamus

regions = c("ATN","ILM","LAT","MED","MTN","VENT","VP") # thalamic mask
thalamus= as.numeric(metadata$CCFparentname %in% regions)
genes   = c("Lgr6","Adra2b","Ret")
df      = data[, colnames(data) %in% genes]
df      = cbind(df, thalamus); df = as.matrix(df)

# get nearest neighbors
locations  = metadata[,c("pos_x","pos_y")]
df[df > 0] = 1
neighbor   = 2
distances  = RANN::nn2(locations, locations, k = neighbor, searchtype = "priority")
distances  = distances$nn.idx[,neighbor]
data_temp  = df[distances,]

# cell pairs belong to thalamus when both cells are thalamic
thalamus = as.logical(df[,4] * data_temp[,4])
non_thalamus = as.logical((1 - df[,4]) * (1 - data_temp[,4]))

df_thalamus = df[thalamus,]
df_non_thamalus = df[non_thalamus,]

data_temp_thalamus = data_temp[thalamus,]
data_temp_non_thalamus = data_temp[non_thalamus,]

# co-expression of Adra2b and Ret (two genes cross-expressing with Lgr6)
tha = t(df_thalamus) %*% df_thalamus
non = t(df_non_thamalus) %*% df_non_thamalus

# total cells expressing each gene in thalamus vs non-thalamus
total_thalamus = colSums(df[df[,4]==1, ]); total_thalamus = total_thalamus[-4]
total_non_thalamus = colSums(df[df[,4]!=1, ]); total_non_thalamus = total_non_thalamus[-4]

# total cross-expressing cells in thalamus vs non-thalamus
cross_thalamus = t(df_thalamus * (1 - data_temp_thalamus)) %*% ((1 - df_thalamus) * data_temp_thalamus)
cross_thalamus = cross_thalamus[-4,-4]; cross_thalamus = cross_thalamus[2,-2]

cross_non_thalamus = t(df_non_thamalus * (1 - data_temp_non_thalamus)) %*% ((1 - df_non_thamalus) * data_temp_non_thalamus)
cross_non_thalamus = cross_non_thalamus[-4,-4]; cross_non_thalamus = cross_non_thalamus[2,-2]

# convert to percentages within thalamus
total = total_thalamus / (total_thalamus + total_non_thalamus); total = total[-2]
cross = cross_thalamus / (cross_thalamus + cross_non_thalamus)

# make bar plot
df = data.frame(total, cross)
df = cbind(df, genes = rownames(df))
df = pivot_longer(df, 1:2, names_to = "type", values_to = "values") %>% as.data.frame()
df$type = factor(df$type, levels = (unique(df$type)))


ggplot(df) + aes(x = genes, y = values, fill = type) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  labs(x = "Gene expression within cells or cross-expression with Lgr6",
       fill = "Expression type", y = "Proportion of cells or cell pairs in the thalamus") +
  scale_fill_manual(values = c("cross" = "#3498db", "total" = "#e74c3c"),
                    labels = c("cross" = "Between cells", "total" = "Within cells")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "top")

ggsave("/Users/AmeerSarwar/Downloads/within_between.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# we find thalamus vs non-thalamus co-expression difference
# we determined this using cross-expression
# check that this also appears in single-cell data
# this ensures that the results are robust
# NOTE: the key finding then is differential co-expression in neighbors vs non-neighbors of Lgr6-positive cells in non-thalamic regions
# naturally, we can only figure this out due to cross-expression, not due to anything else

# cross-expression
cross = cross_expression(data = data[,2:ncol(data)], locations = metadata[,c("pos_x","pos_y")], alpha_co = 0.01, alpha_cross = 0.01)
cross = cross[cross$gene1 %in% "Lgr6" | cross$gene2 %in% "Lgr6",]
cross = cross[cross$cross_sig == 1, ]

# non-thalamic regions
spatial_regions = sort(unique(metadata$CCFparentname))
single_cell_regions = c("AId-AIv-AIp","CTXsp-CLA-EP-LA-BLA-BMA-PA","HIP-CA","STR-sAMY","STR-LSX","STR-STRd","STR-STRv","OLF-COA-PAA-NLOT-TR","HIP","CTXsp-CLA-EP-LA-BLA-BMA-PA","HIP",
                        "HY-MEZ-PVZ-PVR","OLF-AON-TT-DP-PIR-COA-NLOT-PAA-TR","OLF-MOB-AOB","HY-MEZ-PVZ-PVR","RSP","SSp","SSs-GU-VISC-AIp","STR-STRd","VISa","SS-GU-VISC")
single_cell_regions = unique(single_cell_regions)
single_cell_regions = str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.",
      single_cell_regions,".h5ad")

# load scRNA-seq data
df <- NULL

for (i in 1:length(single_cell_regions)){
  
  sc <- zellkonverter::readH5AD(single_cell_regions[i])
  sc <- t(assay(sc, "X"))
  df <- rbind(df, sc)
  print(i/length(single_cell_regions))
}

# subset by common genes
common = intersect(colnames(df), colnames(data))
df = df[, colnames(df) %in% common]
data = data[, colnames(data) %in% common]
df = df[, colnames(data)]
df = df[, colnames(df) %in% unique(c(cross$gene1, cross$gene2))]
df = df[, !colnames(df) %in% "Lgr6"]
non_thalamic = df

# thalamic regions
single_cell_regions = c("TH-AD-AV-AM-IAD-LD","TH-AD-AV-AM-IAD-LD","TH-LGd-IGL-LGv","TH-MD-IMD-PCN-CL","TH-MH-LH-LP","TH-PO-Eth","TH-PVT-PT","TH-RE-RH-CM-IAM-SMT-PR-Xi",
                        "TH-MH-LH-LP","TH-PO-Eth","TH-PVT-PT","TH-VAL-VPM-VPL-VM","TH-PF-SPA-SPFm-VPMpc-VPLpc-RT","TH-VAL-VPM-VPL-VM")
single_cell_regions = unique(single_cell_regions)
single_cell_regions = str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.",
      single_cell_regions,".h5ad")

# load scRNA-seq data
df <- NULL

for (i in 1:length(single_cell_regions)){
  
  sc <- zellkonverter::readH5AD(single_cell_regions[i])
  sc <- t(assay(sc, "X"))
  df <- rbind(df, sc)
}

# subset by common genes
common = intersect(colnames(df), colnames(data))
df = df[, colnames(df) %in% common]
data = data[, colnames(data) %in% common]
df = df[, colnames(data)]
df = df[, colnames(df) %in% unique(c(cross$gene1, cross$gene2))]
df = df[, !colnames(df) %in% "Lgr6"]
thalamic = df

# co-expression in thalamic vs non-thalamic regions
non_thalamic_corr = cor(as.matrix(non_thalamic), method = "spearman")
non_thalamic_corr = upper_tri(non_thalamic_corr)
thalamic_corr     = cor(as.matrix(thalamic), method = "spearman")
thalamic_corr     = upper_tri(thalamic_corr)

# process co-expression for plotting
coexp = c(non_thalamic_corr, thalamic_corr)
coexp = data.frame(corr = coexp, class = c(rep("Non-thalamic", length(non_thalamic_corr)), rep("Thalamic", length(thalamic_corr))))
coexp = data.frame(id = c(1:length(non_thalamic_corr), 1:length(thalamic_corr)), coexp)
coexp$class = factor(coexp$class, levels = c("Thalamic","Non-thalamic"))
pval  = wilcox.test(non_thalamic_corr, thalamic_corr)

# show boxplots
ggplot(coexp) + aes(y = class, x = corr, group = id) +
  geom_violin(aes(group = class), fill = "grey", alpha = 0.5) +
  geom_point(aes(color = as.factor(id)), size = 2) +
  geom_line(aes(color = as.factor(id)), linewidth = 0.2) +
  labs(title = "scRNA-seq dataset",
       x = "Co-expression",
       y = "Tissue",
       subtitle = str_c("P-value = ", signif(pval$p.value, digits = 2))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 90, margin = margin(l = 5, t = 10), size = 10),
        axis.text.x = element_text(margin = margin(b = 2)))
ggsave("/Users/AmeerSarwar/Downloads/violin.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# brain mask
regions = c("ATN","ILM","LAT","MED","MTN","VENT","VP") # missing for mask DORpm
regions = metadata$CCFparentname %in% regions
regions = ifelse(regions, "Thalamic", "Non-thalamic")
regions = data.frame(metadata, regions)
regions$regions = factor(regions$regions, levels = c("Thalamic", "Non-thalamic"))

ggplot(regions) + aes(x = pos_x, y = pos_y, color = regions) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("Non-thalamic" = "gray88", "Thalamic" = "#696969")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(legend.position.inside = c(0.8,0.95),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggsave("/Users/AmeerSarwar/Downloads/mask.png", device = "png", width = 6, height = 4, dpi = 300)

```
