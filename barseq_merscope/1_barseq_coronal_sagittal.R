```{r}

# load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate","pryr",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo","ggrepel","igraph","ggraph")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))

# load coronal data
metadata_coronal = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
data_coronal     = t(as.matrix(readMM(file="/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt")))
colnames(data_coronal) = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv", header = TRUE))$Genes; rownames(data_coronal) = metadata_coronal$sample_id
data_coronal     = data_coronal[rowSums(data_coronal > 0) >= 5 & rowSums(data_coronal) >= 20,]
metadata_coronal = metadata_coronal[metadata_coronal$sample_id %in% rownames(data_coronal),]

# load sagittal data
metadata_sagittal  = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/metadata.csv"))
colnames(metadata_sagittal) = c("sample_id","slice","pos_x","pos_y","CCFx","CCFy","CCFz","CCFparentname","subclass_H2","cluster_H3")
CCFparentname = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/area_names.csv")); CCFparentname_vec = CCFparentname[[1]]
subclass_H2   = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/H2_code.csv", sep = NULL)); subclass_H2_vec = subclass_H2[[1]]
cluster_H3    = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/H3_code.csv", sep = NULL)); cluster_H3_vec  = cluster_H3[[1]]

metadata_sagittal <- metadata_sagittal %>% mutate(CCFparentname = CCFparentname_vec[CCFparentname],
                                subclass_H2   = subclass_H2_vec[subclass_H2],
                                cluster_H3    = cluster_H3_vec[cluster_H3])

data_sagittal = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/expression_matrix.csv")); rownames(data_sagittal) = metadata_sagittal$sample_id
genes         = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/gene_names.csv")); genes = genes[[1]]; colnames(data_sagittal) = genes
data_sagittal = data_sagittal[rowSums(data_sagittal > 0) >= 5 & rowSums(data_sagittal) >= 20,]
metadata_sagittal = metadata_sagittal[metadata_sagittal$sample_id %in% rownames(data_sagittal),]

# align both data
data_coronal  = data_coronal[, intersect(colnames(data_coronal), colnames(data_sagittal))]
data_sagittal = data_sagittal[, intersect(colnames(data_coronal), colnames(data_sagittal))]
data_sagittal = data_sagittal[, colnames(data_coronal)]

# each slice must have a certain min. no. of cells
if (TRUE){
  min_cells = 20000
  
  coronal  = table(metadata_coronal$slice)
  coronal  = as.numeric(names(coronal[coronal >= min_cells]))
  metadata_coronal = metadata_coronal[metadata_coronal$slice %in% coronal,]
  data_coronal     = data_coronal[rownames(data_coronal) %in% metadata_coronal$sample_id,]

  sagittal = table(metadata_sagittal$slice)
  sagittal = as.numeric(names(sagittal[sagittal >= min_cells]))
  metadata_sagittal = metadata_sagittal[metadata_sagittal$slice %in% sagittal,]
  data_sagittal     = data_sagittal[rownames(data_sagittal) %in% metadata_sagittal$sample_id,]
}

```


```{r}

# slice by slice cross-expression correlation within vs between brains
# slices are more similar within the brain than between the brain
# this reflects the fact that coronal/sagittal sections show diff. sampling areas
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# coronal slices
coronal  = as.data.frame(matrix(data = 0, nrow = choose(ncol(data_coronal),  2), ncol = length(unique(metadata_coronal$slice))));  colnames(coronal)  = unique(metadata_coronal$slice)

for (i in 1:ncol(coronal)){
  meta = metadata_coronal[metadata_coronal$slice %in% as.numeric(colnames(coronal))[i],]
  df   = data_coronal[rownames(data_coronal) %in% meta$sample_id,]
  meta = meta[,c("CCFx","CCFy","CCFz")]
  coronal[,i] = cross_expression_correlation(data = df, locations = meta)$correlation
  print(i/ncol(coronal))
}

# sagittal slices
sagittal = as.data.frame(matrix(data = 0, nrow = choose(ncol(data_sagittal), 2), ncol = length(unique(metadata_sagittal$slice)))); colnames(sagittal) = unique(metadata_sagittal$slice)

for (i in 1:ncol(sagittal)){
  meta = metadata_sagittal[metadata_sagittal$slice %in% as.numeric(colnames(sagittal))[i],]
  df   = data_sagittal[rownames(data_sagittal) %in% meta$sample_id,]
  meta = meta[,c("CCFx","CCFy","CCFz")]
  sagittal[,i] = cross_expression_correlation(data = df, locations = meta)$correlation
  print(i/ncol(sagittal))
}

# similarities in slice-by-slice cross-expression
corr_coronal  = cor(coronal,  method = "spearman"); corr_coronal  = upper_tri(corr_coronal)
corr_sagittal = cor(sagittal, method = "spearman"); corr_sagittal = upper_tri(corr_sagittal)
corr_between  = cor(coronal, sagittal, method = "spearman"); dist_corr_mixed = corr_between; corr_between = as.numeric(corr_between)

simm = data.frame(vals = c(corr_coronal, corr_sagittal, corr_between), type = rep(c("Coronal","Sagittal","Mixed"), times = c(length(corr_coronal), length(corr_sagittal), length(corr_between))))
simm$type = factor(simm$type, levels = rev(c("Coronal","Sagittal","Mixed")))

# p-values
pvals = signif(p.adjust(c(wilcox.test(corr_coronal,  corr_sagittal)$p.value, wilcox.test(corr_coronal,  corr_between)$p.value, wilcox.test(corr_sagittal, corr_between)$p.value), method = "BH"), digits = 2)

# plot
ggplot(simm) + aes(x = vals, y = type, fill = type) +
  geom_density_ridges(scale = 1) +
  labs(x = "Cross-expression similarity (correlation) between brain slices", y = "", fill = "",
       subtitle = str_c("Coronal vs Sagittal, p-value = ", pvals[1], "\nCoronal vs Mixed, p-value = ", pvals[2], "\nSagittal vs Mixed, p-value = ", pvals[3])) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("Coronal" = "brown3", "Sagittal" = "steelblue", "Mixed" = "gray88")) +
  theme_classic2() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        legend.position = "none",
        plot.subtitle = element_text(size = 10))
ggsave(filename = "/Users/AmeerSarwar/Downloads/mixed_coronal_sagittal.svg", device = "svg", dpi = 600, width = 6, height = 7)

# slice similarities as a function of distance for mixed brain
# reverse both coronal and sagittal slice names
pp = dist_corr_mixed
rownames(pp) = rev(rownames(pp))
colnames(pp) = rev(colnames(pp))
pp = data.frame(x = rep(colnames(pp), each = nrow(pp)), y = rep(rownames(pp), times = ncol(pp)), corr = as.numeric(pp))
pp$x = as.numeric(pp$x); pp$y = as.numeric(pp$y)
pp = data.frame(x = abs(pp$y - pp$x), y = pp$corr)

ggplot(pp) + aes(x = x, y = y) + geom_point(size = 0.5) + geom_smooth(method = "lm") +
  labs(x = "Distance (in slice ID's of sagittal versus coronal brains)", y = "Correlation",
       subtitle = str_c("Spearman's rho = ", signif(cor(pp$x, pp$y, method = "spearman"), digits = 2))) +
  theme_classic2() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/mixed_coronal_sagittal_distance.svg", device = "svg", dpi = 600, width = 6, height = 4)

```


```{r}

# density of cross-expressing cells in superior/dorsal to inferior/ventral direction in coronal vs sagittal brains

# cells + neighbors
coronal  = metadata_coronal[,c("CCFx","CCFy","CCFz")]
cor_neg  = RANN::nn2(data = coronal, query = coronal, k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2]
cor_df   = data_coronal
cor_df[cor_df > 0] = 1
cor_df_neg = cor_df[cor_neg,]

sagittal = metadata_sagittal[,c("CCFx","CCFy","CCFz")]
sag_neg  = RANN::nn2(data = sagittal, query = sagittal, k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2]
sag_df   = data_sagittal
sag_df[sag_df > 0] = 1 
sag_df_neg = sag_df[sag_neg,]

# densities
genes = data.frame(matrix(data = 0, nrow = ncol(data_coronal), ncol = ncol(data_coronal))); dimnames(genes) = list(colnames(data_coronal), colnames(data_coronal)); diag(genes) = 1

for (i in 1:nrow(genes)){
  
  for (j in 1:ncol(genes)){
    if (i==j){next}
    
    cor_exp = as.numeric((cor_df[,i] > 0 & cor_df_neg[,j] > 0) & (cor_df[,j] == 0 & cor_df_neg[,i] == 0))
    sag_exp = as.numeric((sag_df[,i] > 0 & sag_df_neg[,j] > 0) & (sag_df[,j] == 0 & sag_df_neg[,i] == 0))
    
    cor_dat = data.frame(continuous = coronal$CCFy,  binary = cor_exp)
    sag_dat = data.frame(continuous = sagittal$CCFy, binary = sag_exp)
    
    cor_dat = scale(cor_dat[cor_dat$binary == 1, ]$continuous); cor_dat = as.numeric(cor_dat)
    sag_dat = scale(sag_dat[sag_dat$binary == 1, ]$continuous); sag_dat = as.numeric(sag_dat)
    
    genes[i,j] = wilcox.test(cor_dat, sag_dat)$p.value
  }
  print(i/nrow(genes))
}

```


```{r}

#output = genes
#write.csv(x = output, file = "/inkwell01/ameer/Gene co-expression patterns to define cellular networks/Analysis1 - Cross-Expression/Validation 6/cross_exp_cell_density_dorsal_ventral.csv")

output = read.csv(file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Analysis1 - Cross-Expression/Validation 6/cross_exp_cell_density_dorsal_ventral.csv")
rownames(output) = output$X; output = output[,2:ncol(output)]; output = as.matrix(output)
output = data.frame(x = c(upper_tri(output), upper_tri(t(output))))
output$x = -log10(output$x)

ggplot(output) + aes(x = x) +
  geom_density(fill = "lightblue", color = "darkblue") +
  geom_vline(xintercept = -log10(0.05), color = "darkred", linetype = "dashed") +
  labs(x = "-log10 p-value", y = "Density",
       subtitle = str_c("Proportion without significant difference = ", signif(1 - sum(output$x >= -log10(0.05)) / nrow(output), digits = 2))) +
  theme_classic2() +
  scale_x_continuous(limits = c(0,2)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/density.svg", device = "svg", dpi = 600, width = 6, height = 4)

```


```{r}

# cross-expression per brain region across sagittal and coronal data
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# must have at least a min. number of cells per region
min_cells = 1000
region_coronal  = sort(table(metadata_coronal$CCFparentname))
region_coronal  = names(region_coronal[region_coronal >= min_cells])

region_sagittal = sort(table(metadata_sagittal$CCFparentname))
region_sagittal = names(region_sagittal[region_sagittal >= min_cells])

regions = intersect(region_coronal, region_sagittal)

region_coronal  = matrix(data = 0, ncol = length(regions), nrow = choose(ncol(data_coronal), 2))
colnames(region_coronal) = regions
region_sagittal = region_coronal

for (i in 1:length(regions)){
  
  # coronal
  temp_meta = metadata_coronal[metadata_coronal$CCFparentname %in% regions[i],]
  temp_data = data_coronal[rownames(data_coronal) %in% temp_meta$sample_id,]
  locations = temp_meta[,c("CCFx","CCFy","CCFz")]
  region_coronal[,i] = cross_expression_correlation(data = temp_data, locations = locations)$correlation
  
  # sagittal
  temp_meta = metadata_sagittal[metadata_sagittal$CCFparentname %in% regions[i],]
  temp_data = data_sagittal[rownames(data_sagittal) %in% temp_meta$sample_id,]
  locations = temp_meta[,c("CCFx","CCFy","CCFz")]
  region_sagittal[,i] = cross_expression_correlation(data = temp_data, locations = locations)$correlation
  
  print(i/length(regions))
}

# pre-process
corr = region_coronal;  corr[is.na(corr)] = 0
sagg = region_sagittal; sagg[is.na(sagg)] = 0
simm = cor(corr, sagg)

same = diag(simm)
diff = c(upper_tri(simm), lower_tri(simm))
diff = c(diff, upper_tri(cor(corr)), upper_tri(cor(sagg))) # diff. regions within brains

simm = data.frame(vals = c(same, diff), class = rep(c("Same (between brains)","Different (within or between brains)"), times = c(length(same), length(diff))))

# draw plot
ggplot(simm) + aes(x = vals, y = after_stat(..scaled..), fill = class) +
  geom_density(alpha = 0.8) +
  scale_fill_manual(values  = c("Same (between brains)" = "lightblue", "Different (within or between brains)" = "gray88")) +
  labs(x = "Cross-expression similarity (correlation) between brain regions", y = "Density",
       fill = "", subtitle = str_c("P-value = ", signif(wilcox.test(same,diff)$p.value, digits = 2))) +
  theme_classic2() +
  theme(legend.position = "top",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/density.svg", device = "svg", dpi = 600, width = 6, height = 4)

```