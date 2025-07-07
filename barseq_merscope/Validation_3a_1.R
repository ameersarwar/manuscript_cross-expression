```{r}

# BAR-seq coronal data
# data is shrunk by removing image stitching-related artefacts (cf. Xiaoyin's email)
# data is quality controlled by keeping cells with genes/cell >= 5 and reads/cell >= 20
# data alongside CCF and slide coordinates are saved and can be used for analysis

# load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo","ggrepel")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))

# specify coverage
coverage = "cortex"
#coverage = "whole-brain"

# load data
if (coverage=="whole-brain"){
  coord    = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
  data     = as.matrix(readMM(file="/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt"))
  rownames(data) = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv", header = TRUE))$Genes; colnames(data) = coord$sample_id
  coord_ref = coord
  
} else if (coverage=="cortex"){
  coord    = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
  data     = as.matrix(readMM(file="/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt"))
  rownames(data) = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv", header = TRUE))$Genes; colnames(data) = coord$sample_id
  coord_ref = coord
  
  # subset for cortex
  coord    = coord[which(coord$sample_id %in% as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_cortex.csv"))$sample_id),]
  data     = data[,which(colnames(data)  %in% as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_cortex.csv"))$sample_id)]
}

# choose slice
slice_current = round(runif(1, min = min(unique(coord$slice)), max = max(unique(coord$slice))))
slice_current = 17

coordinates <- coord %>% filter(slice==slice_current); coord_ref <- coord_ref %>% filter(slice==slice_current)
data = data[,which(colnames(data) %in% coordinates$sample_id)]
ggplot(coordinates) + aes(x=pos_x, y=pos_y, color=CCFparentname) + geom_point(size=0.5)

```


```{r}

# overall purpose is to find gene pairs showing cell-cell relations between two individual cells
metadata <- coordinates

# remove cell types appearing less than 30 times within each brain region
metadata <- metadata %>% group_by(CCFparentname, subclass_H2) %>% filter(n() >= 30) %>% group_by(CCFparentname) %>% filter(n_distinct(subclass_H2) > 1) %>% ungroup() %>% as.data.frame()

# remove brain regions with fewer than 200 cells
metadata <- metadata[which(metadata$CCFparentname %in% sort(unique(metadata$CCFparentname))[table(metadata$CCFparentname) >= 200]),] # remove regions with fewer than 200 cells

# remove brain regions with poor shapes
to_remove = c("SSp-un","CTXsp","SSp-tr","ECT","EP","TEa")
metadata  = metadata[which(!metadata$CCFparentname %in% to_remove),]

region = sort(unique(metadata$CCFparentname)); region = c("AUDv","AUDd","AUDp")

# regions of interest
for (i in 1:length(region)){
  metadata_subset1 <- metadata[which(metadata$CCFparentname %in% region[i]),]
  p = ggplot(metadata_subset1) + aes(x=pos_x, y=pos_y, color=subclass_H2) + geom_point(size=0.5) + labs(title = region[i], color="Cell type", x="x coordinates", y="y coordinates") +
    theme_bw() + theme(axis.ticks = element_blank(), axis.text = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10), legend.key.size = unit(0.4, "cm"), legend.justification = c(1,0.5))
  print(p); print(str_c("Percent complete = ",signif(i/length(region)*100, digits = 3)))
}

# reference brain map with regions
coord_ref$map_regions <- ifelse(coord_ref$CCFparentname %in% region, coord_ref$CCFparentname, NA)
unique_names   <- setdiff(sort(unique(coord_ref$map_regions)), "NA")
default_colors <- scales::hue_pal()(length(unique_names))
palette_colors <- brewer.pal(length(unique_names), "Paired")
color_map      <- c(setNames(palette_colors, unique_names), "NA" = "gray")

ggplot(coord_ref) + aes(x=pos_x, y=pos_y, color=map_regions) + geom_point(size=0.2) +
  theme_bw() + labs(x="x coordinates",y="y coordinates", color="Regions") + scale_color_manual(values = color_map, breaks = unique_names) +
  guides(color = guide_legend(override.aes = list(size = 0.8))) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10), legend.key.size = unit(0.4, "cm"), legend.justification = c(1,0.5))
ggsave("/Users/AmeerSarwar/Downloads/brainmap.png")

```


```{r}

# volcano plot showing effect size and p-value for a region of interest
# make volcano plots of two different regions and show most interesting gene pairs
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")
region = c("AUDv","AUDd","AUDp"); region = c("RSPd",  "RSPv") # these regions give different patterns of volcano plots
metadata_subset1 = metadata %>% filter(CCFparentname %in% region)
data_subset1     = t(data); data_subset1 = data_subset1[which(rownames(data_subset1) %in% metadata_subset1$sample_id),]

# p-values and effect sizes
params = celltocell2(data = data_subset1,
                     locations = metadata_subset1[,c("pos_x","pos_y")])

effect = params$effect_size
pvals  = params$pvalue_cross_expression

# account for co-expression
#contrl = params$pvalue_co_expression
#contrl = ifelse(contrl <= 0.05, -1, 1)

#pvals  = pvals * contrl
#pvals  = ifelse(pvals >= 0, pvals, 1)

# combine p-values, effect sizes, and gene names
df = data.frame(pvalues = as.vector(pvals), effect = as.vector(effect))
gene1 = colnames(effect)
gene1 = rep(gene1, each  = length(gene1))
gene2 = colnames(effect)
gene2 = rep(gene2, times = length(gene2))
df    = data.frame(df, gene1, gene2)
df    = df %>% filter(pvalues != 1)

# replace 0 p-values with the minimum p-value
df$pvalues = ifelse(df$pvalues == 0, min(df$pvalues[df$pvalues > 0]), df$pvalues)
df$pvalues = -log10(df$pvalues)

# effect size log transformed, ensuring correct sign after transformation
#df$sign   = sign(df$effect)
#df$effect = df$effect * df$sign
#df$effect = log2(df$effect)
#df$effect[is.na(df$effect) | is.nan(df$effect) | is.infinite(df$effect)] = 0
#df$effect = df$effect * df$sign

# add colors for enriched and depleted gene pairs
class = vector(mode = "character", length = nrow(df))
class[which(df$pvalues >= -log10(0.05) & df$effect >= 1)]  = "High"
class[which(df$pvalues >= -log10(0.05) & df$effect <= -1)] = "Low"
class[which(!class %in% c("High", "Low"))] = "Neither"
df$class = class

# construct labels for annotation
df$label = str_c(df$gene1, df$gene2, sep = "_")
df$label[which(!df$class %in% "High")] = ""

labels   = df[which(!df$label %in% ""),]
labels   = labels[order(labels$pvalues, decreasing = TRUE),]
top_n    = 2
labels   = labels[1:top_n,"label"]

df[which(!df$label %in% labels),"label"] = ""

# make volcano plot
ggplot(df) + aes(y = pvalues, x = effect, color = class) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = label), size = 2, color = "black", box.padding = 0.25, point.padding = 0.25, nudge_y = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotdash") +
  geom_vline(xintercept = c(1,-1), linetype = "dotdash") +
  labs(x = bquote("Effect size, cross-expression standard deviation from co-expression"),
       y = bquote("Cross-expression, -log"[10]*" p-value\n(Co-expression not controlled)"), color = "Cross-expression") +
  scale_color_manual(values = c("High" = "#80CC6D", "Low" = "#FFB347", "Neither" = "gray88"),
                     breaks = c("High", "Low")) +
  scale_x_continuous(limits = c(min(df$effect), -1 * min(df$effect))) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(legend.position = c(0.15,0.85),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/example.svg", device = "svg", height = 4, width = 6)

```