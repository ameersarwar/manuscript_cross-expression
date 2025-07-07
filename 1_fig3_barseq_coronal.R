```{r}

# BAR-seq coronal data
# data is shrunk by removing image stitching-related artefacts (cf. Xiaoyin's email)
# data is quality controlled by keeping cells with genes/cell >= 5 and reads/cell >= 20
# data alongside CCF and slide coordinates are saved and can be used for analysis

# load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo","ggrepel","igraph")
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
metadata <- coordinates

data = t(data)
coordinates = coord[coord$slice %in% slice_current,]
data = data[rownames(data) %in% coordinates$sample_id, ]

```


```{r}

# example of smoothing gene's expression values
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# choose region
region = c("AUDv","AUDp","AUDd")
metadata  = coordinates[coordinates$CCFparentname %in% region, ]
data_temp = data[rownames(data) %in% metadata$sample_id, ]

# layer 2/3 marker gene
gene = "Rorb"
data_temp = data_temp[,gene]

neighbor_to_smooth = c(0, 1, 2, 5)
smoothed = vector(mode = "list", length = length(neighbor_to_smooth))
names(smoothed) = neighbor_to_smooth

for (i in 1:length(smoothed)){
  smooth_data = smooth_cells(data = data_temp, locations = metadata[,c("pos_x","pos_y")], neighbors_smoothed = neighbor_to_smooth[i])
  smoothed[[i]] = smooth_data$smooth_expression
}

# process data for plotting
df = data.frame()
for (i in 1:length(smoothed)){
  df = rbind(df, data.frame(count = smoothed[[i]], class = str_c("Neighbors smoothed = ", neighbor_to_smooth[i])))
}

df[df$class %in% "Neighbors smoothed = 0","class"] = "Raw Counts"
df = data.frame(do.call(rbind, replicate(length(neighbor_to_smooth), metadata, simplify = FALSE)), df)
df = df[order(df$count),]
df$class = factor(df$class, levels = unique(df$class))

# plot these together
ggplot(df) + aes(x = pos_x, y = pos_y, color = count) +
  geom_point(size = 0.5) +
  facet_wrap(~class, nrow = 2, ncol = 2) +
  scale_color_gradient(low = "gray88", high = "red") +
  labs(x = "", y = "", color = "Rorb\nCounts") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#E0F7FA"))
ggsave(filename = "/Users/AmeerSarwar/Downloads/smoothed.png", device = "png", width = 7, height = 5, dpi = 300)

```
