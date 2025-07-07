```{r}

# load packages
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate","scales",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))

# load data
files = list.files(path = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/")
files = str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/", files)
files_data  = files[str_detect(files, "cell_by_gene")]
files_meta  = files[str_detect(files, "cell_metadata")]

data = data.frame()
metadata = data.frame()

for (i in 1:length(files_data)){
  
  # data
  df = as.data.frame(fread(files_data[i])); rownames(df) = df$V1; df = df[,2:ncol(df)]
  df = df[,str_detect(colnames(df), "Blank-", negate = TRUE)]
  df = df[rowSums(df) >= 50,]
  data = rbind(data, df)
  
  # metadata
  meta = as.data.frame(fread(files_meta[i])); meta = data.frame(sample_id = meta$V1, pos_x = meta$center_x, pos_y = meta$center_y)
  meta = meta[meta$sample_id %in% rownames(df),]
  type = str_extract(files_data[i], "Slice\\d+_Replicate\\d+")
  meta$type = rep(type, nrow(meta))
  metadata = rbind(metadata, meta)
  
  print(i/length(files_data))
}

remove(df, meta)

# view plots
ggplot(metadata) + aes(x = pos_x, y = pos_y) + geom_point(size=0) + facet_wrap(~type, scales = "free")

```


```{r}

# similarities in cross-expression b/w slices/replicates
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
type = unique(metadata$type)
simm_corr = data.frame(matrix(data = 0, nrow = choose(ncol(data), 2), ncol = length(type))); colnames(simm_corr) = type

for (i in 1:length(type)){
  simm_corr[,i] = cross_expression_correlation(data = data[metadata$type %in% type[i],], locations = metadata[metadata$type %in% type[i], c("pos_x","pos_y")])$correlation
  print(i/length(type))
}

```


```{r}

slices = unique(str_extract(type, "Slice\\d+"))

for (i in 1:length(slices)){
  
  df = simm_corr[, str_detect(colnames(simm_corr), slices[i])]; colnames(df) = c("x","y","z")
  
  if (i==1){app = str_c("Slice ", i, " (anterior)")}; if (i==2){app = str_c("Slice ", i, " (middle)")}; if (i==3){app = str_c("Slice ", i, " (posterior)")}
  
  p = ggplot(df) + aes(x = x, y = y, color = z) +
    geom_point(size = 0) + geom_abline(intercept = 0, slope = 1) +
    scale_color_viridis_c(option = "viridis") + theme_classic() +
    labs(title = str_c(app, ", cross-expression correlation between replicates"), x = "Replicate 1", y = "Replicate 2", color = "Replicate 3",
         subtitle = str_c("Average Spearman's rho between replicates = ", signif(mean(upper_tri(cor(df, method = "spearman"))), digits = 2))) +
    theme(legend.position = "inside", legend.position.inside = c(0.1,0.7))
  
  print(p)
  #ggsave(filename = str_c("/Users/AmeerSarwar/Downloads/replicate_corrs_",i,".png"), dpi = 600, device = "png", width = 6, height = 4)
}

```