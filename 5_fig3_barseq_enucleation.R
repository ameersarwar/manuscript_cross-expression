```{r}

# enucleation data for controls and experiments
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo","igraph")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))
registerDoParallel(detectCores() - 1)
plan("multisession", workers = 8)

# controls
control1 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Control_1.mat.h5ad.rds")
control2 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Control_2.mat.h5ad.rds")
control3 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Control_3.mat.h5ad.rds")
control4 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Control_4.mat.h5ad.rds")

# experiments
enucleated1 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Enucleated_1.mat.h5ad.rds")
enucleated2 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Enucleated_2.mat.h5ad.rds")
enucleated3 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Enucleated_3.mat.h5ad.rds")
enucleated4 = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/Enucleation/pre_processed/Enucleated_4.mat.h5ad.rds")

# process the datasets
datasets  = ls()[str_detect(ls(), "control|enucleated")]
meta_data = vector(mode = "list", length = length(datasets))
df_data   = vector(mode = "list", length = length(datasets))
names(meta_data) = datasets
names(df_data)   = datasets

for (i in 1:length(datasets)){
  
  # extract dataset
  df       = get(datasets[i])
  
  # data processing
  data     = assays(df)
  data     = t(data$counts)
  data     = data[rowSums(data) >= 20 & rowSums(data > 0) >= 5, str_detect(colnames(data), "unused", negate = TRUE)]
  
  # metadata processing
  metadata = colData(df)
  metadata = as.data.frame(metadata)
  metadata = metadata[rownames(metadata) %in% rownames(data),]
  metadata = metadata %>% dplyr::select(pos_0, pos_1, slice, CCFname, clustid, clustname, CCF_AP_axis, CCF_ML_axis, CCF_DV_axis)
  colnames(metadata)[1:2] = c("pos_x","pos_y")
  
  # store datasets in list objects
  df_data[[i]]   = data
  meta_data[[i]] = metadata
  print(i/length(datasets))
}

# change names
metadata = meta_data
data     = df_data

# store backups
data_subset = data; metadata_subset = metadata

# remove heavy variables
remove(meta_data, df_data, control1, control2, control3, control4, enucleated1, enucleated2, enucleated3, enucleated4, df, datasets)

```


```{r}

# CONTROLS: predict network structure of slices in other datasets

# slices common between all datasets
metadata_controls = metadata_subset[1:(length(metadata_subset) / 2)]
data_controls = data_subset[1:(length(data_subset) / 2)]
controls = names(data_controls)

slices = vector(mode = "list", length = length(metadata_controls))
names(slices) = names(metadata_controls)

for (i in 1:length(metadata_controls)){
  slices[[i]] = sort(unique(metadata_controls[[i]]$slice))
}

slices = Reduce(intersect, slices)

# subset data by common slices
for (i in 1:length(metadata_controls)){
  metadata_controls[[i]] = metadata_controls[[i]][metadata_controls[[i]]$slice %in% slices,]
  data_controls[[i]] = data_controls[[i]][rownames(data_controls[[i]]) %in% rownames(metadata_controls[[i]]),]
}

# leave one out cross-validation
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
auroc    = c()

for (i in 1:length(controls)){
  
  # leave one out and remaining datasets
  left_out = controls[i]
  rest     = controls[!controls %in% left_out]
  
  # cross-expression left out
  cross_slice = matrix(data = 0, nrow = choose(ncol(data[[1]]), 2), ncol = length(slices))
  colnames(cross_slice) = slices
  
  # cross-expression per slice
  for (j in 1:length(slices)){
    
    locations = metadata_controls[[i]][metadata_controls[[i]]$slice %in% slices[j],c("pos_x","pos_y")]
    expr_data = data_controls[[i]][rownames(data_controls[[i]]) %in% rownames(locations),]
    
    expr_mat  = cross_expression(data = expr_data, locations = locations)
    expr_mat  = expr_mat$cross_pvalue
    cross_slice[,j] = expr_mat
  }
  
  cross_slice_left = cross_slice
  
  # cross-expression rest
  cross_slice = matrix(data = 0, nrow = choose(ncol(data[[1]]), 2), ncol = length(slices))
  colnames(cross_slice) = slices
  cross_slice_rest = cross_slice
  
  for (k in 1:length(rest)){
    
    # cross-expression per brain
    rest_temp = rest[k]
    
    # cross-expression per slice
    for (j in 1:length(slices)){
      
      locations = metadata_controls[[rest_temp]][metadata_controls[[rest_temp]]$slice %in% slices[j], c("pos_x","pos_y")]
      expr_data = data_controls[[rest_temp]][rownames(data_controls[[rest_temp]]) %in% rownames(locations),]
      
      expr_mat  = cross_expression(data = expr_data, locations = locations)
      expr_mat  = expr_mat$cross_pvalue
      cross_slice[,j] = expr_mat
    }
    
    # combine cross-expression per brain
    cross_slice_rest = cross_slice_rest + cross_slice
  }
  
  # average cross-expression across brains
  cross_slice_rest = cross_slice_rest / length(rest); remove(cross_slice)
  
  # left out vs rest slice by slice correlation
  corr = correlation(cross_slice_left, cross_slice_rest)
  
  # replace row values by ranks
  corr = apply(corr, 1, order)
  
  # AUROC
  auroc = c(auroc, (diag(corr) - 1) / (length(diag(corr)) - 1))
  print(i/length(controls))
}

auroc_ctrl = auroc

```


```{r}

# ENUCLEATED: predict network structure of slices in other datasets

# slices common between all datasets
metadata_controls = metadata_subset[((length(metadata_subset) / 2) + 1):length(metadata_subset)]
data_controls = data_subset[((length(metadata_subset) / 2) + 1):length(metadata_subset)]
controls = names(data_controls)

slices = vector(mode = "list", length = length(metadata_controls))
names(slices) = names(metadata_controls)

for (i in 1:length(metadata_controls)){
  slices[[i]] = sort(unique(metadata_controls[[i]]$slice))
}

slices = Reduce(intersect, slices)

# subset data by common slices
for (i in 1:length(metadata_controls)){
  metadata_controls[[i]] = metadata_controls[[i]][metadata_controls[[i]]$slice %in% slices,]
  data_controls[[i]] = data_controls[[i]][rownames(data_controls[[i]]) %in% rownames(metadata_controls[[i]]),]
}

# leave one out cross-validation
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
auroc    = c()

for (i in 1:length(controls)){
  
  # leave one out and remaining datasets
  left_out = controls[i]
  rest     = controls[!controls %in% left_out]
  
  # cross-expression left out
  cross_slice = matrix(data = 0, nrow = choose(ncol(data[[1]]), 2), ncol = length(slices))
  colnames(cross_slice) = slices
  
  # cross-expression per slice
  for (j in 1:length(slices)){
    
    locations = metadata_controls[[i]][metadata_controls[[i]]$slice %in% slices[j],c("pos_x","pos_y")]
    expr_data = data_controls[[i]][rownames(data_controls[[i]]) %in% rownames(locations),]
    
    expr_mat  = cross_expression(data = expr_data, locations = locations)
    expr_mat  = expr_mat$cross_pvalue
    cross_slice[,j] = expr_mat
  }
  
  cross_slice_left = cross_slice
  
  # cross-expression rest
  cross_slice = matrix(data = 0, nrow = choose(ncol(data[[1]]), 2), ncol = length(slices))
  colnames(cross_slice) = slices
  cross_slice_rest = cross_slice
  
  for (k in 1:length(rest)){
    
    # cross-expression per brain
    rest_temp = rest[k]
    
    # cross-expression per slice
    for (j in 1:length(slices)){
      
      locations = metadata_controls[[rest_temp]][metadata_controls[[rest_temp]]$slice %in% slices[j], c("pos_x","pos_y")]
      expr_data = data_controls[[rest_temp]][rownames(data_controls[[rest_temp]]) %in% rownames(locations),]
      
      expr_mat  = cross_expression(data = expr_data, locations = locations)
      expr_mat  = expr_mat$cross_pvalue
      cross_slice[,j] = expr_mat
    }
    
    # combine cross-expression per brain
    cross_slice_rest = cross_slice_rest + cross_slice
  }
  
  # average cross-expression across brains
  cross_slice_rest = cross_slice_rest / length(rest); remove(cross_slice)
  
  # left out vs rest slice by slice correlation
  corr = correlation(cross_slice_left, cross_slice_rest)
  
  # replace row values by ranks
  corr = apply(corr, 1, order)
  
  # AUROC
  auroc = c(auroc, (diag(corr) - 1) / (length(diag(corr)) - 1))
  print(i)
}

auroc_expr = auroc

```


```{r}

# plot AUROC distributions
df = rbind(data.frame(auroc = auroc_ctrl, group = rep("Control", length(auroc_ctrl))),
           data.frame(auroc = auroc_expr, group = rep("Enucleated", length(auroc_expr))))
df$group = factor(df$group, levels = c("Control","Enucleated"))

ggplot(df) + aes(x = auroc, fill = group, color = group) +
  geom_density(alpha = 0.5) +
  facet_wrap(~group, ncol = 1) +
  labs(x = "AUROC", y = "Density") +
  scale_y_continuous(breaks = c(0,1,2), labels = c(0,1,2)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        strip.text = element_text(size = 12))

```
