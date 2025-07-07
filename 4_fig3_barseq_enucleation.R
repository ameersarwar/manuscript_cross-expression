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
}

# change names
metadata = meta_data
data     = df_data

# store backups
DATA = data; METADATA = metadata

# remove heavy variables
remove(meta_data, df_data, control1, control2, control3, control4, enucleated1, enucleated2, enucleated3, enucleated4, df, datasets)

```


```{r}

# gene pairs with increased, decreased, unchanged, or non-significant cross-expression b/w controls vs enucleated
# comparison across lateral and non-lateral visual regions
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

anterior          = c("VISa1","VISa2/3","VISa4","VISa5","VISa6a","VISa6b")       # anterior visual area
anterior_medial   = c("VISam1","VISam2/3","VISam4","VISam5","VISam6a","VISam6b") # anterior medial visual area
primary           = c("VISp1","VISp2/3","VISp4","VISp5","VISp6a","VISp6b")       # primary visual area
posterior_medial  = c("VISpm1","VISpm2/3","VISpm4","VISpm5","VISpm6a","VISpm6b") # posterior medial visual area
anterior_lateral  = c("VISal1","VISal2/3","VISal4","VISal5","VISal6a","VISal6b") # anterior lateral visual area
lateral           = c("VISl1","VISl2/3","VISl4","VISl5","VISl6a","VISl6b")       # lateral visual area
posterior_lateral = c("VISpl1","VISpl2/3","VISpl4","VISpl5","VISpl6a","VISpl6b") # posterior lateral visual area
rostral_lateral   = c("VISrl1","VISrl2/3","VISrl4","VISrl5","VISrl6a","VISrl6b") # rostral lateral visual area

regions = list(anterior, anterior_medial, primary, posterior_medial,
               anterior_lateral, lateral, posterior_lateral, rostral_lateral)

names(regions) = c("anterior", "anterior_medial", "primary", "posterior_medial",
                   "anterior_lateral", "lateral", "posterior_lateral", "rostral_lateral")

cross_total = vector(mode = "list", length = length(regions))
names(cross_total) = names(regions)

for (i in 1:length(regions)){
  
  # re-extract data for each region
  metadata_subset = metadata
  data_subset     = data
  region_current  = regions[[i]]
  pvals           = matrix(data = 0, nrow = choose(ncol(data$control1), 2), ncol = length(data))
  colnames(pvals) = names(data)
  
  # cross-expression for controls and experiments per region
  for (j in 1:length(data)){
    
    # subset re-extracted data by regions
    metadata_subset[[j]] = metadata_subset[[j]][metadata_subset[[j]]$CCFname %in% region_current, ]
    data_subset[[j]]     = data_subset[[j]][rownames(data_subset[[j]]) %in% rownames(metadata_subset[[j]]), ]
    
    # cross-expression
    cross = cross_expression(data = data_subset[[j]], locations = metadata_subset[[j]][,c("pos_x","pos_y")])
    
    # save gene names
    if (i == 1 & j == i){
      gene1 = cross$gene1
      gene2 = cross$gene2
    }
    
    # transform cross-expression
    cross = cross$cross_pvalue
    cross = -log10(cross)
    pvals[,j] = cross
  }
  
  # average cross-expression within controls and experiments
  ctrl = apply(pvals[, 1:(length(data) / 2)], MARGIN = 1, FUN = mean)
  expr = apply(pvals[, ((length(data)  / 2) + 1):length(data)], MARGIN = 1, FUN = mean)
  
  # classify gene pair based on control vs enucleated cross-expression
  class = vector(mode = "character", length = length(ctrl))

  class[ctrl >= -log10(0.05) & expr >= -log10(0.05)] = "Robust to enucleation"
  class[ctrl >= -log10(0.05) & expr  < -log10(0.05)] = "Reduced by enucleation"
  class[ctrl  < -log10(0.05) & expr >= -log10(0.05)] = "Enhanced by enucleation"
  class[ctrl  < -log10(0.05) & expr  < -log10(0.05)] = "Non-significant"
  
  # combine results
  cross_total[[i]] = data.frame(gene1, gene2, ctrl, expr, class, region = names(regions[i]))
  print(i/length(regions))
  
}

```


```{r}

# plot controls vs enucleation per gene pair per region
df = cross_total

for (i in 1:length(df)){
  
  # region
  df_temp = df[[i]]
  df_temp$class = factor(df_temp$class, levels = c("Enhanced by enucleation", "Reduced by enucleation", "Robust to enucleation", "Non-significant"))
  region  = unique(df_temp$region)
  
  # control vs enucleated per gene pair
  p = ggplot(df_temp) + aes(x = ctrl, y = expr, color = class) +
    geom_vline(xintercept = -log10(0.05), linetype = "dotdash", alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotdash", alpha = 0.8) +
    geom_point() +
    labs(x = "Controls average p-values (-log10)",
         y = "Enucleated average p-values (-log10)",
         color = "Cross-expression", title = region) +
    scale_color_manual(values = c("Non-significant" = "gray88", "Enhanced by enucleation" = "brown3",
                                  "Reduced by enucleation" = "deepskyblue4", "Robust to enucleation" = "chartreuse3")) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))
  print(p)
}

```


```{r}

# gene pairs with increased, decreased, unchanged, or non-significant cross-expression b/w controls vs enucleated
# comparison across lateral and non-lateral visual regions
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

anterior          = c("VISa1","VISa2/3","VISa4","VISa5","VISa6a","VISa6b")       # anterior visual area
anterior_medial   = c("VISam1","VISam2/3","VISam4","VISam5","VISam6a","VISam6b") # anterior medial visual area
primary           = c("VISp1","VISp2/3","VISp4","VISp5","VISp6a","VISp6b")       # primary visual area
posterior_medial  = c("VISpm1","VISpm2/3","VISpm4","VISpm5","VISpm6a","VISpm6b") # posterior medial visual area
anterior_lateral  = c("VISal1","VISal2/3","VISal4","VISal5","VISal6a","VISal6b") # anterior lateral visual area
lateral           = c("VISl1","VISl2/3","VISl4","VISl5","VISl6a","VISl6b")       # lateral visual area
posterior_lateral = c("VISpl1","VISpl2/3","VISpl4","VISpl5","VISpl6a","VISpl6b") # posterior lateral visual area
rostral_lateral   = c("VISrl1","VISrl2/3","VISrl4","VISrl5","VISrl6a","VISrl6b") # rostral lateral visual area

regions = list(anterior, anterior_medial, primary, posterior_medial,
               anterior_lateral, lateral, posterior_lateral, rostral_lateral)

names(regions) = c("anterior", "anterior_medial", "primary", "posterior_medial",
                   "anterior_lateral", "lateral", "posterior_lateral", "rostral_lateral")

cross_total = matrix(data = 0, nrow = choose(ncol(data$control1), 2), ncol = length(regions))
colnames(cross_total) = names(regions)

# cross-expression across controls and experiments per region
for (i in 1:length(regions)){
  
  # re-extract data for each region
  metadata_subset = metadata
  data_subset     = data
  region_current  = regions[[i]]
  pvals           = matrix(data = 0, nrow = choose(ncol(data$control1), 2), ncol = length(data))
  colnames(pvals) = names(data)
  
  # cross-expression for controls and experiments per region
  for (j in 1:length(data)){
    
    # subset re-extracted data by regions
    metadata_subset[[j]] = metadata_subset[[j]][metadata_subset[[j]]$CCFname %in% region_current, ]
    data_subset[[j]]     = data_subset[[j]][rownames(data_subset[[j]]) %in% rownames(metadata_subset[[j]]), ]
    
    # cross-expression
    cross = cross_expression(data = data_subset[[j]], locations = metadata_subset[[j]][,c("pos_x","pos_y")])
    
    # save gene names
    if (i == 1 & j == i){
      gene1 = cross$gene1
      gene2 = cross$gene2
    }
    
    # transform cross-expression
    cross = cross$cross_pvalue
    cross = -log10(cross)
    pvals[,j] = cross
  }
  
  # compare enucleated vs controls for each gene
  ctrl = pvals[, 1:(length(data) / 2)]
  expr = pvals[, ((length(data)  / 2) + 1):length(data)]
  pval = vector(mode = "numeric", length = nrow(ctrl))
  
  # for each gene pair, enucleated brains have less cross-expression than control brains
  for (k in 1:nrow(ctrl)){
    pval[k] = suppressWarnings(wilcox.test(x = expr[k,], y = ctrl[k,])$p.value)
  }
  
  # combine results
  pval = p.adjust(pval, method = "BH")
  pval = -log10(pval)
  cross_total[,i] = pval
  
  print(i/length(regions))
}

```


```{r}

# pre-process dataset
df = as.data.frame(cross_total)
df = df %>% pivot_longer(1:ncol(df), names_to = "region", values_to = "pvals") %>% as.data.frame()

# plot regions and p-values
ggplot(df) + aes(x = region, y = pvals, fill = region) +
  geom_boxplot()

```


```{r}

# gene pairs with increased, decreased, unchanged, or non-significant cross-expression b/w controls vs enucleated
# comparison across lateral and non-lateral visual regions

anterior          = c("VISa1","VISa2/3","VISa4","VISa5","VISa6a","VISa6b")       # anterior visual area
anterior_medial   = c("VISam1","VISam2/3","VISam4","VISam5","VISam6a","VISam6b") # anterior medial visual area
primary           = c("VISp1","VISp2/3","VISp4","VISp5","VISp6a","VISp6b")       # primary visual area
posterior_medial  = c("VISpm1","VISpm2/3","VISpm4","VISpm5","VISpm6a","VISpm6b") # posterior medial visual area
anterior_lateral  = c("VISal1","VISal2/3","VISal4","VISal5","VISal6a","VISal6b") # anterior lateral visual area
lateral           = c("VISl1","VISl2/3","VISl4","VISl5","VISl6a","VISl6b")       # lateral visual area
posterior_lateral = c("VISpl1","VISpl2/3","VISpl4","VISpl5","VISpl6a","VISpl6b") # posterior lateral visual area
rostral_lateral   = c("VISrl1","VISrl2/3","VISrl4","VISrl5","VISrl6a","VISrl6b") # rostral lateral visual area

regions = list(anterior, anterior_medial, primary, posterior_medial,
               anterior_lateral, lateral, posterior_lateral, rostral_lateral)

names(regions) = c("anterior", "anterior_medial", "primary", "posterior_medial",
                   "anterior_lateral", "lateral", "posterior_lateral", "rostral_lateral")

cross_total = vector(mode = "list", length = length(regions))
names(cross_total) = names(regions)

# cross-expression across controls and experiments per region
for (i in 1:length(regions)){
  
  # re-extract data for each region
  metadata_subset = metadata
  data_subset     = data
  region_current  = regions[[i]]
  pvals           = matrix(data = 0, nrow = choose(ncol(data$control1), 2), ncol = length(data))
  colnames(pvals) = names(data)
  
  # cross-expression for controls and experiments per region
  for (j in 1:length(data)){
    
    # subset re-extracted data by regions
    metadata_subset[[j]] = metadata_subset[[j]][metadata_subset[[j]]$CCFname %in% region_current, ]
    data_subset[[j]]     = data_subset[[j]][rownames(data_subset[[j]]) %in% rownames(metadata_subset[[j]]), ]
    
    # cross-expression
    cross = cross_expression(data = data_subset[[j]], locations = metadata_subset[[j]][,c("pos_x","pos_y")])
    
    # save gene names
    if (i == 1 & j == i){
      gene1 = cross$gene1
      gene2 = cross$gene2
    }
    
    # transform cross-expression
    cross = cross$cross_pvalue
    cross = -log10(cross)
    pvals[,j] = cross
  }
  
  # average cross-expression within controls and experiments
  ctrl = apply(pvals[, 1:(length(data) / 2)], MARGIN = 1, FUN = mean)
  expr = apply(pvals[, ((length(data)  / 2) + 1):length(data)], MARGIN = 1, FUN = mean)

  # classify gene pair based on control vs enucleated cross-expression
  class = vector(mode = "character", length = length(ctrl))

  class[ctrl >= -log10(0.05) & expr >= -log10(0.05)] = "Robust to enucleation"
  class[ctrl >= -log10(0.05) & expr  < -log10(0.05)] = "Reduced by enucleation"
  class[ctrl  < -log10(0.05) & expr >= -log10(0.05)] = "Enhanced by enucleation"
  class[ctrl  < -log10(0.05) & expr  < -log10(0.05)] = "Non-significant"
  
  # combine results
  cross_total[[i]] = data.frame(gene1, gene2, ctrl, expr, class, region = names(regions[i]))
  print(i/length(regions))
}

```


```{r}

# enhancement, reduction, non-change, and non-significance across non-lateral and lateral regions

# compute frequencies
class  = vector(mode = "character")
freq   = vector(mode = "numeric")
region = vector(mode = "character")

for (i in 1:length(cross_total)){
  outcome  = table(cross_total[[i]]$class)
  class    = c(class,  names(outcome))
  freq     = c(freq,   as.numeric(outcome))
  region   = c(region, rep(names(cross_total)[i], length(outcome)))
}

# combine frequencies
df = data.frame(region, class, freq)
df = df[!df$class %in% "Non-significant",]
df$class  = factor(df$class,  levels = c("Enhanced by enucleation","Robust to enucleation","Reduced by enucleation"))
df$region = factor(df$region, levels = c("anterior","anterior_medial","primary","posterior_medial","anterior_lateral","lateral","posterior_lateral","rostral_lateral"))

# average frequencies
non_lateral = df[1:(nrow(df) / 2), ]
lateral     = df[((nrow(df) / 2) + 1):nrow(df), ]

non_lateral_enhance = mean(non_lateral[non_lateral$class %in% "Enhanced by enucleation", "freq"])
non_lateral_reduce  = mean(non_lateral[non_lateral$class %in% "Reduced by enucleation", "freq"])
non_lateral_robust  = mean(non_lateral[non_lateral$class %in% "Robust to enucleation", "freq"])

lateral_enhance = mean(lateral[lateral$class %in% "Enhanced by enucleation", "freq"])
lateral_reduce  = mean(lateral[lateral$class %in% "Reduced by enucleation", "freq"])
lateral_robust  = mean(lateral[lateral$class %in% "Robust to enucleation", "freq"])

# pre-process data
non_lateral = data.frame(avg_freq = c(non_lateral_enhance, non_lateral_robust, non_lateral_reduce),
                     class    = c("Increased","Unchanged","Decreased"), region = "Non-lateral")

lateral = data.frame(avg_freq = c(lateral_enhance, lateral_robust, lateral_reduce),
                     class    = c("Increased","Unchanged","Decreased"), region = "Lateral")

df = rbind(lateral, non_lateral)
df$class  = factor(df$class,  levels = c("Increased","Unchanged","Decreased"))
df$region = factor(df$region, levels = c("Non-lateral","Lateral"))

# plot results
ggplot(df) + aes(x = class, y = avg_freq, fill = region) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "Number of gene pairs (averaged)", fill = "Visual cortex",
       x = "Cross-expression in enucleated versus in control groups") +
  scale_fill_manual(values = c("Non-lateral" = "deepskyblue4", "Lateral" = "brown3")) +
  theme_bw() +
  theme(legend.position.inside = c(0.12,0.85),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)))

```
