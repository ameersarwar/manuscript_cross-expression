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

# process the datasets
datasets  = ls()[str_detect(ls(), "control")]
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
remove(meta_data, df_data, control1, control2, control3, control4, df, datasets)

# we viewed each slice in each brain and chose reasonable looking slices from the middle
brain1 = c(19,21,22); brain2 = c(15,16); brain3 = c(18,19); brain4 = c(15,18)
slices = list(brain1, brain2, brain3, brain4)

for (i in 1:length(metadata)){
  metadata[[i]] = metadata[[i]][metadata[[i]]$slice %in% slices[[i]], ]
  data[[i]] = data[[i]][rownames(data[[i]]) %in% rownames(metadata[[i]]), ]
}

```


```{r}

# single-cell RNA-seq data to assess cell segmentation errors in spatial data
# this is because the cells in scRNA-seq data are intact and no segmentation errors exist
# we do this separately for each brain

# brain 1
meta  = metadata$control1
df_sp = data$control1

# subset scRNA-seq by regions present in the spatial data
spatial_regions = sort(table(meta$CCFname), decreasing = TRUE)
spatial_regions = sort(names(spatial_regions[spatial_regions >= 50]))
single_cell_regions = c("STR-sAMY","ACA","TH-AD-AV-AM-IAD-LD","HY-MEZ-PVZ-PVR","AId-AIv-AIp","CTXsp-CLA-EP-LA-BLA-BMA-PA","HIP-CA","OLF-COA-PAA-NLOT-TR","STR-STRd","HIP","TH-MD-IMD-PCN-CL","TH-MH-LH-LP","HY-LZ","MOp","MOs-FRP",
                        "OLF-PIR","TH-PVT-PT","TH-RE-RH-CM-IAM-SMT-PR-Xi","RSP","TH-PF-SPA-SPFm-VPMpc-VPLpc-RT","SSp","SSs-GU-VISC-AIp","STR-STRd","STR-STRv","TH","TH-VAL-VPM-VPL-VM","SS-GU-VISC")

# remove regions w/o corresponding single-cell regions
# also remove fiber tracts, which are denoted with all lowercase letters
remove_regions = c("CL","CM","GPe","GPi","LA","OLF","PAL","PO","PR","RCH","SI","TRS","V3","VL")
spatial_regions = spatial_regions[!str_detect(spatial_regions, "^[a-z]+$")]
spatial_regions = spatial_regions[!spatial_regions %in% remove_regions]

# subset data by these spatial regions
meta  = meta[meta$CCFname %in% spatial_regions, ]
df_sp = df_sp[rownames(df_sp) %in% rownames(meta), ]

# upload single-cell data
single_cell_regions = unique(single_cell_regions)
single_cell_regions = str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.",
      single_cell_regions,".h5ad")

# load scRNA-seq data
df <- NULL

for (i in 1:length(single_cell_regions)){
  
  sc <- zellkonverter::readH5AD(single_cell_regions[i])
  sc <- t(assay(sc,"X"))
  df <- rbind(df, sc)
  print(i/length(single_cell_regions))
}

# subset by common genes
common = intersect(colnames(df), colnames(df_sp))
df = df[, colnames(df) %in% common]
df_sp = df_sp[, colnames(df_sp) %in% common]
df = df[, colnames(df_sp)]

```


```{r}

# compare spatial and single-cell co-expression
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
spatial = correlation(df_sp)
single  = correlation(df)

# plot spatial and single-cell co-expression against each other
spatial = upper_tri(spatial)
single  = upper_tri(single)

co_exp  = data.frame(spatial, single)

ggplot(co_exp) + aes(x = single, y = spatial) +
  geom_point(size = 0.1, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(subtitle = str_c("Correlation = ", signif(cor(single, spatial, use = "complete.obs"), digits = 2)),
       x = "Co-expression in region matched scRNA-seq data from Zeng", y = "Co-expression in BAR-seq spatial data") +
  theme_bw()

```


```{r}

# we see low co-expression in spatial data than in scRNA-seq data
# this could be because of low counts in spatial data, so genes are less likely to be co-present

# plot average counts of each gene in spatial and scRNA-seq data
spatial = colmeans(as.matrix(df_sp))
single  = colmeans(as.matrix(df))

avg_counts = data.frame(spatial, single)

ggplot(avg_counts) + aes(x = single, y = spatial) +
  geom_point(color = "steelblue", size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(subtitle = str_c("Correlation = ", signif(cor(single, spatial), digits = 2)),
       x = "Average counts in region matched scRNA-seq data from Zeng", y = "Average counts in BAR-seq spatial data") +
  theme_bw()

```

