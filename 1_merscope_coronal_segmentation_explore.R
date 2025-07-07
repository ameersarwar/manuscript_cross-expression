```{r}

# MERFISH brain receptor map
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate","scales",
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
data <- data[,2:ncol(data)]; data = as(as.matrix(data), "sparseMatrix")

```


```{r}

# single-cell RNA-seq data to assess cell segmentation errors in spatial data
# this is because the cells in scRNA-seq data are intact and no segmentation errors exist

# subset scRNA-seq by regions present in the spatial data
spatial_regions = sort(unique(metadata$CCFparentname))
single_cell_regions = c("AId-AIv-AIp","TH-AD-AV-AM-IAD-LD","CTXsp-CLA-EP-LA-BLA-BMA-PA","HIP-CA","STR-sAMY","STR-LSX","STR-STRd","STR-STRv","OLF-COA-PAA-NLOT-TR","HIP",
                        "TH-AD-AV-AM-IAD-LD","TH-LGd-IGL-LGv","TH-MD-IMD-PCN-CL","TH-MH-LH-LP","TH-PO-Eth","TH-PVT-PT","TH-RE-RH-CM-IAM-SMT-PR-Xi","CTXsp-CLA-EP-LA-BLA-BMA-PA",
                        "TH-MH-LH-LP","HIP","TH-PO-Eth","HY-MEZ-PVZ-PVR","TH-PVT-PT","OLF-AON-TT-DP-PIR-COA-NLOT-PAA-TR","OLF-MOB-AOB","HY-MEZ-PVZ-PVR","RSP","SSp","SSs-GU-VISC-AIp",
                        "STR-STRd","TH-VAL-VPM-VPL-VM","VISa","SS-GU-VISC","TH-PF-SPA-SPFm-VPMpc-VPLpc-RT","TH-VAL-VPM-VPL-VM")
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
common = intersect(colnames(df), colnames(data))
df = df[, colnames(df) %in% common]
data = data[, colnames(data) %in% common]
df = df[, colnames(data)]

```


```{r}

# average gene counts across single-cell vs spatial data
spatial = colmeans(as.matrix(data))
single  = colmeans(as.matrix(df))

counts = data.frame(gene = colnames(data), spatial, single)

# plot counts for each gene across single-cell vs spatial data
ggplot(counts) + aes(x = single, y = spatial) +
  geom_point(color = "steelblue") +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Average counts per gene", subtitle = str_c("Spearman's rho = ", signif(cor(counts$spatial, counts$single, method = "spearman"), digits = 2))) +
  theme_classic()

```


```{r}

# some genes have higher spatial than single-cell expression (false positives), increases co-expression
# other genes have lower spatial than single-cell expression (false negatives), decreases co-expression

# normalize by max of each dataset
spatial = spatial / max(spatial)
single  = single  / max(single)

# take difference b/w spatial and single-cell
diff = data.frame(gene = colnames(data), diff = spatial - single)
ggplot(diff) + aes(x = diff) + geom_density() + labs(x = "Spatial - single-cell counts, normalized")

# make groupings based on a tunable parameter
# parameter value can be chosen diagnostically using the distribution
delta_pos = 0.01 # delta for detection/difference, positive/negative
delta_neg = 0.01
delta_neg = delta_neg * -1

pos = diff[diff$diff >= delta_pos,]
neg = diff[diff$diff <= delta_neg,]
btw = diff[diff$diff > delta_neg & diff$diff < delta_pos, ]

# average counts per gene for each category across datasets

# false positives
cc_spatial = data[, colnames(data) %in% pos$gene]
cc_single  = df[, colnames(df) %in% pos$gene]

cc_average = data.frame(spatial = colMeans(cc_spatial), single = colMeans(cc_single))
ggplot(cc_average) + aes(x = single, y = spatial) + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(title = "False positives") +
  scale_x_continuous(limits = c(0, max(cc_average$single))) + scale_y_continuous(limits = c(0, max(cc_average$spatial)))

# false negatives
cc_spatial = data[, colnames(data) %in% neg$gene]
cc_single  = df[, colnames(df) %in% neg$gene]

cc_average = data.frame(spatial = colMeans(cc_spatial), single = colMeans(cc_single))
ggplot(cc_average) + aes(x = single, y = spatial) + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(title = "False negatives") +
  scale_x_continuous(limits = c(0, max(cc_average$single))) + scale_y_continuous(limits = c(0, max(cc_average$spatial)))

# false negatives
cc_spatial = data[, colnames(data) %in% btw$gene]
cc_single  = df[, colnames(df) %in% btw$gene]

cc_average = data.frame(spatial = colMeans(cc_spatial), single = colMeans(cc_single))
ggplot(cc_average) + aes(x = single, y = spatial) + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(title = "In-betweens") +
  scale_x_continuous(limits = c(0, max(cc_average$single))) + scale_y_continuous(limits = c(0, max(cc_average$spatial)))

```


```{r}

# co-expression for each category in spatial vs single-cell data

# false positives
pos_spatial = data[, colnames(data) %in% pos$gene]
pos_single  = df[, colnames(df)   %in% pos$gene]

pos_spatial = cor(as.matrix(pos_spatial), method = "spearman")
pos_single  = cor(as.matrix(pos_single), method = "spearman")

pos = data.frame(spatial = upper_tri(pos_spatial), single = upper_tri(pos_single))

ggplot(pos) + aes(x = single, y = spatial) + geom_point() + geom_abline(intercept = 0, slope = 1)

# false negatives
neg_spatial = data[, colnames(data) %in% neg$gene]
neg_single  = df[, colnames(df)   %in% neg$gene]

neg_spatial = cor(as.matrix(neg_spatial), method = "spearman")
neg_single  = cor(as.matrix(neg_single), method = "spearman")

neg = data.frame(spatial = upper_tri(neg_spatial), single = upper_tri(neg_single))

ggplot(neg) + aes(x = single, y = spatial) + geom_point() + geom_abline(intercept = 0, slope = 1)

# in-between
btw_spatial = data[, colnames(data) %in% btw$gene]
btw_single  = df[, colnames(df)   %in% btw$gene]

btw_spatial = cor(as.matrix(btw_spatial), method = "spearman")
btw_single  = cor(as.matrix(btw_single), method = "spearman")

btw = data.frame(spatial = upper_tri(btw_spatial), single = upper_tri(btw_single))

ggplot(btw) + aes(x = single, y = spatial) + geom_point() + geom_abline(intercept = 0, slope = 1)

```


```{r}

# high spatial vs single-cell co-expression means larger cell boundaries
# low spatial vs single-cell co-expression means smaller cell boundaries
# this is in relation to the in-between gene pairs, i.e., w/o false positives and false negatives

btw_corrs = btw
btw_genes = diff[diff$diff > delta_neg & diff$diff < delta_pos, ]

# !!TO DO!!

# identify gene pairs with very high spatial to single-cell co-expression (above diagonal)
# find total counts of these genes across cells
# cells with very high counts are considered as too big

# identify gene pairs with very low spatial to single-cell co-expression (below diagonal)
# find total counts of these genes across cells
# cells with very low counts are considered too small

```