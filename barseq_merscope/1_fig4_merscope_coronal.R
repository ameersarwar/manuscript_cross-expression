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

source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
coord = rotate_coordinates(x = metadata$pos_x, y = metadata$pos_y, n_degrees = 40, origin = TRUE)
metadata = data.frame(sample_id = metadata$sample_id, coord, CCFparentname = metadata$CCFparentname)

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
  sc <- t(assay(sc, "X"))
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

# compare spatial and single-cell co-expression
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
spatial = correlation(data)
single  = correlation(df)

# plot spatial and single-cell co-expression against each other
spatial = upper_tri(spatial)
single  = upper_tri(single)

co_exp  = data.frame(spatial, single)

ggplot(co_exp) + aes(x = single, y = spatial) +
  geom_point(size = 0.1, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(subtitle = str_c("Pearson's R = ", signif(cor(single, spatial, use = "complete.obs", method = "pearson"), digits = 2)),
       x = "Co-expression in region matched scRNA-seq data", y = "Co-expression in MERSCOPE spatial data") +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave(filename = "/Users/AmeerSarwar/Downloads/cell_segment.png", device = "png", width = 6, height = 4, dpi = 600)

```


```{r}

# we see similar co-expression in spatial and scRNA-seq data
# double check that average expression of each gene is similar across datasets

# plot average counts of each gene in spatial and scRNA-seq data
spatial = colmeans(as.matrix(data))
single  = colmeans(as.matrix(df))

avg_counts = data.frame(spatial, single)

ggplot(avg_counts) + aes(x = single, y = spatial) +
  geom_point(color = "steelblue", size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10(labels = label_number(accuracy = 0.0001)) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  labs(subtitle = str_c("Spearman's Rho = ", signif(cor(single, spatial, method = "spearman"), digits = 2)),
       x = "Average counts in region matched scRNA-seq data", y = "Average counts in MERSCOPE spatial data") +
  theme_classic() +
  theme(plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ggsave(filename = "/Users/AmeerSarwar/Downloads/sc_spatial_average_merscope.svg", device = "svg", width = 6, height = 4, dpi = 600)  

# plot ranks
zz = data.frame(x = rank(single), y = rank(spatial))
ggplot(zz) + aes(x = x, y = y) + geom_point(color = "steelblue") + geom_abline(intercept = 0, slope = 1) + theme_classic() +
  labs(subtitle = str_c("Spearman's Rho = ", signif(cor(rank(single), rank(spatial), method = "spearman"), digits = 2)),
       x = "Rank of average gene counts in scRNA-seq",
       y = "Rank of average gene counts in MERSCOPE spatial data")

```


```{r}

# compare single-cell vs single-nucleus co-expression
# this is to check the upper limit of correlation of between technology co-expression
nucleus = c("ACA","AUD","ENT","MOp","RSP","VIS","VISp","STRd","CTX")
cell    = c("ACA","AUD","ENT","MOp","RSP","VIS","VISp","STR-STRd","CTXsp-CLA-EP-LA-BLA-BMA-PA")

# load scRNA-seq data
single_cell_regions = str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.", cell,".h5ad")
cell <- NULL

for (i in 1:length(single_cell_regions)){
  
  sc <- zellkonverter::readH5AD(single_cell_regions[i])
  sc <- t(assay(sc,"X"))
  cell <- rbind(cell, sc)
  print(i/length(single_cell_regions))
}

# load snRNA-seq data
single_cell_regions = str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/macosko/20221007_atlas_500umis_mt-1pct.h5ad.barseq_and_Vizgen_genes.", nucleus,".h5ad")
nucleus <- NULL

for (i in 1:length(single_cell_regions)){
  
  sc <- zellkonverter::readH5AD(single_cell_regions[i])
  sc <- t(assay(sc,"X"))
  nucleus <- rbind(nucleus, sc)
  print(i/length(single_cell_regions))
}

# single-nucleus data is ensmblID so change to gene names
if (FALSE){
  
  library(httr)
  library(jsonlite)
  
  ensembl_id <- colnames(nucleus)
  gene_name  <- c()

  for (i in 1:length(ensembl_id)){
    query_url <- paste0("https://rest.ensembl.org/lookup/id/", ensembl_id[i], "?content-type=application/json")
    response  <- GET(query_url)
    parsed_response <- fromJSON(rawToChar(response$content))
    gene_name[i] <- parsed_response$display_name
    print(i/length(ensembl_id))
  }
  write.csv(data.frame(gene = gene_name), "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/macosko/gene_names.csv")
}

# subset and order both datasets by genes in spatial data
gene_name = read.csv("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/macosko/gene_names.csv"); gene_name = gene_name$gene
colnames(nucleus) = gene_name

common = intersect(colnames(cell), colnames(data))
cell   = cell[, colnames(cell) %in% common]

common = intersect(colnames(nucleus), colnames(data))
nucleus= nucleus[, colnames(nucleus) %in% common]

cell   = cell[, colnames(cell) %in% colnames(nucleus)]
nucleus= nucleus[, colnames(nucleus) %in% colnames(cell)]

nucleus= nucleus[, colnames(cell)]

# prepare spatial data for the next chunk
temp_df = data
temp_df = temp_df[, colnames(temp_df) %in% colnames(nucleus)]
temp_df = temp_df[, colnames(nucleus)]
n_clear = nucleus

# compute co-expression for cell and nucleus datasets
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
cell    = correlation(cell)
nucleus = correlation(nucleus)

cell    = upper_tri(cell)
nucleus = upper_tri(nucleus)

corr = data.frame(cell, nucleus)

ggplot(corr) + aes(x = cell, y = nucleus) +
  geom_point(size = 0.1, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(subtitle = str_c("Pearson's R = ", signif(cor(cell, nucleus, use = "complete.obs", method = "pearson"), digits = 2)),
       x = "Co-expression in single-cell RNA-seq data", y = "Co-expression in single-nucleus RNA-seq data") +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave(filename = "/Users/AmeerSarwar/Downloads/cell_nucleus_coexpression.png", device = "png", dpi = 600, width = 6, height = 4)

```


```{r}

# compare single-nucleus and spatial co-expression
space = correlation(temp_df)
space = upper_tri(space)

# plot co-expression
corr = data.frame(x = nucleus, y = space)

ggplot(corr) + aes(x = x, y = y) +
  geom_point(size = 0.1, color = "steelblue") + geom_abline(intercept = 0, slope = 1) +
  labs(subtitle = str_c("Spearman's rho = ", signif(cor(space, nucleus, use = "complete.obs", method = "spearman"), digits = 2)),
       x = "Co-expression in snRNA-seq",
       y = "Co-expression in MERSCOPE spatial data") +
  theme_classic()

# compare their average gene expression
spa = colmeans(as.matrix(temp_df))
ncl = colmeans(n_clear)

# plot average gene counts
avg = data.frame(x = ncl, y = spa)

ggplot(avg) + aes(x = x, y = y) +
  geom_point(color = "steelblue") + geom_abline(intercept = 0, slope = 1) +
  scale_x_log10(labels = label_number(accuracy = 0.0001)) + scale_y_log10(labels = label_number(accuracy = 0.0001)) +
  labs(subtitle = str_c("Spearman's Rho = ", signif(cor(ncl, spa, method = "spearman"), digits = 2)),
       x = "Average gene counts in snRNA-seq",
       y = "Average gene counts in MERSCOPE spatial data") +
  theme_classic()

# plot ranks
zz = data.frame(x = rank(ncl), y = rank(spa))
ggplot(zz) + aes(x = x, y = y) + geom_point(color = "steelblue") + geom_abline(intercept = 0, slope = 1) + theme_classic() +
  labs(subtitle = str_c("Spearman's Rho = ", signif(cor(rank(ncl), rank(spa), method = "spearman"), digits = 2)),
       x = "Rank of average gene counts in snRNA-seq",
       y = "Rank of average gene counts in MERSCOPE spatial data")

```


```{r}

# !! skip this section, as the next one has data saved from running this part !!

# transcriptional bursting can be evaluated using another approach
# first, find gene pairs that are cross-expressed
# then find cell-neighbor pairs for each gene pair and do differential expression
# check how many genes are differentially expressed b/w cell-neighbor pairs
# if a lot, then those are cells of different types
# this can be augmented with single-cell data
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

cross = cross_expression(data = data, locations = metadata[,c("pos_x","pos_y")], alpha_co = 0.05)
cross = cross[as.logical(cross$cross_sig), ]

total_DE = rep(NaN, nrow(cross))

for (i in 1:nrow(cross)){
  
  # add gene counts to metadata
  meta  = metadata
  genes = data[, colnames(data) %in% c(cross$gene1[i], cross$gene2[i])]
  colnames(genes) = c("gene1", "gene2")
  meta  = data.frame(meta, genes)
  
  # metadata for neighbors
  distances <- RANN::nn2(data = meta[,c("pos_x","pos_y")], query = meta[,c("pos_x","pos_y")], k = 2, searchtype = "priority")
  distances <- distances$nn.idx[,2]
  meta_neig <- meta[distances,]
  
  # cross-expressing cell-neighbor pairs
  cell1 = meta$gene1 >  0 & meta$gene2 == 0
  cell2 = meta$gene1 == 0 & meta$gene2 >  0
  
  neig1 = meta_neig$gene1 == 0 & meta_neig$gene2 >  0
  neig2 = meta_neig$gene1 >  0 & meta_neig$gene2 == 0
  
  # gene A to B cell-neighbor pairs
  pairs = as.logical(cell1 * neig1)
  cell1 = meta[pairs, ]
  neig1 = meta_neig[pairs, ]
  
  # add sample ID to expression matrices
  dat   = data.frame(sample_id = meta$sample_id, data)
  
  # compile expression matrices
  cell1 = dat[dat$sample_id %in% unique(cell1$sample_id), ]
  neig1 = dat[dat$sample_id %in% unique(neig1$sample_id), ]
  
  # gene B to A cell-neighbor pairs
  pairs = as.logical(cell2 * neig2)
  cell2 = meta[pairs, ]
  neig2 = meta_neig[pairs, ]
  
  # compile expression matrices
  cell2 = dat[dat$sample_id %in% unique(cell2$sample_id), ]
  neig2 = dat[dat$sample_id %in% unique(neig2$sample_id), ]
  
  # differential analyses
  min_cells = 2   # ignore if number of cells/neighbors is less than a threshold (power)
  
  # A to B DE
  AB_pvals  = rep(NaN, ncol(cell1)-1)
  
  for (j in 2:ncol(cell1)){
    
    # ignore if number of cells/neighbors is less than a threshold (power)
    if (nrow(cell1) < min_cells | nrow(neig1) < min_cells) {break}
    
    # ignore gene if standard deviation is 0 in cell/neighbor data
    if (sd(cell1[,j]) == 0 | sd(neig1[,j]) == 0) {next}
    
    # perform DE
    pp = t.test(cell1[,j], neig1[,j])
    AB_pvals[j-1] = pp$p.value
  }
  
  # B to A DE
  BA_pvals  = rep(NaN, ncol(cell2)-1)
  
  for (j in 2:ncol(cell2)){
    
    # ignore if number of cells/neighbors is less than a threshold (power)
    if (nrow(cell2) < min_cells | nrow(neig2) < min_cells) {break}
    
    # ignore gene if standard deviation is 0 in cell/neighbor data
    if (sd(cell2[,j]) == 0 | sd(neig2[,j]) == 0) {next}
    
    # perform DE
    pp = t.test(cell2[,j], neig2[,j])
    BA_pvals[j-1] = pp$p.value
  }
  
  # combine DE results
  pvalues = c(AB_pvals, BA_pvals)
  pvalues = pvalues[!is.nan(pvalues)]
  pvalues = p.adjust(pvalues, method = "BH")
  total_DE[i] = sum(pvalues <= 0.05)
  
  # show proportion complete
  print(signif(i/nrow(cross), digits = 2))
}

# save file
write.csv(total_DE, file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Analysis/Validation 4/merscope_DE.csv")

```


```{r}

# load file
DE = read.csv(file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Analysis/Validation 4/merscope_DE.csv")
DE = DE$x
#DE = DE[DE != 0]
DE = data.frame(DE)

ggplot(DE) + aes(x = DE, after_stat(..scaled..)) +
  geom_density(alpha = 0.5, color = "darkblue", fill = "steelblue") +
  scale_x_continuous(breaks = seq(0,600,100), labels = seq(0,600,100)) +
  labs(x = "DE genes between cell-neighbor pairs",
       y = "Density",
       title = "MERSCOPE data") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        title = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ggsave(filename = "/Users/AmeerSarwar/Downloads/DE_genes.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# develop two algorithms: cross-expression as correlation and differential cross-expression
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

corr = cross_expression_correlation(data, locations = metadata[,c("pos_x","pos_y")], output_matrix = TRUE)

N    = 20
vals = c(upper_tri(corr))
top  = sort(vals, decreasing = TRUE )[1:N]
top  = data.frame(top, row = 0, col = 0)

for (i in 1:N){
  
  # top
  id = which(corr == top$top[i], arr.ind = TRUE)
  top$row[i] = id[1]
  top$col[i] = id[2]
}

# tissue expression plots
for (i in 1:N){
  
  # top
  p = tissue_expression_plot(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = colnames(data)[top$row[i]], gene2 = colnames(data)[top$col[i]], point_size = 0.1, cross_expression = FALSE); print(p)
  p = tissue_expression_plot(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = colnames(data)[top$row[i]], gene2 = colnames(data)[top$col[i]], point_size = 0.1, cross_expression = TRUE) ; print(p)
  
  p = spatial_enrichment(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = colnames(data)[top$row[i]], gene2 = colnames(data)[top$col[i]])
  print(p$pvalue)
  plot(density(p$target))
  plot(density(p$null))
  
  #print(i/N)
}

```


```{r}

# find all anatomical marker gene pairs
cross   = cross_expression(data = data, locations = metadata[,c("pos_x","pos_y")])
cross   = cross[cross$cross_sig == 1,]
markers = data.frame(gene1 = cross$gene1, gene2 = cross$gene2, pvalue = 1)

for (i in 1:nrow(markers)){
  pp = spatial_enrichment(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = markers$gene1[i], gene2 = markers$gene2[i])
  markers[i,"pvalue"] = pp$pvalue
  print(i/nrow(markers))
}

mark = markers
mark$pvalue = p.adjust(mark$pvalue, method = "BH")

mark = mark[mark$pvalue <= 0.05,]

for (i in 1:nrow(mark)){
  p = tissue_expression_plot(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = mark$gene1[i], gene2 = mark$gene2[i], point_size = 0, scale_bar = 0)
  print(p)
  print(i/nrow(mark))
}

# !!! SORT BY SPATIAL ENRICHMENT P-VALUES AND PLOT AGAIN !!!

```



```{r}

# !!! DIFFERENTIAL CROSS-EXPRESSION ALGORITHM !!!




```