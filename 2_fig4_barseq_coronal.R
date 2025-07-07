```{r}

# BAR-seq coronal data
# data is shrunk by removing image stitching-related artefacts (cf. Xiaoyin's email)
# data is quality controlled by keeping cells with genes/cell >= 5 and reads/cell >= 20
# data alongside CCF and slide coordinates are saved and can be used for analysis

# load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate","scales",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo","ggrepel","igraph","ggraph")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))

# specify coverage
coverage = "cortex"
coverage = "whole-brain"

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

coordinates <- coord %>% filter(slice==slice_current); coord_ref <- coord_ref %>% filter(slice==slice_current)
data = data[,which(colnames(data) %in% coordinates$sample_id)]; data = t(data); data = as(data, "sparseMatrix")

```


```{r}

# single-cell RNA-seq data to assess cell segmentation errors in spatial data
# this is because the cells in scRNA-seq data are intact and no segmentation errors exist

# subset scRNA-seq by cortical regions present in the spatial data
spatial_regions = sort(unique(coordinates$CCFparentname))
single_cell_regions = c("ACA","AI-CLA","AId-AIv-AIp","AUD","CTXsp-CLA-EP-LA-BLA-BMA-PA","ENT","HIP-CA","PTLp","RSP","SSp","TEa-PERI-ECT")
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
df = df[, colnames(df) %in% colnames(data)]
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
  labs(subtitle = str_c("Spearman's rho = ", signif(cor(single, spatial, method = "spearman"), digits = 2)),
       x = "Co-expression in scRNA-seq data",
       y = "Co-expression in BARseq spatial data") +
  theme_classic() +
  theme(plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ggsave(filename = "/Users/AmeerSarwar/Downloads/cell_segment_barseq.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# we see low co-expression in spatial data than in scRNA-seq data
# this could be because of low counts in spatial data, so genes are less likely to be co-present

# plot average counts of each gene in spatial and scRNA-seq data
spatial = colmeans(as.matrix(data))
single  = colmeans(as.matrix(df))

# normalize the average counts
#spatial = rank(spatial)
#single  = rank(single)

spatial = spatial / max(spatial)
single  = single  / max(single)

avg_counts = data.frame(spatial, single)

ggplot(avg_counts) + aes(x = single, y = spatial) +
  geom_point(color = "steelblue", size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(subtitle = str_c("Spearman's Rho = ", signif(cor(single, spatial, method = "spearman"), digits = 2)),
       x = "Normalized average counts in scRNA-seq data",
       y = "Normalized average counts in BARseq spatial data") +
  theme_classic() +
  #scale_x_log10() + scale_y_log10() +
  theme(plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ggsave(filename = "/Users/AmeerSarwar/Downloads/sc_spatial_average_barseq.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# we see less co-expression and fewer average counts in spatial data
# this might be because the boundaries are too small
# compare spatial data with single-nucleus RNA-seq data (smaller boundaries)

# nucleus regions
nucleus = c("ACA","AUD","CTX","HPF","HY","LSX","MB","MOp","OLF","PALd","PALm","PALv","RSP","STRd","SUB","TH")

# load snRNA-seq data
single_cell_regions = str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/macosko/20221007_atlas_500umis_mt-1pct.h5ad.barseq_and_Vizgen_genes.", nucleus,".h5ad")
nucleus <- NULL

for (i in 1:length(single_cell_regions)){
  
  sc <- zellkonverter::readH5AD(single_cell_regions[i])
  sc <- t(assay(sc,"X"))
  nucleus <- rbind(nucleus, sc)
  print(i/length(single_cell_regions))
}

# subset and order both datasets by genes in spatial data
gene_name = read.csv("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/macosko/gene_names.csv"); gene_name = gene_name$gene
colnames(nucleus) = gene_name

nucleus = nucleus[, colnames(nucleus) %in% colnames(data)]
nucleus = nucleus[, -55] # removing an all 0's columns with label Ngr1, a gene otherwise present in the data

# load spatial data
coord = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
space = as.matrix(readMM(file="/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt"))
rownames(space) = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv", header = TRUE))$Genes; colnames(space) = coord$sample_id

# choose few slices from the middle
slices = 7:24
space  = space[, colnames(space) %in% coord[coord$slice %in% slices, "sample_id"]]
space  = t(space)

# re-order gene names
nucleus = nucleus[, colnames(space)]

# compute co-expression
spatial_data = correlation(space)
nucleus_data = correlation(nucleus)

spatial_data = upper_tri(spatial_data)
nucleus_data = upper_tri(nucleus_data)

# plot co-expression against each other
corr = data.frame(space = spatial_data, nucleus = nucleus_data)

ggplot(corr) + aes(x = nucleus, y = space) +
  geom_point(size = 0.1, color = "steelblue") + geom_abline(intercept = 0, slope = 1) +
  labs(subtitle = str_c("Spearman's Rho = ", signif(cor(nucleus_data, spatial_data, method = "spearman"), digits = 2)),
       x = "Co-expression in snRNA-seq data",
       y = "Co-expression in BARseq spatial data",) +
  theme_classic()

```


```{r}

# compare average counts of each gene across single-nucleus and spatial data
avg_space   = colmeans(space)
avg_nucleus = colmeans(nucleus)

avg_space   = rank(avg_space)
avg_nucleus = rank(avg_nucleus)

# plot average gene counts across datasets
avg = data.frame(x = avg_nucleus, y = avg_space)

ggplot(avg) + aes(x = x, y = y) +
  geom_point(color = "steelblue") + geom_abline(intercept = 0, slope = 1) +
  labs(subtitle = str_c("Spearman's Rho = ", signif(cor(avg_space, avg_nucleus, method = "spearman"), digits = 2)),
       x = "Rank of average gene counts in snRNA-seq",
       y = "Rank of average gene counts in BARseq spatial data") +
  theme_classic()

```


```{r}

# we want to assess whether cross-expression is drive by cell type differences
# find cross-expressing genes between glutamate cell type
# compare cell-neighbor matches at H2 or H3 subtype within these cells
# if cell-neighbor pairs are mostly of different types, then cell type differences explain cross expression

# load data
coord = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
data  = as.matrix(readMM(file="/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt"))
rownames(data) = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv", header = TRUE))$Genes; colnames(data) = coord$sample_id
data  = t(data)

# cross-expression between glutamate cell-neighbor and then lower cell type levels explored
type = "GABAergic"; type = "Glutamatergic"

# remove unclear and low quality cells at H2 and H3 levels
to_remove = c("Unclear","Low Quality","","\\?")
coord = coord[!coord$subclass_H2 %in% to_remove,]
data  = data[rownames(data) %in% coord$sample_id,]

out_cells = vector(mode = "numeric", length = nrow(coord))
for (i in 1:length(to_remove)){if (to_remove[i] == ""){next}; out_cells = as.numeric(str_detect(string = coord$cluster_H3, pattern = to_remove[i])) + out_cells}

coord = coord[!as.logical(out_cells),]
data  = data[rownames(data) %in% coord$sample_id,]

# cross-expression
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
cross = cross_expression(data = data, locations = coord[,c("CCFx","CCFy","CCFz")])
cross = cross[cross$cross_pvalue <= 0.05, ]

# fraction of cell-neighbor pairs that are different
fraction    = data.frame(cross[,c("gene1","gene2")], diff = rep(0, nrow(cross)), same = rep(0, nrow(cross)), pairs = rep(0, nrow(cross)))
other_fraction = fraction

# neighbor distances
distances = RANN::nn2(data = coord[,c("CCFx","CCFy","CCFz")], query = coord[,c("CCFx","CCFy","CCFz")], k = 2, searchtype = "priority")$nn.idx[,2]

for (i in 1:nrow(cross)){

    # add gene counts to metadata
    meta  = coord
    genes = data[, colnames(data) %in% c(cross$gene1[i], cross$gene2[i])]
    colnames(genes) = c("gene1", "gene2")
    meta  = data.frame(meta, genes)
    
    # metadata for neighbors
    meta_neig = meta[distances,]
    
    # cross-expressing cell-neighbor pairs
    cell1 = meta$gene1 >  0 & meta$gene2 == 0
    cell2 = meta$gene1 == 0 & meta$gene2 >  0
  
    neig1 = meta_neig$gene1 == 0 & meta_neig$gene2 >  0
    neig2 = meta_neig$gene1 >  0 & meta_neig$gene2 == 0
    
    # subset datasets by cross-expressing cell-neighbor pairs, A to B
    pairs = as.logical(cell1 * neig1)
    cell1 = meta[pairs,]
    neig1 = meta_neig[pairs,]
    
    # subset datasets by cross-expressing cell-neighbor pairs, B to A
    pairs = as.logical(cell2 * neig2)
    cell2 = meta[pairs,]
    neig2 = meta_neig[pairs,]
    
    # combine datasets
    cells = rbind(cell1, cell2)
    neigs = rbind(neig1, neig2)
    
    # subset at class H1 level (both cells and neighbors are glutamate)
    H1_only = cells$class_H1 %in% type & neigs$class_H1 %in% type
    cells   = cells[H1_only,]
    neigs   = neigs[H1_only,]
    
    # characterize cell types of cell-neighbor pairs
    fraction[i,"diff"]  = sum(cells$cluster_H3 != neigs$cluster_H3)
    fraction[i,"same"]  = sum(cells$cluster_H3 == neigs$cluster_H3)
    fraction[i,"pairs"] = nrow(cells)
    
    print(i/nrow(cross))
    next
    
    # null for cell-neighbor pairs not involved in cross-expression
    # simply remove cell-neighbor pairs involved in cross-expression
    to_remove = str_c(cells$sample_id, neigs$sample_id, sep = "___")
    all_pairs = str_c(meta$sample_id, meta_neig$sample_id, sep = "___")
    to_keep   = !all_pairs %in% to_remove
    
    other_cells = meta[to_keep,]
    other_neigs = meta_neig[to_keep,]
    
    # subset at class H1 level (both cells and neighbors are glutamate)
    H1_only     = other_cells$class_H1 %in% type & other_neigs$class_H1 %in% type
    other_cells = other_cells[H1_only,]
    other_neigs = other_neigs[H1_only,]
    
    # remove pairs where any gene is expressed
    to_remove = other_cells$gene1 > 0 | other_cells$gene2 > 0 | other_neigs$gene1 > 0 | other_neigs$gene2 > 0
    other_cells = other_cells[!to_remove,]
    other_neigs = other_neigs[!to_remove,]
    
    # randomly sample to match size of cross-expressing cell-neighbor pairs
    # else the global properties get reflected since most cell-neighbor pairs do not cross-express
    rand_ids = sample(1:nrow(other_cells), size = nrow(cells))
    
    other_cells = other_cells[rand_ids,]
    other_neigs = other_neigs[rand_ids,]
    
    # characterize cell types of cell-neighbor pairs
    other_fraction[i,"diff"]  = sum(other_cells$cluster_H3 != other_neigs$cluster_H3)
    other_fraction[i,"same"]  = sum(other_cells$cluster_H3 == other_neigs$cluster_H3)
    other_fraction[i,"pairs"] = nrow(other_cells)
    
    print(i/nrow(cross))
  }

# plot 1
ggplot(fraction) + aes(x = same, y = diff) + geom_point(size = 2, color = "steelblue") + geom_abline(intercept = 0, slope = 1) + theme_classic2() +
  labs(x = "Cell-neighbor pairs with the same cell type labels", y = "Cell-neighbor pairs with different cell type labels", title = str_c(type, " cells"),
      subtitle = str_c("Right-tailed Wilcoxon signed-rank test,\nDifferent ≥ Same, p-value = ", signif(wilcox.test(x = fraction$diff, y = fraction$same, alternative = "greater", paired = TRUE)$p.value, digits = 2))) +
  scale_x_log10() + scale_y_log10() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/diff_vs_same_labels1.svg", device = "svg", width = 6, height = 4, dpi = 600)

# plot 2
prop = data.frame(prop = fraction$diff / fraction$pairs)

ggplot(prop) + aes(x = prop, y = after_stat(scaled)) + geom_density(fill = "steelblue", color = "darkblue") + theme_bw() +
  labs(x = "Proportion of cell-neighbor pairs with different cell type labels", y = "Density",
       title = str_c(type, " cells"), subtitle = str_c("Median = ", signif(median(prop$prop), digits = 2))) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/diff_vs_same_labels2.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# this section can be skipped

# transcriptional bursting/ stochastic expression
# choose a cell type (H2 level) and find cross-expressing genes
# show that at a lower level (H3) this is cell type variation
# aka, genes are consistently cross-expressed in different cell types
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# choose cell type (use frequency at H2 and H3 types to make your choice)
cell_types = table(coordinates$subclass_H2)
cell_types = names(cell_types[cell_types > 200])

# storage dataframe
cell_neig = data.frame()

for (k in 1:length(cell_types)){
  
  # choose H2 cell type
  type = cell_types[k]
  coordinates_temp = coordinates[coordinates$subclass_H2 %in% type,]
  data_temp = data[rownames(data) %in% coordinates_temp$sample_id,]
  
  # ignore H2 cell types if H3 are not defined
  if (length(unique(coordinates_temp$cluster_H3)) <= 2){next}
  
  # cross-expressing genes
  cross = cross_expression(data = data_temp, locations = coordinates_temp[,c("pos_x","pos_y")])
  cross = cross[cross$cross_pvalue <= 0.05, ]
  
  # fraction of cell-neighbor pairs that are different
  fraction = data.frame(cross[,c("gene1","gene2")],
                        diff  = rep(0, nrow(cross)),
                        pairs = rep(0, nrow(cross)))

  for (i in 1:nrow(cross)){
    
    # add gene counts to metadata
    meta  = coordinates_temp
    genes = data_temp[, colnames(data_temp) %in% c(cross$gene1[i], cross$gene2[i])]
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
    
    # subset datasets by cross-expressing cell-neighbor pairs, A to B
    pairs = as.logical(cell1 * neig1)
    cell1 = meta[pairs,]
    neig1 = meta_neig[pairs,]
    
    # subset datasets by cross-expressing cell-neighbor pairs, B to A
    pairs = as.logical(cell2 * neig2)
    cell2 = meta[pairs,]
    neig2 = meta_neig[pairs,]
    
    # combine datasets
    cells = rbind(cell1, cell2)
    neigs = rbind(neig1, neig2)
    
    # characterize cell types of cell-neighbor pairs
    fraction[i,"diff"]  = sum(cells$cluster_H3 != neigs$cluster_H3)
    fraction[i,"pairs"] = nrow(cells)
    
  }
  
  # add H2 cell type and compute fraction diff. type in cell-neighbor pairs
  fraction = fraction[fraction$pairs > 10, ]
  fraction$prop = fraction$diff / fraction$pairs
  fraction = data.frame(fraction, type)
  
  # append to master dataframe
  cell_neig = rbind(cell_neig, fraction)
  total = length(unique(coordinates_temp$cluster_H3))
  print(str_c(type, ", subtypes = ", total))
}

# plot cell-neighbor pairs vs pairs w/ different cell type annotations
ggplot(cell_neig) + aes(x = pairs, y = diff, color = factor(type)) +
  geom_point(size = 0.8) +
  labs(title = "H2 to H3 cell types",
       x = "Cell-neighbor pairs involved in cross-expression",
       y = "Cell-neighbor pairs between different cell types",
       color = "H2 cell types") +
  geom_abline(slope = 1, intercept = 0) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw()


# generate null distribution by shuffling the diff/pairs values
# this is to test if the proportion difference is significant
null = matrix(0, nrow = 1000, ncol = length(unique(cell_neig$type)))
colnames(null) = unique(cell_neig$type)
rownames(null) = 1:1000

for (i in 1:1000){
  
  # shuffle diff/pairs within each cell type
  type = unique(cell_neig$type)
  
  for (j in 1:length(type)){
    
    temp = cell_neig[cell_neig$type %in% type[j], ]
    temp$diff = sample(temp$diff)
    
    # compute proportions of diff cell-neighbor annotations
    temp$prop1 = temp$diff / temp$pairs
    prop1 = temp$prop1
    prop1 = prop1[prop1 <= 1]
    
    null[i,j] = median(prop1)
  }
}

# plot density of proportion alongside p-value
ggplot(cell_neig) + aes(x = prop) +
  geom_density(color = "steelblue", fill = "lightblue", alpha = 0.5) +
  labs(subtitle = str_c("p-value = ", "a"),
       x = "Proportion of cell-neighbor pairs with different cell type annotation",
       y = "Density") +
  theme_bw()

ggplot(cell_neig) + aes(x = prop, y = factor(type)) +
  geom_boxplot() +
  scale_x_continuous(breaks = seq(0,1,0.1), labels = seq(0,1,0.1), limits = c(0,1))


nullx = as.data.frame(null)
nullx = pivot_longer(nullx, cols = everything(), names_to = "type", values_to = "prop") %>% as.data.frame()

ggplot(nullx) + aes(x = prop, y = factor(type)) +
  geom_boxplot() +
  scale_x_continuous(breaks = seq(0,1,0.1), labels = seq(0,1,0.1), limits = c(0,1))


# combined null distribution
null = c()

for (i in 1:1000){
  
  temp = cell_neig
  temp$diff = sample(temp$diff)
    
  # compute proportions of diff cell-neighbor annotations
  temp$prop1 = temp$diff / temp$pairs
  prop1 = temp$prop1
  prop1 = prop1[prop1 <= 1]
  null[i] = mean(prop1)
}


1 - (sum(mean(cell_neig$prop) >= null) / length(null))

nullx = data.frame(value = null, type = "null")
prop  = data.frame(value = cell_neig$prop, type = "empirical")
nullx = rbind(nullx, prop)

ggplot(nullx) + aes(x = value, after_stat(scaled), color = factor(type)) +
  geom_density()

```


```{r}

# co-expression within cells vs co-expression vs cross-expression with neighbors
corr_cell  = cor(data)
corr_neig  = cor(data, data[RANN::nn2(data = coord[,c("CCFx","CCFy","CCFz")], query = coord[,c("CCFx","CCFy","CCFz")], k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2],])
corr_cross = cross_expression_correlation(data = data, locations = coord[,c("CCFx","CCFy","CCFz")], output_matrix = TRUE)

df = data.frame(correlation_within_cells = upper_tri(corr_cell), correlation_between_neighbors = upper_tri(corr_neig), correlation_between_cross_expressing_neighbors = upper_tri(corr_cross))
df = as.data.frame(pivot_longer(data = df, cols = 2:3, names_to = "type", values_to = "vals"))

# plot
ggplot(df) + aes(x = correlation_within_cells, y = vals, color = type) + geom_point(size = 0) + facet_wrap(~type) +
  theme(legend.position = "none") + labs(y = "Correlation")



```