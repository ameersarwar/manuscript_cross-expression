```{r}

# gene panel includes cell type and neuromodulator genes on sagittal slices
# whole-brain BAR-seq data registered to the Allen Common Coordinate Framework version 3 (CCFv3)
# data is quality controlled by keeping cells with genes/cell >= 5 and reads/cell >= 20

# load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))

# load metadata and change column names
metadata = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/metadata.csv"))
colnames(metadata) = c("sample_id","slice","pos_x","pos_y","CCFx","CCFy","CCFz","CCFparentname","subclass_H2","cluster_H3")

# load and replace names of CCFparentname, subclass_H2, and cluster_H3
CCFparentname = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/area_names.csv")); CCFparentname_vec = CCFparentname[[1]]
subclass_H2   = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/H2_code.csv", sep = NULL)); subclass_H2_vec = subclass_H2[[1]]
cluster_H3    = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/H3_code.csv", sep = NULL)); cluster_H3_vec  = cluster_H3[[1]]

metadata <- metadata %>% mutate(CCFparentname = CCFparentname_vec[CCFparentname],
                                subclass_H2   = subclass_H2_vec[subclass_H2],
                                cluster_H3    = cluster_H3_vec[cluster_H3])

# load expression matrix, gene names, and replace columns with genes
data  = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/expression_matrix.csv"))
genes = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuromodulator_sagittal_barseq/gene_names.csv")); genes = genes[[1]]

colnames(data) = genes
rownames(data) = metadata$sample_id

# remove cells without genes/cell >= 5, reads/cell >= 20, and unused genes
data = data[rowSums(data) >= 20, ]
data = data[rowSums(data > 0) >= 5, ]
data = data[,which(!str_detect(colnames(data), "unused"))]

metadata = metadata[which(metadata$sample_id %in% rownames(data)),]

# remove cells not assigned to brain regions
metadata = metadata[which(!is.na(metadata$CCFparentname)),]
data     = data[which(rownames(data) %in% metadata$sample_id),]

# choose slice
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
#slice_current = 6:16 # other slices are too small
#metadata = metadata[metadata$slice %in% slice_current,]

metadata[,c("CCFx","CCFy")] = rotate_coordinates(x = metadata$CCFx, y = metadata$CCFy, flip_y = TRUE)
data = data[rownames(data) %in% metadata$sample_id,]; data = as(as.matrix(data), "sparseMatrix")

# remove slices w/ fewer than a min. no. of cells
if (TRUE){
  min_cells = 20000
  slice = table(metadata$slice)
  slice = as.numeric(names(slice[slice >= min_cells]))
  metadata = metadata[metadata$slice %in% slice,]
  data = data[rownames(data) %in% metadata$sample_id,]
}

```



```{r}

# example ligands and receptors tissue expression map + network
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
genes = c("Adcyap1","Adcyap1r1","Adra1a","Adrb1","Cck","Cckbr","Chrm1","Chrm2","Chrm3","Chrna7","Drd1","Htr1f","Htr2c","Htr3a","Npy","Npy1r","Oprd1","Oprk1","Oprm1","Pdyn","Penk","Sst","Sstr2","Vip","Vipr1")
slice = unique(metadata$slice)
pvalues = as.data.frame(matrix(data = 0, nrow = choose(ncol(data[, ]),2), ncol = length(slice))); colnames(pvalues) = slice

for (i in 1:length(slice)){
  temp_meta = metadata[metadata$slice %in% slice[i], ]
  temp_data = data[rownames(data) %in% temp_meta$sample_id, ]
  locations = temp_meta[,c("CCFx","CCFy","CCFz")]
  cross     = cross_expression(data = temp_data, locations = locations)
  nam_gene  = data.frame(gene1 = cross$gene1, gene2 = cross$gene2)
  pvalues[,i] = cross$cross_padj
  
  # plot tissue expression and spatial enrichment
  cross = cross[cross$cross_padj <= 0.05 & cross$cross_padj > 0,]; if (nrow(cross) == 0){next}
  print(i/length(slice))
  next
  
  for (j in 1:nrow(cross)){
    
    # tissue expression
    p = tissue_expression_plot(data = temp_data, locations = locations, gene1 = cross$gene1[j], gene2 = cross$gene2[j], cross_expression = FALSE) + theme(legend.position = "top"); print(p)
    ggsave(filename = str_c("/Users/AmeerSarwar/Downloads/LR/",i,"_",j,"cross_FALSE.png"), device = "png", dpi = 300, width = 5, height = 4)
    
    p = tissue_expression_plot(data = temp_data, locations = locations, gene1 = cross$gene1[j], gene2 = cross$gene2[j], cross_expression = TRUE) + theme(legend.position = "top"); print(p)
    ggsave(filename = str_c("/Users/AmeerSarwar/Downloads/LR/",i,"_",j,"cross_TRUE.png"), device = "png", dpi = 300, width = 5, height = 4)

    # spatial enrichment
    enrich = spatial_enrichment(data = temp_data, locations = locations, gene1 = cross$gene1[j], gene2 = cross$gene2[j])
    pvalue = enrich$pvalue
    enrich = data.frame(x = c(enrich$target, enrich$null), type = rep(c("Cross-expressing","Random"), times = c(length(enrich$target),length(enrich$null))))
    enrich$type = factor(enrich$type, levels = c("Random","Cross-expressing"))
    
    p = ggplot(enrich) + aes(x = x, fill = type, y = after_stat(scaled)) + geom_density(alpha = 0.8) +
      theme_classic2() + guides(fill = guide_legend(reverse = TRUE)) +
      labs(x = "Distance to cells", y = "Density", fill = "") +
      scale_fill_manual(values = c("Random" = "gray88","Cross-expressing" = "lightblue")) +
      theme(legend.position = "top", legend.text = element_text(size = 12), axis.title.x = element_text(size = 12))
    print(p)
    ggsave(filename = str_c("/Users/AmeerSarwar/Downloads/LR/",i,"_",j,"space_enrich.svg"), device = "svg", dpi = 600, width = 5, height = 4)
  }
}

pvals = pvalues
pvals[pvals <= 0.05] = -1
pvals[pvals != -1] = 0
pvals[pvals == -1] = 1
pvals = rowSums(pvals)
pvals = pvals[pvals != 0]
pvals = table(pvals)
pvals = data.frame(slices = as.numeric(names(pvals)), freq = as.numeric(pvals))

ggplot(pvals) + aes(x = slices, y = freq) +
  geom_bar(stat = "identity", fill = "lightblue", color = "darkblue") +
  theme_classic2() +
  labs(x = "Number of slices where genes are cross-expressed", y = "Number of gene pairs") +
  scale_x_continuous(breaks = pvals$slices, labels = pvals$slices) +
  coord_flip() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave("/Users/AmeerSarwar/Downloads/plot_freq.svg", device = "svg", dpi = 600, width = 6, height = 4)

#pvalues = data.frame(gene1 = nam_gene$gene1, gene2 = nam_gene$gene2, min_pval = rowMins(as.matrix(pvalues), value = TRUE))

```


```{r}

# single-cell RNA-seq data to assess cell segmentation errors in spatial data
# this is because the cells in scRNA-seq data are intact and no segmentation errors exist

# subset scRNA-seq by regions present in the spatial data
spatial_regions = sort(unique(metadata$CCFparentname))
single_cell_regions = c("ACA","CB-CBN","CB-VERM-CBN","CB-VERM","CNU-PAL","HB-MY","HIP-CA","HIP","HY-LZ","HY-MEZ-PVZ-PVR","MB-IC","MB-PRT","MB-SC-IC","MB-SC","MO-FRP","MOs-FRP",
                        "MY-MY-sat-MY-mot-Anterior","MY-MY-sat-MY-mot-Posterior","MY-MY-sat-MY-mot-Venteral-Anterior","MY-MY-sat-MY-mot-Venteral-Posterior","MY-MY-sat-MY-mot",
                        "OLF-AON-TT-DP-PIR-COA-NLOT-PAA-TR","OLF-AON-TT-DP","OLF-COA-PAA-NLOT-TR","OLF-MOB-AOB","OLF-PIR","PL-ILA-ORB","PONS-Pmot-Psat-Anterior","PONS-Pmot-Psat-Posterior",
                        "RSP","STR-STRv")
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
  labs(subtitle = str_c("Correlation = ", signif(cor(single,spatial, use = "complete.obs"), digits = 2)),
       x = "Co-expression in region matched scRNA-seq data from Zeng", y = "Co-expression in BAR-seq spatial data") +
  theme_bw()


```


```{r}

# we see low co-expression in spatial data than in scRNA-seq data
# this could be because of low counts in spatial data, so genes are less likely to be co-present

# plot average counts of each gene in spatial and scRNA-seq data
spatial = colmeans(as.matrix(data))
single  = colmeans(as.matrix(df))

avg_counts = data.frame(spatial, single)

ggplot(avg_counts) + aes(x = single, y = spatial) +
  geom_point(color = "steelblue", size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(subtitle = str_c("Correlation = ", signif(cor(single,spatial), digits = 2)),
       x = "Average counts in region matched scRNA-seq data from Zeng", y = "Average counts in MERSCOPE spatial data") +
  theme_bw()

```


```{r}

# analysis for ligand-receptor pairs
# NOTE: reload the data by commenting out the slice information
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# ligands and receptors
genes = c("Adcyap1","Adcyap1r1","Adra1a","Adrb1","Cck","Cckbr","Chrm1","Chrm2","Chrm3","Chrna7","Drd1","Htr1f","Htr2c","Htr3a","Npy","Npy1r","Oprd1","Oprk1","Oprm1","Pdyn","Penk","Sst","Sstr2","Vip","Vipr1")
dat = data[,colnames(data) %in% genes]
slice_current = unique(metadata$slice)
pvalues = data.frame()

for (i in 1:length(slice_current)){
  
  temp_meta = metadata[metadata$slice %in% slice_current[i],]
  temp_dat  = dat[rownames(dat) %in% temp_meta$sample_id,]
  locations = temp_meta[,c("CCFx","CCFy","CCFz")]

  temp_dat  = cross_expression(data = temp_dat, locations = locations)
  temp_dat  = temp_dat[temp_dat$cross_sig == 1,]
  temp_dat  = temp_dat[temp_dat$cross_padj != 0,]
  temp_dat  = data.frame(temp_dat, slice = rep(slice_current[i], nrow(temp_dat)))

  pvalues   = rbind(pvalues, temp_dat) 
  print(i/length(slice_current))
}

# subset by frequent gene pairs
freq = str_c(pvalues$gene1, pvalues$gene2, sep = "_")
freq = table(freq)
freq = freq[freq > 2]
freq = names(freq)
targ = str_c(pvalues$gene1, pvalues$gene2, sep = "_")

pvalues = pvalues[targ %in% freq,]
pvalues = pvalues[order(pvalues$gene1),]

# plot each tissue cross-expression
for (i in 1:nrow(pvalues)){
  
  if (i %in% c(4,5,6)){next}
  
  temp_meta = metadata[metadata$slice %in% pvalues[i,]$slice,]
  temp_dat  = dat[rownames(dat) %in% temp_meta$sample_id,]
  locations = temp_meta[,c("CCFx","CCFy")]
  
  p1 = tissue_expression_plot(data = temp_dat, locations = locations, gene1 = pvalues$gene1[i], gene2 = pvalues$gene2[i], cross_expression = FALSE, point_size = 0, scale_bar = 0); print(p1)
  p2 = tissue_expression_plot(data = temp_dat, locations = locations, gene1 = pvalues$gene1[i], gene2 = pvalues$gene2[i], cross_expression = TRUE , point_size = 0, scale_bar = 0); print(p2)
  
  # save plots
  if (!TRUE){
    ggsave(plot = p1, filename = str_c("/Users/AmeerSarwar/Downloads/FALSE", i,".png"), width = 7, height = 4, device = "png", dpi = 300)
    ggsave(plot = p2, filename = str_c("/Users/AmeerSarwar/Downloads/TRUE" , i,".png"), width = 7, height = 4, device = "png", dpi = 300)
  }
  print(i/nrow(pvalues))
}

```


```{r}

# genes can cross-express until (say) 25 nearest neighbors
# compute cross-expression until this number per slice
# take the proportion with significant cross-expression
# this means that the gene expression programs are robust
# take the max and plot proportions for three pair types
# both ligand/receptor (2), one ligand/receptor (1), no ligand/receptor (0)

slice_current = unique(metadata$slice)
neigs = 1:10
cross = data.frame(matrix(data = 0, nrow = choose(ncol(data),2), ncol = length(slice_current)))
colnames(cross) = slice_current

for (i in 1:length(slice_current)){
  
  # extract slice
  temp_meta  = metadata[metadata$slice %in% slice_current[i],]
  temp_data  = data[rownames(data) %in% temp_meta$sample_id,]
  temp_meta  = temp_meta[,c("CCFx","CCFy","CCFz")]
  
  temp_cross = data.frame(matrix(data = 0, nrow = choose(ncol(data),2), ncol = length(neigs)))
  colnames(temp_cross) = neigs
  
  # significant cross-expression until specified nearest neighbor
  for (j in 1:length(neigs)){
    xx = cross_expression(data = temp_data, locations = temp_meta, neighbor = neigs[j])
    temp_cross[,j] = xx$cross_padj
  }
  
  # proportion of nearest neighbors with significant cross-expression
  temp_cross = t(apply(X = temp_cross, MARGIN = 1, FUN = function(x) p.adjust(x, method = "BH")))
  temp_cross[temp_cross <= 0.05] = -1
  temp_cross[temp_cross != -1] = 0
  temp_cross[temp_cross == -1] = 1
  temp_cross = rowSums(temp_cross) / max(neigs)
  
  # store proportion for the slice
  cross[,i] = temp_cross
  print(i/length(slice_current))
}

# genes active at least once
pp = data.frame(gene1 = xx$gene1, gene2 = xx$gene2, exp = as.numeric(rowSums(cross) >= 1))
genes = c("Adcyap1","Adcyap1r1","Adra1a","Adrb1","Cck","Cckbr","Chrm1","Chrm2","Chrm3","Chrna7","Drd1","Htr1f","Htr2c","Htr3a","Npy","Npy1r","Oprd1","Oprk1","Oprm1","Pdyn","Penk","Sst","Sstr2","Vip","Vipr1")

LR = sum(pp[  (pp$gene1 %in% genes & pp$gene2 %in% genes)  ,]$exp) / choose(length(genes),2)
LR_not = sum(pp[  !(pp$gene1 %in% genes | pp$gene2 %in% genes)  ,]$exp) / choose(ncol(data) - length(genes),2)

# network edges
genes = c("Adcyap1","Adcyap1r1","Adra1a","Adrb1","Cck","Cckbr","Chrm1","Chrm2","Chrm3","Chrna7","Drd1","Htr1f","Htr2c","Htr3a","Npy","Npy1r","Oprd1","Oprk1","Oprm1","Pdyn","Penk","Sst","Sstr2","Vip","Vipr1")
express = cross
express = apply(X = express, MARGIN = 1, FUN = function(x) median(x[x != 0])) # median computed w/o 0's
express[is.na(express)] = 0

nam = matrix(data = 0, nrow = ncol(data), ncol = ncol(data))
nam = which(upper.tri(nam), arr.ind = TRUE)
nam = data.frame(gene1 = colnames(data)[nam[,1]], gene2 = colnames(data)[nam[,2]])

express = data.frame(nam, prop = express)

type = vector(mode = "numeric", length = nrow(express))
type[express$gene1 %in% genes & express$gene2 %in% genes] = 2
type[xor(express$gene1 %in% genes, express$gene2 %in% genes)] = 1

express = data.frame(express, type)
express = data.frame(express, edge = as.numeric(express$prop > 0))

# node characteristics
node = table(c(express[express$edge == 1,]$gene1, express[express$edge == 1,]$gene2))
node = data.frame(gene = names(node), degree = as.numeric(node), type = as.numeric(names(node) %in% genes))

ggplot(node) + aes(x = reorder(gene, degree), y = degree, fill = factor(type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 1) +
  theme_classic() + guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "darkred"), labels = c("0" = "No", "1" = "Yes")) +
  labs(x = "Cross-expressing genes", y = "Node degree", fill = "Ligand/Receptor",
       subtitle = str_c("P-value = ", signif(wilcox.test(node[node$type == 1,]$degree, node[node$type == 0,]$degree, alternative = "greater")$p.value, digits = 2))) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.9),
        legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12))

edge = express[express$edge == 1,]
trend_node = data.frame()

for (i in 0:(length(unique(edge$type)) - 1)){
  
  trend = table(c(edge[edge$type == i,]$gene1, edge[edge$type == i,]$gene2))
  trend = data.frame(genes = names(trend), degree = as.numeric(trend), type = rep(i,length(trend)))
  avg_prop = vector(mode = "numeric", length = nrow(trend))
  
  for (j in 1:nrow(trend)){
    tt = edge[edge$type == i,]
    avg_prop[j] = median(tt[tt$gene1 %in% trend$genes[j] | tt$gene2 %in% trend$genes[j],]$prop)
  }
  trend = data.frame(trend, avg_prop)
  trend_node = rbind(trend_node, trend)
}

for (i in 0:(length(unique(edge$type)) - 1)){
  
  trend = trend_node[trend_node$type == i,c("degree","avg_prop")]
  p = ggplot(trend) + aes(x = degree, y = avg_prop) +
    geom_point(size = 2) + geom_smooth(method = "lm") +
    theme_classic2() +
    labs(x = "Node degree", y = "Proportion of cross-expressing neighbors (n = 25)",
       subtitle = str_c("Spearman's rho = ", signif(cor(trend$degree, trend$avg_prop, method = "spearman"), digits = 2))) +
    theme(plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
  print(p)
}

# genes with higher node degree show stronger connections
# ligands/receptors with higher node degree show weaker conncetions
# so, trade-off between strong and specific vs weak and widespread connections

# global statistics
edge  = express[express$edge == 1,]
fact1 = as.numeric((table(edge$type)[2] + table(edge$type)[3]) / sum(table(edge$type))) # 53% of cross-expression links contain at least one ligand/receptor
fact2 = length(genes) / ncol(data) # only 19% of the panel contains ligand/receptor genes

# save network
write.csv(x = node, file = "/Users/AmeerSarwar/Downloads/node_sagittal_barseq.csv", row.names = FALSE, quote = FALSE)
write.csv(x = express[express$edge == 1,], file = "/Users/AmeerSarwar/Downloads/edge_sagittal_barseq.csv", row.names = FALSE, quote = FALSE)

# consider ligand/receptor only network
edge = express[express$type == 2,]
edge = edge[edge$edge == 1,]
node = table(c(edge$gene1, edge$gene2))
node = data.frame(gene = names(node), degree = as.numeric(node))

# save network
write.csv(x = node, file = "/Users/AmeerSarwar/Downloads/LR_node_sagittal_barseq.csv", row.names = FALSE, quote = FALSE)
write.csv(x = edge, file = "/Users/AmeerSarwar/Downloads/LR_edge_sagittal_barseq.csv", row.names = FALSE, quote = FALSE)

```


```{r}

# similarity between slice-specific networks as a function of distance between slices
slice_current = unique(metadata$slice)
network       = data.frame(matrix(data = 0, nrow = choose(ncol(data),2), ncol = length(slice_current)))
colnames(network) = slice_current

for (i in 1:length(slice_current)){
  
  temp_meta = metadata[metadata$slice %in% slice_current[i],]
  temp_dat  = data[rownames(data) %in% temp_meta$sample_id,]
  
  xx = cross_expression_correlation(data = temp_dat, locations = temp_meta[,c("CCFx","CCFy","CCFz")])
  network[,i] = xx$correlation
  print(i/length(slice_current))
}

corr = cor(network, method = "spearman")
corr = data.frame(which(upper.tri(corr), arr.ind = TRUE), corr = upper_tri(corr))
corr = data.frame(dist = corr[,2] - corr[,1], corr = corr$corr)

ggplot(corr) + aes(x = dist, y = corr) +
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Distance (in slices, 300µm between slices)", y = "Correlation",
       subtitle = str_c("Spearman's rho = ", signif(cor(corr$dist, corr$corr, method = "spearman"), digits = 2))) +
  scale_x_continuous(breaks = min(corr$dist):max(corr$dist), labels = min(corr$dist):max(corr$dist)) +
  scale_y_continuous(breaks = seq(from = signif(min(corr$corr), digits = 1), to = signif(max(corr$corr), digits = 1), by = 0.1),
                     labels = seq(from = signif(min(corr$corr), digits = 1), to = signif(max(corr$corr), digits = 1), by = 0.1)) +
  theme_classic() + theme(axis.text = element_text(size = 14),
                          axis.title = element_text(size = 16),
                          plot.subtitle = element_text(size = 16),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ggsave(filename = "/Users/AmeerSarwar/Downloads/corr_plot.svg", device = "svg", dpi = 600, width = 8, height = 5)

```
