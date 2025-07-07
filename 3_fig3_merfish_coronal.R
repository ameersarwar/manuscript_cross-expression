```{r}

# MERFISH brain receptor map
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo","igraph","viridis")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))
registerDoParallel(detectCores() - 1)
plan("multisession", workers = 8)

# choose slice and replicate (mouse)
replicate = 2; slice = 2; cluster = FALSE

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

# rotate slice (counter-clockwise) and center to origin
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
degs <- 40
rot  <- rotate_coordinates(x = metadata$pos_x, y = metadata$pos_y, n_degrees = degs, center = TRUE)

metadata$pos_x = rot$pos_x
metadata$pos_y = rot$pos_y

```


```{r}

xx = data

# testing no. of cells chunk
df = data[,c("Gabbr1","Ltb4r2")]
colnames(df) = c("gene1","gene2")
locations = metadata[,c("pos_x","pos_y")]

cross_expression(data = df, locations = locations)

# manual
df_neig = df[RANN::nn2(data = locations, query = locations, k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2],]

k = sum(df$gene1 > 0 & df$gene2 == 0)                                           # sample size: number of cells with gene A (but not B); aka "number of balls drawn"
q = sum(df$gene1 > 0 & df$gene2 == 0 & df_neig$gene1 == 0 & df_neig$gene2 > 0)  # observed successes: number of neighbors with gene B (but not A); aka "white balls seen"
m = sum(df_neig$gene1 == 0 & df_neig$gene2 > 0)                                 # success states: number of possible neighbors with gene B (but not A); aka "white balls in population"
n = nrow(df) - m                                                                # failure states: total population - success states; aka "black/other/irrelevant balls in population"

phyper(q = q - 1, m = m, n = n, k = k, lower.tail = FALSE)

```


```{r}

# how many neighboring cells are neighbored by source cells?
# note that a cell cannot have more than one nearest neighbor (obviously)
locations = metadata[,c("pos_x","pos_y")]
loc_cells = 1:nrow(data)
loc_neigs = loc_cells[RANN::nn2(data = locations, query = locations, k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2]]

freq = as.numeric(table(loc_neigs))
freq = as.numeric(table(freq)) / sum(as.numeric(table(freq)))
freq = data.frame(x = 1:length(freq), y = freq)

ggplot(freq) + aes(x = x, y = y) + geom_point(size = 5, color = "steelblue") +
  geom_segment(aes(x = x, xend = x, y = 0, yend = y), linewidth = 2, color = "steelblue") +
  labs(x = "Neighbor targeted as nearest by a given number of cells", y = "Proportion") +
  theme_bw() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
ggsave("/Users/AmeerSarwar/Downloads/prop.svg", device = svg, dpi = 600, width = 6, height = 4)

# see if results are consistent when one-cell one-neighbor, many-cells one-neighbor, and full dataset?

# one-cell one-neighbor
freq = table(loc_neigs)
freq = as.numeric(names(freq[freq == 1]))
filt = loc_cells %in% freq
filt = locations[filt,]
filt1 = cross_expression_correlation(data = data[loc_cells %in% freq, 2:ncol(data)], locations = filt)$correlation

# many-cell one-neighbor
freq = table(loc_neigs)
freq = as.numeric(names(freq[freq != 1]))
filt = loc_cells %in% freq
filt = locations[filt,]
filt2 = cross_expression_correlation(data = data[loc_cells %in% freq, 2:ncol(data)], locations = filt)$correlation

# full dataset
unfilt = cross_expression_correlation(data = data[,2:ncol(data)], locations = locations)$correlation

simm = data.frame(x = filt1, y = filt2, full = unfilt)

ggplot(simm) + aes(x = x, y = y, color = unfilt) +
  geom_point(size = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "One-cell one-neighbor",
       y = "Many-cells one-neighbor",
       color = "Full dataset",
       title = "Cross-expression correlations",
       subtitle = str_c("Average Pearson's R between replicates = ", signif(mean(upper_tri(cor(simm))), digits = 2))) +
  scale_color_viridis_c() +
  theme_classic2() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.75),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("/Users/AmeerSarwar/Downloads/simm.png", device = png, dpi = 600, width = 6, height = 4)

```


```{r}

# network lacking co-expression and w/o sparse connections (p-value = 0)
# idea is to only find interesting cross-expression signal w/o co-expression
cross = cross_expression(data = data[,2:ncol(data)], locations = metadata[,c("pos_x","pos_y")], alpha_co = 0.05)
cross = cross[as.logical(cross$cross_sig),]
cross = cross[order(cross$gene1),]

# prune the network by removing Lgr6 and genes not in the main network
#to_remove = c("Lgr6","Ptgdr2","Ptafr","Abcc9","Mrgprb4","Vmn1r50","S1pr4","Ntrk2","Tas2r7","Taar7e")
#cross = cross[!(cross$gene1 %in% to_remove | cross$gene2 %in% to_remove),]

if (FALSE){
  for (i in 1:nrow(cross)){
    p = tissue_expression_plot(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = cross$gene1[i], gene2 = cross$gene2[i], point_size = 0.1, scale_bar = 0, cross_expression = TRUE ); print(p)
    print(i/nrow(cross))
  }
}

# nodes and edges
node = table(c(cross$gene1, cross$gene2))
node = data.frame(gene = names(node), degree = as.numeric(node))

# remove genes with 1 node degree
to_remove = c(node[node$degree == 1,]$gene)
cross = cross[!(cross$gene1 %in% to_remove | cross$gene2 %in% to_remove),]
node = table(c(cross$gene1, cross$gene2))
node = data.frame(gene = names(node), degree = as.numeric(node))

# co-expression patterns
co   = cor(data[,node$gene])
ht = Heatmap(co, name = " ", heatmap_legend_param = list(title = "Correlation", title_gp = gpar(fontsize = 12), title_position = "leftcenter-rot"))
draw(ht)

svg("/Users/AmeerSarwar/Downloads/heatmap.svg", width = 8, height = 6)
draw(ht)
dev.off()

# save files
# manually assign community based on co-expression dendrogram heatmap
node$comm = c(1,2,1,4,4,2,3,3,3,4,1,2,4,2,3,3,3,3)

write.csv(node,  file = "/Users/AmeerSarwar/Downloads/merfish_small_net_node.csv", row.names = FALSE, quote = FALSE)
write.csv(cross, file = "/Users/AmeerSarwar/Downloads/merfish_small_net_edge.csv", row.names = FALSE, quote = FALSE)

# broad cell type markers correlation with our network genes
astrocytes = c("Aldh1l1","Aqp4","Gfap","Sox9")
microglia  = c("Cx3cr1","Adgre1","Csf1r","P2ry12")
oligo      = c("Olig1","Pdgfra")
neuron     = c("Chat","Th","Gabbr1","Gabbr2","Gad1","Slc32a1","Grin2b")
endo       = "Igf1r"
cells      = c(astrocytes, microglia, oligo, neuron, endo)

corr = cor(data[,node$gene], data[,cells])
Heatmap(corr)

# create heatmap

# row annotations
group_labels = rep("Other", nrow(corr))
names(group_labels) = rownames(corr)

group_labels[c("C3ar1","Cd300c2","Gpr183","P2yr13")] = "Microglia"
group_labels[c("Adora2b","Egfr","Lgr6","Ppp1r3g")] = "Astrocytes"
group_labels[c("Adgrg5","Agtr1b","Gpr4")] = "Neurons"
group_labels[c("Emcn","Epha2","Flt4","Pth1r","Slco1a4","Tek","Tie1")] = "Endothelial"
group_colors = c("Microglia" = "#FFD92F", "Astrocytes" = "#FC8D62", "Neurons" = "#66C2A5", "Endothelial" = "#E78AC3")

row_ha = rowAnnotation(Group = anno_simple(group_labels, col = group_colors), annotation_label = "")
annotation_legend = Legend(labels = names(group_colors), legend_gp = gpar(fill = group_colors), title = "")

ht = Heatmap(corr, name = " ", right_annotation = row_ha, heatmap_legend_param = list(title = "Correlation", title_gp = gpar(fontsize = 12), title_position = "leftcenter-rot"))
draw(ht, annotation_legend_list = list(annotation_legend))

# save heatmap
svg("/Users/AmeerSarwar/Downloads/heatmap.svg", width = 8, height = 6)
draw(ht, annotation_legend_list = list(annotation_legend))
dev.off()

# GO enrichment
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/GO_enrichment.R")
test_genes = node$gene
go = GO_enrichment(test_genes = test_genes, background_genes = colnames(data)[2:ncol(data)])

pvals = go$p_values
pvals = pvals[pvals$pvals_orig != 1,] # most GO groups are useless but we don't know which ones
pvals$pvals_FDR = p.adjust(pvals$pvals_orig, method = "BH")
pvals$sig_GO = as.numeric(pvals$pvals_FDR <= 0.05)
pvals$pvals_orig_log10 = -log10(pvals$pvals_orig)
pvals$pvals_FDR_log10  = -log10(pvals$pvals_FDR)

pvals = pvals[as.logical(pvals$sig_GO),]
pvals = pvals[order(pvals$pvals_FDR_log10),]

pvals[4,]$group = "positive regulation of vascular\nendothelial growth factor production"

ggplot(pvals) + aes(x = pvals_FDR_log10, y = fct_reorder(group, pvals_FDR_log10)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
  theme_classic() +
  labs(x = "-log10 p-value", y = "") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="darkred") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("/Users/AmeerSarwar/Downloads/GO.svg", device = svg, dpi = 600, width = 6, height = 4)

```


```{r}

# network analysis using cross-expression
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# p-values
pvals = cross_expression(data = data[,2:ncol(data)], locations = metadata[,c("pos_x","pos_y")], output_matrix = TRUE)
pvals = pvals$Cross_with_FDR
pvals = ifelse(pvals <= 0.05, 1, 0)

# edge list
df = data.frame(which(upper.tri(pvals), arr.ind = TRUE))
df = data.frame(df,
                gene1 = colnames(pvals)[df$row],
                gene2 = colnames(pvals)[df$col],
                edge  = upper_tri(pvals))

# two-step network
n_network = pvals %*% pvals
diag(n_network) = 0
df = data.frame(df, two_step = upper_tri(n_network)); second_order = df

two_step_threshold = 4
df$two_step = as.integer(df$two_step >= two_step_threshold)

# node degree
node = df[df$edge == 1,]
node = table(c(node$gene1, node$gene2))
node = data.frame(gene = names(node), degree = as.vector(node))

# 2nd order community
coms = df[df$two_step == 1, c("gene1","gene2")]
coms = graph_from_data_frame(coms, directed = FALSE)
coms = membership(cluster_louvain(coms))

# degree + community
node$coms = vector(mode = "numeric", length = nrow(node))
node[node$gene %in% names(coms),"coms"] = as.vector(coms)

# 1st and 2nd order edge combinations
type = vector(mode = "numeric", length = nrow(df))
type[df$edge == 1 & df$two_step == 0] = 1
type[df$edge == 0 & df$two_step == 1] = 2
type[df$edge == 1 & df$two_step == 1] = 3

df = data.frame(df, edge_type = type)

# save results
df = df[df$edge == 1 | df$two_step == 1,!colnames(df) %in% c("row","col")]
write.csv(df,   file = "/Users/AmeerSarwar/Downloads/edges.csv", row.names = FALSE, quote = FALSE)
write.csv(node, file = "/Users/AmeerSarwar/Downloads/nodes.csv", row.names = FALSE, quote = FALSE)

# node degree distribution
nn = node[order(node$degree),]

ggplot(nn) + aes(x = degree) + geom_histogram(binwidth = 1, fill = "lightblue", color = "black", width = 0.1) +
  theme_classic2() + labs(x = "Node degree", y = "Frequency") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  annotate("segment", x = 38, xend = 39.5, y = 10, yend = 3, arrow = arrow(length = unit(0.3, "cm")), color = "black", linewidth = 0.8) +
  annotate("text", x = 37.5, y = 10, label = "Gpr20", color = "black", vjust = -1)
ggsave("/Users/AmeerSarwar/Downloads/node_degree.svg", device = svg, dpi = 600, width = 6, height = 4)

```


```{r}

# co-expression of genes cross-expressed with Gpr20
target = "Gpr20"
genes = df[df$gene1 %in% target | df$gene2 %in% target,]

go = genes[genes$edge_type %in% c(1,3),]
go = unique(c(go$gene1, go$gene2))
go = go[!go %in% target]

co = cor(data[,go])
Heatmap(co, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))

# co-expression of genes cross-expressed with Gpr20 against marker genes
astrocytes = c("Aldh1l1","Aqp4","Gfap","Sox9")
microglia  = c("Cx3cr1","Adgre1","Csf1r","P2ry12")
oligo      = c("Olig1","Pdgfra")
neuron     = c("Chat","Th","Gabbr1","Gabbr2","Gad1","Slc32a1","Grin2b")
endo       = "Igf1r"
cells      = c(astrocytes, microglia, oligo, neuron, endo)

co = cor(data[,go], data[,cells])
Heatmap(co, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8))

# co-expression of genes cross-expressed with Gpr20 against marker genes using cells involved in Gpr20 cross-expression
target_cells = data
target_neigs = data[RANN::nn2(data = metadata[,c("pos_x","pos_y")], query = metadata[,c("pos_x","pos_y")], k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2],]
coexpr = matrix(data = 0, nrow = length(go), ncol = length(cells), dimnames = list(go, cells))

for (i in 1:length(go)){
  
  # cell-neighbor relations
  locus   = target_cells[,c(target, go[i])] # chosen cell
  project = target_neigs[,c(target, go[i])] # its neighbor
  
  # target in cells and other genes in neighbors
  hits1 = (locus[,target] > 0 & locus[,go[i]] == 0) & (project[,target] == 0 & project[,go[i]] > 0)
  hits1 = target_neigs[hits1,]
  
  # target in neighbors and other genes in cells
  hits2 = (locus[,target] == 0 & locus[,go[i]] > 0) & (project[,target] > 0 & project[,go[i]] == 0)
  hits2 = target_cells[hits2,]
  
  # combine both directions and do not remove duplicates (because of multi-mapping of cell-neighbor relations)
  hits  = rbind(hits1, hits2)
  
  # compute co-expression of cross-expressing gene with cell type markers
  coexpr[i,] = suppressWarnings(cor(hits[,go[i]], hits[,cells]))
}

# plot heatmap of co-expression

# pre-process
corr = coexpr
corr[is.na(corr)] = 0

# annotations
group_labels = rep("Other", ncol(corr))
names(group_labels) = colnames(corr)

group_labels[c("Aldh1l1","Aqp4","Gfap","Sox9")] = "Astrocytes"
group_labels[c("Cx3cr1","Adgre1","Csf1r","P2ry12")] = "Microglia"
group_labels[c("Olig1","Pdgfra")] = "Oligodendrocytes"
group_labels[c("Chat","Th","Gabbr1","Gabbr2","Gad1","Slc32a1","Grin2b")] = "Neurons"
group_labels[c("Igf1r")] = "Endothelial"
group_colors = c("Astrocytes" = "#FC8D62", "Microglia" = "#FFD92F", "Oligodendrocytes" = "lightblue", "Neurons" = "#66C2A5", "Endothelial" = "#E78AC3")

col_ha = columnAnnotation(Group = anno_simple(group_labels, col = group_colors), annotation_label = "")
annotation_legend = Legend(labels = names(group_colors), legend_gp = gpar(fill = group_colors), title = "")

# draw heatmap
ht = Heatmap(corr, name = " ", cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8), bottom_annotation = col_ha,
        heatmap_legend_param = list(title = "Correlation", title_gp = gpar(fontsize = 12), title_position = "leftcenter-rot"))
draw(ht, annotation_legend_list = list(annotation_legend))

# save heatmap
svg("/Users/AmeerSarwar/Downloads/heatmap.svg", width = 8, height = 6)
draw(ht, annotation_legend_list = list(annotation_legend))
dev.off()

# plot distributions for each marker gene
type_astro = as.vector(corr[,colnames(corr) %in% astrocytes])
type_micro = as.vector(corr[,colnames(corr) %in% microglia])
type_oligo = as.vector(corr[,colnames(corr) %in% oligo])
type_neuro = as.vector(corr[,colnames(corr) %in% neuron])
type_endo  = as.vector(corr[,colnames(corr) %in% endo])

types = data.frame(vals = c(type_astro, type_micro, type_oligo, type_neuro, type_endo),
                   type = rep(c("Astrocytes", "Microglia", "Oligodendrocytes", "Neurons", "Endothelial"),
                              times = c(length(type_astro), length(type_micro), length(type_oligo), length(type_neuro), length(type_endo))))
types = types[types$vals != 0,]
types$type = factor(types$type, levels = c("Endothelial", "Neurons", "Oligodendrocytes", "Microglia", "Astrocytes"))

# make plot
ggplot(types) + aes(y = type, x = vals, fill = type) + geom_violin(trim = FALSE, scale = "width") +
  scale_fill_manual(values = c("Astrocytes" = "#FC8D62", "Microglia" = "#FFD92F", "Oligodendrocytes" = "lightblue", "Neurons" = "#66C2A5", "Endothelial" = "#E78AC3")) +
  labs(y = "", x = "Co-expression with genes cross-expressed with Gpr20", subtitle = "Co-expression in cells involved in\nGpr20-specific cross-expression") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12))
ggsave("/Users/AmeerSarwar/Downloads/type_coexpr.svg", device = svg, dpi = 600, width = 6, height = 4)

# pairwise p-values
pvalues = list(type_astro[type_astro != 0], type_micro[type_micro != 0], type_oligo[type_oligo != 0], type_neuro[type_neuro != 0], type_endo[type_endo != 0])
names(pvalues) = c("Astrocytes", "Microglia", "Oligodendrocytes", "Neurons", "Endothelial")
pvals   = matrix(data = 0, nrow = length(pvalues), ncol = length(pvalues), dimnames = list(names(pvalues), names(pvalues)))

for (i in 1:nrow(pvals)){
  for (j in 1:ncol(pvals)){
    # row types greater than column types
    pvals[i,j] = wilcox.test(pvalues[[i]], pvalues[[j]], exact = FALSE, alternative = "greater")$p.value
  }
}

pvals[upper.tri(pvals)] = 1; diag(pvals) = 1
pvals[lower.tri(pvals)] = p.adjust(lower_tri(pvals), method = "BH")

```


```{r}

# highlight cells expressing a certain gene, e.g., Gpr20
target = "Gpr20"

# slice 1, replicate 1 as example
replicate = 1
slice     = 1

# load data
expr <- fread(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice",slice,"_Replicate",replicate,"_cell_by_gene.csv")); expr  = as.data.frame(expr)
post <- fread(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice",slice,"_Replicate",replicate,"_cell_metadata.csv")); post = as.data.frame(post)
  
# pre-process data
post = post[,4:5]; colnames(post) = c("pos_x","pos_y")
post = data.frame(post, light = expr[,target] > 0)
post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 175)

# further pre-process
post$pos_x = scale_01(post$pos_x)
post$pos_y = scale_01(post$pos_y)
  
post$light = factor(post$light, levels = c("FALSE","TRUE"))
post = post[order(post$light),]

# make plot
ggplot(post) + aes(x = pos_x, y = pos_y, color = light) + geom_point(size = 0) +
  scale_color_manual(values = c("FALSE" = "gray88", "TRUE" = "brown3"), labels = c("FALSE" = "No", "TRUE" = "Yes")) + theme_classic() +
  labs(color = "Gpr20", x = "x coordinates", y = "y coordinates") +
  guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 2))) +
  
  # draw rectangles
  geom_rect(aes(xmin = 0.68, xmax = 0.77, ymin = 0.24, ymax = 0.43), alpha = 0, color = "#444444", linewidth = 0.1) +
  geom_rect(aes(xmin = 0.60, xmax = 0.70, ymin = 0.10, ymax = 0.20), alpha = 0, color = "#444444", linewidth = 0.1) +
  geom_rect(aes(xmin = 0.58, xmax = 0.65, ymin = 0.35, ymax = 0.45), alpha = 0, color = "#444444", linewidth = 0.1) +
  geom_rect(aes(xmin = 0.25, xmax = 0.32, ymin = 0.35, ymax = 0.58), alpha = 0, color = "#444444", linewidth = 0.1) +
  geom_rect(aes(xmin = 0.68, xmax = 0.77, ymin = 0.45, ymax = 0.60), alpha = 0, color = "#444444", linewidth = 0.1) +
  geom_rect(aes(xmin = 0.29, xmax = 0.37, ymin = 0.59, ymax = 0.70), alpha = 0, color = "#444444", linewidth = 0.1) +
  geom_rect(aes(xmin = 0.57, xmax = 0.63, ymin = 0.66, ymax = 0.73), alpha = 0, color = "#444444", linewidth = 0.1) +
  geom_rect(aes(xmin = 0.56, xmax = 0.65, ymin = 0.76, ymax = 0.87), alpha = 0, color = "#444444", linewidth = 0.1)

ggsave(filename = "/Users/AmeerSarwar/Downloads/blood_full.png", device = png, dpi = 300, width = 6, height = 4)

  # draw rectangles (slice 2, replicate 2) -- keep this just in case
  #geom_rect(aes(xmin = 0.86, xmax = 0.96, ymin = 0.60, ymax = 0.76), alpha = 0, color = "#444444", linewidth = 0.1) +
  #geom_rect(aes(xmin = 0.28, xmax = 0.32, ymin = 0.13, ymax = 0.20), alpha = 0, color = "#444444", linewidth = 0.1) +
  #geom_rect(aes(xmin = 0.00, xmax = 0.06, ymin = 0.64, ymax = 0.68), alpha = 0, color = "#444444", linewidth = 0.1) +
  #geom_rect(aes(xmin = 0.30, xmax = 0.33, ymin = 0.47, ymax = 0.55), alpha = 0, color = "#444444", linewidth = 0.1)

# areas within rectangles
rectangles_df = data.frame(xmin = c(0.68, 0.60, 0.58, 0.25, 0.68, 0.29, 0.57, 0.56),
                           xmax = c(0.77, 0.70, 0.65, 0.32, 0.77, 0.37, 0.63, 0.65),
                           ymin = c(0.24, 0.10, 0.35, 0.35, 0.45, 0.59, 0.66, 0.76),
                           ymax = c(0.43, 0.20, 0.45, 0.58, 0.60, 0.70, 0.73, 0.87))

for (i in 1:nrow(rectangles_df)){
  
  # pre-process coordinates
  rect = post
 
  xmin = rectangles_df$xmin[i]
  xmax = rectangles_df$xmax[i]
  ymin = rectangles_df$ymin[i]
  ymax = rectangles_df$ymax[i]
  
  rect = rect[rect$pos_x >= xmin & rect$pos_x <= xmax & rect$pos_y >= ymin & rect$pos_y <= ymax,]
  
  # make plot
  p = ggplot(rect) + aes(x = pos_x, y = pos_y, color = light) + geom_point(size = 3) +
    scale_color_manual(values = c("FALSE" = "gray88", "TRUE" = "brown3"), labels = c("FALSE" = "No", "TRUE" = "Yes")) + theme_bw() +
    labs(color = "Gpr20", x = "x coordinates", y = "y coordinates") +
    guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 2))) +
    theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  if (FALSE){ggsave(str_c("/Users/AmeerSarwar/Downloads/zoom_blood",i,".svg"), device = svg, dpi = 600, width = 6, height = 4)}
  print(p)
}

# spatial enrichment: 1:10 neighbors (proportion with Gpr20) vs 1:10 random cells (proportion with Gpr20)
# performed on original dataset
neigs = 500
neigs = RANN::nn2(data = metadata[,c("pos_x", "pos_y")], query = metadata[,c("pos_x", "pos_y")], k = neigs + 1, treetype = "kd", searchtype = "standard")$nn.idx
tar_data = data[,target]
store = matrix(data = 0, ncol = 2, nrow = ncol(neigs) - 1); colnames(store) = c("Neighbors","Random")
store = as.data.frame(store)

for (i in 2:ncol(neigs)){
  
  # total cells with Gpr20 at n-th neighbor
  store[i-1,"Neighbors"] = t(tar_data) %*% tar_data[neigs[,i]]
  
  # total cells with Gpr20 at random neighbor positions
  rand_ids = sample(x = 1:length(tar_data), size = length(tar_data), replace = FALSE)
  store[i-1,"Random"] = t(tar_data) %*% tar_data[rand_ids]
}

# pre-process for plotting
store = data.frame(neig_id = 1:nrow(store), store)
store$cum_neig = cumsum(store$Neighbors / sum(store$Neighbors))
store$cum_rand = cumsum(store$Random / sum(store$Random))

store1 = as.data.frame(pivot_longer(data = store, cols = 2:3, names_to = "type", values_to = "vals"))
store2 = as.data.frame(pivot_longer(data = store, cols = 4:5, names_to = "type", values_to = "vals"))

# make plot 1
ggplot(store1) + aes(x = neig_id, y = vals, color = type) + geom_point(size = 1) + theme_classic2() +
  scale_color_manual(values = c("Neighbors" = "brown3", "Random" = "gray"), labels = c("Random" = "Random cells")) +
  scale_x_continuous(breaks = c(1,100,200,300,400,500), labels = c(1,100,200,300,400,500)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Neighbor", color = "", subtitle = "Spatial autocorrelation of Gpr20-positive cells",
       y = "Number of Gpr20-expressing neighors of\nGpr20-positive source cells")
ggsave("/Users/AmeerSarwar/Downloads/auto_corr.svg", device = svg, dpi = 600, width = 6, height = 4)

# make plot 2
area_neig = signif(pracma::trapz(x = scale_01(store$neig_id), y = store$cum_neig), digits = 2)
area_rand = signif(pracma::trapz(x = scale_01(store$neig_id), y = store$cum_rand), digits = 2)
pp = signif(wilcox.test(x = store$cum_neig, y = store$cum_rand, alternative = "greater", paired = TRUE)$p.value, digits = 2)

ggplot(store2) + aes(x = neig_id, y = vals, color = type) + geom_line(linewidth = 1) + theme_classic2() +
  geom_segment(aes(x = 1, xend = nrow(store), y = 0, yend = 1), linetype = "dashed", color = "#444444", linewidth = 0.5) +
  scale_color_manual(values = c("cum_neig" = "brown3",    "cum_rand" = "gray"),
                     labels = c("cum_neig" = "Neighbors", "cum_rand" = "Random cells")) +
  scale_x_continuous(breaks = c(1,100,200,300,400,500), labels = c(1,100,200,300,400,500)) +
  labs(x = "Neighbor", y = "Cumulative sum (normalized)", color = "",
       subtitle = str_c("Neighbors AUC = ", area_neig, ", Random cells AUC = ",
                        area_rand, "\nRight-tailed Wilcoxon signed-rank test,\nNeighbors ≥ Random cells, p-value = ", pp))
ggsave("/Users/AmeerSarwar/Downloads/auto_corr_cum.svg", device = svg, dpi = 600, width = 6, height = 4)

```


```{r}

# replicability of blood-like pattern in other slices and replicates
blood = data.frame(replicate = rep(1:3, each = 3), slice = rep(1:3, times = 3))

for (i in 1:nrow(blood)){
  
  # choose slices and replicates
  replicate = blood$replicate[i]
  slice     = blood$slice[i]
  
  # load data
  expr <- fread(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice",slice,"_Replicate",replicate,"_cell_by_gene.csv")); expr  = as.data.frame(expr)
  post <- fread(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice",slice,"_Replicate",replicate,"_cell_metadata.csv")); post = as.data.frame(post)
  
  # pre-process data
  post = post[,4:5]; colnames(post) = c("pos_x","pos_y")
  post = data.frame(post, light = expr[,target] > 0)
  
  # rotate coordinates
  if (i == 1){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 175)}
  if (i == 2){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 165)}
  if (i == 3){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 233)}
  if (i == 4){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 235)}
  if (i == 5){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 40)}
  if (i == 6){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 165)}
  if (i == 7){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 223)}
  if (i == 8){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 353)}
  if (i == 9){post[,c("pos_x","pos_y")] = rotate_coordinates(x = post$pos_x, y = post$pos_y, n_degrees = 65)}
  
  # further pre-process
  post$pos_x = scale_01(post$pos_x)
  post$pos_y = scale_01(post$pos_y)
  
  post$light = factor(post$light, levels = c("FALSE","TRUE"))
  post = post[order(post$light),]

  # make plot
  p = ggplot(post) + aes(x = pos_x, y = pos_y, color = light) + geom_point(size = 0) +
    scale_color_manual(values = c("FALSE" = "gray88", "TRUE" = "brown3"), labels = c("FALSE" = "No", "TRUE" = "Yes")) + theme_classic() +
    labs(color = "Gpr20", x = "x coordinates", y = "y coordinates") +
    guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 2)))
  if (FALSE){ggsave(filename = str_c("/Users/AmeerSarwar/Downloads/blood_full_",i,"_brain.png"), device = png, dpi = 300, width = 6, height = 4)}
  print(p)
  print(i/nrow(blood))
}

```


```{r}

# target gene, e.g., Gpr20, GO enrichment
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/GO_enrichment.R")
target = "Gpr20"

genes = df[df$gene1 %in% target | df$gene2 %in% target,]

# first order (or third order: first + second)
go = genes[genes$edge_type %in% c(1,3),]
go = unique(c(go$gene1, go$gene2))
go = go[!go %in% target]

background_genes = colnames(data)
background_genes = background_genes[!background_genes %in% target]

go = GO_enrichment(test_genes = go, background_genes = background_genes)
go = go$p_values
first = go[as.logical(go$sig_GO),]

# plot GO enrichment
first = first[order(first$pvals_FDR_log10),]

ggplot(first) + aes(x = pvals_FDR_log10, y = fct_reorder(group, pvals_FDR_log10)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
  theme_classic() +
  labs(x = "-log10 p-value", y = "") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="darkred") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("/Users/AmeerSarwar/Downloads/GO.svg", device = svg, dpi = 600, width = 6, height = 4)

```


```{r}

# tissue expression plots for genes cross-expressed with selected genes
selected_gene = "Gpr20"
gaba = df[df$gene1 %in% selected_gene | df$gene2 %in% selected_gene,]
gaba = unique(c(gaba$gene1, gaba$gene2))
gaba = sort(gaba[!gaba %in% selected_gene])

for (i in 1:length(gaba)){
  p = tissue_expression_plot(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = selected_gene, gene2 = gaba[i], cross_expression = TRUE , neighbor = 1, point_size = 0.1, scale_bar = 0); print(p)
  p = tissue_expression_plot(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = selected_gene, gene2 = gaba[i], cross_expression = FALSE , neighbor = 1, point_size = 0.1, scale_bar = 0); print(p)
  print(i/length(gaba))
}

```


```{r}

# what happens to network as diff. count threshold is applied?
min_mol = 1:10
thresh  = matrix(data = 0, ncol = length(min_mol), nrow = choose(ncol(data)-1, 2))
colnames(thresh) = min_mol

for (i in 1:length(min_mol)){
  temp_data   = data[,2:ncol(data)]
  temp_data[temp_data < min_mol[i]] = 0
  thresh[,i]  = cross_expression_correlation(data = temp_data, locations = metadata[,c("pos_x","pos_y")])$correlation
  print(i/length(min_mol))
}

simm = cor(thresh, use = "complete.obs")
map  = Heatmap(simm, cluster_rows = FALSE, cluster_columns = FALSE, name = "Cross-expression\nCorrelation",
               column_title = "Minimum number of molecules needed to threshold", column_title_side = "bottom",
               row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 15),
               column_title_gp = gpar(fontsize = 18))

svg("/Users/AmeerSarwar/Downloads/heatmap.svg", width = 10, height = 8)
map
dev.off()

map = data.frame(x = upper_tri(simm))

ggplot(map) + aes(x = x, y = after_stat(..scaled..)) +
  geom_density(fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = median(map$x), color = "brown3", linetype = "dashed") +
  labs(x = "Cross-expression similarity (correlation) for different noise thresholds",
       y = "Density", subtitle = str_c("Median Pearson's R = ", signif(median(map$x), digits = 2))) +
  theme_classic2() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/density.svg", device = "svg", dpi = 600, width = 6, height = 4)

```


```{r}

# anatomical markers Lgr6 vs Adra2b and Lgr6 vs Ret as cross-expression spatial enrichment
pp = spatial_enrichment(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = "Adra2b", gene2 = "Ret")
pp = data.frame(x = c(pp$target, pp$null), type = rep(c("Cross-expressing", "Random"), times = c(length(pp$target), length(pp$null))))
pp$type = factor(pp$type, levels = c("Random", "Cross-expressing"))

ggplot(pp) + aes(x = x, fill = type, y = after_stat(scaled)) + geom_density(alpha = 0.8) +
  theme_classic2() + labs(x = "Distance to cells", y = "Density", fill = "") +
  scale_fill_manual(values = c("Random" = "gray88","Cross-expressing" = "lightblue")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.position = "top",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/spatial_enrichment.svg", device = "svg", dpi = 600, width = 5, height = 4)

p = tissue_expression_plot(data = data, locations = metadata[,c("pos_x","pos_y")], gene1 = "Adra2b", gene2 = "Ret", cross_expression = FALSE)
p + theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          legend.position = "none")

ggsave("/Users/AmeerSarwar/Downloads/brain.png", device = png, dpi = 600, width = 4, height = 2)

```


```{r}

# co-expression within cells vs co-expression vs cross-expression with neighbors
corr_cell  = cor(data)
corr_neig  = cor(data, data[RANN::nn2(data = metadata[,c("pos_x","pos_y")], query = metadata[,c("pos_x","pos_y")], k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2],])
corr_cross = cross_expression_correlation(data = data, locations = metadata[,c("pos_x","pos_y")], output_matrix = TRUE)

df = data.frame(correlation_within_cells = upper_tri(corr_cell), correlation_between_neighbors = upper_tri(corr_neig), correlation_between_cross_expressing_neighbors = upper_tri(corr_cross))
df = as.data.frame(pivot_longer(data = df, cols = 2:3, names_to = "type", values_to = "vals"))

# plot
ggplot(df) + aes(x = correlation_within_cells, y = vals, color = type) + geom_point(size = 0) + facet_wrap(~type) +
  theme(legend.position = "none") + labs(y = "Correlation")




df = data[,sample(1:109, size = 500, replace = TRUE)]
id = sample(1:nrow(data), size = 100000, replace = FALSE)
df = df[id,]
cr = coord[id,c("CCFx","CCFy","CCFz")]

start = Sys.time()
xx = cross_expression(data = df, locations = cr)
end = Sys.time() - start

```


```{r}

library(arrow)
df = read_parquet("/Users/AmeerSarwar/Downloads/transcripts.parquet")
df = as.data.frame(df)

fwrite(df, "/Users/AmeerSarwar/Downloads/transcripts.csv")
adata = AnnData(X = df)


```