```{r}

# BAR-seq coronal data
# data is shrunk by removing image stitching-related artefacts (cf. Xiaoyin's email)
# data is quality controlled by keeping cells with genes/cell >= 5 and reads/cell >= 20
# data alongside CCF and slide coordinates are saved and can be used for analysis

# load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus","Seurat","Rcpp","RcppArmadillo")
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

coordinates <- coord %>% filter(slice==slice_current); coord_ref <- coord_ref %>% filter(slice==slice_current)
data = data[,which(colnames(data) %in% coordinates$sample_id)]
ggplot(coordinates) + aes(x=pos_x, y=pos_y, color=CCFparentname) + geom_point(size=0.5)

```


```{r}

# cross-expression is orthogonal to cell types
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
region = "SSp-bfd";  region = c("AUDv","AUDd","AUDp")
metadata = coordinates
metadata_subset1 <- metadata[which(metadata$CCFparentname %in% region),]
metadata_subset1 <- metadata_subset1 %>% filter(pos_x <= 8000)
metadata_subset1 <- metadata_subset1[!metadata_subset1$subclass_H2 %in% "Unclear",]
data_subset1     <- t(data); data_subset1 <- data_subset1[which(rownames(data_subset1) %in% metadata_subset1$sample_id),]

# example
tissue.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x","pos_y")],
            gene_pair_1 = "Gfra1",
            gene_pair_2 = "Foxp2",
            cross_expression_highlight = TRUE,
            cell_type_labels = metadata_subset1$subclass_H2,
            legend_title = "") +
  theme(legend.background = element_blank(),
        legend.key = element_blank())

```


```{r}

# cross-expression is orthogonal to cell types
# for each gene pair check proportions of cell-neighbor pairs w/ same vs diff. labels
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# combine multiple slices, 10-30
coord    = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
data     = as.matrix(readMM(file="/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt"))
rownames(data) = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv", header = TRUE))$Genes; colnames(data) = coord$sample_id
coord    = coord[which(coord$sample_id %in% as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_cortex.csv"))$sample_id),]
data     = data[,which(colnames(data)  %in% as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_cortex.csv"))$sample_id)]

cc_total = coord[coord$slice %in% 10:30, ]
df_total = t(data[, colnames(data) %in% cc_total$sample_id])

# cross-expression
cross = cross_expression(data = df_total, locations = cc_total[,c("pos_x","pos_y")])
#cross = cross[as.logical(cross$cross_sig),]
cross = cross[cross$cross_pvalue <= 0.05, ]
locations = cc_total[,c("pos_x","pos_y")]

prop = data.frame(matrix(data = 0, nrow = nrow(cross), ncol = 2))
prop = data.frame(data.frame(cross$gene1, cross$gene2), prop)
colnames(prop) = c("gene1","gene2","same","num")

for (i in 1:nrow(cross)){
  
  # gene pairs
  df = df_total[, c(cross$gene1[i], cross$gene2[i])]
  
  # neighbors
  distances <- RANN::nn2(locations, locations, k = 2, searchtype = "priority")
  distances <- distances$nn.idx[,2]
  df_temp   <- df[distances,]
  
  # pairs
  pair1 = df[,1] > 0 & df[,2] == 0 & df_temp[,2] > 0 & df_temp[,1] == 0
  pair2 = df[,1] == 0 & df[,2] > 0 & df_temp[,2] == 0 & df_temp[,1] > 0
  
  cells = c(names(pair1[pair1]), names(pair2[pair2]))
  neigs = c(rownames(df_temp[pair1, ]), rownames(df_temp[pair2, ]))
  
  # compare cell-neighbor labels
  cells = cc_total[match(cells, cc_total$sample_id), "subclass_H2"]
  neigs = cc_total[match(neigs, cc_total$sample_id), "subclass_H2"]
  
  # save proportions
  prop[i,"same"] = sum(cells == neigs) / length(cells)
  prop[i,"num"] = length(cells)
  print(i/nrow(cross))
}

```


```{r}

# dot plot
ggplot(prop) + aes(x = gene1, y = gene2, color = same, size = num) +
  geom_point() +
  scale_color_gradient(low = "#80CC6D", high = "#FFB347") +
  labs(x = "", y = "", color = "Cell type\npurity", size = "Cell-neighbor\npairs") +
  guides(size = guide_legend(order = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7))
ggsave(filename = "/Users/AmeerSarwar/Downloads/cell_type_orthogonal1.svg", device = "svg", height = 4, width = 6, dpi = 600)

# purity vs distance
ggplot(prop) + aes(x = num, y = same, color = same) +
  geom_point(size = 2) +
  scale_color_gradient(low = "#80CC6D", high = "#FFB347") +
  geom_smooth(method = "lm", se = TRUE, color = "darkgray") +
  labs(x = "Number of cross-expressing cell-neighbor pairs",
       y = "Cell type purity",
       title = str_c("Spearman's rho = ", signif(cor(prop$same, prop$num, method = "spearman"), digits = 2))) +
  theme_classic() + theme(legend.position = "none")
ggsave(filename = "/Users/AmeerSarwar/Downloads/cell_type_orthogonal2.svg", device = "svg", height = 4, width = 5, dpi = 600)

```


```{r}

# densities of cell type purity and number of neighbors
same = data.frame(x = prop$same)
num  = data.frame(x = prop$num)

ggplot(same) + aes(x = x, y = ..scaled..) +
  geom_density(color = "darkgreen", fill = "#80CC6D") +
  labs(x = "Proportion of cross-expressing cell pairs with the same cell type label", y = "Density") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/purity.svg", device = "svg", width = 6, height = 4, dpi = 600)

ggplot(num) + aes(x = x, y = ..scaled..) +
  geom_density(color = "#737000", fill = "#FFB347") +
  labs(x = "Number of cross-expressing cell pairs", y = "Density") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/number.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# assess whether cells cross-express w/ other cells of the same vs. different types
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
coord    = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
data     = as.matrix(readMM(file="/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt"))
rownames(data) = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv", header = TRUE))$Genes; colnames(data) = coord$sample_id
coord    = coord[which(coord$sample_id %in% as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_cortex.csv"))$sample_id),]
data     = data[,which(colnames(data)  %in% as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_cortex.csv"))$sample_id)]

slices   = table(coord$slice)
slices   = as.numeric(names(slices[slices >= 10000]))

cc_total = coord[coord$slice %in% slices, ]
df_total = t(data[, colnames(data) %in% cc_total$sample_id])

to_remove= c("Low Quality","Unclear")

cc_total = cc_total[!cc_total$subclass_H2 %in% to_remove,]
df_total = df_total[rownames(df_total) %in% cc_total$sample_id,]

cell_types = unique(cc_total$subclass_H2)
cell_types = matrix(data = 0, nrow = length(cell_types), ncol = length(cell_types), dimnames = list(cell_types, cell_types))

slices  = unique(cc_total$slice)

for (i in 1:length(slices)){
  
  # cell and neighbor info per brain slice
  temp_coord = cc_total[cc_total$slice %in% slices[i],]
  temp_data  = df_total[rownames(df_total) %in% temp_coord$sample_id,]
  
  neig_coord = RANN::nn2(data = temp_coord[,c("CCFx","CCFy","CCFz")], query = temp_coord[,c("CCFx","CCFy","CCFz")], k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2]
  neig_data  = temp_data[neig_coord,]
  neig_coord = temp_coord[neig_coord,]
  
  # cell type frequencies per cross-expressing gene pair in the slice
  exp_cross = cross_expression(data = temp_data, locations = temp_coord[,c("CCFx","CCFy","CCFz")])
  exp_cross = exp_cross[as.logical(exp_cross$cross_sig), c("gene1","gene2")]
  exp_cross = rbind(exp_cross, data.frame(gene1 = exp_cross$gene2, gene2 = exp_cross$gene1))
  if (nrow(exp_cross)==0){next}
  
  for (j in 1:nrow(exp_cross)){
    
    gene1 = exp_cross$gene1[j]
    gene2 = exp_cross$gene2[j]
    
    cross = (temp_data[,gene1] > 0) & (temp_data[,gene2] == 0) & (neig_data[,gene1] == 0) & (neig_data[,gene2] > 0)
    row_cell = temp_coord$subclass_H2[cross]
    col_cell = neig_coord$subclass_H2[cross]
    
    # add counter to cell types
    for (id in 1:length(row_cell)){cell_types[which(rownames(cell_types) %in% row_cell[id]), which(colnames(cell_types) %in% col_cell[id])] = cell_types[which(rownames(cell_types) %in% row_cell[id]), which(colnames(cell_types) %in% col_cell[id])] + 1}
  }
  print(i/length(slices))
}

# normalize by cell type pair frequency
output = cell_types

slices = unique(cc_total$slice)
freq_back = output
freq_back[freq_back != 0] = 0

genes = 0

for (i in 1:length(slices)){
  
  temp_coord = cc_total[cc_total$slice %in% slices[i],]
  neig_coord = RANN::nn2(data = temp_coord[,c("CCFx","CCFy","CCFz")], query = temp_coord[,c("CCFx","CCFy","CCFz")], k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2]
  neig_coord = temp_coord[neig_coord,]
  
  exp_cross  = cross_expression(data = df_total[rownames(df_total) %in% temp_coord$sample_id,], locations = temp_coord[,c("pos_x","pos_y")])
  genes      = genes + sum(exp_cross$cross_sig)
  
  slice_freq = freq_back
  slice_freq[slice_freq != 0] = 0
  
  for (type1 in 1:nrow(freq_back)){
    
    ids = temp_coord$subclass_H2 %in% rownames(freq_back)[type1]
    ids = neig_coord$subclass_H2[ids]
    
    for (type2 in 1:nrow(freq_back)){slice_freq[type1,type2] = sum(ids %in% colnames(freq_back)[type2])}
  }
  freq_back = freq_back + slice_freq
  print(i/length(slices))
}

# normalize
freq_norm  = freq_back*genes
freq_norm[freq_norm == 0] = 1
cell_norm  = cell_types
cell_norm  = cell_norm / freq_norm

# view relationships of frequent cross-expressing cell types
ht = Heatmap(cell_norm, name = "Normalized Frequencies", row_names_gp = gpar(fontsize = 6.5), column_names_gp = gpar(fontsize = 6.5), show_heatmap_legend = TRUE,
             heatmap_legend_param = list(legend_direction = "horizontal"))
draw(ht)

# save heatmap
svg(filename = "/Users/AmeerSarwar/Downloads/heatmap_cell_neighbor_types.svg", width = 5, height = 4)
draw(ht)
dev.off()

# distribution
distribution = cell_norm
diag(distribution) = 0
distribution = c(upper_tri(distribution), upper_tri(t(distribution)))
distribution = data.frame(x = distribution)

ggplot(distribution) + aes(x = x * 100, y = ..scaled..) +
  geom_density(color = "steelblue", fill = "lightblue") +
  labs(x = "Norm. freq. of neighbor types with which cell types cross-expresses (x100)", y = "Density",
       subtitle = str_c("Median (x100) = ", signif(median(distribution$x)*100, digits = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.subtitle = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/number.svg", device = "svg", width = 6, height = 4, dpi = 600)

```


```{r}

# cross-expression b/w cell types
# pick cells of one type and ask if neighbors are consistently of another type
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/Distance Metrics/jaccard.R")

coord = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv"))
coord = coord[which(coord$sample_id %in% as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_cortex.csv"))$sample_id),]
coord = coord[!coord$subclass_H2 %in% c("Unclear","Low Quality"),]

cells = coord$subclass_H2
neigs = cells[RANN::nn2(data = coord[,c("CCFx","CCFy","CCFz")], query = coord[,c("CCFx","CCFy","CCFz")], k = 2, treetype = "kd", searchtype = "standard")$nn.idx[,2]]

# p-value
types = sort(unique(cells))
pvals = matrix(data = 0, nrow = length(types), ncol = length(types), dimnames = list(types, types))

for (i in 1:length(types)){
  
  for (j in 1:length(types)){
    
    k = cells %in% types[i]
    m = neigs %in% types[j]
    q = as.numeric(k) %*% as.numeric(m)
    pvals[i,j] = phyper(q = q - 1, m = sum(m), n = length(cells) - sum(m), k = sum(k), lower.tail = FALSE)
  }
  print(i/length(types))
}

outcome = pvals
diag(outcome) = 1

# heatmap
p_value_vector <- as.vector(outcome)
fdr_corrected_vector <- p.adjust(p_value_vector, method = "BH")
fdr_corrected_matrix <- matrix(fdr_corrected_vector, nrow = nrow(outcome), byrow = TRUE, dimnames = list(rownames(outcome), colnames(outcome)))
matt = fdr_corrected_matrix
matt[matt == 0] = min(matt[matt != 0])
matt = -log10(matt)

Heatmap(matt)

# further analysis
val_ids = which(upper.tri(outcome), arr.ind = TRUE)
vals    = c(upper_tri(outcome), upper_tri(t(outcome)))
cell_n  = rownames(outcome)[val_ids[,1]]
neig_n  = colnames(outcome)[val_ids[,2]]
cell_nn = c(cell_n, neig_n)
neig_nn = c(neig_n, cell_n)

pvalues = data.frame(cells = cell_nn, neigs = neig_nn, pvals = vals)
pvalues$pvals = p.adjust(pvalues$pvals, method = "BH")

pvalues = cbind(pvalues[1:(nrow(pvalues)/2),], data.frame(pvals2 = pvalues[((nrow(pvalues)/2)+1):nrow(pvalues),"pvals"]))
pvalues = data.frame(cells = pvalues$cells, neigs = pvalues$neigs, pvals = Pmin(pvalues$pvals, pvalues$pvals2))

# network analysis
net = pvalues[pvalues$pvals <= 0.05, c("cells","neigs")]
net$edge = rep(1,nrow(net))
deg = table(c(net$cells, net$neigs))
deg = data.frame(cell = names(deg), deg = as.numeric(deg))
deg = deg[order(deg$deg),]

write.csv(deg, file = "/Users/AmeerSarwar/Downloads/cell_type_deg.csv",  quote = FALSE, row.names = FALSE)
write.csv(net, file = "/Users/AmeerSarwar/Downloads/cell_type_edge.csv", quote = FALSE, row.names = FALSE)

```


```{r}

# cross-expression recovers layer boundaries
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")
region = unique(metadata$CCFparentname)
region = region[!(region %in% c("PIR", "RSPd", "RSPv"))]
metadata_subset1 <- metadata[which(metadata$CCFparentname %in% region),]
cells  = unique(metadata_subset1$subclass_H2)
cells  = cells[!(cells %in% c("RSP DL", "RSP UL"))]
metadata_subset1 <- metadata_subset1[which(metadata_subset1$subclass_H2 %in% cells),]
data_subset1     <- t(data); data_subset1 <- data_subset1[which(rownames(data_subset1) %in% metadata_subset1$sample_id), ]

gene1 = "Cdh13"
gene2 = "Foxp2"

tissue.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x", "pos_y")],
            gene_pair_1 = gene1,
            gene_pair_2 = gene2,
            point_size = 0.1,
            scale_bar_length = 0) +
    theme(legend.position = "none",
          panel.grid = element_blank())
ggsave(filename = "/Users/AmeerSarwar/Downloads/layer_recovery_full.png", device = "png", height = 4, width = 6, dpi = 600)

tissue.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x", "pos_y")],
            gene_pair_1 = gene1,
            gene_pair_2 = gene2,
            cross_expression_highlight = TRUE,
            point_size = 0.1,
            scale_bar_length = 0) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"))
ggsave(filename = "/Users/AmeerSarwar/Downloads/layer_recovery_cross.png", device = "png", height = 4, width = 6, dpi = 600)

# cortical slice with cell types
df = metadata_subset1

ggplot(df) + aes(x = pos_x, y = pos_y, color = subclass_H2) +
  geom_point(size = 0.1) +
  labs(x="", y="", color="") +
  scale_color_brewer(palette = "Paired") +
  geom_segment(aes(x = (0.95 * max(pos_x)) - 500, xend = 0.95 * max(pos_x),
                       y = min(pos_y) + (0.05 * min(pos_y)), yend = min(pos_y) + (0.05 * min(pos_y))), color = "black", linewidth = 1) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(legend.position.inside = c(0.15,0.4),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
ggsave(filename = "/Users/AmeerSarwar/Downloads/cell_type_brainmap.png", device = "png", height = 4, width = 6, dpi = 600)

```