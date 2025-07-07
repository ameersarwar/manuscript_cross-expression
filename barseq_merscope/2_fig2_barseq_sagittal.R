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
slices   = c(12,13,14,15) # nice looking slices

# communicating genes
ligands_peptide     = c("Cck",  "Sst",  "Vip",  "Penk", "Penk", "Pdyn", "Npy",  "Adcyap1")
receptors_peptide   = c("Cckbr","Sstr2","Vipr1","Oprm1","Oprd1","Oprk1","Npy1r","Adcyap1r1")
receptors_monoamine = c("Chrm1","Chrm2","Chrm3","Chrna7", "Adra1a","Adrb1","Htr1f","Htr2c","Htr3a","Drd1") # Adra2a not found
genes = sort(c(ligands_peptide, receptors_peptide, receptors_monoamine))

# comment out if considering full dataset
#data  = data[,which(colnames(data) %in% genes)]

```


```{r}

# two good looking examples of cross-expression
# best results
# slice 3: Sst to Sstr2 in VISC
# slice 5: Sstr2 to Sst in SSp-n and Penk to Oprd1 in VISC
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

examples = data.frame(slice = c(3,5), region = c("VISC","SSp-n"), gene1 = c("Sst","Sstr2"), gene2 = c("Sstr2","Sst"))

for (i in 1:nrow(examples)){
  
  meta_temp = metadata[metadata$slice %in% examples$slice[i] & metadata$CCFparentname %in% examples$region[i], ]
  df_temp   = data[rownames(data) %in% meta_temp$sample_id, c(examples$gene1[i], examples$gene2[i])]
  
  print(cross_expression(data = df_temp, locations = meta_temp[,c("pos_x","pos_y")]))
  
  # full tissue
  p = tissue_expression_plot(data = df_temp, locations = meta_temp[,c("pos_x","pos_y")], gene1 = examples$gene1[i], gene2 = examples$gene2[i], point_size = 1, cross_expression = FALSE)
  print(p)
  
  p = tissue_expression_plot(data = df_temp, locations = meta_temp[,c("pos_x","pos_y")], gene1 = examples$gene1[i], gene2 = examples$gene2[i], point_size = 1, cross_expression = TRUE)
  print(p)
  
  p = tissue_expression_plot(data = df_temp, locations = meta_temp[,c("pos_x","pos_y")], gene1 = examples$gene1[i], gene2 = examples$gene2[i], point_size = 1, cross_expression = FALSE) +
    theme_bw() + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.border = element_blank(),
                       legend.text = element_text(size = 12), legend.position = "none", legend.background = element_blank(), legend.box = element_blank(), legend.key = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  print(p)
  
  
  # cross-expression
  p = tissue_expression_plot(data = df_temp, locations = meta_temp[,c("pos_x","pos_y")], gene1 = examples$gene1[i], gene2 = examples$gene2[i], point_size = 1, cross_expression = TRUE) +
    theme_bw() + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.border = element_blank(),
                       legend.text = element_text(size = 12), legend.position = "none", legend.background = element_blank(), legend.box = element_blank(), legend.key = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
  print(p)
  
  # bullseye
  score = bullseye_scores(data = df_temp, locations = meta_temp[,c("pos_x","pos_y")], ratio_to_co = TRUE, log_2 = FALSE)
  p = bullseye_plot(scores = score[1,3:ncol(score)])
  print(p)
  ggsave(str_c("/Users/AmeerSarwar/downloads/bull_plot_",i,".svg"), device = "svg", dpi = 600, width = 5, height = 5)
}


```


```{r}

# two good looking examples of cross-expression
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")

# best results
# slice 3: Sst to Sstr2 in VISC
# slice 5: Sstr2 to Sst in SSp-n and Penk to Oprd1 in VISC

# example 1
ligand   = "Sst"
receptor = "Sstr2"
region   = "VISC"
slice_cr = 3

metadata_subset1 = metadata %>% filter(slice %in% slice_cr, CCFparentname %in% region)
data_subset1     = data[which(rownames(data) %in% metadata_subset1$sample_id), ]
genes_total1     = celltocell(data = data_subset1,
                              locations = metadata_subset1[,c("pos_x","pos_y")],
                              gene_pair_1 = ligand,
                              gene_pair_2 = receptor)

tissue.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x","pos_y")],
            gene_pair_1 = ligand,
            gene_pair_2 = receptor,
            legend_title = "")
ggsave(filename = "/Users/AmeerSarwar/Downloads/example1_full.svg", device = "svg", width = 8, height = 6, dpi = 600)

tissue.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x","pos_y")],
            gene_pair_1 = ligand,
            gene_pair_2 = receptor,
            cross_expression_highlight = TRUE)
ggsave(filename = "/Users/AmeerSarwar/Downloads/example1_cross.svg", device = "svg", width = 8, height = 6, dpi = 600)

bull_stats = bullseye.statistics(data_subset1, metadata_subset1[,c("pos_x","pos_y")], ligand, receptor)
bullseye.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x","pos_y")],
            gene_pair_1 = ligand,
            gene_pair_2 = receptor)
ggsave(filename = "/Users/AmeerSarwar/Downloads/example1_bull.svg", device = "svg", width = 6, height = 6, dpi = 600)

metadata %>% filter(slice %in% slice_cr) %>% mutate(CCFparentname = if_else(CCFparentname %in% region, CCFparentname, NA_character_)) %>%
ggplot() + aes(x=pos_x, y=pos_y, color=CCFparentname) + geom_point(size=0.1) + labs(x="", y="", color="") +
  geom_segment(aes(x = (0.95 * max(pos_x)) - 1000, xend = 0.95 * max(pos_x),
                       y = min(pos_y) + (0.05 * min(pos_y)), yend = min(pos_y) + (0.05 * min(pos_y))), color = "black", linewidth = 1) +
  scale_color_manual(values = c("VISC" = "brown3"), na.value = "gray88") +
  theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none",
                     panel.border = element_blank(), panel.grid = element_blank())
ggsave(filename = "/Users/AmeerSarwar/Downloads/example1_map.png", device = "png", width = 6, height = 4, dpi = 300)


# example 2
ligand   = "Sstr2"
receptor = "Sst"
region   = "SSp-n"
slice_cr = 5

metadata_subset1 = metadata %>% filter(slice %in% slice_cr, CCFparentname %in% region)
data_subset1     = data[which(rownames(data) %in% metadata_subset1$sample_id), ]
genes_total2     = celltocell(data = data_subset1,
                              locations = metadata_subset1[,c("pos_x","pos_y")],
                              gene_pair_1 = ligand,
                              gene_pair_2 = receptor)

tissue.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x","pos_y")],
            gene_pair_1 = ligand,
            gene_pair_2 = receptor,
            legend_title = "")
ggsave(filename = "/Users/AmeerSarwar/Downloads/example2_full.svg", device = "svg", width = 8, height = 6, dpi = 600)

tissue.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x","pos_y")],
            gene_pair_1 = ligand,
            gene_pair_2 = receptor,
            cross_expression_highlight = TRUE)
ggsave(filename = "/Users/AmeerSarwar/Downloads/example2_cross.svg", device = "svg", width = 8, height = 6, dpi = 600)

bull_stats = bullseye.statistics(data_subset1, metadata_subset1[,c("pos_x","pos_y")], ligand, receptor)
bullseye.plot(data = data_subset1,
            locations = metadata_subset1[,c("pos_x","pos_y")],
            gene_pair_1 = ligand,
            gene_pair_2 = receptor)
ggsave(filename = "/Users/AmeerSarwar/Downloads/example2_bull.svg", device = "svg", width = 6, height = 6, dpi = 600)

metadata %>% filter(slice %in% slice_cr) %>% mutate(CCFparentname = if_else(CCFparentname %in% region, CCFparentname, NA_character_)) %>%
ggplot() + aes(x=pos_x, y=pos_y, color=CCFparentname) + geom_point(size=0.1) + labs(x="", y="", color="") +
  geom_segment(aes(x = (0.95 * max(pos_x)) - 1000, xend = 0.95 * max(pos_x),
                       y = min(pos_y) + (0.05 * min(pos_y)), yend = min(pos_y) + (0.05 * min(pos_y))), color = "black", linewidth = 1) +
  scale_color_manual(values = c("SSp-n" = "brown3"), na.value = "gray88") +
  theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none",
                     panel.border = element_blank(), panel.grid = element_blank())
ggsave(filename = "/Users/AmeerSarwar/Downloads/example2_map.png", device = "png", width = 6, height = 4, dpi = 300)

```


```{r}

# bullseye densities for significant and non-significant gene pairs
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")

region   = "SSp-n"
slice_cr = 5

metadata_subset1 = metadata %>% filter(slice %in% slice_cr, CCFparentname %in% region)
data_subset1     = data[which(rownames(data) %in% metadata_subset1$sample_id), ]

# p-values for cross-expression
outcome = celltocell2(data = data_subset1,
                 locations = metadata_subset1[,c("pos_x","pos_y")])

alpha = 0.05

# significant pairs cross-expression w/o co-expression
outcome_sig = outcome$pvalue_cross_expression
outcome_sig[outcome_sig == 0] = 1
outcome_ctr = outcome$pvalue_co_expression
outcome_ctr = ifelse(outcome_ctr <= alpha, -1, 1)

outcome_sig = outcome_sig * outcome_ctr
outcome_sig = ifelse(outcome_sig < 0, 1, outcome_sig)
outcome_sig = Diag.fill(outcome_sig, v = 1)
outcome_sig = which(outcome_sig <= alpha, arr.ind = TRUE)
size_comp   = nrow(outcome_sig)

outcome_sig = bullseye.statistics(data = data_subset1,
                    locations = metadata_subset1[,c("pos_x","pos_y")],
                    gene_pair_1 = outcome_sig[,1],
                    gene_pair_2 = outcome_sig[,2],
                    window_sizes = 1:10)

# non-significant pairs
outcome_nsig = outcome$pvalue_cross_expression
outcome_nsig[outcome_nsig == 0] = 1
outcome_nsig = Diag.fill(outcome_nsig, v = 1)
outcome_nsig = which(outcome_nsig > alpha, arr.ind = TRUE)

size_comp    = sample(1:nrow(outcome_nsig), size = size_comp, replace = FALSE)
outcome_nsig = outcome_nsig[size_comp,]

outcome_nsig = bullseye.statistics(data = data_subset1,
                    locations = metadata_subset1[,c("pos_x","pos_y")],
                    gene_pair_1 = outcome_nsig[,1],
                    gene_pair_2 = outcome_nsig[,2],
                    window_sizes = 1:10)

# combine significant and non-significant pairs' bullseye densities
outcome_sig$class  = rep("sig",  nrow(outcome_sig))
outcome_nsig$class = rep("nsig", nrow(outcome_nsig))

df = rbind(outcome_sig, outcome_nsig)
df$class = factor(df$class, levels = c("sig","nsig"))

# compare significant and non-significant pairs' bullseye densities
ggplot(df) + aes(x = as.factor(neighbors.in.window), y = log10(bullseye.score), fill = class) +
  geom_boxplot(outlier.size = 0.1) +
  labs(x = "Neighbors",
       y = expression("Bullseye score, log"[10]*""),
       fill = "Cross-expression") +
  scale_x_discrete(labels = c("0" = "Cell")) +
  scale_fill_manual(values = c("sig" = "#4A82D9", "nsig" = "#C3E8FF"),
                    labels = c("sig" = "Significant", "nsig" = "Not significant")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85,0.85),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave(filename = "/Users/AmeerSarwar/Downloads/global_bull.svg", device = "svg", height = 4, width = 6, dpi = 600)

# cross-expression to co-expression ratio (effect size) for significant and non-significant groups
self_sig  = df %>% filter(neighbors.in.window %in% 0, class %in% "sig")  %>% dplyr::select(bullseye.score)
self_nsig = df %>% filter(neighbors.in.window %in% 0, class %in% "nsig") %>% dplyr::select(bullseye.score)

neig_sig  = df %>% filter(neighbors.in.window %in% 1, class %in% "sig")  %>% dplyr::select(bullseye.score)
neig_nsig = df %>% filter(neighbors.in.window %in% 1, class %in% "nsig") %>% dplyr::select(bullseye.score)

sig  = neig_sig$bullseye.score  / self_sig$bullseye.score
nsig = neig_nsig$bullseye.score / self_nsig$bullseye.score

sig   = data.frame(bull = sig,  class = rep("sig",  length(sig)))
nsig  = data.frame(bull = nsig, class = rep("nsig", length(nsig)))

df2   = rbind(sig, nsig)

# subplot of these effect size distributions
ggplot(df2) + aes(y = class, fill = class, x = bull) +
  geom_boxplot(outlier.size = 0.1) +
  labs(y = "", x = "Cross-expression to co-expression") +
  scale_fill_manual(values =  c("sig" = "#4A82D9", "nsig" = "#C3E8FF")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.background = element_blank(),
        panel.background = element_blank())
ggsave(filename = "/Users/AmeerSarwar/Downloads/ratio_bull.svg", device = "svg", height = 2, width = 3, dpi = 600)

```
