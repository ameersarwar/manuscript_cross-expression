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

# comment out if considering full dataset
#cortex = c("CTX","VISam","AOB","AUDpo","COA","TT","ILA","VISli","PERI","VISpl","VISpm","FRP","ORBm","VISal","AUDd","VISrl","VISpor","AIv","SSp-un","BMA","BLA","CLA","PA","LA","EP","SSp-tr","VISl",
 #          "ORBvl","VISa","ECT","GU","ACAv","AUDv","AUDp","DG","VISC","COAp","SSp-ll","PL","AIp","ORBl","CTXsp","ACAd","AId","RSPagl","TEa","ENTm","CTXpl","RSPv","SSp-n","SSp-ul","ENTl","RSPd","SSp-m",
  #         "SSp-bfd","CA","RHP","MOp","SSs","VISp","MOs","OLF")

#metadata = metadata[metadata$CCFparentname %in% cortex,]
#data     = data[rownames(data) %in% metadata$sample_id,]

```


```{r}

# change in network topology for adjacent slices
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

slices = unique(metadata$slice)
result = matrix(data = 0, nrow = choose(ncol(data), 2), ncol = length(slices))

for (i in 1:length(slices)){
  
  # extract data and metadata for each slice
  temp_meta = metadata[metadata$slice %in% slices[i],]
  temp_data = data[rownames(data) %in% temp_meta$sample_id,]
  
  # cross-expression between gene pairs
  cross = cross_expression(data = temp_data, locations = temp_meta[,c("pos_x","pos_y")])
  
  # store results
  result[,i] = cross$cross_pvalue
  print(i/length(slices))
}

```


```{r}

# slice by slice correlation
df = cor(result, method = "spearman")

# distance vs correlation
df = data.frame(slice_locus  = rep(slices, each  = nrow(df)),
                slice_target = rep(slices, times = nrow(df)),
                values = as.vector(df))

# process data set
df$slice_locus = df$slice_locus - 1
df$slice_target= df$slice_target - 1
df$slice_locus = str_c("Slice ", df$slice_locus)
df$slice_locus = factor(df$slice_locus, levels = str_c("Slice ", slices - 1))

# plot
ggplot(df) + aes(x = slice_target, y = values, color = slice_locus) +
  geom_point(size = 0.7) +
  geom_line(linewidth = 0.2) +
  facet_wrap(~slice_locus, ncol = 5, nrow = 3) +
  labs(x = "Adjacent slice", y = "Correlation") +
  scale_x_continuous(breaks = seq(from = min(slices - 1), to = max(slices - 1), by = 2)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "#BCD4E6"),
        strip.text = element_text(size = 8))

```


```{r}

# plot correlations as a function of distance/slice
df = cor(result, method = "spearman")
n  = nrow(df)
correlations = c()
distances    = c()

# correlation versus distance/slice
for (i in 1:(n-1)) {
    for (j in (i+1):n) {
        correlations <- c(correlations, df[i, j])
        distances    <- c(distances, abs(j - i))
    }
}

# pre-processing and plotting
plot_data <- data.frame(correlation = correlations, distance = distances)

ggplot(plot_data) + aes(x = distance, y = correlation) +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Distance (in slices)", y = "Correlation",
       subtitle = str_c("Spearman's rho = ", signif(cor(plot_data$correlation, plot_data$distance, method = "spearman"), digits = 2))) +
  scale_x_continuous(breaks = 1:(n-1), labels = 1:(n-1)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15, margin = margin(t = 10, unit = "pt")),
        axis.title.y = element_text(size = 15, margin = margin(r = 10, unit = "pt")),
        axis.text = element_text(size = 12),
        plot.subtitle = element_text(size = 15),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))

ggsave(filename = "/Users/AmeerSarwar/Downloads/correlation_vs_distance.svg", bg = "transparent", device = "svg", dpi = 600, height = 4, width = 6)

```