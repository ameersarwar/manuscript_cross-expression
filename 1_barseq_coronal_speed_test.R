```{r}

# BAR-seq coronal data
# data is shrunk by removing image stitching-related artefacts (cf. Xiaoyin's email)
# data is quality controlled by keeping cells with genes/cell >= 5 and reads/cell >= 20
# data alongside CCF and slide coordinates are saved and can be used for analysis

# load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate","pryr",
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
#slice_current = round(runif(1, min = min(unique(coord$slice)), max = max(unique(coord$slice))))
#slice_current = 17
#metadata <- coordinates

#coordinates <- coord %>% filter(slice==slice_current); coord_ref <- coord_ref %>% filter(slice==slice_current)
#data = data[,which(colnames(data) %in% coordinates$sample_id)]
data = t(data); data = as(data, "sparseMatrix")

# remove slices with less than min. number of cells
if (TRUE){
  min_cells = 20000
  slices    = table(coord$slice)
  slices    = as.numeric(names(slices[slices >= 20000]))
  coord     = coord[coord$slice %in% slices,]
  data      = data[rownames(data) %in% coord$sample_id,]
}

```


```{r}

# slice-specific cross-expression networks vs slice distance
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")
slice_current = unique(coord$slice)
corr = as.data.frame(matrix(data = 0, nrow = choose(ncol(data),2), ncol = length(slice_current)))
colnames(corr) = slice_current

for (i in 1:length(slice_current)){
  
  temp_meta = coord[coord$slice %in% slice_current[i],]
  temp_data = data[rownames(data) %in% temp_meta$sample_id,]
  locations = temp_meta[,c("CCFx","CCFy","CCFz")]
  
  xx = cross_expression_correlation(data = temp_data, locations = locations)
  corr[,i] = xx$correlation
  print(i/length(slice_current))
}

nets = corr
nets = cor(nets, method = "spearman")
nams = which(upper.tri(nets), arr.ind = TRUE)
nets = data.frame(x = nams[,2] - nams[,1], y = upper_tri(nets))

ggplot(nets) + aes(x = x, y = y) +
  geom_point(size = 0.4) + geom_smooth(method = "lm") +
  labs(x = "Distance (in slices, 200µm between slices)", y = "Correlation", subtitle = str_c("Spearman's rho = ", signif(cor(nets$x, nets$y, method = "spearman"), digits = 2))) +
  theme_classic2()
ggsave(filename = "/Users/AmeerSarwar/Downloads/corr_vs_dist.svg", device = "svg", dpi = 600, width = 6, height = 4)

```


```{r}

# determine speed and memory usage at which software runs
# do this for many cells across number of gene pairs
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CrossExpression.R")

# subset data
data = data[1:5000,]
coordinates = coordinates[1:5000, c("pos_x","pos_y")]

# storage matrices
cells = c(50000, 100000, 200000)
genes = c(1000, 5000, 10000, 20000)

time  = as.data.frame(matrix(data = 0, nrow = length(cells), ncol = length(genes)))
rownames(time) = cells; colnames(time) = genes

# time and memory benchmarks
for (i in 1:length(genes)){
  
  # sample genes
  data_temp = data[, sample(1:ncol(data), size = genes[i], replace = TRUE)]
  
  # sample cells
  for (j in 1:length(cells)){
    
    n = nrow(data)
    repeat_times = cells[j] / n
    indices = rep(1:n, each = repeat_times)
    data_temp = data_temp[indices,]
    locations = coordinates[indices,]
    
    # benchmark
    start_time = Sys.time()
    
    # cross-expression
    xx = cross_expression(data = data_temp, locations = locations)
    
    # compute and store time and memory usage
    diff_time = Sys.time() - start_time
    diff_time = str_c(as.numeric(diff_time), " ", attr(diff_time, "units"))
    time[j,i] = diff_time
    
    # clear memory and show counter
    gc()
    print(str_c("genes = ", genes[i], ", cells = ", cells[j]))
  }
}

# save copy
#write.csv(time, file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Analysis/Validation 5/time.csv")

```


```{r}

# investigate software speed
library(tidyverse)
speed <- read.csv(file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Analysis/Validation 5/time.csv")
genes <- c(500, 1000, 2000, 5000, 8000)
colnames(speed) <- c("cells", genes)
speed <- speed[3:nrow(speed),c(1,4:6)]

# pivot data to long format
long_speed <- speed %>% pivot_longer(cols = -cells, names_to = "genes", values_to = "time") %>% as.data.frame()
long_speed <- data.frame(genes = long_speed$genes, cells = long_speed$cells, time = long_speed$time)

# process time and units
times <- str_split(long_speed$time, pattern = " ", simplify = TRUE)
value <- as.numeric(times[,1])
units <- times[,2]

# convert time on the same scale
value[units == "secs"] = value[units == "secs"] / 60

long_speed$time <- value

# convert cells to factors
long_speed$cells <- long_speed$cells / 1000
long_speed$cells <- factor(long_speed$cells)

# convert genes to gene pairs
long_speed$genes <- as.numeric(long_speed$genes)
#long_speed$genes <- choose(long_speed$genes, 2)

# plot graph
ggplot(long_speed) + aes(x = genes, y = time, color = cells) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "Genes", y = "Time (mins)", color = "Cells\n(thousands)") +
  theme_classic() +
  theme(legend.position = c(0.2,0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

ggsave("/Users/AmeerSarwar/Downloads/speed.svg", width = 4, height = 3, dpi = 600)

```


```{r}

# supplementary tables showing MERFISH spatail, scRNA-seq, and snRNA-seq brain regions used to compare co-expression

# single-cell and spatial
single_cell_regions = c("AId-AIv-AIp","TH-AD-AV-AM-IAD-LD","CTXsp-CLA-EP-LA-BLA-BMA-PA","HIP-CA","STR-sAMY","STR-LSX","STR-STRd","STR-STRv","OLF-COA-PAA-NLOT-TR","HIP", "TH-AD-AV-AM-IAD-LD","TH-LGd-IGL-LGv","TH-MD-IMD-PCN-CL","TH-MH-LH-LP","TH-PO-Eth","TH-PVT-PT","TH-RE-RH-CM-IAM-SMT-PR-Xi","CTXsp-CLA-EP-LA-BLA-BMA-PA", "TH-MH-LH-LP","HIP","TH-PO-Eth","HY-MEZ-PVZ-PVR","TH-PVT-PT","OLF-AON-TT-DP-PIR-COA-NLOT-PAA-TR","OLF-MOB-AOB","HY-MEZ-PVZ-PVR","RSP","SSp","SSs-GU-VISC-AIp","STR-STRd","TH-VAL-VPM-VPL-VM","VISa","SS-GU-VISC","TH-PF-SPA-SPFm-VPMpc-VPLpc-RT","TH-VAL-VPM-VPL-VM")
single_cell_regions = unique(unlist(str_split(single_cell_regions, "-")))

spatial = c("AIp","ATN","BLA","BMA","CA","cc","CEA","CNU","COA","cst","CTX","CTXpl","CTXsp","DG","DORpm","EP","EPI","epsc","fa","fiber tracts","fxpo","fxs","hc","HIP","IB","IIn","ILM","LAT","lfbst","LZ","MED","MEZ","mfbc","mfbse","mfbsma","mfsbshy","MTN","OLF","PALd","PALv","PVR","PVZ","RSPagl","RSPd","RSPv","sAMY","SSp-bfd","SSp-ll","SSp-tr","SSp-ul","SSp-un","SSs","st","STRd","VENT","VISa","VISC","VL","VP","VS")
spatial = unique(spatial)

supp1 = data.frame(regions = c(single_cell_regions, spatial), data = rep(c("MERFISH", "scRNA-seq"), times = c(length(spatial), length(single_cell_regions))))
write.csv(supp1, file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Analysis1 - Cross-Expression/Validation 5/Supplementary Table 1.csv", row.names = FALSE)

# single-cell and single-nucleus
nucleus = c("ACA","AUD","ENT","MOp","RSP","VIS","VISp","STRd","CTX")
cell    = c("ACA","AUD","ENT","MOp","RSP","VIS","VISp","STR-STRd","CTXsp-CLA-EP-LA-BLA-BMA-PA")
cell    = unique(unlist(str_split(cell, "-")))

supp2   = data.frame(regions = c(cell, nucleus), data = rep(c("scRNA-seq","snRNA-seq"), times = c(length(cell),length(nucleus))))
write.csv(supp2, file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Analysis1 - Cross-Expression/Validation 5/Supplementary Table 2.csv", row.names = FALSE)

```