```{r}

# MERFISH brain receptor map
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel","ape","RColorBrewer","proxy","Rfast","prabclus")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))

# load dataset
data     <- fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice2_Replicate2_cell_by_gene.csv")
metadata <- fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/MERFISH_receptor_map/Slice2_Replicate2_cell_metadata.csv")
colnames(data)[1] <- "sample_id"; colnames(metadata)[1] <- "sample_id"; metadata <- metadata[,-c("fov","volume","min_x","max_x","min_y","max_y")]
data <- as.data.frame(data); data <- data[which(rowSums(data[,2:ncol(data)]) >= 50),]; metadata <- metadata[which(metadata$sample_id %in% data$sample_id),]
data <- data[,which(!str_detect(colnames(data), "^Blank-"))]; colnames(metadata)[2:3] <- c("pos_x","pos_y")

# plot on brain slice
ggplot() + aes(x=metadata$pos_x, y=metadata$pos_y) + geom_point(size=0.1)

```