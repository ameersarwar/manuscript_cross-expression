```{r}

# MERFISH brain receptor map
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix","grid","spdep","diptest","ggbeeswarm","Signac","metafor","ggforce","anndata","reticulate",
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

```


```{r}

## compare co-expression between scRNA-seq and spatial data for cell segmentation artefacts

# brain regions and corresponding atlas regions are saved for each slice and replicate
matches <- read.csv(str_c("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/region_to_atlas_correspondance/Validation 2/slice",slice,"_replicate",replicate,".csv"))

# remove non-matching (NA) and poorly matched regions
matches <- matches %>% filter(!is.na(corresponding_atlas_regions), match_quality != "Poor"); matches <- matches[,c("data_regions","corresponding_atlas_regions")]

# merge data_regions that form an atlas region and merge atlas regions that form a data region (separation designated by "/")
matches <- matches %>% group_by(corresponding_atlas_regions) %>% summarise(data_regions = paste(data_regions, collapse = "/")) %>% as.data.frame(); matches <- matches[,c("data_regions","corresponding_atlas_regions")]

# complete mapping between data and atlas genes
df   <- zellkonverter::readH5AD("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.ACA.h5ad")
df   <- t(assay(df, "X"))
data <- data[,c(1,which(colnames(data) %in% colnames(df)))]

```


```{r}

# compare merged data vs atlas gene co-expression for significant and non-significant gene pairs
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/Distance Metrics/correlation.R")

results = vector(mode = "list", length = nrow(matches)); names(results) = matches$data_regions

for (i in 1:nrow(matches)){
  
  # spatial data
  target = matches[i,"data_regions"]
  target = str_split(target,"/")[[1]]
  spatial_meta = metadata[which(metadata$CCFparentname %in% target),]
  spatial_data = data[which(data$sample_id %in% spatial_meta$sample_id),]; spatial_data <- spatial_data[,2:ncol(spatial_data)]
  
  # significant genes
  quantile_threshold = 0.05
  genes_total <- celltocell(spatial_data, locations = spatial_meta[,c("pos_x","pos_y")], neighbor = 1, alpha_neighbor = 0.05, alpha_cell = 0, bidirectional = FALSE, show_one_direction = FALSE)
  
  # spatial co-expression
  spatial_data   = apply(spatial_data, 2, function(x) {q <- quantile(x[x > 0], probs = quantile_threshold); ifelse(x > 0 & x < q, 0, x) })
  spatial_co_exp = correlation(spatial_data)
  
  # distribution of significant gene pairs
  mask = matrix(FALSE, nrow = nrow(spatial_co_exp), ncol = ncol(spatial_co_exp)); mask[cbind(genes_total$row.index, genes_total$col.index)] = TRUE
  significant_spatial = spatial_co_exp[mask]; significant_spatial = significant_spatial[!(is.na(significant_spatial) | is.nan(significant_spatial) | is.infinite(significant_spatial))]
  
  # distribution of non-significant gene pairs
  non_significant_spatial = spatial_co_exp[!mask]; non_significant_spatial = non_significant_spatial[!(is.na(non_significant_spatial) | is.nan(non_significant_spatial) | is.infinite(non_significant_spatial))]
  
  
  # atlas data
  target = matches[i,"corresponding_atlas_regions"]
  target = str_split(target,"/")[[1]]
  prefix = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes."
  target = str_c(prefix,target,".h5ad")
  
  # load, concatenate, and re-order to match column order of spatial data
  df = zellkonverter::readH5AD("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.ACA.h5ad")
  df = t(assay(df, "X"))
  df <- df[1,1:ncol(df)]
  for (k in 1:length(target)){dats = zellkonverter::readH5AD(target[k]); dats = t(assay(dats, "X")); df = rbind(df,dats)}
  df = df[,which(colnames(df) %in% colnames(spatial_data))]
  df = df[, colnames(spatial_data)]
  
  # atlas co-expression
  df = apply(df, 2, function(x) {q <- quantile(x[x > 0], probs = quantile_threshold); ifelse(x > 0 & x < q, 0, x) })
  atlas_co_expr = correlation(df)
  
  # distribution of significant gene pairs
  significant_atlas = atlas_co_expr[mask]; significant_atlas = significant_atlas[!(is.na(significant_atlas) | is.nan(significant_atlas) | is.infinite(significant_atlas))]
  
  # distribution of non-significant gene pairs
  non_significant_atlas = atlas_co_expr[!mask]; non_significant_atlas = non_significant_atlas[!(is.na(non_significant_atlas) | is.nan(non_significant_atlas) | is.infinite(non_significant_atlas))]
  
  
  # convert results to dataframe
  significant_spatial <- data.frame(values = significant_spatial)
  significant_spatial <- significant_spatial %>% mutate(class = "significant_spatial", spatial_regions = matches[i,"data_regions"], atlas_regions = matches[i,"corresponding_atlas_regions"])
  
  non_significant_spatial <- data.frame(values = non_significant_spatial)
  non_significant_spatial <- non_significant_spatial %>% mutate(class = "non_significant_spatial", spatial_regions = matches[i,"data_regions"], atlas_regions = matches[i,"corresponding_atlas_regions"])
  
  significant_atlas <- data.frame(values = significant_atlas)
  significant_atlas <- significant_atlas %>% mutate(class = "significant_atlas", spatial_regions = matches[i,"data_regions"], atlas_regions = matches[i,"corresponding_atlas_regions"])
  
  non_significant_atlas <- data.frame(values = non_significant_atlas)
  non_significant_atlas <- non_significant_atlas %>% mutate(class = "non_significant_atlas", spatial_regions = matches[i,"data_regions"], atlas_regions = matches[i,"corresponding_atlas_regions"])
  
  # combine results
  output = rbind(significant_spatial, significant_atlas, non_significant_spatial, non_significant_atlas)
  output$class = factor(output$class, levels = c("significant_spatial", "significant_atlas", "non_significant_spatial", "non_significant_atlas"))
  
  results[[i]] = output
  
  print(i/nrow(matches))
}

```


```{r}

## plot significant and non-significant gene pairs' co-expression across spatial and atlas datasets

for (i in 1:length(results)){
  
  output = results[[i]]
  
  p <- ggplot(output) + aes(x = values, y = class, fill = class) +
  geom_violin(scale = "width", adjust = 1.5) +
  labs(subtitle=str_c("Data (spatial), ", unique(output$spatial_regions), "\nAtlas (non-spatial), ", unique(output$atlas_regions), "\nQuantile, ", quantile_threshold),
       x="Co-expression", y="", fill="") +
  guides(fill = guide_legend(reverse = TRUE)) + theme_bw() +
  theme(axis.text.y = element_blank())
  
  print(p)
  ggsave(filename = str_c("/Users/AmeerSarwar/Downloads/picture_",i,".png"))
  
  print(i/length(results))
}

```


```{r}

# sample cells with gene A and ask if they also contain gene B for spatial and atlas datasets
# cells with gene A whose neighbors contain gene B are colored/shaped differently
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/CellComms.R")
source("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Functions/Distance Metrics/correlation.R")

results_celltocell = vector(mode = "list", length = nrow(matches)); names(results_celltocell) = matches$data_regions

for (i in 1:nrow(matches)){
  
  # spatial data
  target = matches[i,"data_regions"]
  target = str_split(target,"/")[[1]]
  spatial_meta = metadata[which(metadata$CCFparentname %in% target),]
  spatial_data = data[which(data$sample_id %in% spatial_meta$sample_id),]; spatial_data <- spatial_data[,2:ncol(spatial_data)]
  
  # atlas data
  target = matches[i,"corresponding_atlas_regions"]
  target = str_split(target,"/")[[1]]
  prefix = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes."
  target = str_c(prefix,target,".h5ad")
  
  # load, concatenate, and re-order to match column order of spatial data
  df = zellkonverter::readH5AD("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.ACA.h5ad")
  df = t(assay(df, "X"))
  df <- df[1,1:ncol(df)]
  for (k in 1:length(target)){dats = zellkonverter::readH5AD(target[k]); dats = t(assay(dats, "X")); df = rbind(df,dats)}
  df = df[,which(colnames(df) %in% colnames(spatial_data))]
  df = df[, colnames(spatial_data)]
  
  
  # number of co-expressing cells per gene pair for spatial and atlas datasets
  spatial      = fisher_similarity(spatial_data, spatial_data); spatial = spatial$cells
  atlas        = fisher_similarity(as.matrix(df), as.matrix(df)); atlas = atlas$cells
  
  # co-expression p-values using the celltocelletric test
  pval_spatial = matrix(0, nrow = nrow(spatial), ncol = ncol(spatial)); dimnames(pval_spatial) = dimnames(spatial)
  pval_atlas   = matrix(0, nrow = nrow(atlas),   ncol = ncol(atlas));   dimnames(pval_atlas)   = dimnames(atlas)
  
  for (j in 1:nrow(spatial)){
    for (k in 1:ncol(spatial)){
      pval_spatial[j,k] = phyper(q=spatial[j,k] - 1, m=spatial[k,k], n=nrow(spatial_data) - spatial[k,k], k=spatial[j,j], lower.tail = FALSE)
      pval_atlas[j,k]   = phyper(q=atlas[j,k]   - 1, m=atlas[k,k],   n=nrow(df)           - atlas[k,k],   k=atlas[j,j],   lower.tail = FALSE)
    }
  }
  
  pval_spatial = -log10(pval_spatial); pval_spatial[is.infinite(pval_spatial)] = 0
  pval_atlas   = -log10(pval_atlas);   pval_atlas[is.infinite(pval_atlas)]     = 0
  
  # genes with significant cell-cell relations
  genes_total  = celltocell(spatial_data, locations = spatial_meta[,c("pos_x","pos_y")], neighbor = 1, alpha_neighbor = 0.05, alpha_cell = 0, bidirectional = FALSE, show_one_direction = FALSE)
  genes_total$p.value = -log10(genes_total$p.value.adj)
  genes_matrix = matrix(0, nrow = nrow(spatial), ncol = ncol(spatial)); dimnames(genes_matrix) = dimnames(spatial)
  
  for (k in 1:nrow(genes_total)){genes_matrix[genes_total$row.index[k], genes_total$col.index[k]] = genes_total$p.value[k]}
  
  p_values = data.frame(coexp_spatial=as.vector(pval_spatial), coexp_atlas=as.vector(pval_atlas), neighbor_sig=as.vector(genes_matrix))
  p_values$neighbor_sig = ifelse(p_values$neighbor_sig != 0, 1, 0)
  
  results_celltocell[[i]] = p_values
  print(i/nrow(matches))
}

```


```{r}

# plot atlas and spatial co-expression p-values (log), with neighbor significant highlighted
for (i in 1:nrow(matches)){

p_values = results_celltocell[[i]] #%>% filter(coexp_atlas <= -log10(1-0.05))
pvals_sig_neighbor  = p_values %>% filter(neighbor_sig != 0)
pvals_nsig_neighbor = p_values %>% filter(neighbor_sig == 0)

p <- ggplot() + 
  geom_point(data=pvals_nsig_neighbor, aes(x=coexp_atlas, y=coexp_spatial), color="lightblue", alpha=0.5, size=0.8) +
  geom_rug(data  =pvals_nsig_neighbor, aes(x=coexp_atlas, y=coexp_spatial), color="lightblue", alpha=0.3, sides="tr", linewidth=0.2) +
  
  geom_point(data=pvals_sig_neighbor, aes(x=coexp_atlas, y=coexp_spatial), color="coral", alpha=0.8, size=0.8) +
  geom_rug(data  =pvals_sig_neighbor, aes(x=coexp_atlas, y=coexp_spatial), color="coral", alpha=0.3, sides="bl", linewidth=0.2) +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x="Atlas co-expression -log10 p-values", y="Spatial data co-expression -log10 p-values",
       title=str_c("celltocelletric test for co-expression\nData (spatial), ",matches$data_regions[i],"; Atlas (non-spatial), ",matches$corresponding_atlas_regions[i],
       "\nCell-cell relation significant (ORANGE) and non-significant (BLUE)"),
       subtitle=str_c("ORANGE: spatial > atlas: ",sum(pvals_sig_neighbor$coexp_spatial > pvals_sig_neighbor$coexp_atlas),
                      "\nORANGE: spatial < atlas: ",sum(pvals_sig_neighbor$coexp_spatial < pvals_sig_neighbor$coexp_atlas),
                      "\nBLUE: spatial > atlas: ",sum(pvals_nsig_neighbor$coexp_spatial > pvals_nsig_neighbor$coexp_atlas),
                      "\nBLUE: spatial < atlas: ",sum(pvals_nsig_neighbor$coexp_spatial < pvals_nsig_neighbor$coexp_atlas))) + 
  theme_classic() + coord_fixed(ratio = 0.5)

ggsave(str_c("/Users/AmeerSarwar/Downloads/image",i,".png"))
print(p)

print(i/nrow(matches))
}

```



```{r}

# save ratios of the total number of gene pairs with higher co-expression in spatial > atlas to spatial < atlas for different atlas (and spatial) p-values
atlas_only_reference = TRUE # FALSE means that both spatial and atlas are considered as 'reference' (TRUE is more accurate)
ratio = data.frame(class="X",orange=0, blue=0); ratio = ratio[-1,]

for (i in 1:nrow(matches)){
  
# compute and save ratios for number of gene pairs with spatial > atlas to spatial < atlas for different atlas p-values
# significant lack of co-expression (left)
if (atlas_only_reference) {p_values = results_celltocell[[i]] %>% filter(coexp_atlas <= -log10(1-0.05))}
if (!atlas_only_reference){p_values = results_celltocell[[i]] %>% filter(coexp_atlas <= -log10(1-0.05) & coexp_spatial <= -log10(1-0.05))}
pvals_sig_neighbor  = p_values %>% filter(neighbor_sig != 0)
pvals_nsig_neighbor = p_values %>% filter(neighbor_sig == 0)

ratio_temp = data.frame(class="significant lack of co-expression",
                        orange=sum(pvals_sig_neighbor$coexp_spatial > pvals_sig_neighbor$coexp_atlas)/
                          sum(pvals_sig_neighbor$coexp_spatial < pvals_sig_neighbor$coexp_atlas),
                        blue=sum(pvals_nsig_neighbor$coexp_spatial > pvals_nsig_neighbor$coexp_atlas)/
                          sum(pvals_nsig_neighbor$coexp_spatial < pvals_nsig_neighbor$coexp_atlas))

ratio = rbind(ratio,ratio_temp)

# non-significant co-expression (middle)
if (atlas_only_reference) {p_values = results_celltocell[[i]] %>% filter(coexp_atlas > -log10(1-0.05) & coexp_atlas < -log10(0.05))}
if (!atlas_only_reference){p_values = results_celltocell[[i]] %>% filter((coexp_atlas > -log10(1-0.05) & coexp_atlas < -log10(0.05)) & (coexp_spatial > -log10(1-0.05) & coexp_spatial < -log10(0.05)))}
pvals_sig_neighbor  = p_values %>% filter(neighbor_sig != 0)
pvals_nsig_neighbor = p_values %>% filter(neighbor_sig == 0)

ratio_temp = data.frame(class="non-significant co-expression",
                        orange=sum(pvals_sig_neighbor$coexp_spatial > pvals_sig_neighbor$coexp_atlas)/
                          sum(pvals_sig_neighbor$coexp_spatial < pvals_sig_neighbor$coexp_atlas),
                        blue=sum(pvals_nsig_neighbor$coexp_spatial > pvals_nsig_neighbor$coexp_atlas)/
                          sum(pvals_nsig_neighbor$coexp_spatial < pvals_nsig_neighbor$coexp_atlas))

ratio = rbind(ratio,ratio_temp)

# significant co-expression (right)
if (atlas_only_reference) {p_values = results_celltocell[[i]] %>% filter(coexp_atlas >= -log10(0.05))}
if (!atlas_only_reference){p_values = results_celltocell[[i]] %>% filter(coexp_atlas >= -log10(0.05) & coexp_spatial >= -log10(0.05))}
pvals_sig_neighbor  = p_values %>% filter(neighbor_sig != 0)
pvals_nsig_neighbor = p_values %>% filter(neighbor_sig == 0)

ratio_temp = data.frame(class="significant co-expression",
                        orange=sum(pvals_sig_neighbor$coexp_spatial > pvals_sig_neighbor$coexp_atlas)/
                          sum(pvals_sig_neighbor$coexp_spatial < pvals_sig_neighbor$coexp_atlas),
                        blue=sum(pvals_nsig_neighbor$coexp_spatial > pvals_nsig_neighbor$coexp_atlas)/
                          sum(pvals_nsig_neighbor$coexp_spatial < pvals_nsig_neighbor$coexp_atlas))

ratio = rbind(ratio,ratio_temp)
}

# clean up the data
ratio = pivot_longer(ratio, cols = 2:3, values_to = "ratio", names_to = "orange_blue") %>% as.data.frame()
ratio = ratio %>% mutate(orange_blue = recode(orange_blue, "orange" = "significant", "blue" = "non-significant"))
ratio$class = str_replace(ratio$class," co-expression", "")
ratio$class = factor(ratio$class, levels = c("significant", "non-significant", "significant lack of"))

```



```{r}

# boxplots of ratios of spatial > atlas to spatial < atlas for different co-expression levels in atlas
if (atlas_only_reference){yy="Atlas co-expression"}; if (!atlas_only_reference){yy="Atlas and spatial co-expression"}

ggplot(ratio) + aes(y=class, x=ratio, fill=orange_blue) +
  geom_boxplot() + labs(x="Ratio of the number of gene pairs with higher co-expression\nin the spatial data to the number of gene pairs with\nhigher co-expression in the atlas data",
                        y=yy, subtitle=str_c("Spatial to atlas paired brain regions: ",nrow(matches))) +
  scale_fill_manual(values = c("significant"="coral","non-significant"="lightblue")) + theme_bw() + scale_x_log10() +
  guides(fill = guide_legend(title = "Cell-cell\nrelations", labels = c("Non-significant", "Significant"), reverse = TRUE))
ggsave("/Users/AmeerSarwar/Downloads/picture.png")

```


```{r}

# find gene pairs with significant cell-cell relations and lack of significant co-expression and plot their co-expression across spatial and atlas datasets

bursting_results = vector(mode = "list", length = nrow(matches)); names(bursting_results) = matches$data_regions

for (i in 1:nrow(matches)){
  
  # spatial data
  target = matches[i,"data_regions"]
  target = str_split(target,"/")[[1]]
  spatial_meta = metadata[which(metadata$CCFparentname %in% target),]
  spatial_data = data[which(data$sample_id %in% spatial_meta$sample_id),]; spatial_data <- spatial_data[,2:ncol(spatial_data)]
  
  # atlas data
  target = matches[i,"corresponding_atlas_regions"]
  target = str_split(target,"/")[[1]]
  prefix = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes."
  target = str_c(prefix,target,".h5ad")
  
  # load, concatenate, and re-order to match column order of spatial data
  df = read_h5ad("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/single_cell/zeng/AIT21.0.merged.with_multiome.barseq_and_Vizgen_genes.ACA.h5ad")
  df = df$to_df(); df <- df[1,1:ncol(df)]; df <- df[-1,]
  for (k in 1:length(target)){dats = read_h5ad(target[k]); dats = dats$to_df(); df = rbind(df,dats)}
  df = df[,which(colnames(df) %in% colnames(spatial_data))]
  df = df[, colnames(spatial_data)]
  
  
  # number of co-expressing cells per gene pair for spatial and atlas datasets
  spatial      = fisher_similarity(spatial_data, spatial_data); spatial = spatial$cells
  atlas        = fisher_similarity(df, df); atlas = atlas$cells
  
  # co-expression p-values using the celltocelletric test
  pval_spatial = matrix(0, nrow = nrow(spatial), ncol = ncol(spatial)); dimnames(pval_spatial) = dimnames(spatial)
  pval_atlas   = matrix(0, nrow = nrow(atlas),   ncol = ncol(atlas));   dimnames(pval_atlas)   = dimnames(atlas)
  
  for (j in 1:nrow(spatial)){
    for (k in 1:ncol(spatial)){
      pval_spatial[j,k] = phyper(q=spatial[j,k] - 1, m=spatial[k,k], n=nrow(spatial_data) - spatial[k,k], k=spatial[j,j], lower.tail = FALSE)
      pval_atlas[j,k]   = phyper(q=atlas[j,k]   - 1, m=atlas[k,k],   n=nrow(df)           - atlas[k,k],   k=atlas[j,j],   lower.tail = FALSE)
    }
  }
  
  # cell-cell relations gene pairs w/o significant co-expression
  genes_total = celltocell(spatial_data, locations = spatial_meta[,c("pos_x","pos_y")], neighbor = 1, alpha_neighbor = 0.05, alpha_cell = 0.05, bidirectional = FALSE, show_one_direction = FALSE)
  mask = matrix(FALSE, nrow = nrow(pval_spatial), ncol = ncol(pval_spatial), dimnames = dimnames(pval_spatial))
  mask[cbind(genes_total$row.index, genes_total$col.index)] = TRUE
  
  pval_spatial = pval_spatial[mask]
  pval_atlas   = pval_atlas[mask]
  
  bursting_results[[i]] = data.frame(pval_spatial,pval_atlas)
  print(i/nrow(matches))
}

```



```{r}

# plot spatial against atlas co-expression p-values (for cell-cell relations genes w/o co-expression)

for (i in 1:nrow(matches)){
  
  target = bursting_results[[i]]
  target = -log10(target)
  
  # atlas significant co-expression
  outcome = target %>% filter(pval_atlas >= -log10(0.05))
  outcome = outcome[order(outcome$pval_atlas),]
  n = 0.1; total = nrow(outcome)
  outcome = outcome[1:(nrow(outcome) - n*nrow(outcome)),]
  
  p1 = ggplot(outcome) + aes(x=pval_atlas,y=pval_spatial) + geom_point(size=0.5) +
    labs(subtitle=str_c("Gene pairs = ",total,", trimmed (for visualization) = atlas top ",n*100,"%, correlation = ", signif(cor(outcome$pval_spatial,outcome$pval_atlas), digits = 2)),
         x="Atlas co-expression -log10 p-values",y="Spatial co-expression -log10 p-values",
         title=str_c("Significant co-expression in atlas\nSpatial, ",matches$data_regions[i],"; Atlas, ",matches$corresponding_atlas_regions[i])) +
    geom_rug(sides="bl") + theme_bw()
  ggsave(str_c("/Users/AmeerSarwar/Downloads/picture",i,".png"))
  
  # atlas non-significant co-expression
  outcome = target %>% filter(pval_atlas < -log10(0.05))
  
  p2 = ggplot(outcome) + aes(x=pval_atlas,y=pval_spatial) + geom_point(size=0.5) +
    labs(subtitle=str_c("Gene pairs = ",nrow(outcome),", correlation = ", signif(cor(outcome$pval_spatial,outcome$pval_atlas), digits = 2)),
         x="Atlas co-expression -log10 p-values",y="Spatial co-expression -log10 p-values",
         title=str_c("Non-significant co-expression in atlas\nSpatial, ",matches$data_regions[i],"; Atlas, ",matches$corresponding_atlas_regions[i])) +
    geom_rug(sides="bl") + theme_bw()
  ggsave(str_c("/Users/AmeerSarwar/Downloads/picture",i,".png"))
  
  print(p2); print(p1)
  print(i/nrow(matches))
}

```
