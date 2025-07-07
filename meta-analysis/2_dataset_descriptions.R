```{r}

# load libraries
library(tidyverse)
library(data.table)
library(Matrix)
library(Rfast)
library(matrixStats)
library(ggridges)
library(reticulate)
library(anndata)
library(scales)
library(ComplexHeatmap)
library(forcats)
library(igraph)
library(mclust)
library(future.apply)
library(UpSetR)
library(gtools)

# source functions
source("/inkwell05/ameer/functions/0_source_functions.R")

```


```{r}

# total number of cells and genes in each dataset
data = data.frame(data = names(df), cells = 0, genes = 0)

for (i in 1:length(df)){
  data[i,"cells"] = nrow(df[[i]])
  data[i,"genes"] = ncol(df[[i]])
}

# save data
#fwrite(data, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/total_cells_genes.csv")

```


```{r}

# extract gene panels and number of samples per dataset
dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

nam  = str_remove(dirr, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
nam  = str_remove(nam, ".csv")
nam  = str_remove(nam, "_Mouse")

genes = vector(mode = "list", length = length(dirr))
names(genes) = nam

for (i in 1:length(genes)){
  
  data = as.data.frame(fread(dirr[i]))
  genes[[i]] = unique(c(data$gene1_ensembl, data$gene2_ensembl))
  
  print(i)
}

# overlap in gene panels
overlap = matrix(data = 0, nrow = length(genes), ncol = length(genes), dimnames = list(names(genes), names(genes)))
for (i in 1:length(genes)){for (j in 1:length(genes)){overlap[i,j] = length(intersect(genes[[i]], genes[[j]]))}}

# save dataset
#fwrite(as.data.frame(overlap), "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_gene_panel_overlap.csv")

```


```{r}

# cross-expression profile similarity between slices within the datasets
# we aim to remove slices whose cross-expression profiles diverge from those of other slices in the same dataset
# the results are saved and subsequent analyses (scripts) ignore slices identified in the saved results
dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

targ = list.files("/inkwell05/ameer/databases/spatial/")
to_remove = c(NA,0,0,1,NA,NA,NA,0,1,3,4,1,1) # worst NN slices to remove per dataset (NA means 1 slice in dataset)

df_total  = data.frame()

for (i in 1:length(targ)){
  
  # skip when data contains 1 slice or slices are highly similar to each other
  temp_dirr = dirr[str_detect(dirr, targ[i])]
  temp_dirr = mixedsort(temp_dirr)
  if (length(temp_dirr)==1){next}
  
  # extract slice-specific profiles
  profiles  = matrix(data = 0, ncol = length(temp_dirr), nrow = length(as.data.frame(fread(temp_dirr[1]))$correlation))
  nam = str_remove(temp_dirr, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
  nam = str_remove(nam, "_cross_expression.csv")
  colnames(profiles) = nam; rownames(profiles) = 1:nrow(profiles)
  
  for (j in 1:length(temp_dirr)){profiles[,j] = as.data.frame(fread(temp_dirr[j]))$correlation; print(j/length(temp_dirr))}
  
  # correlation b/w cross-expression profiles
  profiles = cor(profiles)
  
  # find slices whose profiles are very different from those of other slices
  y = sort(colMedians(profiles))
  x = scale_01(rank(y))

  # plot median correlation of each slice against others vs the ranks of correlations
  df = data.frame(data = targ[i], slice = names(y), x = x, y = y); rownames(df) = 1:nrow(df)
  df$retained = c(rep("No", to_remove[i]), rep("Yes", nrow(df) - to_remove[i])) # whether slices meet threshold for retention
  
  # store data
  df_total = rbind(df_total, df)
}

# save dataset
#fwrite(df_total, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_slices_with_divergent_cross-expression.csv")

# make plot
df_total$data = str_remove(df_total$data, "_Mouse")

shape_map = c("10x_AD-controls_2023_Xenium" = 0, "10x_replicates4_2023_Xenium" = 1, "Gillis_Unpublished_BARseq" = 2, "Vizgen_2022_MERSCOPE" = 5,
              "Zador_2024_BARseq" = 6, "Zeng_2023_MERSCOPE" = 3, "Zhuang_2023_MERFISH" = 4, "Zhuang_2024_MERFISH" = 8, "Wang_2023_STARmap" = 9)

ggplot(df_total) + aes(x = x, y = y, color = retained, shape = data) + geom_point(size = 5) +
  scale_y_continuous(limits = c(0, 1)) + theme_classic() + scale_shape_manual(values = shape_map) +
  labs(x = "Normalized rank of average similarity (median correlation) between\nslices' cross-expression profiles",
       y = "Average similarity (median correlation) between\nslices' cross-expression profiles",
       color = "Slice retained", shape = "Data (≥ 2 slices)") +
  theme(legend.position = "right", axis.text = element_text(size = 20), axis.title = element_text(size = 25),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25))

```


```{r}

# we will compute unique and total number of significant edges as a function of number of datasets used

# load datasets
dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

nam  = str_remove(dirr, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
nam  = str_remove(nam, ".csv")
nam  = str_remove(nam, "_Mouse")

data = vector(mode = "list", length = length(dirr))
names(data) = nam

samples = c()

for (i in 1:length(data)){
  
  data[[i]] = as.data.frame(fread(dirr[i]))
  samples = c(samples, length(unique(data[[i]]$sample_id)))
  
  print(i)
}

```


```{r}

# no. of unique significant connections vs no. of slices used
# no. of total significant connections vs no. of slices used
output = vector(mode = "list", length = length(data))
names(output) = nam

for (i in 1:length(data)){
  
  # process each dataset
  df = data[[i]]
  
  # initialize variables for sampling
  iter = 100
  result_total = as.data.frame(matrix(data = 0, nrow = samples[i], ncol = iter + 1))
  colnames(result_total) = c("samples", str_c("total_iter_", 1:iter))
  result_total$samples = 1:samples[i]
  
  result_unique = result_total
  colnames(result_unique) = c("samples", str_c("unique_iter_", 1:iter))
  
  for (j in 1:samples[i]){
    
    # total and unique cross-expression links averaged over `iter` iterations
    for (k in 1:iter){
      
      # select random sample
      id_sample = unique(df$sample_id)
      id_sample = sample(x = id_sample, size = result_total$samples[j], replace = TRUE)
      
      # total cross-expression links
      result_total[j,k+1] = sum(df[df$sample_id %in% id_sample,]$cross_sig)
      
      # unique cross-expression links
      x = df[df$sample_id %in% id_sample,c("gene1_ensembl","gene2_ensembl","cross_sig")]
      x = x[x$cross_sig == 1,]
      result_unique[j,k+1] = nrow(unique(x))
    }
    print(str_c("i = ", i, "; j = ", signif(j/samples[i], digits = 2)))
  }
  
  # append dataset results to output list
  result = cbind(result_total, result_unique[,-c(iter + 2)])
  output[[i]] = result
}

# attach dataset info to each matrix in the list
result = output

outcome = data.frame()

for (i in 1:length(result)){
  
  temp_result = result[[i]]
  temp_result = data.frame(data = nam[i], sample_id = temp_result$samples, temp_result[,2:ncol(temp_result)])
  outcome = rbind(outcome, temp_result)
}

# save results
#fwrite(outcome, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_total_unique_cross_expressing_genes.csv")

```


```{r}

# cross-expression reproducibility within and between datasets
output = data.frame()

for (i in 1:length(data)){
  
  for (j in 1:length(data)){
    
    # extract data
    df1 = data[[i]]
    df2 = data[[j]]
    
    # subset to common gene pairs
    genes_df1 = df1[df1$sample_id == unique(df1$sample_id)[1],c("gene1_ensembl","gene2_ensembl")]
    genes_df2 = df2[df2$sample_id == unique(df2$sample_id)[1],c("gene1_ensembl","gene2_ensembl")]
  
    genes_df1_posX = str_c(genes_df1$gene1_ensembl, genes_df1$gene2_ensembl)
    genes_df1_posY = str_c(genes_df1$gene2_ensembl, genes_df1$gene1_ensembl)
    
    genes_df2_posX = str_c(genes_df2$gene1_ensembl, genes_df2$gene2_ensembl)
    genes_df2_posY = str_c(genes_df2$gene2_ensembl, genes_df2$gene1_ensembl)
    
    posN = c(intersect(genes_df1_posX, genes_df2_posX),
             intersect(genes_df1_posX, genes_df2_posY),
             intersect(genes_df1_posY, genes_df2_posX),
             intersect(genes_df1_posY, genes_df2_posY))
    
    genes_df1 = genes_df1[(genes_df1_posX %in% posN) | (genes_df1_posY %in% posN),]
    genes_df2 = genes_df2[(genes_df2_posX %in% posN) | (genes_df2_posY %in% posN),]

    # re-order them similarly
    genes_df1_pairs = str_c(pmin(genes_df1$gene1_ensembl, genes_df1$gene2_ensembl), 
                             pmax(genes_df1$gene1_ensembl, genes_df1$gene2_ensembl))

    genes_df2_pairs = str_c(pmin(genes_df2$gene1_ensembl, genes_df2$gene2_ensembl), 
                             pmax(genes_df2$gene1_ensembl, genes_df2$gene2_ensembl))

    matching_indices = match(genes_df1_pairs, genes_df2_pairs)
    genes_df2 = genes_df2[matching_indices, ]
    
    # apply common and ordered references to each sample in each dataset
    result_df1 = data.frame()
    
    for (m in 1:length(unique(df1$sample_id))){
      
      temp_df = df1[df1$sample_id == unique(df1$sample_id)[m],]
      
      genes1 = str_c(pmin(genes_df1$gene1_ensembl, genes_df1$gene2_ensembl),
                     pmax(genes_df1$gene1_ensembl, genes_df1$gene2_ensembl))
      
      genes2 = str_c(pmin(temp_df$gene1_ensembl, temp_df$gene2_ensembl),
                     pmax(temp_df$gene1_ensembl, temp_df$gene2_ensembl))

      temp_df = temp_df[genes2 %in% genes1, ]
      
      genes2 = str_c(pmin(temp_df$gene1_ensembl, temp_df$gene2_ensembl),
                     pmax(temp_df$gene1_ensembl, temp_df$gene2_ensembl))
      
      matching_indices = match(genes1, genes2)
      temp_df = temp_df[matching_indices, ]
      result_df1 = rbind(result_df1, temp_df)
    }
    
    
    result_df2 = data.frame()
    
    for (m in 1:length(unique(df2$sample_id))){
      
      temp_df = df2[df2$sample_id == unique(df2$sample_id)[m],]
      
      genes1 = str_c(pmin(genes_df2$gene1_ensembl, genes_df2$gene2_ensembl),
                     pmax(genes_df2$gene1_ensembl, genes_df2$gene2_ensembl))
      
      genes2 = str_c(pmin(temp_df$gene1_ensembl, temp_df$gene2_ensembl),
                     pmax(temp_df$gene1_ensembl, temp_df$gene2_ensembl))

      temp_df = temp_df[genes2 %in% genes1, ]
      
      genes2 = str_c(pmin(temp_df$gene1_ensembl, temp_df$gene2_ensembl),
                     pmax(temp_df$gene1_ensembl, temp_df$gene2_ensembl))
      
      matching_indices = match(genes1, genes2)
      temp_df = temp_df[matching_indices, ]
      result_df2 = rbind(result_df2, temp_df)
    }
    
    df1 = result_df1
    df2 = result_df2
   
    # extract sample-specific information
    sample_id = unique(df1$sample_id)
    corr_df1  = matrix(data = 0, nrow = as.numeric(table(df1$sample_id)[1]), ncol = length(sample_id))
    colnames(corr_df1) = sample_id
    
    roc_labels1_df1 = corr_df1
    roc_labels2_df1 = corr_df1
    
    for (k in 1:length(sample_id)){
    
      corr_df1[,k]        = df1[df1$sample_id %in% sample_id[k],]$cross_correlation
      roc_labels1_df1[,k] = df1[df1$sample_id %in% sample_id[k],]$cross_sig
      roc_labels2_df1[,k] = as.integer(df1[df1$sample_id %in% sample_id[k],]$cross_pvalue <= 0.05)
    }
    
    roc_values_df1 = corr_df1
    
    sample_id = unique(df2$sample_id)
    corr_df2  = matrix(data = 0, nrow = as.numeric(table(df2$sample_id)[1]), ncol = length(sample_id))
    colnames(corr_df2) = sample_id
    
    roc_labels1_df2 = corr_df2
    roc_labels2_df2 = corr_df2
    
    for (k in 1:length(sample_id)){
    
      corr_df2[,k]        = df2[df2$sample_id %in% sample_id[k],]$cross_correlation
      roc_labels1_df2[,k] = df2[df2$sample_id %in% sample_id[k],]$cross_sig
      roc_labels2_df2[,k] = as.integer(df2[df2$sample_id %in% sample_id[k],]$cross_pvalue <= 0.05)
    }
    
    roc_values_df2 = corr_df2

    # compute reproducibility statistics
    correlations  = cor(corr_df1, corr_df2)
    
    auroc_postFDR_1 = auroc(labels = roc_labels1_df2, scores = roc_values_df1)
    auroc_preFDR_1  = auroc(labels = roc_labels2_df2, scores = roc_values_df1)
    
    auroc_postFDR_2 = auroc(labels = roc_labels1_df1, scores = roc_values_df2)
    auroc_preFDR_2  = auroc(labels = roc_labels2_df1, scores = roc_values_df2)
    
    # combine correlation and AUROC with data and sample ID's 
    result = data.frame(data1 = names(data)[i],
                        data2 = names(data)[j],
                        sample1 = rep(rownames(correlations), times = ncol(correlations)),
                        sample2 = rep(colnames(correlations), each  = nrow(correlations)),
                        correlations = as.vector(correlations),
                        auroc_postFDR_1 = as.vector(auroc_postFDR_1),
                        auroc_preFDR_1  = as.vector(auroc_preFDR_1),
                        auroc_postFDR_2 = as.vector(auroc_postFDR_2),
                        auroc_preFDR_2  = as.vector(auroc_preFDR_2))

  # append results together
  output = rbind(output, result)
  
  print(str_c("i = ", i, "; j = ", j))
  gc()
  }
}

# save results
#fwrite(output, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_replicability_across_samples_and_datasets.csv")

```


```{r}

# !! Results already saved !! Do NOT run this code block !!

# save cross-expression matrices per slice
for (i in 1:length(data)){
  
  temp_data = data[[i]]
  
  
  
  for (j in 1:length(unique(temp_data$sample_id))){
    
    temp_sample = temp_data[temp_data$sample_id %in% unique(temp_data$sample_id)[j],]
    fwrite(temp_sample, str_c("/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/slice_level/",
                              names(data)[i], "_sample_", unique(temp_data$sample_id)[j], ".csv"))
    
    print(str_c("i = ", i, "; j = ", signif(j/length(unique(temp_data$sample_id)), digits = 2)))
    gc()
  }
}

```