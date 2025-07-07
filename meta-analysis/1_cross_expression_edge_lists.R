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

# source functions
source("/inkwell05/ameer/functions/0_source_functions.R")

```


```{r}

# cross-expression profiles

dirr  = "/inkwell05/ameer/databases/spatial/"
dirr  = str_c(dirr, list.files(dirr))
genes = as.data.frame(fread("/inkwell05/ameer/databases/scRNAseq/Zeng_2023_Mouse_10x/metadata/gene_master_list.csv"))

for (i in 1:length(dirr)){
  
  # extract datasets
  temp_dirr = dirr[i]
  temp_dirr = str_c(temp_dirr, "/", list.files(temp_dirr))
  
  metadata  = temp_dirr[str_detect(temp_dirr, "metadata")]
  metadata  = str_c(metadata, "/", list.files(metadata))
  
  df = temp_dirr[str_detect(temp_dirr, "expression")]
  df = str_c(df, "/", list.files(df))
  
  # cross-expression profile per slice
  for (j in 1:length(df)){
    
    # load datasets and perform QC
    temp_meta = spatial_QC(path_to_metadata = metadata[j], path_to_expression = df[j])$metadata
    temp_data = spatial_QC(path_to_metadata = metadata[j], path_to_expression = df[j])$data
    
    # cross-expression
    data = data.frame(cross_expression(data = temp_data, locations = temp_meta[,c("x","y")]),
                      correlation = cross_expression_correlation(data = temp_data, locations = temp_meta[,c("x","y")])$correlation)[,c(-1,-2)]
    
    # append gene information: matches for genes
    gene1_hits1 = match(data$gene1, genes$gene_identifier) # ensembl id hits
    gene1_hits2 = match(data$gene1, genes$gene_symbol)     # gene name hits
  
    gene2_hits1 = match(data$gene2, genes$gene_identifier)
    gene2_hits2 = match(data$gene2, genes$gene_symbol)

    if (!(sum(is.na(gene1_hits1)) == nrow(data))){
      gene1_ensembl = genes$gene_identifier[gene1_hits1]
      gene1_symbol  = genes$gene_symbol[gene1_hits1]
      gene1_description = genes$name[gene1_hits1]
    
      gene2_ensembl = genes$gene_identifier[gene2_hits1]
      gene2_symbol  = genes$gene_symbol[gene2_hits1]
      gene2_description = genes$name[gene2_hits1]
    }
  
    if (!(sum(is.na(gene1_hits2)) == nrow(data))){
      gene1_ensembl = genes$gene_identifier[gene1_hits2]
      gene1_symbol  = genes$gene_symbol[gene1_hits2]
      gene1_description = genes$name[gene1_hits2]
    
      gene2_ensembl = genes$gene_identifier[gene2_hits2]
      gene2_symbol  = genes$gene_symbol[gene2_hits2]
      gene2_description = genes$name[gene2_hits2]
    }
  
    # combine gene information with cross-expression
    data = data.frame(gene1_ensembl, gene2_ensembl, gene1_symbol, gene2_symbol, gene1_description, gene2_description, data[,c(-1,-2)])

    # append name
    nam = str_remove(metadata[j], dirr[i])
    nam = str_remove(nam, "/metadata/")
    nam = str_remove(nam, "_metadata.csv")
    
    # save results
    fwrite(data, str_c("/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/", nam, "_cross_expression.csv"))
    
    print(str_c("i = ", signif(i/length(dirr), digits = 2), "; j = ", signif(j/length(df), digits = 2) ))
  }
}

```


```{r}

# cross-expression reproducibility

dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

datasets  = list.files("/inkwell05/ameer/databases/spatial/")
sig_genes = data.frame()
auroc_result = data.frame()

for (i in 1:length(datasets)){
  
  # dataset 1
  temp_dirr1 = dirr[str_detect(dirr, datasets[i])]
  data1 = vector(mode = "list", length = length(temp_dirr1))
  nam   = str_remove(temp_dirr1, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
  names(data1) = str_remove(nam, "_cross_expression.csv")
  for (k in 1:length(temp_dirr1)){data1[[k]] = as.data.frame(fread(temp_dirr1[k]))}
  
  # create gene pairs data 1
  df1 = data1[[1]]
  genes_df1 = df1[,c("gene1_ensembl","gene2_ensembl")]
  genes_df1_posX = str_c(genes_df1$gene1_ensembl, genes_df1$gene2_ensembl)
  genes_df1_posY = str_c(genes_df1$gene2_ensembl, genes_df1$gene1_ensembl)
  
  for (j in 1:length(datasets)){
    
    # dataset 2
    temp_dirr2 = dirr[str_detect(dirr, datasets[j])]
    data2 = vector(mode = "list", length = length(temp_dirr2))
    nam   = str_remove(temp_dirr2, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
    names(data2) = str_remove(nam, "_cross_expression.csv")
    for (k in 1:length(temp_dirr2)){data2[[k]] = as.data.frame(fread(temp_dirr2[k]))}

    # create gene pairs data 2
    df2 = data2[[1]]
    genes_df2 = df2[,c("gene1_ensembl","gene2_ensembl")]
    genes_df2_posX = str_c(genes_df2$gene1_ensembl, genes_df2$gene2_ensembl)
    genes_df2_posY = str_c(genes_df2$gene2_ensembl, genes_df2$gene1_ensembl)
    
    # subset to common genes
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
    
    # initialize AUROC
    result1 = matrix(data = 0, nrow = nrow(genes_df1), ncol = length(temp_dirr1)); colnames(result1) = names(data1) # correlation
    result2 = matrix(data = 0, nrow = nrow(genes_df1), ncol = length(temp_dirr1)); colnames(result2) = names(data1) # p-value
    result3 = matrix(data = 0, nrow = nrow(genes_df2), ncol = length(temp_dirr2)); colnames(result3) = names(data2) # significance (labels)
    
    for (k in 1:length(temp_dirr1)){result1[,k] = data1[[k]][as.numeric(rownames(genes_df1)), "correlation"]}
    for (k in 1:length(temp_dirr1)){result2[,k] = data1[[k]][as.numeric(rownames(genes_df1)), "cross_padj"]}
    for (k in 1:length(temp_dirr2)){result3[,k] = data2[[k]][as.numeric(rownames(genes_df2)), "cross_sig"]}
    
    # AUROC using correlation and p-value
    corr_auroc = auroc(scores = result1, labels = result3)
    pval_auroc = auroc(scores = result2, labels = result3)
    
    result_auroc = data.frame(data1_scores  = datasets[i],
                              data2_labels  = datasets[j],
                              slice1_scores = rep(colnames(corr_auroc), each  = nrow(corr_auroc)),
                              slice2_labels = rep(rownames(corr_auroc), times = ncol(corr_auroc)),
                              auroc_corr = as.vector(corr_auroc),
                              auroc_pval = as.vector(pval_auroc))
    
    result_auroc$slice1_scores = str_remove(result_auroc$slice1, str_c(result_auroc$data1, "_"))
    result_auroc$slice2_labels = str_remove(result_auroc$slice2, str_c(result_auroc$data2, "_"))
    
    # append AUROC results together
    auroc_result = rbind(auroc_result, result_auroc)

    # significant genes in dataset 2 (used as prediction labels)
    hits = rowSums(result3) > 0
    result3 = data2[[k]][as.numeric(rownames(genes_df2)), ]
    result3 = data.frame(data1_scores = datasets[i], data2_labels = datasets[j], result3[hits, 1:6])
    
    # append significant genes together
    sig_genes = rbind(sig_genes, result3)

    print(str_c("i = ", signif(i/length(datasets), digits = 2) , "; j = ", signif(j/length(datasets), digits = 2) ))
  }
}

# save auroc scores and significant genes
fwrite(auroc_result, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_replicability_across_samples_and_datasets.csv")

sig_genes = sig_genes[order(match(sig_genes$data, datasets)), ]
fwrite(sig_genes, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_significant_genes_across_samples_and_datasets.csv")

```


```{r}

# !! DO NOT RUN THIS CODE BLOCK !!
# !! ALREADY RAN AND OVERWROTE THE RESULTS !!

# append info. about gene names, ensembl ID's, and gene description from master list

dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

for (i in 1:length(dirr)){

  genes = as.data.frame(fread("/inkwell05/ameer/databases/scRNAseq/Zeng_2023_Mouse_10x/metadata/gene_master_list.csv"))  
  data  = as.data.frame(fread(dirr[i]))
  
  # matches for genes
  gene1_hits1 = match(data$gene1, genes$gene_identifier) # ensembl id hits
  gene1_hits2 = match(data$gene1, genes$gene_symbol)     # gene name hits
  
  gene2_hits1 = match(data$gene2, genes$gene_identifier)
  gene2_hits2 = match(data$gene2, genes$gene_symbol)

  if (!(sum(is.na(gene1_hits1)) == nrow(data))){
    gene1_ensembl = genes$gene_identifier[gene1_hits1]
    gene1_symbol  = genes$gene_symbol[gene1_hits1]
    gene1_description = genes$name[gene1_hits1]
    
    gene2_ensembl = genes$gene_identifier[gene2_hits1]
    gene2_symbol  = genes$gene_symbol[gene2_hits1]
    gene2_description = genes$name[gene2_hits1]
  }
  
  if (!(sum(is.na(gene1_hits2)) == nrow(data))){
    gene1_ensembl = genes$gene_identifier[gene1_hits2]
    gene1_symbol  = genes$gene_symbol[gene1_hits2]
    gene1_description = genes$name[gene1_hits2]
    
    gene2_ensembl = genes$gene_identifier[gene2_hits2]
    gene2_symbol  = genes$gene_symbol[gene2_hits2]
    gene2_description = genes$name[gene2_hits2]
  }
  
  # combine gene information
  info_gene = data.frame(gene1_ensembl, gene2_ensembl, gene1_symbol, gene2_symbol, gene1_description, gene2_description)

  # append gene information to data
  data = data.frame(info_gene, data[,c(-1,-2)])
  
  # save datasets
  nam = str_replace_all(dirr[i], "cross_expression_edge_lists", "edge_list_temp")
  fwrite(data, nam)

  print(i)
}

```
