# load libraries
library(tidyverse)
library(data.table)
library(Matrix)
library(Rfast)
library(matrixStats)
library(ggridges)
library(reticulate)
library(anndata)
library(gtools)

# source functions
source("/inkwell05/ameer/functions/0_source_functions.R")

# we will predict the significance (alpha = 0.05) of one slice using p-values or correlation of another slice
# cross-expression reproducibility

dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

datasets  = list.files("/inkwell05/ameer/databases/spatial/")
sig_genes = data.frame()
auroc_result = data.frame()

for (i in 1:length(datasets)){
  
  # dataset 1
  temp_dirr1 = mixedsort(dirr[str_detect(dirr, datasets[i])])
  
  data1 = vector(mode = "list", length = length(temp_dirr1))
  nam   = str_remove(temp_dirr1, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
  names(data1) = str_remove(nam, "_cross_expression.csv")
  for (k in 1:length(temp_dirr1)){data1[[k]] = as.data.frame(fread(temp_dirr1[k]))}
  
  # storage variables
  temp_auroc = data.frame()
  temp_genes = data.frame()
  
  for (j in 1:length(datasets)){
    
    # dataset 2
    temp_dirr2 = mixedsort(dirr[str_detect(dirr, datasets[j])])
    
    data2 = vector(mode = "list", length = length(temp_dirr2))
    nam   = str_remove(temp_dirr2, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
    names(data2) = str_remove(nam, "_cross_expression.csv")
    for (k in 1:length(temp_dirr2)){data2[[k]] = as.data.frame(fread(temp_dirr2[k]))}
    
    # create gene pairs data 1
    df1 = data1[[1]]
    genes_df1 = df1[,c("gene1_ensembl","gene2_ensembl")]
    genes_df1_posX = str_c(genes_df1$gene1_ensembl, genes_df1$gene2_ensembl)
    genes_df1_posY = str_c(genes_df1$gene2_ensembl, genes_df1$gene1_ensembl)
    
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
    
    result_auroc = data.frame(data1  = datasets[i],
                              data2  = datasets[j],
                              slice1 = rep(colnames(corr_auroc), each  = nrow(corr_auroc)),
                              slice2 = rep(rownames(corr_auroc), times = ncol(corr_auroc)),
                              auroc_corr = as.vector(corr_auroc),
                              auroc_pval = as.vector(pval_auroc))
    
    result_auroc$slice1 = str_remove(result_auroc$slice1, str_c(result_auroc$data1, "_"))
    result_auroc$slice2 = str_remove(result_auroc$slice2, str_c(result_auroc$data2, "_"))
    
    # append AUROC results together
    temp_auroc = rbind(temp_auroc, result_auroc)
    
    # significant genes in dataset 2 (used as prediction labels)
    hits = rowSums(result3) > 0
    if (sum(hits) == 0){next} # skip if no significant genes are found
    result3 = data2[[k]][as.numeric(rownames(genes_df2)), ]
    result3 = data.frame(data1 = datasets[i], data2 = datasets[j], result3[hits, 1:6])
    
    # append significant genes together
    temp_genes = rbind(temp_genes, result3)
    
    print(str_c("i = ", signif(i/length(datasets), digits = 2) , "; j = ", signif(j/length(datasets), digits = 2) ))
  }
  
  # collate datasets
  auroc_result = rbind(auroc_result, temp_auroc)
  sig_genes    = rbind(sig_genes, temp_genes)
}

# save auroc scores and significant genes
fwrite(auroc_result, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_replicability_across_samples_and_datasets.csv")
fwrite(sig_genes, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_significant_genes_across_samples_and_datasets.csv")

print("Successful completion")