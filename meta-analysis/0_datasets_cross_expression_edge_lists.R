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

# our purpose is to compute and save cross-expression profiles (p-values and correlations) per slice per dataset

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

print("Successful completion")