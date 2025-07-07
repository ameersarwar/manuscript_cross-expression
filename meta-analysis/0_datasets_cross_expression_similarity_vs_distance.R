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
library(patchwork)
library(circlize)
library(imager)

# source functions
source("/inkwell05/ameer/functions/0_source_functions.R")

# we will compute the similarity (Jacquard, raw overlap, and correlation) between slices' cross-expression profiles as a function of their distance

target = c("Gillis_Unpublished_Mouse_BARseq", "Vizgen_2022_Mouse_MERSCOPE", "Zador_2024_Mouse_BARseq", "Zeng_2023_Mouse_MERSCOPE", "Zhuang_2023_Mouse_MERFISH")
dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

outcome = data.frame()

for (i in 1:length(target)){
  
  # load datasets
  to_load  = mixedsort(dirr[str_detect(dirr, target[i])])
  nam = str_remove(to_load, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
  nam = str_remove(nam, "_cross_expression.csv")
  profiles_sig  = matrix(data = 0, ncol = length(to_load), nrow = nrow(as.data.frame(fread(to_load[1]))))  
  dimnames(profiles_sig) = list(1:nrow(profiles_sig), nam)
  profiles_corr = profiles_sig
  
  for (j in 1:length(to_load)){
    profiles_sig[,j]  = as.data.frame(fread(to_load[j]))$cross_sig
    profiles_corr[,j] = as.data.frame(fread(to_load[j]))$correlation
    print(j/length(to_load))
  }
  
  # process slices
  # note: Gillis_Unpublished_Mouse_BARseq and Zeng_2023_Mouse_MERSCOPE are already ordered
  # note: Vizgen_2022_Mouse_MERSCOPE needs slice ordering between datasets
  # note: Zador_2024_Mouse_BARseq and Zhuang_2023_Mouse_MERFISH require separating slices into brains (already ordered)
  
  # Vizgen_2022_Mouse_MERSCOPE (order slices across brains by CCF coordinates)
  if (target[i] == "Vizgen_2022_Mouse_MERSCOPE"){
    
    dirr_temp = "/inkwell05/ameer/databases/spatial/"
    dirr_temp = str_c(dirr_temp, list.files(dirr_temp), "/")
    to_load   = dirr_temp[str_detect(dirr_temp, target[i])]
    to_load   = str_c(to_load, list.files(to_load))
    to_load   = str_c(to_load[str_detect(to_load, "metadata")], "/")
    to_load   = mixedsort(str_c(to_load, list.files(to_load)))
    
    to_sort = c()
    for (q in 1:length(to_load)){to_sort = c(to_sort, median(as.data.frame(fread(to_load[q]))[,"CCFx"], na.rm = TRUE))}
    
    profiles_sig  = profiles_sig[, order(to_sort)]
    profiles_corr = profiles_corr[, order(to_sort)]
  }
  
  # Zador_2024_Mouse_BARseq (separate dataset into brains, which are already ordered by CCF coordinates)
  if (target[i] == "Zador_2024_Mouse_BARseq"){
    
    to_extract = c("control1","control2","control3","control4","enucleated1","enucleated2","enucleated3","enucleated4","pilot")
    results = data.frame()
    
    for (q in 1:length(to_extract)){
      
      xx = profiles_sig[, str_detect(colnames(profiles_sig), to_extract[q])]
      
      distance    = abs(outer(1:ncol(xx), 1:ncol(xx), "-"))
      jaccard_idx = jaccard(xx)
      raw_totals  = t(xx) %*% xx
      corr        = cor(profiles_corr[, str_detect(colnames(profiles_corr), to_extract[q])])
      
      diag(distance) = NA; diag(jaccard_idx) = NA; diag(raw_totals) = NA; diag(corr) = NA
      xx = data.frame(brain = str_c(target[i], "_", to_extract[q]), distance = as.vector(distance), jaccard = as.vector(jaccard_idx), same_genes = as.vector(raw_totals), correlation = as.vector(corr))
      results = rbind(results, xx)
    }
  }
  
  # Zhuang_2023_Mouse_MERFISH (separate dataset into brains, which are already ordered by CCF coordinates)
  if (target[i] == "Zhuang_2023_Mouse_MERFISH"){
    
    to_extract = c("sagittal_3","sagittal_23","coronal_66","coronal_147")
    results = data.frame()
    
    for (q in 1:length(to_extract)){
      
      xx = profiles_sig[, str_detect(colnames(profiles_sig), to_extract[q])]
      
      distance    = abs(outer(1:ncol(xx), 1:ncol(xx), "-"))
      jaccard_idx = jaccard(xx)
      raw_totals  = t(xx) %*% xx
      corr        = cor(profiles_corr[, str_detect(colnames(profiles_corr), to_extract[q])])
      
      diag(distance) = NA; diag(jaccard_idx) = NA; diag(raw_totals) = NA; diag(corr) = NA
      xx = data.frame(brain = str_c(target[i], "_", to_extract[q]), distance = as.vector(distance), jaccard = as.vector(jaccard_idx), same_genes = as.vector(raw_totals), correlation = as.vector(corr))
      results = rbind(results, xx)
    }
  }
  
  # if brains within a dataset do not need to be split, compute Jacquard index, shared genes, and correlations b/w slices' cross-expression profiles vs distance between them
  if (!(target[i] == "Zador_2024_Mouse_BARseq" | target[i] == "Zhuang_2023_Mouse_MERFISH")){
    
    distance    = abs(outer(1:ncol(profiles_sig), 1:ncol(profiles_sig), "-"))
    jaccard_idx = jaccard(profiles_sig)
    raw_totals  = t(profiles_sig) %*% profiles_sig
    corr        = cor(profiles_corr)
    
    diag(distance) = NA; diag(jaccard_idx) = NA; diag(raw_totals) = NA; diag(corr) = NA
    results = data.frame(brain = target[i], distance = as.vector(distance), jaccard = as.vector(jaccard_idx), same_genes = as.vector(raw_totals), correlation = as.vector(corr))
  }
  
  # collate and save data
  outcome = rbind(outcome, results)
  fwrite(outcome, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_similarity_vs_distance.csv")
  print(str_c("Dataset no. ", signif(i/length(target), digits = 2), " is completed and saved"))
}