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

# we will predict the significance (alpha = 0.05) of one slice using correlation of other slices within a brain
# we will do this for each held-out slice as a function of including more slices
# the output is AUROC, labels are the p-values (1/0) of held-out slice, and prediction scores are average cross-expression correlation of slices included

target = c("Zador_2024_Mouse_BARseq", "Zhuang_2023_Mouse_MERFISH") # datasets with multiple brains needing splitting
dirr = "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/"
dirr = str_c(dirr, list.files(dirr))

datasets = "/inkwell05/ameer/databases/spatial/"
datasets = str_c(datasets, list.files(datasets))

output = data.frame()

# parallel processing setup
memory_per_worker_gb = 10             # in GB -- note: this is the only parameter you can change (estimate memory based on dataset size and how often it is repeated in the loops)
total_cores = parallel::detectCores() # number of cores
total_memory_gb = as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)) / 1024^2 # memory in GB

# calculate maximum number of possible workers given memory allocation per worker
max_workers = floor((total_memory_gb * 0.8) / memory_per_worker_gb) # maximum number of workers
num_workers = min(total_cores - 2, max_workers)                     # maximum workers possible given available cores (and leave 2 cores)

# set memory allocation and parallel processing plan
plan(multisession, workers = num_workers)
options(future.globals.maxSize = memory_per_worker_gb * 1024^3)

# process each dataset in turn
for (i in 1:length(datasets)){
  
  # extract relevant dataset
  temp_dirr = mixedsort(dirr[str_detect(dirr, str_remove(datasets[i], "/inkwell05/ameer/databases/spatial/"))])
  if (length(temp_dirr)==1){next} # skip if data contains only one sample
  
  print(str_c("Loading dataset no. = ", signif(i / length(datasets), 2)))
  
  edges = matrix(data = 0, nrow = nrow(as.data.frame(fread(temp_dirr[1]))), ncol = length(temp_dirr))
  nam   = str_remove(temp_dirr, "/inkwell05/ameer/cross_expression_brain_meta_analysis/cross_expression_edge_lists/")
  nam   = str_remove(nam, "_cross_expression.csv")
  nam   = str_remove(nam, str_c(datasets[i], "_"))
  colnames(edges) = nam
  rownames(edges) = 1:nrow(edges)
  score = edges
  for (k in 1:length(temp_dirr)){
    x = as.data.frame(fread(temp_dirr[k]))
    edges[,k] = x$cross_sig
    score[,k] = x$correlation
    print(k/length(temp_dirr))
  }
  
  print(str_c("Processing dataset no. = ", signif(i/length(datasets), 2)))
  
  
  # cross-validation: how well do samples (avg. correlation) predict a held-out sample's p-values (1/0)
  # predict the held-out sample as a function of the number of samples used for prediction
  
  if (sum(str_detect(datasets[i], target)) == 0){
    
    iter = 100
    outcome = matrix(0, nrow = ncol(score) - 1, ncol = ncol(score) + 1)
    dimnames(outcome) = list(1:nrow(outcome), 1:ncol(outcome))
    outcome[, 1] = 1:nrow(outcome)
    
    # cross-validation by leaving out one slice and predicting using other slices
    res_list = future_lapply(1:ncol(edges), function(p) {
      
      to_predict = as.numeric(edges[, p])
      scores = score[, -p]
      
      # repeat the procedure N times and average the data
      sapply(1:ncol(scores), function(q) {
        temp_auroc = 0
        for (k in 1:iter) {
          temp_scores = rowMedians(as.matrix(scores[, sample(1:ncol(scores), size = q)]))
          temp_auroc = temp_auroc + auroc(scores = temp_scores, labels = to_predict)
        }
        
        # print out log file and collect garbage
        print(str_c("Dataset no. = ", signif(i / length(datasets), 2), "; slice no. = ", signif(p / ncol(edges), 2), "; iteration no. = ", signif(q / iter, 2)))
        gc()
        
        # return results
        temp_auroc / iter
        
      })
    }, future.seed = TRUE)
    
    # fill in the outcome matrix
    for (p in seq_along(res_list)) {
      outcome[, p + 1] = res_list[[p]]
    }
    
    # process and annotate dataset
    result = rowMedians(outcome[,2:ncol(outcome)], na.rm = TRUE)
    result = data.frame(data = str_remove(datasets[i], "/inkwell05/ameer/databases/spatial/"), slices_included = 1:length(result), auroc = as.numeric(result))
    
    # append results across datasets
    output = rbind(output, result)
    
    # save results (and keep overwriting the saved file)
    fwrite(output, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_predictions_vs_slices_cross_validation.csv")
    print(str_c("Dataset no. ", signif(i/length(datasets), digits = 2), " is completed and saved"))
    gc()
  }
  
  
  # process Zador_2024_Mouse_BARseq
  if (str_detect(datasets[i], "Zador_2024_Mouse_BARseq")){
    
    to_extract = c("control1","control2","control3","control4","enucleated1","enucleated2","enucleated3","enucleated4","pilot")
    results = data.frame()
    
    for (k in 1:length(to_extract)){
      
      temp_edges = edges[, str_detect(colnames(edges), to_extract[k])]
      temp_score = score[, str_detect(colnames(edges), to_extract[k])]
      
      iter = 100
      outcome = matrix(0, nrow = ncol(temp_score) - 1, ncol = ncol(temp_score) + 1)
      dimnames(outcome) = list(1:nrow(outcome), 1:ncol(outcome))
      outcome[, 1] = 1:nrow(outcome)
      
      # cross-validation by leaving out one slice and predicting using other slices
      res_list = future_lapply(1:ncol(temp_edges), function(p) {
        
        to_predict = as.numeric(temp_edges[, p])
        scores     = temp_score[, -p]
        
        # repeat the procedure N times and average the data
        sapply(1:ncol(scores), function(q) {
          temp_auroc = 0
          for (mm in 1:iter) {
            temp_scores = rowMedians(as.matrix(scores[, sample(1:ncol(scores), size = q)]))
            temp_auroc = temp_auroc + auroc(scores = temp_scores, labels = to_predict)
          }
          
          # print out progress and collect garbage
          print(str_c("Dataset no. = ", signif(i / length(datasets), 2), "; brain no. = ", signif(k / length(to_extract), 2), "; slice no. = ", signif(p / ncol(temp_edges), 2), "; iteration no. = ", signif(mm / iter, 2)))
          gc()
          
          # return results
          temp_auroc / iter
          
        })
      }, future.seed = TRUE)
      
      # fill in the outcome matrix
      for (p in seq_along(res_list)) {
        outcome[, p + 1] = res_list[[p]]
      }
      
      # process and annotate dataset
      result = rowMedians(outcome[,2:ncol(outcome)], na.rm = TRUE)
      result = data.frame(data = str_c(str_remove(datasets[i], "/inkwell05/ameer/databases/spatial/"), "_", to_extract[k]), slices_included = 1:length(result), auroc = as.numeric(result))
      
      # append results across datasets
      output = rbind(output, result)
      
      # save results (and keep overwriting the saved file)
      fwrite(output, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_predictions_vs_slices_cross_validation.csv")
      print(str_c("Dataset no. ", signif(i/length(datasets), digits = 2), " and brain no. = ", signif(k / length(to_extract), 2), " is completed and saved"))
      gc()
    }
  }
  
  
  # process Zhuang_2023_Mouse_MERFISH
  if (str_detect(datasets[i], "Zhuang_2023_Mouse_MERFISH")){
    
    to_extract = c("sagittal","coronal_66","coronal_147")
    results = data.frame()
    
    for (k in 1:length(to_extract)){
      
      temp_edges = edges[, str_detect(colnames(edges), to_extract[k])]
      temp_score = score[, str_detect(colnames(edges), to_extract[k])]
      
      iter = 100
      outcome = matrix(0, nrow = ncol(temp_score) - 1, ncol = ncol(temp_score) + 1)
      dimnames(outcome) = list(1:nrow(outcome), 1:ncol(outcome))
      outcome[, 1] = 1:nrow(outcome)
      
      # cross-validation by leaving out one slice and predicting using other slices
      res_list = future_lapply(1:ncol(temp_edges), function(p) {
        
        to_predict = as.numeric(temp_edges[, p])
        scores     = temp_score[, -p]
        
        # repeat the procedure N times and average the data
        sapply(1:ncol(scores), function(q) {
          temp_auroc = 0
          for (mm in 1:iter) {
            temp_scores = rowMedians(as.matrix(scores[, sample(1:ncol(scores), size = q)]))
            temp_auroc = temp_auroc + auroc(scores = temp_scores, labels = to_predict)
          }
          
          # print out progress and collect garbage
          print(str_c("Dataset no. = ", signif(i / length(datasets), 2), "; brain no. = ", signif(k / length(to_extract), 2), "; slice no. = ", signif(p / ncol(temp_edges), 2), "; iteration no. = ", signif(mm / iter, 2)))
          gc()
          
          # return results
          temp_auroc / iter
          
        })
      }, future.seed = TRUE)
      
      # fill in the outcome matrix
      for (p in seq_along(res_list)) {
        outcome[, p + 1] = res_list[[p]]
      }
      
      # process and annotate dataset
      result = rowMedians(outcome[,2:ncol(outcome)], na.rm = TRUE)
      result = data.frame(data = str_c(str_remove(datasets[i], "/inkwell05/ameer/databases/spatial/"), "_", to_extract[k]), slices_included = 1:length(result), auroc = as.numeric(result))
      
      # append results across datasets
      output = rbind(output, result)
      
      # save results (and keep overwriting the saved file)
      fwrite(output, "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/datasets_cross_expression_predictions_vs_slices_cross_validation.csv")
      print(str_c("Dataset no. ", signif(i/length(datasets), digits = 2), " and brain no. = ", signif(k / length(to_extract), 2), " is completed and saved"))
      gc()
    }
  }
}

print("Successful completion")