# load functions and libraries
source("/inkwell05/ameer/functions/0_source_functions.R")
library(peakRAM)

# load data
data = spatial_QC(path_to_expression = "/inkwell05/ameer/databases/spatial/Vizgen_2022_Mouse_MERSCOPE/expression/Vizgen_2022_Mouse_MERSCOPE_brain_2_slice_2_expression.h5ad",
                  path_to_metadata   = "/inkwell05/ameer/databases/spatial/Vizgen_2022_Mouse_MERSCOPE/metadata/Vizgen_2022_Mouse_MERSCOPE_brain_2_slice_2_metadata.csv")

locations = data$metadata[,c("x","y")]
data = data$data

print("Data and locations loaded")

# increase the size of the data to 500,000 cells and 20,000 genes
print("Increasing dataset size...")

data = data[1:50000, 1:100]
locations = locations[1:50000,]

data = data[, rep(seq_len(ncol(data)), times = 200)]
data = data[rep(seq_len(nrow(data)), times = 10), ]

locations = locations[rep(seq_len(nrow(locations)), times = 10), ]

print("Data size increased")
gc()

# specify how many genes and cells to include
cells = c(10000, 20000, 50000, 100000, 200000, 300000, 400000, 500000)
genes = c(1000, 2000, 5000, 10000, 15000, 20000)
combo = expand.grid(cells = cells, genes = genes)
combo = data.frame(combo, ram_MB = 0, time_sec = 0)

# benchmark memory (RAM in MB) and time (in seconds) usage
for (i in 1:nrow(combo)) {
  
  # run benchmark
  benchmark = tryCatch(
  {
    peakRAM({X = cross_expression(data = data[1:combo$cells[i], 1:combo$genes[i]], locations = locations[1:combo$cells[i], ]); if (exists("X")){rm(X)}})
  }, error = function(e) {
    list(Peak_RAM_Used_MiB = -1, Elapsed_Time_sec = -1)
  }
  )
  
  # store results
  combo[i, "ram_MB"]   = benchmark$Peak_RAM_Used_MiB
  combo[i, "time_sec"] = benchmark$Elapsed_Time_sec
  gc()
  
  # save results
  fwrite(x = combo, file = "/inkwell05/ameer/cross_expression_brain_meta_analysis/plot_datasets/speed_memory_benchmark.csv")
  print(signif(i/nrow(combo), digits = 2))
}

print("Successful completion")