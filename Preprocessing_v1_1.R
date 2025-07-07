```{r}

## whole-brain BAR-seq data is registered to the Allen Common Coordinate Framework version 3 (CCFv3)
## data is shrunk by removing image stitching related artefacts (cf. Xiaoyin's email)
## data is quality controlled by keeping cells with genes/cell >= 5 and reads/cell >= 20
## data alongside CCF and slide coordinates are saved and can be used for analysis

## load libraries
suppressPackageStartupMessages(library(xfun))
pkgs = c("SingleCellExperiment","tidyverse","data.table","dendextend","fossil","gridExtra","gplots","metaSEM","foreach","Matrix",
         "matrixStats","remotes","ggpubr","ggridges","patchwork","dineq","ComplexHeatmap","fastcluster","doParallel")
suppressPackageStartupMessages(xfun::pkg_attach(pkgs))

## load data
coord    = as.data.frame(fread(file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/neuron_position_whole-brain.csv"))
data     = readRDS("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/barseq2_whole-brain.rds")
data     = SingleCellExperiment::counts(data)

CCF      = data.frame(as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/cell_id.csv")),
                      as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/CCF_info.csv")))
colnames(CCF) = c("sample_id","CCFx","CCFy","CCFz","CCFparentname","CCFname","CCFnano")

## keep coordinates for plotting whole-brain images
fwrite(x = data.frame(coord[which(coord$sample_id %in% CCF$sample_id),],
                      CCF[,2:ncol(CCF)]), file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/all_cell_coordinates.csv")

## subset data and coordinates by CCF (registered and truncated data after removing image stitching artefacts; cf. Xiaoyin's email)
data     = data[,which(colnames(data) %in% as.character(CCF$sample_id))]
coord    = coord[which(coord$sample_id %in% CCF$sample_id),]
coord    = data.frame(coord, CCF[,2:ncol(CCF)])

## quality control to keep cells with genes/cell >= 5 and reads/cell >= 20
data     = suppressWarnings(data[,which(apply(X = data > 0, MARGIN = 2, FUN = sum) >= 5)])
data     = suppressWarnings(data[,which(apply(X = data,     MARGIN = 2, FUN = sum) >= 20)])

#  subset data by relevant genes
gens     = read.csv("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/genes.csv")
data     = data[which(rownames(data) %in% gens$Genes),]

## subset coordinates by quality controlled cells and add labels
coord    = coord[which(as.character(coord$sample_id) %in% colnames(data)),]
labels   = as.data.frame(fread("/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/labels.csv", header = TRUE))
labels   = data.frame(labels[,c("sample_id","slice","class")], labels[,"subclass"], labels[,c("cluster","notes")])
colnames(labels)[c(3,4,5)] = c("class_H1","subclass_H2","cluster_H3")
coord    = data.frame(coord,labels[,3:ncol(labels)])

```


```{r}

## re-arrange the slices to brain images in the correct orientation
final_coord <- coord

# target slice = n
slice_target = 5

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

final_coord[which(coord$slice==slice_target),"pos_x"] = -coord[which(coord$slice==slice_target),]$pos_x

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 6

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

final_coord[which(coord$slice==slice_target),"pos_x"] = -coord[which(coord$slice==slice_target),]$pos_x

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 7

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

final_coord[which(coord$slice==slice_target),"pos_x"] = -coord[which(coord$slice==slice_target),]$pos_x

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 8

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

final_coord[which(coord$slice==slice_target),"pos_x"] = -coord[which(coord$slice==slice_target),]$pos_x

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 9

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

final_coord[which(coord$slice==slice_target),"pos_y"] = -coord[which(coord$slice==slice_target),]$pos_y

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 10

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

final_coord[which(coord$slice==slice_target),"pos_y"] = -coord[which(coord$slice==slice_target),]$pos_y

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 11

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 12

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 13

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 14

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 15

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 16

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 17

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 18

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 19

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 20

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 21

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 22

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 23

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 24

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 25

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 26

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_x = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord$pos_x
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 27

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

slice_target = 27

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_y, y = -pos_x) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 28

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 29

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 30

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 31

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_y = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 32

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_y = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 33

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_y = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 34

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y) %>% mutate(pos_y = -pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 35

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 36

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 37

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 38

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 39

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)


# target slice = n
slice_target = 40

coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

x_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_y)
y_coord <- final_coord %>% filter(slice==slice_target) %>% select(pos_x) %>% mutate(pos_x = -pos_x)

final_coord[which(coord$slice==slice_target),"pos_x"] = x_coord
final_coord[which(coord$slice==slice_target),"pos_y"] = y_coord

final_coord %>% filter(slice==slice_target) %>%
ggplot(.) +
  aes(x = pos_x, y = pos_y) +
  geom_point(size=0.1, alpha=0.3)

```


```{r}
## save processed data
coord <- final_coord
writeMM(as(data, "dgCMatrix"), file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_expression_matrix.txt")
fwrite(x = coord, file = "/Users/AmeerSarwar/Desktop/Gene co-expression patterns to define cellular networks/Raw Data/BARseq_processed_coordinates_plus_CCF.csv")

```
