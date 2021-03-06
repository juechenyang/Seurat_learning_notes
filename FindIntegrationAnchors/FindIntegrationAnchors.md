Learn FindIntegrationAnchors
================
Juechen Yang
5/18/2022

## load libraries

``` r
library(Seurat)
library(SeuratData)
library(patchwork)
#InstallData("ifnb")
LoadData("ifnb")
```

    ## An object of class Seurat 
    ## 14053 features across 13999 samples within 1 assay 
    ## Active assay: RNA (14053 features, 0 variable features)

## preprocess

``` r
ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
```

## initailiza parameters

``` r
object.list = ifnb.list
assay = NULL
reference = NULL 
anchor.features = features
scale = TRUE
normalization.method = c("LogNormalize", "SCT")
sct.clip.range = NULL
reduction = c("cca", "rpca", "rlsi")
l2.norm = TRUE
dims = 1:30
k.anchor = 5 
k.filter = 200
k.score = 30
max.features = 200
nn.method = "annoy" 
n.trees = 50
eps = 0
verbose = FALSE
```

## finalized methods

``` r
normalization.method <- normalization.method[1]
reduction <- reduction[1]
```

## define a future function

``` r
library("future.apply")
```

    ## Loading required package: future

``` r
my.lapply <- ifelse(test = verbose && nbrOfWorkers() == 
    1, yes = pblapply, no = future_lapply)
```

## get number of cells for each dataset

``` r
object.ncells <- sapply(X = object.list, FUN = function(x) dim(x = x)[2])
```

## get default assay for each dataset

``` r
assay <- sapply(X = object.list, FUN = DefaultAssay)
```

## check duplicate cell names

``` r
# initialize slot tools for Integration
object.list <- lapply(X = object.list, FUN = function(obj) {
  slot(object = obj, name = "tools")$Integration <- NULL
  return(obj)
})
# check duplicate cell names
CheckDuplicateCellNames <- function(object.list, verbose = TRUE, stop = FALSE) {
  cell.names <- unlist(x = lapply(X = object.list, FUN = colnames))
  if (any(duplicated(x = cell.names))) {
    if (stop) {
      stop("Duplicate cell names present across objects provided.")
    }
    if (verbose) {
      warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    }
    object.list <- lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        return(RenameCells(
          object = object.list[[x]],
          new.names = paste0(Cells(x = object.list[[x]]), "_", x)
        ))
      }
    )
  }
  return(object.list)
}
object.list <- CheckDuplicateCellNames(object.list = object.list)
```

## scale data for only those selected features for integration

``` r
slot <- "data"
object.list <- my.lapply(X = object.list, FUN = function(object){
  ScaleData(object = object, features = anchor.features, 
    verbose = FALSE)
})
nn.reduction <- reduction
internal.neighbors <- list()
```

## build combinations

``` r
combinations <- expand.grid(1:length(x = object.list), 1:length(x = object.list))
combinations <- combinations[combinations$Var1 < combinations$Var2, 
    , drop = FALSE]
```

## get cumulative number of cells without last dataset

``` r
objects.ncell <- sapply(X = object.list, FUN = ncol)
offsets <- as.vector(x = cumsum(x = c(0, objects.ncell)))[1:length(x = object.list)]
```

## create diet seurat obj
