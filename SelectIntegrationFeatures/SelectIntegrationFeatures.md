Learn SelectIntegrationFeatures
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
```

## initializa parameters and preprocess

``` r
nfeatures = 2000
assay = NULL
verbose = TRUE
fvf.nfeatures = 2000
object.list = ifnb.list
```

## assign assays to each obj in the list

``` r
if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    for (ii in length(x = object.list)) {
      DefaultAssay(object = object.list[[ii]]) <- assay[ii]
    }
}else {
  assay <- sapply(X = object.list, FUN = DefaultAssay)
}
```

## check if variable features exist

``` r
for (ii in 1:length(x = object.list)) {
  if (length(x = VariableFeatures(object = object.list[[ii]])) == 
    0) {
    if (verbose) {
      message(paste0("No variable features found for object", 
        ii, " in the object.list. Running FindVariableFeatures ..."))
    }
    object.list[[ii]] <- FindVariableFeatures(object = object.list[[ii]], 
      nfeatures = fvf.nfeatures, verbose = verbose, 
      ...)
  }
}
```

## fetch highly variable feature names from all seurat obj in the list

``` r
var.features <- unname(obj = unlist(x = lapply(X = 1:length(x = object.list), 
    FUN = function(x) VariableFeatures(object = object.list[[x]], 
      assay = assay[x]))))
```

## find shared highly variable genes

``` r
var.features <- sort(x = table(var.features), decreasing = TRUE)
#make sure highly variable genes should present in all seurat obj 
#(even may not highly variable in every)
for (i in 1:length(x = object.list)) {
  var.features <- var.features[names(x = var.features) %in% 
    rownames(x = object.list[[i]][[assay[i]]])]
}
```

## select features that have highly variable freq larger than tie value

``` r
tie.val <- var.features[min(nfeatures, length(x = var.features))]
features <- names(x = var.features[which(x = var.features > 
    tie.val)])
length(features)
```

    ## [1] 787

``` r
features[1:20]
```

    ##  [1] "A2M"         "ABHD12"      "ABI3"        "AC002456.2"  "AC006129.4" 
    ##  [6] "AC007228.11" "ACP5"        "ACRBP"       "ACSL1"       "ACTB"       
    ## [11] "ACTG1"       "ADA"         "ADK"         "ADM"         "ADORA3"     
    ## [16] "ADPRM"       "ADTRP"       "AGPAT9"      "AHSP"        "AIF1"

## fetch variable features for each dataset

``` r
vf.list <- lapply(X = object.list, FUN = VariableFeatures)
vf.list[[1]][1:10]
```

    ##  [1] "HBB"    "HBA2"   "HBA1"   "CCL3"   "CCL4"   "CXCL10" "TXN"    "CCL7"  
    ##  [9] "CCL2"   "GNLY"

## re-rank features based on their median rank for all vairable gene set from each dataset

``` r
#features are genes that shared by more than tie.val datasets
if (length(x = features) > 0) {
  feature.ranks <- sapply(X = features, FUN = function(x) {
    ranks <- sapply(X = vf.list, FUN = function(vf) {
      if (x %in% vf) {
        return(which(x = x == vf))
      }
      return(NULL)
    })
    median(x = unlist(x = ranks))
  })
  features <- names(x = sort(x = feature.ranks))
}
```

## grab the features that equal to tie.val

``` r
features.tie <- var.features[which(x = var.features == tie.val)]
```

## do the same re-rank for tie features

``` r
tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
  ranks <- sapply(X = vf.list, FUN = function(vf) {
    if (x %in% vf) {
      return(which(x = x == vf))
    }
    return(NULL)
  })
  median(x = unlist(x = ranks))
})
```

## return a vector with features in the front and tie.features in the back

``` r
features <- c(features, names(x = head(x = sort(x = tie.ranks), 
    nfeatures - length(x = features))))
```
