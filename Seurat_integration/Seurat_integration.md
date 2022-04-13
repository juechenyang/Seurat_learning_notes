Integration pipeline
================
Juechen Yang
4/13/2022

# load libraries

``` r
library(Seurat)
library(SeuratData)
library(patchwork)
LoadData("ifnb")
```

    ## An object of class Seurat 
    ## 14053 features across 13999 samples within 1 assay 
    ## Active assay: RNA (14053 features, 0 variable features)

# learn FindVariableFeatures

## create object data for convenience

``` r
ifnb.list <- SplitObject(ifnb, split.by = "stim")
object = ifnb.list[[1]]
#normalization
object = NormalizeData(object = object, verbose = F)
```

## get assay data for object

Because FindVariableFeatures needs the count matrix for input, so it is
required to fetch the assay data from the Seurat object for downstream
calculation

``` r
object = GetAssayData(object)
#initializa verbose
verbose = FALSE
```

## computes the mean and variance for each row(gene)

``` r
SparseRowVar2 <- function(mat, mu, display_progress) {
    .Call('_Seurat_SparseRowVar2', PACKAGE = 'Seurat', mat, mu, display_progress)
}
hvf.info <- data.frame(mean = rowMeans(x = object))
hvf.info$variance <- SparseRowVar2(mat = object, mu = hvf.info$mean, 
                      display_progress = verbose)
hvf.info[1:10,]
```

    ##                      mean    variance
    ## AL627309.1    0.001290017 0.002743297
    ## RP11-206L10.2 0.001149103 0.002199545
    ## LINC00115     0.010434879 0.020462062
    ## NOC2L         0.148424623 0.274238307
    ## KLHL17        0.003280740 0.006510261
    ## PLEKHN1       0.004709366 0.007916569
    ## HES4          0.030695027 0.064723658
    ## ISG15         0.648610052 1.218974327
    ## AGRN          0.003782200 0.005732422
    ## C1orf159      0.005994965 0.011509557
