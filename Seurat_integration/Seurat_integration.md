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
clip.max = "auto"
if (clip.max == "auto") {
      clip.max <- sqrt(x = ncol(x = object))
}
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

## assign 0 to some metrics (I don’t know why doing this)

``` r
hvf.info$variance.expected <- 0
hvf.info$variance.standardized <- 0
not.const <- hvf.info$variance > 0
```

## produce a fit using loess

This is try to regress and smooth the curve of log10(mean) and
log10(var). Loess Regression is the most common method used to smooth a
volatile time series. It is a non-parametric methods where least squares
regression is performed in localized subsets, which makes it a suitable
candidate for smoothing any numerical vector.

``` r
loess.span = 0.3
fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
      data = hvf.info[not.const, ], span = loess.span)
```

## get the expected variance for gene that has larger than 0 variance

``` r
hvf.info$variance.expected[not.const] <- 10^fit$fitted
```

## calculate variance of zscore which is captured by mean and expected variance

``` r
SparseRowVarStd <- function(mat, mu, sd, vmax, display_progress) {
    .Call('_Seurat_SparseRowVarStd', PACKAGE = 'Seurat', mat, mu, sd, vmax, display_progress)
}
hvf.info$variance.standardized <- SparseRowVarStd(mat = object, 
      mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected), 
      vmax = clip.max, display_progress = verbose)
```

## we selected top n genes that have highest standardized variance

``` r
hvf.info = hvf.info[order(hvf.info$variance.standardized, decreasing = T),]
top_n = 2000
top.genes = rownames(hvf.info)[1:top_n]
top.genes[1:20]
```

    ##  [1] "TIMP1"    "FTL"      "HLA-DRA"  "C15orf48" "IL8"      "SOD2"    
    ##  [7] "TYROBP"   "CCL2"     "FCER1G"   "FTH1"     "GNLY"     "CD74"    
    ## [13] "HBG1"     "HBA1"     "S100A8"   "CD63"     "S100A4"   "HLA-DRB1"
    ## [19] "HBA2"     "S100A11"

# learn SelectIntegrationFeatures

## initializa parameters and preprocess

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
