# scNCA 

scNCA is an R package for integrating multiple sources of temporal single cell RNA-seq data and correct the batch effects

## 1. Installation

The recommended installation method for `scNCA` is using `install_github` command from `devtools` library.  You will first have to have [devtools](https://github.com/hadley/devtools) package installed.

```r
library(devtools)
install_github('gongx030/scNCA')
```

A number of needed packages are installed in this process.

## 2. Quick Start

We first load the scNCA and related packages:
```r
library(scNCA)
library(SummarizedExperiment)
library(Matrix)
```

### 3 Using scNCA to integrate simulated balanced and imbalanced temporal scRNA-seq dataset.

Let us first simulate a balanced temporal scRNA-seq data with 2,000 genes, three lineages, three batches and five time points.  In the balanced data, each batch covers the cells from all lineages and time points.  We assume that there are 100 cells for each lineage/batch combination so there are in total 900 cells:
```r
ncells <- cbind(
  c(100, 100, 100),   # batch i
  c(100, 100, 100),   # batch ii
  c(100, 100, 100)    # batch iii
)
time.points <- cbind(
  c(1, 1, 1, 1, 1),   # batch i
  c(1, 1, 1, 1, 1),   # batch ii
  c(1, 1, 1, 1, 1)    # batch iii
)
```

```r
set.seed(1)
se <- simulateTemporalRNAseq(ncells = ncells, time.points = time.points)
```
The returned value `se` is an object of SummarizedExperiment class, which is commonly used to store the large-scale data in biology:
```shell
> se
class: SummarizedExperiment
dim: 2000 900
metadata(0):
assays(1): counts
rownames: NULL
rowData names(0):
colnames: NULL
colData names(5): batch group time time.table cell_label
```
The simulated scRNA-seq data can be accessed at `assays(se)$counts`.  In the `colData` section, `colData(se)$batch` indicates the batch label of each cell (1,2 and 3), `colData(se)$time` indicates the time index of each cell (1,2,3,4 and 5), `colData(se)$group` represents the simulated lineages (1,2 and 3), and `colData(se)$time.table` represents a 900 by 5 binary matrix indicates the time index of each cell.  

The input data will be scaled through cosine normalization, as suggested by [Haghverdi et al.](https://www.nature.com/articles/nbt.4091) to better capture the similar cells between batches. 
```r
X <- assays(se)$counts
X <- as.matrix(X %*% Diagonal(x = 1 / sqrt(colSums(X^2))))
```

The `scNCA` function integrates the temporal scRNA-seq data from three batches: 
```r
mf <- scNCA(X, batch = colData(se)$batch, time.table = colData(se)$time.table)
```
```shell
[2018-12-11 18:00:56] iter=    1 | loss=  -249.916
[2018-12-11 18:00:59] iter=  500 | loss=  -865.967
[2018-12-11 18:01:02] iter= 1000 | loss=  -887.694
[2018-12-11 18:01:05] iter= 1500 | loss=  -906.931
[2018-12-11 18:01:08] iter= 2000 | loss=  -924.331
[2018-12-11 18:01:11] iter= 2500 | loss=  -940.254
[2018-12-11 18:01:14] iter= 3000 | loss=  -954.904
[2018-12-11 18:01:17] iter= 3500 | loss=  -968.450
[2018-12-11 18:01:19] iter= 4000 | loss=  -981.072
[2018-12-11 18:01:22] iter= 4500 | loss=  -992.871
[2018-12-11 18:01:25] iter= 5000 | loss= -1003.848
```
The batch-free low dimensional representation of the input data can be accessed by 
```shell
> head(t(mf$Y))
          [,1]       [,2]      [,3]       [,4]        [,5]
[1,] 1.0050181 -1.9978689 1.3551587  1.9200974  0.16518551
[2,] 0.9516414 -2.1473274 1.5028186  2.3225472  0.08014079
[3,] 1.3695000 -1.8931294 1.2898606  1.8284671  0.09509533
[4,] 2.4596052 -0.9768379 0.1806223  2.2318383  0.10208578
[5,] 0.7063559  0.1022112 1.3341723 -0.2666387 -0.74208700
[6,] 1.1895738 -2.6235787 1.5147607  2.8488985 -0.05536205
```

Next, we will use [DDRTree](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5764547/) to visualize the integrated data. Qiu et al. has  shown that DDRTree excels at revealing the structure of complex trajectories with multiple branches from temporal scRNA-seq data.  We first define the shape and the color of each cells:
```r
n2bg <- c(
  '1_1' = 'red', '1_2' = 'lightblue', '1_3' = 'green',  # batch i
  '2_1' = 'red', '2_2' = 'lightblue', '2_3' = 'green',  # batch ii
  '3_1' = 'black', '3_2' = 'black', '3_3' = 'black'     # batch iii
)
n2pch <- c(
  '1_1' = 21, '1_2' = 21, '1_3' = 21, # batch i
  '2_1' = 24, '2_2' = 24, '2_3' = 24, # batch ii
  '3_1' = 3, '3_2' = 3, '3_3' = 3     # batch iii
)
n2col <- c(
  '1_1' = 'black', '1_2' = 'black', '1_3' = 'black',  # batch i
  '2_1' = 'black', '2_2' = 'black', '2_3' = 'black',  # batch ii
  '3_1' = 'red', '3_2' = 'lightblue', '3_3' = 'green' # batch iii
)
bg <- n2bg[colData(se)$cell_label]
pch <- n2pch[colData(se)$cell_label]
col2 <- n2col[colData(se)$cell_label]
```

Next, we run the DDRTree on the results from scNCA. 
```r
library(DDRTree)
dt <- DDRTree(mf$Y, verbose = TRUE, sigma = 1e-3, maxIter = 5)
y <- t(dt$Z)
plot(y[, 1], y[, 2], bg = bg, pch = pch, col = col2, cex = 1.5, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
```
![alt text](/docs/images/scnca_balanced_synthetic.png)

We may use tsBET to test the batch effects of the temporal scRNA-seq data. 
```r
tsBET(dt$stree, batch = colData(se)$batch, time.table = colData(se)$time.table)
```
```shell
$p.value
[1] 0.3910734

$neighborhood_size
[1] 100

$resample
[1] 100

$M
[1] 900
```
The first argument for `tsBET` is an adjacency matrix of cells.  Through our study, we use the adjacency matrix returned by DDRTree as the input. 

Similarly, we can generate an imbalanced temporal scRNA-seq data where each batch only covers a subset of lineages and time points: 
```
ncells <- cbind(
  c(100, 100, 0),   # batch i
  c(100, 100, 100),   # batch ii
  c(0, 100, 100)    # batch iii
)
time.points <- cbind(
  c(1, 1, 1, 0, 0),   # batch i
  c(0, 1, 1, 1, 0),   # batch ii
  c(0, 0, 1, 1, 1)    # batch iii
)
```
And we use scNCA to integrate the data and DDRTree to visualize the results:
```r
set.seed(1)
se <- simulateTemporalRNAseq(ncells = ncells, time.points = time.points)
X <- assays(se)$counts
X <- as.matrix(X %*% Diagonal(x = 1 / sqrt(colSums(X^2))))
mf <- scNCA(X, batch = colData(se)$batch, time.table = colData(se)$time.table)

bg <- n2bg[colData(se)$cell_label]
pch <- n2pch[colData(se)$cell_label]
col2 <- n2col[colData(se)$cell_label]

dt <- DDRTree(mf$Y, verbose = TRUE, sigma = 1e-3, maxIter = 5)
y <- t(dt$Z)
plot(y[, 1], y[, 2], bg = bg, pch = pch, col = col2, cex = 1.5, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
```
![alt text](/docs/images/scnca_imbalanced_synthetic.png)

And we test the batch effect:
```r
tsBET(dt$stree, batch = colData(se)$batch, time.table = colData(se)$time.table)
```
```shell
$p.value
[1] 0.4087084

$neighborhood_size
[1] 100

$resample
[1] 100

$M
[1] 700
```

### 4 Using scNCA to integrate scRNA-seq data of mouse early cardiovascular development
we used scNCA to integrate 3,387 single cells from six published temporal scRNA-seq datasets on mouse cardiovascular development from epiblast (E6.5) to the four-chambered heart (E9.5).  These datasets included the single cells from the whole epiblast3-5, Mesp1+ cardiovascular progenitors4, Flk1+ mesodermal progenitors5, Etv2+ hemato-endothelial progenitors6, as well as two scRNA-seq of distinct cardiac regions at E8.5 and E9.57,26.  We also added 264 Nkx2-5-EYFP+ single cells from E7.75 and E8.5.  Previous analysis showed that Nkx2-5 is expressed in multipotent progenitors in the mouse embryo, and targeted disruption of Nkx2-5, in the mouse, results in perturbed heart morphogenesis and perturbed endothelial and hematopoietic development at approximately E9.527.  Thus, the Nkx2-5-EYFP+ scRNA-seq from E7.75 and E8.5 bridged the gap between early cardiovascular progenitors represented by Mesp1+, Flk1+ and Etv2+ cells, and the late whole heart scRNA-seq data from E9.5.  

We first load the dataset as a SummarizedExperiment object:
```r
library(SummarizedExperiment)
se <- local({get(load(url('https://s3.msi.umn.edu/scNCA/data/mouse_heart_development.rda')))})
```
```shell
> se
class: SingleCellExperiment
dim: 23429 3387
metadata(1): log.exprs.offset
assays(2): counts logcounts
rownames(23429): 0610005C13Rik 0610007C21Rik ... a l7Rn6
rowData names(7): symbol mean ... p.value FDR
colnames: NULL
colData names(4): population time batch time.table
reducedDimNames(0):
spikeNames(0):
> table(colData(se)$batch, colData(se)$time)
             E6.5 E7.0-E7.25 E7.5 E7.75 E8.25-E8.5 E9.5
  DeLaughter    0          0    0     0          0  178
  Lescroart   172        341    0     0          0    0
  Li            0          0    0     0        140  556
  Mohammed    250          0    0     0          0    0
  Scialdone   501        138  259   307          0    0
  Gong          0         83    0   200        262    0
```

The raw read counts can be accessed at `assays(se)$counts` and the log-transformed and normalized can be accessed at `assays(se)$logcounts`.  The raw read counts were first normalized by the devolution-based size factors, as implemented by the function computeSumFactors in the [scran package](https://bioconductor.org/packages/release/bioc/html/scran.html). We used the `trendVar` function in the `scran` package to fitted a mean-variance trend model using endogenous genes.  The highly variable genes (HVGs) with a false discovery rate of 5% or less were used for the integration analysis.  

```r
n <- !is.na(rowData(se)$FDR) & rowData(se)$FDR < 0.05 # find HVG genes
X <- assays(se)$logcounts[n, ]
Xcos <- as.matrix(X %*% Diagonal(x = 1 / sqrt(colSums(X^2)))) # cosine scaling
```

We use scNCA to integrate six datasets: 
```r
set.seed(1)
mf <- scNCA(Xcos, batch = colData(se)$batch, time.table = colData(se)$time.table)
```
It would take ~6 mins on a 3.3GHz Intel Core i6 machine.  Since the optimization is performed by TensorFlow, running the optimization in parallel with GPUs would greatly speed up the process for large datasets. 

Then, we will use DDRTree to visualize the integrated dataset:
```r
library(Rtsne)
library(DDRTree)
set.seed(2)
y <- Rtsne(t(as.matrix(mf$Y)), check_duplicates = FALSE)$Y # initialize DDRTree space
dt <- DDRTree(t(y), verbose = TRUE, sigma = 1e-3, maxIter = 5)
```

Prepare the cell colors and make the plots:
```r
par(mfrow = c(1, 2))
study2bg <- c(
  'Mohammed' = 'gold',
  'DeLaughter' = 'red',
  'Gong' = 'green',
  'Lescroart' = 'blue',
  'Li' = 'brown',
  'Scialdone' = 'cyan'
)
bg <- study2bg[colData(se)$batch]
y <- t(dt$Z) + rnorm(prod(dim(dt$Z)), mean = 0, sd = diff(range(dt$Z)) / 100)
plot(y[, 1], y[, 2], bg = bg, pch = 21, cex = 0.4, col = bg, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'By batch')


tc <- colorpanel(100, low = 'blue', mid = 'gray', high = 'red')
bg.time <- sapply(time, function(t) tc[round(as.numeric(t) * 100 / nlevels(time))])
plot(y[, 1], y[, 2], bg = bg.time, pch = 21, cex = 0.4, col = bg.time, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'By time index')
```
![alt text](/docs/images/scnca_heart_development.png)

The integrated scRNA-seq dataset of mouse cardiovascular development can be explored at [https://z.umn.edu/scNCA](https://z.umn.edu/scNCA).  

# 4. Session Information
```shell
> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices datasets  utils     methods   base

other attached packages:
 [1] DDRTree_0.1.5               irlba_2.3.2                 Rtsne_0.13                  scNCA_1.0.1                 scran_1.8.4                 SingleCellExperiment_1.2.0  Matrix_1.2-14               SummarizedExperiment_1.10.1 DelayedArray_0.6.6          BiocParallel_1.14.2         matrixStats_0.54.0
[12] Biobase_2.40.0              GenomicRanges_1.32.7        GenomeInfoDb_1.16.0         IRanges_2.14.12             S4Vectors_0.18.3            BiocGenerics_0.26.0         BiocInstaller_1.30.0

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0         colorspace_1.3-2         rjson_0.2.20             dynamicTreeCut_1.63-1    rprojroot_1.3-2          XVector_0.20.0           base64enc_0.1-3          fs_1.2.6                 rstudioapi_0.8           remotes_2.0.1            DT_0.4                   tximport_1.8.0
 [13] scater_1.8.4             pkgload_1.0.2            spam_2.2-0               jsonlite_1.5             cluster_2.0.7-1          tfruns_1.4               shinydashboard_0.7.1     shiny_1.1.0              compiler_3.5.1           backports_1.1.2          assertthat_0.2.0         lazyeval_0.2.1
 [25] limma_3.36.5             cli_1.0.1                later_0.7.5              htmltools_0.3.6          prettyunits_1.0.2        tools_3.5.1              bindrcpp_0.2.2           igraph_1.2.2             dotCall64_1.0-0          gtable_0.2.0             glue_1.3.0               GenomeInfoDbData_1.1.0
 [37] reshape2_1.4.3           dplyr_0.7.7              maps_3.3.0               Rcpp_1.0.0               gdata_2.18.0             DelayedMatrixStats_1.2.0 stringr_1.3.1            ps_1.2.0                 mime_0.6                 gtools_3.8.1             devtools_2.0.1           statmod_1.4.30
 [49] edgeR_3.22.5             zlibbioc_1.26.0          MASS_7.3-51              scales_1.0.0             promises_1.0.1           rhdf5_2.24.0             fields_9.6               memoise_1.1.0            reticulate_1.10          gridExtra_2.3            ggplot2_3.1.0            stringi_1.2.4
 [61] tensorflow_1.9           desc_1.2.0               caTools_1.17.1.1         pkgbuild_1.0.2           rlang_0.3.0              pkgconfig_2.0.2          bitops_1.0-6             lattice_0.20-35          purrr_0.2.5              Rhdf5lib_1.2.1           bindr_0.1.1              htmlwidgets_1.3
 [73] processx_3.2.0           tidyselect_0.2.5         plyr_1.8.4               magrittr_1.5             R6_2.3.0                 gplots_3.0.1             pillar_1.3.0             whisker_0.3-2            withr_2.1.2              RCurl_1.95-4.11          tibble_1.4.2             crayon_1.3.4
 [85] KernSmooth_2.23-15       viridis_0.5.1            usethis_1.4.0            locfit_1.5-9.1           grid_3.5.1               data.table_1.11.8        FNN_1.1.2.1              callr_3.0.0              digest_0.6.18            xtable_1.8-3             httpuv_1.4.5             munsell_0.5.0
 [97] beeswarm_0.2.3           viridisLite_0.3.0        vipor_0.4.5              sessioninfo_1.1.0
```
