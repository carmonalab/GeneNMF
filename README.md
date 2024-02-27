# GeneNMF: unsupervised discovery of gene programs in single-cell data

<p align="center">
  <img height="80" src="inst/RSticker_GeneNMF.png">
</p>

Non-negative matrix factorization is a method for the analysis of high dimensional data that allows extracting sparse and meaningful features from a set of non-negative data vectors. It is well suited for decomposing scRNA-seq data, effectively reducing large complex matrices ($10^4$ of genes times $10^5$ of cells) into a few interpretable gene programs. It has been especially used to extract recurrent gene programs in cancer cells (see e.g. [Barkely et al. (2022)](https://www.nature.com/articles/s41588-022-01141-9) and [Gavish et al. (2023)](https://www.nature.com/articles/s41586-023-06130-4)), which are otherwise difficult to integrate and analyse jointly.

**GeneNMF** is a package that implements methods for NMF decomposition for single-cell omics data. It can be applied directly on Seurat objects to reduce the dimensionality of the data and to detect robust gene programs across multiple samples.  

## Installation
```r
library(remotes)
remotes::install_github("carmonalab/GeneNMF")
```

## Test your installation
```r
library(GeneNMF)
data(sampleObj)
sampleObj <- runNMF(sampleObj, k=5)
```

## GeneNMF demo
Find a demo of the functionalities of GeneNMF in the following tutorial: [HTML](https://carmonalab.github.io/GeneNMF.demo/NMF_demo_PBMC.html) and [repository](https://github.com/carmonalab/GeneNMF.demo).

## References
Coming soon!
