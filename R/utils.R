#Calculate Jaccard Index
jaccardIndex <- function(a, b) {
  i = length(intersect(a, b))
  u = length(a) + length(b) - i
  return (i/u)
}

#Center matrix
centerData <- function(data, non_negative=T) {
  mean <- apply(data, 1, mean)
  data <- data - mean
  if (non_negative) {
    data[data<0] <- 0
  }
  return(data)
}

#Calculate entropy
getEntropy <- function(tab, pseudo=0.1) {
  p <- (tab+pseudo) / sum(tab+pseudo)
  entropy <- -sum(p * log2(p), na.rm = TRUE)
  return(entropy)
}


#' Find variable features
#'
#' Select highly variable genes (HVG) from an expression matrix. Genes from a blocklist
#' (e.g. cell cycling genes, mitochondrial genes) can be excluded from the list of
#' variable genes, as well as genes with very low or very high average expression
#'
#' @param obj A Seurat object containing an expression matrix
#' @param nfeatures Number of top HVG to be returned
#' @param genesBlocklist Optionally takes a vector or list of vectors of gene names.
#'     These genes will be ignored for HVG detection. This is useful to mitigate effect
#'     of genes associated with technical artifacts or batch effects
#'     (e.g. mitochondrial, heat-shock response). If set to `NULL` no genes will be excluded
#' @param min.exp Minimum average normalized expression for HVG. If lower, the gene will be excluded
#' @param max.exp Maximum average normalized expression for HVG. If higher, the gene will be excluded
#' @return Returns a list of highly variable genes
#' @import Seurat
#' 
FindVariableFeatures_wfilters <- function(
    obj,
    nfeatures=2000,
    genesBlockList="default",
    min.exp=0.01,
    max.exp=3)
{
  
  assay <- DefaultAssay(obj)
  #Calculate a fixed number of HVG, then filtered to nfeat at the end
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 10000, verbose=F)
  
  varfeat <- VariableFeatures(obj)
  
  if (is.vector(genesBlockList)) {
    genes.block <- genesBlockList # user-provided vector
  } else {
    genes.block <- NULL # No excluded genes
  }
  
  varfeat <- setdiff(varfeat, genes.block)
  
  #Also remove genes that are very poorly or always expressed (=not really variable genes)
  means <- apply(GetAssayData(obj, assay=assay, slot="data")[varfeat,], 1, mean)
  removeGenes2 <- names(means[means<min.exp | means>max.exp])
  
  varfeat <- setdiff(varfeat, removeGenes2)
  n <- min(length(varfeat), nfeatures)
  
  VariableFeatures(obj) <- varfeat[1:n]
  
  return(obj)
}  

#Find highly variable genes in a list of Seurat objects
findHVG <- function(obj.list, nfeatures=2000,
                    min.exp=0.01, max.exp=3.0, hvg.blocklist=NULL) {
  obj.list <- lapply(obj.list, function(x){
    
    ncalc <- min(5*nfeatures, nrow(x))
    x <- FindVariableFeatures_wfilters(x, nfeatures=ncalc, min.exp=min.exp, max.exp=max.exp,
                                       genesBlockList=hvg.blocklist)
    x
  })
  hvg <- Seurat::SelectIntegrationFeatures(obj.list,nfeatures = nfeatures,
                                           verbose = FALSE)
  return(hvg)
}