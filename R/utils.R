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

#Find highly variable genes in a list of Seurat objects
findHVG <- function(obj.list, nfeatures=2000,
                    min.exp=0.01, max.exp=3.0) {
  obj.list <- lapply(obj.list, function(x){
    
    ncalc <- min(5*nfeatures, nrow(x))
    x <- Seurat::FindVariableFeatures(x, nfeat=ncalc, verbose=FALSE)
    varfeat <- VariableFeatures(x)
    #Also exclude genes that are very poorly or always expressed
    means <- apply(GetAssayData(x, slot="data")[varfeat, ], 1, mean)
    removeGenes <- names(means[means<min.exp | means>max.exp])
    
    varfeat <- setdiff(varfeat, removeGenes)
    n <- min(length(varfeat), nfeatures)
    VariableFeatures(x) <- varfeat[1:n]
    x
  })
  hvg <- Seurat::SelectIntegrationFeatures(obj.list,nfeatures = nfeatures,
                                           verbose = FALSE)
  return(hvg)
}