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
findHVG <- function(obj.list, nfeatures=2000) {
  obj.list <- lapply(obj.list, function(x){
    FindVariableFeatures.STACAS(x, nfeat=nfeatures)
  })
  hvg <- SelectIntegrationFeatures(obj.list, nfeatures = nfeatures)
  return(hvg)
}