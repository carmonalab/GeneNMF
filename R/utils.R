#Calculate Jaccard Index
jaccardIndex <- function(a, b) {
  i = length(intersect(a, b))
  u = length(a) + length(b) - i
  return (i/u)
}

jaccardSimilarity <- function(gene.vectors) {
  nprogs <- length(gene.vectors)
  J <- matrix(data=0, ncol=nprogs, nrow = nprogs)
  colnames(J) <- names(gene.vectors)
  rownames(J) <- names(gene.vectors)
  for (i in 1:nprogs) {
    for (j in 1:nprogs) {
      J[i,j] <- jaccardIndex(names(gene.vectors[[i]]), names(gene.vectors[[j]]))
    }  
  }
  return(J)
}

#Calculate cosine similarity
cosineSimilarity <- function(gene.vectors) {
  gene.table <- geneList2table(gene.vectors)
  cosine.matrix <- cosine(as.matrix(gene.table))
  return(cosine.matrix)
}

#From list of genes to complete table, with imputed zeros
geneList2table <- function(gene.vectors) {
  allgenes <- unique(unlist(lapply(gene.vectors, names)))
  gene.table <- lapply(gene.vectors, function(x) {
    zeros.names <- setdiff(allgenes, names(x))
    zeros <- rep(0, length(zeros.names))
    names(zeros) <- zeros.names
    x <- c(x, zeros)
    x[allgenes]
  })
  return(as.data.frame(gene.table))
}

#Calculate entropy
getEntropy <- function(tab, pseudo=0.1) {
  p <- (tab+pseudo) / sum(tab+pseudo)
  entropy <- -sum(p * log2(p), na.rm = TRUE)
  return(entropy)
}

#Find highly variable genes in a list of Seurat objects
findHVG <- function(obj.list, nfeatures=2000,
                    min.exp=0.01, max.exp=3.0, hvg.blocklist=NULL) {
  obj.list <- lapply(obj.list, function(x){
    
    ncalc <- min(5*nfeatures, nrow(x))
    x <- findVariableFeatures_wfilters(x, nfeatures=ncalc, min.exp=min.exp, max.exp=max.exp,
                                       genesBlockList=hvg.blocklist)
    x
  })
  hvg <- Seurat::SelectIntegrationFeatures(obj.list,nfeatures = nfeatures,
                                           verbose = FALSE)
  return(hvg)
}

#Calculate metrics for meta-programs
get_metaprogram_consensus <- function(nmf.genes=NULL,
                                      nprograms=10,
                                      min.confidence=0,
                                      max.genes=200,
                                      comb.function=median,
                                      cl_members=NULL) {
  
  markers.consensus <- lapply(seq(1, nprograms), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    gene.vectors <- nmf.genes[which.samples]
    gene.table <- geneList2table(gene.vectors)
    
    genes.avg <- apply(gene.table, 1, comb.function)
    genes.confidence <- apply(gene.table, 1, function(x){sum(x>0)/ length(x)})
    
    genes.avg <- genes.avg[genes.confidence > min.confidence]
    genes.avg <- sort(genes.avg, decreasing = T)
    
    head(genes.avg, min(length(genes.avg), max.genes))
  })
  
  names(markers.consensus) <- paste0("MetaProgram",seq(1,nprograms))
  return(markers.consensus)
}

#Calculate metrics for meta-programs
get_metaprogram_metrics <- function(J=NULL, Jdist=NULL,
                                   markers.consensus=NULL,
                                   cl_members=NULL) {
  nprograms <- length(markers.consensus)
  all.samples <- unique(gsub("\\.k\\d+\\.p\\d+","",colnames(J)))
  sample.coverage <- lapply(seq(1, nprograms), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    ss <- gsub("\\.k\\d+\\.p\\d+","",which.samples)
    ss <- factor(ss, levels=all.samples)
    ss.tab <- table(ss)
    #Percent samples represented
    sum(ss.tab>0)/length(ss.tab)
  })
  names(sample.coverage) <- paste0("MetaProgram",seq(1,nprograms))
  
  #calculate MP silhouettes
  sil <- cluster::silhouette(cl_members, dist=Jdist)
  sil.widths <- summary(sil)$clus.avg.widths
  names(sil.widths) <- paste0("MetaProgram",seq(1,nprograms))
  
  #calculate MP internal average similarity
  clusterSim <- rep(NA,nprograms)
  for(i in seq_len(nprograms)){
    selectMP <- which(cl_members==i)
    if (length(selectMP) > 1) { #needs at least two values
      selectJ <- J[selectMP,selectMP]
      value <- round(mean(selectJ[upper.tri(selectJ)]),3)
    } else {
      value <- 0
    }
    clusterSim[i] <- value
  }
  #number of genes in each meta-program
  metaprograms.length <- unlist(lapply(markers.consensus,length))
  
  #number of programs in meta-program
  metaprograms.size <- as.character(table(cl_members))
  
  metaprograms.metrics <- data.frame(
    sampleCoverage=unlist(sample.coverage),
    silhouette=sil.widths,
    meanSimilarity=clusterSim,
    numberGenes=metaprograms.length,
    numberPrograms=metaprograms.size)
  
  rownames(metaprograms.metrics) <- paste0("MetaProgram",seq(1,nprograms))
  
  return(metaprograms.metrics)
}

#Split positive and negative components of PCA, and reorder by variance
nonNegativePCA <- function(pca, k) {
  
  pcaP <- pca$rotation
  pcaN <- pca$rotation
  colnames(pcaP) <- paste0(colnames(pcaP),"p")
  colnames(pcaN) <- paste0(colnames(pcaN),"n")
  
  pcaP[pcaP<0] <- 0
  pcaN[pcaN>0] <- 0
  pcaN <- abs(pcaN)
  
  sumP <- apply(pcaP, 2, sum)
  sumN <- apply(pcaN, 2, sum)
  
  wP <- pca$sdev * sumP/(sumP+sumN)
  wN <- pca$sdev * sumN/(sumP+sumN)
  wSort <- sort(c(wP,wN), decreasing = T)
  
  #Collate, re-rank components, and rescale coefficients
  pca_abs <- cbind(pcaP, pcaN)
  pca_abs <- pca_abs[,names(wSort)[1:k]]
#  pca_abs <- apply(pca_abs, 2, function(x){x/sum(x)})
  return(pca_abs)
} 

#Weighting factor matrix by feature specificity
weightedLoad <- function(matrix, w) {
  rownorm <- apply(matrix, 1, normVector)
  spec <- apply(rownorm, 2, max)
  spec.w <- spec^w
  matrix <- matrix * spec.w
  #renormalize
  apply(matrix, 2, normVector)
}

normVector <- function(vector) {
  s <- sum(vector)
  if (s>0) {
    vector <- vector/s
  }
  return(vector)
}
