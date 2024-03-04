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
    x = c(x, zeros)
    x[allgenes]
  })
  return(as.data.frame(gene.table))
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
                                      cl_members=NULL) {
  
  markers.consensus <- lapply(seq(1, nprograms), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    gene.vectors <- nmf.genes[which.samples]
    gene.table <- geneList2table(gene.vectors)
    
    score.avg <- sort(apply(gene.table, 1, mean), decreasing=T)
    head(names(score.avg), min(length(score.avg), max.genes))
  })
  
#  markers.consensus <- lapply(seq(1, nprograms), function(c) {
#    which.samples <- names(cl_members)[cl_members == c]
#    genes <- nmf.genes[which.samples]
#    genes.confidence <- sort(table(unlist(genes)), decreasing = T)/(length(which.samples))
#    genes.unique <- names(genes.confidence)[genes.confidence > min.confidence]
#    head(genes.unique, min(length(genes.unique), max.genes))
#  })
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
  
  #calculate MP internal average Jaccard similarity
  clusterJaccard <- rep(NA,nprograms)
  for(i in seq_len(nprograms)){
    selectMP <- which(cl_members==i)
    if (length(selectMP) > 1) { #needs at least two values
      selectJ <- J[selectMP,selectMP]
      value <- round(mean(selectJ[upper.tri(selectJ)]),3)
    } else {
      value <- 0
    }
    clusterJaccard[i] <- value
  }
  #number of genes in each meta-program
  metaprograms.length <- unlist(lapply(markers.consensus,length))
  
  #number of programs in meta-program
  metaprograms.size <- as.character(table(cl_members))
  
  metaprograms.metrics <- data.frame(
    sampleCoverage=unlist(sample.coverage),
    silhouette=sil.widths,
    meanJaccard=clusterJaccard,
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
weightedLoad <- function(matrix) {
  spec <- apply(matrix, 1, function(x){max(x/sum(x))})
  matrix * spec
}
