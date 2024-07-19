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
cosineSimilarity <- function(gene.table) {
  cosine.matrix <- lsa::cosine(as.matrix(gene.table))
  return(cosine.matrix)
}

#From list of genes to complete table, with imputed zeros
geneList2table <- function(gene.vectors) {
  names <- names(gene.vectors)
  gene.vectors <- lapply(names, function(n) {
    g <- gene.vectors[[n]]
    colnames(g) <- paste(n,seq(1,ncol(g)),sep=".")
    g
  })
  gene.table <- Reduce(f=cbind, x=gene.vectors)
  return(gene.table)
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
get_metaprogram_consensus <- function(nmf.wgt,
                                      nMP=10,
                                      min.confidence=0.5,
                                      weight.explained=0.5,
                                      max.genes=200,
                                      cl_members=NULL) {
  
  #calculate genes that explain 80% of weight in individual samples
  #this is used to calculate gene confidence 
  nmf.genes.single <- getNMFgenes(nmf.res=nmf.wgt,
                                     specificity.weight=NULL,
                                     weight.explained=0.8,
                                     max.genes=1000) 
  
  markers.consensus <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    gene.table <- geneList2table(nmf.wgt)[,which.samples]
    
    genes.avg <- apply(as.matrix(gene.table), 1, function(x){
      mean <- mean(x)
      sd <- sd(x)
      x.out <- x[x>mean-2*sd & x<mean+2*sd]  #remove outliers
      mean(x.out)
    })
    genes.avg <- sort(genes.avg, decreasing = T)
    genes.pass <- weightCumul(genes.avg, weight.explained=weight.explained)
    
    this <- nmf.genes.single[which.samples]
    genes.only <- lapply(this, names)
    genes.sum <- sort(table(unlist(genes.only)), decreasing=T)
    genes.confidence <- genes.sum/length(this)
    genes.confidence <- genes.confidence[genes.confidence > min.confidence]
    
    genes.pass <- genes.pass[names(genes.pass) %in% names(genes.confidence)]
    
    head(genes.pass, min(length(genes.pass), max.genes))
  })
  
  names(markers.consensus) <- paste0("MetaProgram",seq(1,nMP))
  return(markers.consensus)
}

#Calculate metrics for meta-programs
get_metaprogram_metrics <- function(J=NULL, Jdist=NULL,
                                   markers.consensus=NULL,
                                   cl_members=NULL) {
  nMP <- length(markers.consensus)
  all.samples <- unique(gsub("\\.k\\d+\\.\\d+","",colnames(J)))
  sample.coverage <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    ss <- gsub("\\.k\\d+\\.\\d+","",which.samples)
    ss <- factor(ss, levels=all.samples)
    ss.tab <- table(ss)
    #Percent samples represented
    sum(ss.tab>0)/length(ss.tab)
  })
  names(sample.coverage) <- paste0("MetaProgram",seq(1,nMP))
  
  #calculate MP silhouettes
  sil <- cluster::silhouette(cl_members, dist=Jdist)
  sil.widths <- summary(sil)$clus.avg.widths
  names(sil.widths) <- paste0("MetaProgram",seq(1,nMP))
  
  #calculate MP internal average similarity
  clusterSim <- rep(NA,nMP)
  for(i in seq_len(nMP)){
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
  
  rownames(metaprograms.metrics) <- paste0("MetaProgram",seq(1,nMP))
  
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
weightedLoadings <- function(nmf.res,
                        specificity.weight=5) {
  
  nmf.wgt <- lapply(nmf.res, function(model) {
    wgtLoad(model$w, w=specificity.weight)
  })
  return(nmf.wgt)
}

wgtLoad <- function(matrix, w) {
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

weightCumul <- function(vector, weight.explained=0.5) {
    x.sorted <- sort(vector, decreasing = T)
    cs <- cumsum(x.sorted)
    norm.cs <- normVector(cs)
    norm.cs <- norm.cs/max(norm.cs)
    x.sorted[norm.cs<weight.explained]
}