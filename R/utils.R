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
                                      max.genes=500,
                                      cl_members=NULL) {
  
    markers.consensus <- lapply(seq(1, nprograms), function(c) {
      which.samples <- names(cl_members)[cl_members == c]
      genes <- nmf.genes[which.samples]
      genes.confidence <- sort(table(unlist(genes)), decreasing = T)/(length(which.samples))
      genes.unique <- names(genes.confidence)[genes.confidence > min.confidence]
      head(genes.unique, min(length(genes.unique), max.genes))
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