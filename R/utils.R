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

#Downsample vector, while keeping at least one element per class
downsampleMin <- function(vector, size=500) {
  if (length(vector) <= size) {
    return(vector)
  }
  vecout <- names(vector)[!duplicated(vector)] #first match
  vector <- names(vector)
  vector <- setdiff(vector, vecout)
  if (length(vecout) >= size) {
    return(vecout)
  }
  ss <- sample(vector, size = size - length(vecout))
  return(c(vecout, ss))
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
    gene.table <- GeneNMF:::geneList2table(nmf.wgt)[,which.samples]
    
    genes.avg <- apply(as.matrix(gene.table), 1, function(x){
      mean <- mean(x)
      if (length(x) >=3) { #remove outliers (SD only with 3 or more points)
        sd <- sd(x)
        x.out <- x[x>mean-2*sd & x<mean+2*sd]  
      } else {
        x.out <- mean
      }
      
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

#In which samples was a MP detected?
get_metaprogram_composition <- function(J=NULL,
                                    markers.consensus=NULL,
                                    cl_members=NULL) {
  nMP <- length(markers.consensus)
  all.samples <- unique(gsub("\\.k\\d+\\.\\d+","",colnames(J)))
  sample.comp <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    ss <- gsub("\\.k\\d+\\.\\d+","",which.samples)
    ss <- factor(ss, levels=all.samples)
    table(ss)
  })
  names(sample.comp) <- paste0("MetaProgram",seq(1,nMP))
  composition <- do.call(rbind, sample.comp)
  return(composition)
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



getCoph <- function(data, k_rank=2:15, n_subsets = 100, subset_size = 0.5) {
  # data：expr mat
  # k_rank：rang of k to serch
  # n_subsets：number of subsets
  ## requirements:
  # get a cpp function to update consensus matrix
source(NMF)
source(Rcpp)
source(Matrix)
source(RcppML)
source(ProjectSVR) 
source(NMF)
sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix update_consensus_matrix_cpp(NumericVector subset_indices, NumericMatrix consensus_matrix, CharacterVector cluster_assignment) {
int n = subset_indices.size();

for (int j = 0; j < n; ++j) {
  int idx_j = subset_indices[j] - 1; // R索引从1开始，C++索引从0开始
  for (int l = j; l < n; ++l) {
    int idx_l = subset_indices[l] - 1;
    if (cluster_assignment[j] == cluster_assignment[l]) {
      consensus_matrix(idx_j, idx_l) += 1;
      consensus_matrix(idx_l, idx_j) += 1;
    }
  }
}

return consensus_matrix;
}
')
  # body:
  coph<-NULL
  for(k in k_rank){
    n <- ncol(data)  
    consensus_matrix <- matrix(0, n, n)  
    
    for (i in 1:n_subsets) {
      set.seed(10+i) 
      subset_indices <- sample(1:n, size = ceiling(n * subset_size))
      subset_data <- data[,subset_indices]
      nmf_result <- RcppML::nmf(subset_data, k=k, seed = 10+i, tol = 1e-4)
      cppversion <- packageVersion("RcppML")
      if(cppversion >= "0.5.6"){
        nmf_result<-list(w=nmf_result@w,d=nmf_result@d,h=nmf_result@h,tol=nmf_result@misc$tol,iter=nmf_result@misc$iter)
      }
      h <- nmf_result$h 
      cluster_assignment <- row.names(h)[apply(h, 2, which.max)]  # confirmed cluster assignment 
      # update consensus matrix
      consensus_matrix <- update_consensus_matrix_cpp(subset_indices, consensus_matrix, cluster_assignment)
    }
    # calculate coph
    consensus_matrix <- consensus_matrix / n_subsets
    diag(consensus_matrix) <- 1
    coph<-c(coph, NMF::cophcor(consensus_matrix))
  }
  return(coph)
}

optimal.rank <- function(coph){
    coph_diff <- NULL
  for (i in 2:length(coph)){
    coph_diff <- c(coph_diff, coph[i-1]-coph[i])
  }
  k.best <- which.max(coph_diff)+1
  return(k.best)
}

