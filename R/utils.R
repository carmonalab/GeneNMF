#Calculate Jaccard Index
jaccardIndex <- function(a, b) {
  i = length(intersect(a, b))
  u = length(a) + length(b) - i
  return (i/u)
}

jaccardSimilarity <- function(gene.vectors) {
  nprogs <- length(gene.vectors)
  all.genes <- unique(unlist(lapply(gene.vectors, names)))
  # binary membership matrix: genes x programs
  mem <- matrix(0L, nrow=length(all.genes), ncol=nprogs,
                dimnames=list(all.genes, names(gene.vectors)))
  for (i in seq_len(nprogs)) {
    mem[names(gene.vectors[[i]]), i] <- 1L
  }
  # intersection counts via cross-product; union = |a| + |b| - intersection
  inter <- crossprod(mem)
  sizes <- colSums(mem)
  union <- outer(sizes, sizes, "+") - inter
  J <- inter / union
  diag(J) <- 1
  return(J)
}

#Calculate cosine similarity
cosineSimilarity <- function(gene.table) {
  cosine.matrix <- lsa::cosine(as.matrix(gene.table))
  return(cosine.matrix)
}

#Downsample vector, while keeping at least one element per class
downsampleMin <- function(vector, size=500, seed=NULL) {
  if (length(vector) <= size) {
    return(vector)
  }
  #reproducible sampling, without altering the global random number stream
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old.seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old.seed, envir = globalenv()))
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()))
    }
    set.seed(seed)
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

#Index of NMF components for the programs (columns) of a loadings matrix
programIds <- function(matrix) {
  ids <- colnames(matrix)
  if (is.null(ids)) {
    ids <- seq_len(ncol(matrix))
  }
  gsub("^pattern", "", ids)
}

#Sample of origin for NMF programs, named as [sample].k[k].[component]
programSample <- function(programs) {
  sub("\\.k\\d+\\.\\d+$", "", programs)
}

#Assign unique names to a list of samples
nameSamples <- function(obj.list) {
  names <- names(obj.list)
  if (is.null(names)) {
    names <- rep("", length(obj.list))
  }
  unnamed <- is.na(names) | names == ""
  names[unnamed] <- paste0("sample", seq_along(obj.list))[unnamed]
  if (any(duplicated(names))) {
    stop("Sample names (names of obj.list) must be unique")
  }
  names(obj.list) <- names
  return(obj.list)
}

#Name of the meta-program each program belongs to. 'programs.clusters' stores
#the index of the MP in the list of meta-programs
clusterNames <- function(mp.res) {
  mp.names <- names(mp.res$metaprograms.genes)
  cl_members <- mp.res$programs.clusters
  if (any(cl_members > length(mp.names), na.rm = TRUE)) {
    #objects from older versions, where the index was not updated on dropping MPs
    cl_names <- ifelse(is.na(cl_members), NA, paste0("MP", cl_members))
  } else {
    cl_names <- mp.names[cl_members]
  }
  names(cl_names) <- names(cl_members)
  return(cl_names)
}

#From list of genes to complete table, with imputed zeros
geneList2table <- function(gene.vectors) {
  names <- names(gene.vectors)
  gene.vectors <- lapply(names, function(n) {
    g <- gene.vectors[[n]]
    colnames(g) <- paste(n,programIds(g),sep=".")
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
findHVG <- function(obj.list, nfeatures=2000, assay=NULL,
                    min.exp=0.01, max.exp=3.0, hvg.blocklist=NULL) {
  obj.list <- lapply(obj.list, function(x){
    
    ncalc <- min(5*nfeatures, nrow(x))
    x <- findVariableFeatures_wfilters(x, nfeatures=ncalc, min.exp=min.exp, max.exp=max.exp,
                                       genesBlockList=hvg.blocklist, assay=assay)
    x
  })
  if (!is.null(assay)) {
    assay <- rep(assay, length(obj.list))
  }
  hvg <- Seurat::SelectIntegrationFeatures(obj.list,nfeatures = nfeatures,
                                           assay = assay, verbose = FALSE)
  return(hvg)
}

#Calculate metrics for meta-programs
get_metaprogram_consensus <- function(nmf.wgt,
                                      nMP=10,
                                      min.confidence=0.5,
                                      weight.explained=0.5,
                                      max.genes=200,
                                      cl_members=NULL,
                                      weight.explained.single=0.9,
                                      max.genes.single=1000) {
  
  #calculate genes that explain 90% of weight in individual programs
  #this is used to calculate gene confidence 
  nmf.genes.single <- getNMFgenes(nmf.res=nmf.wgt,
                                     specificity.weight=NULL,
                                     weight.explained=weight.explained.single,
                                     max.genes=max.genes.single) 
  
  gene.table.all <- geneList2table(nmf.wgt)
  
  markers.consensus <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    gene.table <- gene.table.all[,which.samples]
    
    genes.avg <- apply(as.matrix(gene.table), 1, function(x){
      mean <- mean(x)
      if (length(x) >=3) { #remove outliers (SD only with 3 or more points)
        sd <- sd(x)
        x.out <- x[x>mean-3*sd & x<mean+3*sd]  
      } else {
        x.out <- mean
      }
      
      mean(x.out)
    })
    #first criterion: explain x% of total weight
    genes.avg <- sort(genes.avg, decreasing = T)
    genes.pass <- weightCumul(genes.avg, weight.explained=weight.explained)
    
    #second criterion: consistently detect a given gene across runs
    this <- nmf.genes.single[which.samples]
    genes.only <- lapply(this, names)
    genes.sum <- sort(table(unlist(genes.only)), decreasing=T)
    genes.confidence <- genes.sum/length(this)
    genes.confidence <- genes.confidence[genes.confidence >= min.confidence]
    
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
  all.samples <- unique(programSample(colnames(J)))
  sample.coverage <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    ss <- programSample(which.samples)
    ss <- factor(ss, levels=all.samples)
    ss.tab <- table(ss)
    #Percent samples represented
    sum(ss.tab>0)/length(ss.tab)
  })
  names(sample.coverage) <- paste0("MetaProgram",seq(1,nMP))
  
  #calculate MP silhouettes
  sil <- cluster::silhouette(cl_members, dist=Jdist)
  if (inherits(sil, "silhouette")) {
    sil.widths <- summary(sil)$clus.avg.widths
  } else { #not defined, e.g. for a single cluster
    sil.widths <- rep(NA_real_, nMP)
  }
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
  metaprograms.size <- as.numeric(table(cl_members))
  
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
  all.samples <- unique(programSample(colnames(J)))
  sample.comp <- lapply(seq(1, nMP), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    ss <- programSample(which.samples)
    ss <- factor(ss, levels=all.samples)
    table(ss)
  })
  names(sample.comp) <- paste0("MetaProgram",seq(1,nMP))
  composition <- do.call(rbind, sample.comp)
  return(composition)
}


#Weighting factor matrix by feature specificity
weightedLoadings <- function(nmf.res,
                        specificity.weight=5) {
  
  nmf.wgt <- lapply(nmf.res, function(model) {
    wgt <- wgtLoad(model$w, w=specificity.weight)
    #identify programs by the index of their NMF component
    colnames(wgt) <- seq_len(ncol(wgt))
    wgt
  })
  #drop empty programs (all-zero loadings)
  nempty <- 0
  nmf.wgt <- lapply(nmf.wgt, function(wgt) {
    keep <- colSums(wgt) > 0
    nempty <<- nempty + sum(!keep)
    wgt[, keep, drop=FALSE]
  })
  if (nempty > 0) {
    message(sprintf("Dropped %i empty NMF programs (all-zero loadings)", nempty))
  }
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

#Top elements that together explain a given fraction of total weight. The
#element that crosses the threshold is included, to return at least one element
weightCumul <- function(vector, weight.explained=0.5) {
    x.sorted <- sort(vector, decreasing = T)
    norm.x <- normVector(x.sorted)
    norm.x <- norm.x[norm.x > 0]
    cs <- cumsum(norm.x)
    n <- sum(cs < weight.explained) + 1
    head(norm.x, n)
}

#Input matrix for RcppML::nmf. Factorization of sparse matrices requires the
#Matrix package to be attached; otherwise, use a dense matrix
nmfInput <- function(mat) {
  if (!is.matrix(mat) && !"package:Matrix" %in% search()) {
    mat <- as.matrix(mat)
  }
  return(mat)
}

check_cpp_version <- function(model) {
  cppversion <- packageVersion("RcppML")
  if(cppversion >= "0.5.6"){
    model<-list(w=model@w, 
                d=model@d, 
                h=model@h, 
                tol=model@misc$tol,
                iter=model@misc$iter)
  }
  return(model)
}
