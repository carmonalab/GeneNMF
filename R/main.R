#' RunNMF on a list of Seurat objects
#'
#' Get the gene expression matrix from a Seurat object, optionally centered
#' and/or subset on highly variable genes
#'
#' @param obj.list A list of Seurat objects
#' @param k Number of target components for NMF (can be a vector)
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param calculate_hvg Whether to compute highly variable features
#' @param hvg List of pre-calculated variable genes to subset the matrix.
#'     If hvg=NULL it calculates them using 
#' @param nfeatures Number of HVG, if calculate_hvg=TRUE
#' @param do_centering Whether to center the data matrix
#' @param hvg.blocklist Optionally takes a vector or list of vectors of gene names. These genes will be ignored for HVG detection. This is useful to mitigate effect of genes associated with technical artifacts and batch effects (e.g. mitochondrial), and to exclude TCR and BCR adaptive immune (clone-specific) receptors.
#'  If set to `NULL` no genes will be excluded
#' @param L1 L1 regularization term for NMF
#' @param min.cells.per.sample Minimum numer of cells per sample (smaller samples will be ignored)
#' @param min.exp Minimum average log-expression value for retaining genes
#' @param max.exp Maximum average log-expression value for retaining genes

#' @param seed Random seed     
#'     
#' @return Returns a list of NMF results
#'
#' @examples
#' # NMF_results <- multiNMF(obj.list)
#' 
#' @importFrom RcppML nmf
#' @export  

multiNMF <- function(obj.list, assay="RNA", slot="data", k=5:6,
                   hvg=NULL, nfeatures = 2000, L1=c(0,0),
                   min.exp=0.01, max.exp=3.0, do_centering=TRUE,
                   min.cells.per.sample = 10,
                   hvg.blocklist=NULL, seed=123) {
  
  set.seed(seed)
  
  #exclude small samples
  nc <- lapply(obj.list, ncol)
  obj.list <- obj.list[nc > min.cells.per.sample]
  
  if (is.null(hvg)) {
    hvg <- findHVG(obj.list, nfeatures=nfeatures,
                   min.exp=min.exp, max.exp=max.exp, hvg.blocklist=hvg.blocklist)
  }
  
  #run NMF by sample and k
  nmf.res <- lapply(obj.list, function(this) {
    
    mat <- getDataMatrix(obj=this, assay=assay, slot=slot,
                      hvg=hvg, do_centering=do_centering)
    
    res.k <- lapply(k, function(k.this) {
      
      model <- RcppML::nmf(mat, k = k.this, L1 = L1, verbose=FALSE, seed=seed)
      
      rownames(model$h) <- paste0("pattern",1:nrow(model$h))
      colnames(model$h) <- colnames(mat)
      rownames(model$w) <- rownames(mat)
      colnames(model$w) <- paste0("pattern",1:ncol(model$w))
      model
    })
    names(res.k) <- paste0("k",k)
    res.k
  })
  nmf.res <- unlist(nmf.res, recursive = F)
  
  return(nmf.res)
}  

#' Get list of genes for each NMF program
#'
#' Run it over a list of NMF models obtained using `multiNMF`
#'
#' @param nmf.res A list of NMF models obtained from `multiNMF`
#' @param method Parameter passed to `NMF::extractFeatures` to obtain top genes
#'     for each program. When 'method' is a number between 0 and 1, it indicates
#'     the minimum relative basis contribution above which the feature is
#'     selected, i.e. how specific is a gene for a given program.
#' @param max.genes Max number of genes for each programs     
#' @return Returns a list of top genes for each gene program
#'
#' @examples
#' # nmf_genes <- getNMFgenes(nmf_results.list)
#' 
#' @importFrom NMF extractFeatures
#' @export  

getNMFgenes <- function(nmf.res, method=0.5, max.genes=50) {
  
  nmf.genes <- lapply(nmf.res, function(model) {
    
    emb <- model$h
    load <- model$w
    
    m <- NMF::extractFeatures(load, method=method)
    m <- lapply(m, function(x){
      genes <- rownames(load)[x]
      head(genes, min(length(genes), max.genes))
    })
    
    names(m) <- paste0("p",seq(1,length(m)))
    m
  })
  
  nmf.genes <- unlist(nmf.genes, recursive = F)
  return(nmf.genes)
}

#' Extract consensus gene programs (meta-programs)
#'
#' Run it over a list of NMF models obtained using `multiNMF`; it will determine
#' gene programs that are consistently observed across samples and values of k.
#'
#' @param nmf.res A list of NMF models obtained from `multiNMF`
#' @param method Parameter passed to `NMF::extractFeatures` to obtain top genes
#'     for each program
#' @param max.genes Max number of genes for each programs
#' @param hclust.method Method to build similarity tree between individual programs
#' @param nprograms Total number of consensus programs
#' @param min.confidence Percentage of programs in which a gene is seen (out of programs in the corresponding program tree branch/cluster), to be retained
#'      in the consensus metaprograms
#' @param plot.tree Whether to plot (and return) the Jaccard similarity tree between gene programs
#' @return Returns a list with i) 'metaprograms.genes' top genes for each meta-program; ii) 'metaprograms.metrics' dataframe with meta-programs stats: a) freq. of samples where the MP is present, b) average silhouette width, c) mean Jaccard similarity, and d) number of genes in MP; iii) 'programs.jaccard': matrix of Jaccard similarities between meta-programs; iv) 
#' 'programs.tree': hierarchical clustering of meta-programs (hclust tree); v) 'programs.clusters': meta-program assignment to each program
#'
#' @examples
#' # markers <- getMetaPrograms(nmf_results.list)
#' 
#' @importFrom NMF extractFeatures
#' @importFrom stats cutree
#' @export  

getMetaPrograms <- function(nmf.res, method=0.5, max.genes=50,
                                    hclust.method="ward.D2", nprograms=10,
                                    min.confidence=0.5) {
  
  nmf.genes <- getNMFgenes(nmf.res=nmf.res, method=method, max.genes=max.genes) 
  
  nprogs <- length(nmf.genes)
  J <- matrix(data=0, ncol=nprogs, nrow = nprogs)
  colnames(J) <- names(nmf.genes)
  rownames(J) <- names(nmf.genes)
  
  for (i in 1:nprogs) {
    for (j in 1:nprogs) {
      J[i,j] <- jaccardIndex(nmf.genes[[i]], nmf.genes[[j]])
    }  
  }
  
  Jdist <- as.dist(1-J)
  
  tree <- hclust(Jdist, method=hclust.method)
  
  cl_members <- cutree(tree, k = nprograms)
  
  markers.consensus <- lapply(seq(1, nprograms), function(c) {
    which.samples <- names(cl_members)[cl_members == c]
    genes <- nmf.genes[which.samples]
    genes.confidence <- sort(table(unlist(genes)), decreasing = T)/(length(which.samples))
    genes.unique <- names(genes.confidence)[genes.confidence > min.confidence]
    head(genes.unique, min(length(genes.unique), max.genes))
  })
  names(markers.consensus) <- paste0("MetaProgram",seq(1,nprograms))
  
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
  
  #calcuate metaprogram silhouettes
  sil <- cluster::silhouette(cl_members,dist=Jdist)
  sil.widths <- summary(sil)$clus.avg.widths
  names(sil.widths) <- paste0("MetaProgram",seq(1,nprograms))
  
  #calculate metaprogram internal average Jaccard similarity
  clusterJaccard <- rep(NA,nprograms)
  for(i in seq_len(nprograms)){
    selectMP <- which(cl_members==i)
    selectJ <- J[selectMP,selectMP]
    value <- round(mean(selectJ[upper.tri(selectJ)]),3)
    clusterJaccard[i] <- value
  }
  
  metaprograms.length <- unlist(lapply(markers.consensus,length))
  
  metaprograms.metrics <- data.frame(
             sampleCoverage=unlist(sample.coverage),
             silhouette=sil.widths,
             meanJaccard=clusterJaccard,
             numberGenes=metaprograms.length)
  
  rownames(metaprograms.metrics) <- paste0("MetaProgram",seq(1,nprograms))
  
  output.object <- list()
  output.object[["metaprograms.genes"]] <- markers.consensus
  output.object[["metaprograms.metrics"]] <- metaprograms.metrics
  output.object[["programs.jaccard"]] <- J
  output.object[["programs.tree"]] <- tree
  output.object[["programs.clusters"]] <- cl_members
  return(output.object)
}  

#' Visualizations for meta-programs
#'
#' Generates clustered heatmap and dendrogram for meta-program similarities (Jaccard index)
#'
#' @param mp.res The meta-programs object generated by `getMetaPrograms`
#' @param jaccard.cutoff Min and max values for plotting the Jaccard index
#' @return Returns a list with a clustered heatmap and a dendrogram plots of metaprograms similaritites
#'
#' @examples
#' # nmf_genes <- plotMetaPrograms(mp.res)
#' 
#' @importFrom pheatmap pheatmap
#' @importFrom viridis viridis
#' @export  

plotMetaPrograms <- function(mp.res, method=0.5, max.genes=50,
                            jaccard.cutoff=c(0,0.8), 
                            heatmap.clustering_method = "ward.D2",
                            heatmap.scale = "none",
                            heatmap.color = NULL,
                            heatmap.main = "Clustered Heatmap",
                            ...) {
  
  if (is.null(heatmap.color)) {
    heatmap.color <- viridis(100, option="A", direction=-1)
  }
  
  output.object <- list()
  tree <- mp.res[["programs.tree"]]
  cl_members <- mp.res[["programs.clusters"]]

  suppressPackageStartupMessages(require(dendextend))
  dendro <- as.dendrogram(tree)
  labs.order <- labels(dendro)
  cluster.order <- unique(cl_members[labs.order])
  nprograms <- length(cluster.order)
  
  plot(x = tree, labels =  row.names(tree), cex = 0.3)
  treePlot <- dendextend::rect.dendrogram(tree = dendro, k = nprograms, which = 1:nprograms,
                              border = 1:nprograms, cluster = cl_members, text=cluster.order)
  
  output.object[["tree"]] <- treePlot
  
    
  J <- mp.res[["programs.jaccard"]]
    
  J.plot <- J
  J.plot[J<jaccard.cutoff[1]] <- jaccard.cutoff[1]
  J.plot[J>jaccard.cutoff[2]] <- jaccard.cutoff[2]
  
  ph <- pheatmap(J.plot,
                 clustering_method = heatmap.clustering_method,
                 scale = heatmap.scale,
                 color = heatmap.color,
                 main = heatmap.main,
                 cluster_rows = tree,
                 cluster_cols = tree,
                 cutree_rows = nprograms,
                 cutree_cols = nprograms,
                 ...
  )
  output.object[["heatmap"]] <- ph
  
  return(output.object)

  }  

#' Run Gene set enrichment analysis
#'
#' Utility function to run GSEA over a list of genes and obtain enriched sets.
#'
#' @param genes A vector of genes
#' @param universe Background universe of gene symbols (passed on to `fgsea::fora`)
#' @param category GSEA main category (e.g. "H" or "C5")
#' @param subcategory GSEA subcategory
#' @param pval.thr Min p-value to include results
#' @return Returns a table of enriched gene programs from GSEA
#'
#' @examples
#' # gsea_res <- runGSEA(genevector)
#' 
#' @export  

runGSEA <- function(genes, universe=NULL,
                    category="H", subcategory=NULL,
                    pval.thr=0.05) {
  
  
  if (any(duplicated(genes))) {
    genes <- genes[!duplicated(genes)]
  }
  
  msig_df <- msigdbr::msigdbr(species = "Homo sapiens", category = category, subcategory=subcategory)
  
  #msig_list <- msig_df %>% split(x = .$gene_symbol, f = .$gs_name)
  msig_list <- split(x=msig_df$gene_symbol, f=msig_df$gs_name)
  
  fgRes <- fgsea::fora(pathways = msig_list,
                       genes = genes,
                       universe = universe)
  
  fgRes <- fgRes[fgRes$pval <= pval.thr,]
  return(fgRes)
}


#' Compute NMF as a low-dim embedding for Seurat
#'
#' compute and load NMF embeddings for a single object, and store them in Seurat data
#' structure. They can be used as an alternative to PCA for downstream analyses.
#' 
#' @param obj A seurat object
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param k Number of components for low-dim representation
#' @param hvg Which genes to use for the reduction
#' @param new.reduction Name of new dimensionality reduction
#' @param do_centering Whether to center the data matrix
#' @param exclude_ribo_mito Exclude ribosomal and mitochondrial genes from
#'     data matrix
#' @param L1 L1 regularization term for NMF
#' @param seed Random seed
#' @return Returns a Seurat object with a new dimensionality reduction (NMF)
#'
#' @examples
#' # seurat <- RunNMF(seurat, k=8)
#' @importFrom RcppML nmf
#' @export  
RunNMF <- function(obj, assay="RNA", slot="data", k=10,
                   new.reduction="NMF", seed=123,
                   L1=c(0,0), hvg=NULL,
                   do_centering=TRUE,
                   exclude_ribo_mito=FALSE) {
  
  
  set.seed(seed)
  
  if (is.null(hvg)) {
    hvg <- VariableFeatures(obj, assay=assay)
  }
  if (is.null(hvg) | length(hvg)==0) {
    stop("No variable features found. Please run FindVariableFeatures() or specify genes with 'hvg' parameter")
  }
  
  mat <- getDataMatrix(obj=obj, assay=assay, slot=slot,
                    hvg=hvg, do_centering=do_centering,
                    exclude_ribo_mito=exclude_ribo_mito)
  
  model <- RcppML::nmf(mat, k = k, L1 = L1, verbose=FALSE)
  
  rownames(model$h) <- paste0(new.reduction,"_",1:nrow(model$h))
  colnames(model$h) <- colnames(mat)
  
  rownames(model$w) <- rownames(mat)
  colnames(model$w) <- paste0(new.reduction,"_",1:ncol(model$w))
  
  #New dim reduction
  obj@reductions[[new.reduction]] <- new("DimReduc",
                                         cell.embeddings = t(model$h),
                                         feature.loadings = model$w,
                                         assay.used = assay,
                                         stdev = model$d,
                                         key = paste0(new.reduction,"_"),
                                         global = FALSE)
  return(obj)
}  

#' Extract data matrix from Seurat object
#'
#' Get the gene expression matrix from a Seurat object, optionally centered
#' and/or subset on highly variable genes
#'
#' @param obj Seurat object
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param hvg List of variable genes to subset the matrix. If NULL, uses
#'     all genes
#' @param do_centering Whether to center the data matrix
#' @param exclude_ribo_mito Exclude ribosomal and mitochondrial genes from
#'     data matrix
#'     
#' @return Returns a data matrix (cells per genes), subset according to
#' the given parameters
#'
#' @examples
#' # matrix <- getDataMatrix(seurat_object)
#' 
#' @importFrom Seurat GetAssayData
#' @export  
getDataMatrix <- function(obj, assay="RNA", slot="data", hvg=NULL, do_centering=TRUE) {
  
  mat <- GetAssayData(obj, assay=assay, layer=slot)
  
  #subset on HVG
  if (!is.null(hvg)) mat <- mat[hvg,]
  
  if (do_centering) {
    mat <- centerData(mat)
  }
  return(mat)
}
