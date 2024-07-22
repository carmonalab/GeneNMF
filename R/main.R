#' Run NMF on a list of Seurat objects
#'
#' Given a list of Seurat objects, run non-negative matrix factorization on 
#' each sample individually, over a range of target NMF components (k).
#'
#' @param obj.list A list of Seurat objects
#' @param k Number of target components for NMF (can be a vector)
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param hvg List of pre-calculated variable genes to subset the matrix.
#'     If hvg=NULL it calculates them automatically
#' @param nfeatures Number of HVG, if calculate_hvg=TRUE
#' @param center Whether to center the data matrix
#' @param scale Whether to scale the data matrix
#' @param hvg.blocklist Optionally takes a vector or list of vectors of gene
#'     names. These genes will be ignored for HVG detection. This is useful
#'     to mitigateeffect of genes associated with technical artifacts and
#'     batch effects (e.g. mitochondrial), and to exclude TCR and BCR 
#'     adaptive immune(clone-specific) receptors. If set to `NULL` no genes 
#'     will be excluded
#' @param L1 L1 regularization term for NMF
#' @param min.cells.per.sample Minimum numer of cells per sample (smaller 
#'     samples will be ignored)
#' @param min.exp Minimum average log-expression value for retaining genes
#' @param max.exp Maximum average log-expression value for retaining genes

#' @param seed Random seed     
#'     
#' @return Returns a list of NMF programs, one for each sample and for each
#'     value of 'k'. The format of each program in the list follosw the
#'     structure of \code{\link[RcppML]{nmf}} factorization models.
#'
#' @examples
#' library(Seurat)
#' data(sampleObj)
#' geneNMF_programs <- multiNMF(list(sampleObj), k=5)
#' 
#' @importFrom RcppML nmf
#' @export  
multiNMF <- function(obj.list, assay="RNA", slot="data", k=5:6,
                   hvg=NULL, nfeatures = 2000, L1=c(0,0),
                   min.exp=0.01, max.exp=3.0,
                   center=FALSE, scale=FALSE,
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
                      hvg=hvg, center=center,
                      scale=scale)
    
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
  nmf.res <- unlist(nmf.res, recursive = FALSE)
  
  return(nmf.res)
}  

#' Run PCA on a list of Seurat objects
#'
#' Given a list of Seurat objects, run non-negative PCA factorization on 
#' each sample individually.
#'
#' @param obj.list A list of Seurat objects
#' @param k Number of target components for PCA
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param hvg List of pre-calculated variable genes to subset the matrix.
#'     If hvg=NULL it calculates them automatically
#' @param nfeatures Number of HVG, if calculate_hvg=TRUE
#' @param center Whether to center the data matrix
#' @param scale Whether to scale the data matrix
#' @param hvg.blocklist Optionally takes a vector or list of vectors of gene
#'     names. These genes will be ignored for HVG detection. This is useful
#'     to mitigateeffect of genes associated with technical artifacts and
#'     batch effects (e.g. mitochondrial), and to exclude TCR and BCR 
#'     adaptive immune(clone-specific) receptors. If set to `NULL` no genes 
#'     will be excluded
#' @param min.cells.per.sample Minimum numer of cells per sample (smaller 
#'     samples will be ignored)
#' @param min.exp Minimum average log-expression value for retaining genes
#' @param max.exp Maximum average log-expression value for retaining genes

#' @param seed Random seed     
#'     
#' @return Returns a list of non-negative PCA programs, one for each sample.
#'     The format of each program in the list follows the
#'     structure of \code{\link[RcppML]{nmf}} factorization models.
#'
#' @examples
#' library(Seurat)
#' data(sampleObj)
#' geneNMF_programs <- multiPCA(list(sampleObj), k=5)
#' 
#' @importFrom irlba prcomp_irlba
#' @export  

multiPCA <- function(obj.list, assay="RNA", slot="data", k=4:5,
                     hvg=NULL, nfeatures = 500,
                     min.exp=0.01, max.exp=3.0,
                     min.cells.per.sample = 10,
                     center=FALSE, scale=FALSE,
                     hvg.blocklist=NULL, seed=123) {
  
  set.seed(seed)
  
  #exclude small samples
  nc <- lapply(obj.list, ncol)
  obj.list <- obj.list[nc > min.cells.per.sample]
  
  if (is.null(hvg)) {
    hvg <- findHVG(obj.list, nfeatures=nfeatures,
                             min.exp=min.exp, max.exp=max.exp, hvg.blocklist=hvg.blocklist)
  }
  
  #run PCA by sample
  pca.res <- lapply(obj.list, function(this) {
    
    mat <- getDataMatrix(obj=this, assay=assay, slot=slot,
                                   hvg=hvg, center=center,
                                   scale=scale, non_negative = FALSE)
    res.k <- lapply(k, function(k.this) {
      
      pca <- prcomp_irlba(t(as.matrix(mat)), center=F, scale.=F, n=k.this)
      rownames(pca$rotation) <- rownames(mat)
      
      nn_pca <- nonNegativePCA(pca, k=k.this) 
      
      npca.obj <- list(w = nn_pca, h = NULL)
      npca.obj
    })
    names(res.k) <- paste0("k",k)
    res.k
  })
  pca.res <- unlist(pca.res, recursive = FALSE)
  
  return(pca.res)
}  


#' Get list of genes for each NMF program
#'
#' Run it over a list of NMF models obtained using \code{multiNMF()}
#'
#' @param nmf.res A list of NMF models obtained using \code{multiNMF()}
#' @param specificity.weight A parameter controlling how specific gene
#'     should be for each program. `specificity.weight=0` no constraint on
#'     specificity, and positive values impose increasing specificity.
#' @param weight.explained Fraction of NMF weights explained by selected
#'     genes. For example if weight.explained=0.5, all genes that together
#'     account for 50\% of NMF weights are used to return program signatures.
#' @param max.genes Max number of genes for each program 
#'     
#' @return Returns a list of top genes for each gene program found
#'     by \code{multiNMF()}
#'
#' @examples
#' library(Seurat)
#' data(sampleObj)
#' geneNMF_programs <- multiNMF(list(sampleObj), k=5)
#' geneNMF_genes <- getNMFgenes(geneNMF_programs)
#' 
#' @importFrom utils head
#' @export
  
getNMFgenes <- function(nmf.res,
                        specificity.weight=5,
                        weight.explained=0.5,
                        max.genes=200) {
  
  
  if (!is.null(specificity.weight)) {
    nmf.res <- weightedLoadings(nmf.res, specificity.weight=specificity.weight)
  }
  
  nmf.genes <- lapply(nmf.res, function(model) {
    
    gene.pass <- apply(model, 2, function(x) {
      weightCumul(x, weight.explained = weight.explained)
    })
    
    m <- lapply(gene.pass, function(g) {
      head(g, min(length(g), max.genes))
    })
    
    #drop empty programs
    isna <- lapply(m, function(x) {all(is.na(x))})
    m <- m[!as.numeric(isna)]
    
    names(m) <- seq(1,length(m))
    m
  })
  
  nmf.genes <- unlist(nmf.genes, recursive = FALSE)
  return(nmf.genes)
}

#' Extract consensus gene programs (meta-programs)
#'
#' Run it over a list of NMF models obtained using \code{\link{multiNMF}}; it will 
#' determine gene programs that are consistently observed across samples 
#' and values of k.
#'
#' @param nmf.res A list of NMF models obtained from \code{\link{multiNMF}}
#' @param nMP Total number of meta-programs
#' @param metric Metric to calculate pairwise similarity between programs     
#' @param max.genes Max number of genes for each programs
#' @param hclust.method Method to build similarity tree between individual programs
#' @param specificity.weight A parameter controlling how specific gene
#'     should be for each program. `specificity.weight=0` no constraint on
#'     specificity, and positive values impose increasing specificity.
#' @param weight.explained Fraction of NMF weights explained by selected
#'     genes. For example if weight.explained=0.5, all genes that together
#'     account for 50\% of NMF weights are used to return program signatures.
#' @param min.confidence Percentage of programs in which a gene is seen 
#'      (out of programs in the corresponding program tree branch/cluster), to be 
#'      retained in the consensus metaprograms
#' @param remove.empty Whether to remove meta-programs with no genes above 
#' confidence threshold      
#' @return Returns a list with i) 'metaprograms.genes' top genes for each 
#'     meta-program; ii) 'metaprograms.metrics' dataframe with meta-programs 
#'     statistics: a) freq. of samples where the MP is present, b) average 
#'     silhouette width, c) mean similarity (cosine or Jaccard), d) number of genes in MP, 
#'     e) number of gene programs in MP; iii) 'programs.similarity': matrix of 
#'     similarities (Jaccard or cosine) between meta-programs; iv) 'programs.tree': 
#'     hierarchical clustering of meta-programs (hclust tree); v) 
#'     'programs.clusters': meta-program identity for each program
#'
#' @examples
#' library(Seurat)
#' data(sampleObj)
#' geneNMF_programs <- multiNMF(list(sampleObj), k=5)
#' geneNMF_metaprograms <- getMetaPrograms(geneNMF_programs, nMP=3)
#' 
#' @importFrom stats cutree dist as.dist hclust
#' @importFrom cluster silhouette
#' @importFrom lsa cosine
#' @importFrom utils head
#' @export  
getMetaPrograms <- function(nmf.res,
                            nMP=10,
                            specificity.weight=5,
                            weight.explained=0.5,
                            max.genes=200,
                            metric = c("cosine","jaccard"),
                            hclust.method="ward.D2",
                            min.confidence=0.5,
                            remove.empty=TRUE) {
  
  metric = metric[1]
  
  nmf.wgt <- weightedLoadings(nmf.res=nmf.res,
                           specificity.weight=specificity.weight)
  
  if (metric == "cosine") {
    J <- cosineSimilarity(geneList2table(nmf.wgt))  
  } else if (metric == "jaccard") {
    nmf.genes <- getNMFgenes(nmf.res=nmf.wgt,
                             specificity.weight=NULL, #because it was precalculated
                             weight.explained=weight.explained,
                             max.genes=max.genes) 
    J <- jaccardSimilarity(nmf.genes)
  } else {
    stop("Unknown distance metric.")
  }
  Jdist <- as.dist(1-J)
  #Cluster programs by gene overlap
  tree <- hclust(Jdist, method=hclust.method)
  cl_members <- cutree(tree, k = nMP)
  
  #Get consensus markers for MPs
  markers.consensus <- get_metaprogram_consensus(nmf.wgt=nmf.wgt,
                                                 nMP=nMP,
                                                 min.confidence=min.confidence,
                                                 weight.explained=weight.explained,
                                                 max.genes=max.genes,
                                                 cl_members=cl_members)
  #Get meta-program metrics
  metaprograms.metrics <- get_metaprogram_metrics(J=J, Jdist=Jdist,
                                                  markers.consensus=markers.consensus,
                                                  cl_members=cl_members)
  
  #Remove any empty meta-program
  if (remove.empty) {
    keep <- metaprograms.metrics$numberGenes > 0
    metaprograms.metrics <- metaprograms.metrics[keep,]
    markers.consensus <- markers.consensus[keep]
    if (sum(!keep)>0) {
      message(sprintf("Dropped %i empty meta-programs", sum(!keep)))
    }
  }
  
  #Rename meta-programs after sorting them by metrics
  ord <- order(metaprograms.metrics$sampleCoverage,
        metaprograms.metrics$silhouette,
        decreasing = T)
  old.names <- names(markers.consensus)[ord]
  new.names <- paste0("MP",seq_along(ord))
  
  markers.consensus <- markers.consensus[ord]
  names(markers.consensus) <- new.names
  metaprograms.metrics <- metaprograms.metrics[ord,]
  rownames(metaprograms.metrics) <- new.names
  
  map.index <- seq_along(old.names)
  names(map.index) <- as.numeric(gsub("MetaProgram","",old.names))
  cl_members.new <- map.index[as.character(cl_members)]

  names(cl_members.new) <- names(cl_members)
  
  #weights of individual genes in each MP
  markers.consensus <- lapply(markers.consensus, function(m){m/sum(m)})
  markers.genes <- lapply(markers.consensus, names)
  
  output.object <- list()
  output.object[["metaprograms.genes"]] <- markers.genes
  output.object[["metaprograms.genes.weights"]] <- markers.consensus
  output.object[["metaprograms.metrics"]] <- metaprograms.metrics
  output.object[["programs.similarity"]] <- J
  output.object[["programs.tree"]] <- tree
  output.object[["programs.clusters"]] <- cl_members.new
  return(output.object)
}  

#' Visualizations for meta-programs
#'
#' Generates a clustered heatmap for meta-program similarities (by Jaccard 
#' index or Cosine similarity). This function is intended to be run on the object
#' generated by \code{\link[GeneNMF]{getMetaPrograms}}, which contains a pre-calculated 
#' tree of pairwise similarities between clusters (as a 'hclust' object).
#'
#' @param mp.res The meta-programs object generated by \code{\link{getMetaPrograms}}
#' @param similarity.cutoff Min and max values for similarity metric
#' @param scale Heatmap rescaling (passed to pheatmap as 'scale')
#' @param palette Heatmap color palette (passed to pheatmap as 'color')
#' @param annotation_colors Color palette for MP annotations 
#' @param main Heatmap title 
#' @param show_rownames Whether to display individual program names as rows
#' @param show_colnames Whether to display individual program names as cols
#' @param ... Additional parameters for pheatmap
#' @return Returns a clustered heatmap of MP similarities, in ggplot2 format
#'
#' @examples
#' library(Seurat)
#' data(sampleObj)
#' geneNMF_programs <- multiNMF(list(sampleObj), k=5)
#' geneNMF_metaprograms <- getMetaPrograms(geneNMF_programs, nMP=3)
#' plotMetaPrograms(geneNMF_metaprograms)
#' 
#' @importFrom pheatmap pheatmap
#' @importFrom viridis viridis
#' @importFrom stats as.dendrogram
#' @export  
plotMetaPrograms <- function(mp.res,
                            similarity.cutoff=c(0,1),
                            scale = "none",
                            palette = viridis(100, option="A", direction=-1),
                            annotation_colors = NULL,
                            main = "Clustered Heatmap",
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            ...) {

  J <- mp.res[["programs.similarity"]]
  tree <- mp.res[["programs.tree"]]
  cl_members <- mp.res[["programs.clusters"]]

  cl_names <- names(cl_members)
  cl_members[!is.na(cl_members)] <- paste0("MP",cl_members[!is.na(cl_members)])
  names(cl_members) <- cl_names

  #Recover order of MP clusters
  labs.order <- labels(as.dendrogram(tree))
  cluster.order <- unique(cl_members[labs.order])
  nMP <- length(cluster.order)
  
  #Annotation column
  annotation_col <- as.data.frame(cl_members)
  colnames(annotation_col) <- "Metaprogram"
  annotation_col[["Metaprogram"]] <- factor(cl_members, levels=cluster.order)

  #Apply trimming to similarity for plotting  
  J[J<similarity.cutoff[1]] <- similarity.cutoff[1]
  J[J>similarity.cutoff[2]] <- similarity.cutoff[2]
  
  ph <- pheatmap(J,
                 scale = scale,
                 color = palette,
                 main = main,
                 cluster_rows = tree,
                 cluster_cols = tree,
                 cutree_rows = nMP,
                 cutree_cols = nMP,
                 annotation_col = annotation_col,
                 annotation_row = annotation_col,
                 annotation_colors = annotation_colors,
                 annotation_names_col = FALSE,
                 annotation_names_row = FALSE,
                 show_rownames = show_rownames,
                 show_colnames = show_colnames,
                 ...
  )
  return(ph)
}  

#' Run Gene set enrichment analysis
#'
#' Utility function to run Gene Set Enrichment Analysis (GSEA) against gene 
#' sets from MSigDB. Note: this is an optional function, which is conditional
#' to the installation of suggested packages \code{fgsea} and \code{msigdbr}.
#'
#' @param genes A vector of genes
#' @param universe Background universe of gene symbols (passed on to \code{fgsea::fora})
#' @param category GSEA main category (e.g. "H" or "C5")
#' @param subcategory GSEA subcategory
#' @param species Species for GSEA analysis. For a list of the available species,
#'     type \code{msigdbr::msigdbr_species()}
#' @param pval.thr Min p-value to include results
#' @return Returns a table of enriched gene programs from GSEA
#'
#' @examples
#' data(sampleObj)
#' geneset <- c("BANK1","CD22","CD79A","CD19","IGHD","IGHG3","IGHM")
#' #test is conditional on availability of suggested packages
#' if (requireNamespace("fgsea", quietly=TRUE) &
#'    requireNamespace("msigdbr", quietly=TRUE)) {
#'    gsea_res <- runGSEA(geneset,
#'        universe=rownames(sampleObj),
#'        category = "C8")
#' }
#' 
#' @export  
runGSEA <- function(genes, universe=NULL,
                    category="H", subcategory=NULL,
                    species="Homo sapiens",
                    pval.thr=0.05) {
  
  
  if (!requireNamespace("fgsea", quietly = TRUE) |
      !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Function 'runGSEA' requires the 'fgsea' and 'msigdbr' packages.
            Please install them.", call. = FALSE)
  }  
  
  if (any(duplicated(genes))) {
    genes <- genes[!duplicated(genes)]
  }
  
  msig_df <- msigdbr::msigdbr(species = species, category = category, subcategory=subcategory)
  msig_list <- split(x=msig_df$gene_symbol, f=msig_df$gs_name)
  
  fgRes <- fgsea::fora(pathways = msig_list,
                       genes = genes,
                       universe = universe)
  
  fgRes <- fgRes[fgRes$pval <= pval.thr,]
  return(fgRes)
}


#' Compute NMF as a low-dim embedding for Seurat
#'
#' Compute NMF embeddings for single-cell dataset, and store them in the Seurat data
#' structure. They can be used as an alternative to PCA for downstream analyses.
#' 
#' @param obj A seurat object
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param k Number of components for low-dim representation
#' @param hvg Which genes to use for the reduction
#' @param new.reduction Name of new dimensionality reduction
#' @param center Whether to center the data matrix
#' @param scale Whether to scale the data matrix
#' @param L1 L1 regularization term for NMF
#' @param seed Random seed
#' @return Returns a Seurat object with a new dimensionality reduction (NMF)
#'
#' @examples
#' data(sampleObj)
#' sampleObj <- runNMF(sampleObj, k=8)
#' @importFrom RcppML nmf
#' @importFrom methods new
#' @export  
runNMF <- function(obj, assay="RNA", slot="data", k=10,
                   new.reduction="NMF", seed=123,
                   L1=c(0,0), hvg=NULL,
                   center=FALSE, scale=FALSE) {
  
  
  set.seed(seed)
  
  if (is.null(hvg)) {
    hvg <- VariableFeatures(obj, assay=assay)
  }
  if (is.null(hvg) | length(hvg)==0) {
    stop("No variable features found. Please run FindVariableFeatures() or specify genes with 'hvg' parameter")
  }
  
  mat <- getDataMatrix(obj=obj, assay=assay, slot=slot,
                    hvg=hvg, center=center, scale=scale)
  
  model <- RcppML::nmf(mat, k = k, L1 = L1, verbose=FALSE, seed = seed)
  
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
#' @param center Whether to center the data matrix
#' @param scale Whether to scale the data matrix
#' @param non_negative Enforce non-negative values for NMF
#'     
#' @return Returns a sparse data matrix (cells per genes), subset 
#' according to the given parameters
#'
#' @examples
#' data(sampleObj)
#' matrix <- getDataMatrix(sampleObj)
#' 
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix t
#' @export  
getDataMatrix <- function(obj, assay="RNA", slot="data", hvg=NULL,
                          center=FALSE, scale=FALSE,
                          non_negative=TRUE) {
  
  mat <- GetAssayData(obj, assay=assay, layer=slot)
  
  #subset on HVG
  if (!is.null(hvg)) mat <- mat[hvg,]
  
  #Center and rescale
  mat <- t(scale(Matrix::t(mat), center=center, scale=scale))
  
  #check scaling with rowsum=0 (gives NaN)
  if (scale) {
    mat[is.na(mat)] <- 0
  }
  if (non_negative) {
    mat[mat<0] <- 0
  }
  return(mat)
}

#' Find variable features
#'
#' Select highly variable genes (HVG) from an expression matrix. Genes from a blocklist
#' (e.g. cell cycling genes, mitochondrial genes) can be excluded from the list of
#' variable genes, as well as genes with very low or very high average expression
#'
#' @param obj A Seurat object containing an expression matrix
#' @param nfeatures Number of top HVG to be returned
#' @param genesBlockList Optionally takes a vector or list of vectors of gene names.
#'     These genes will be ignored for HVG detection. This is useful to mitigate effect
#'     of genes associated with technical artifacts or batch effects
#'     (e.g. mitochondrial, heat-shock response). If set to `NULL` no genes will be excluded
#' @param min.exp Minimum average normalized expression for HVG. If lower, the gene will be excluded
#' @param max.exp Maximum average normalized expression for HVG. If higher, the gene will be excluded
#' @return Returns the input Seurat object \code{obj} with the calculated highly 
#'     variable features accessible through \code{VariableFeatures(obj)} 
#' @examples
#' data(sampleObj)
#' sampleObj <- findVariableFeatures_wfilters(sampleObj, nfeatures=100)
#' 
#' @export
#' @import Seurat
#' 
findVariableFeatures_wfilters <- function(
    obj,
    nfeatures=2000,
    genesBlockList=NULL,
    min.exp=0.01,
    max.exp=3)
{
  
  assay <- DefaultAssay(obj)
  #Calculate a fixed number of HVG, then filtered to nfeat at the end
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 10000, verbose=FALSE)
  
  varfeat <- VariableFeatures(obj)
  
  if (is.vector(genesBlockList)) {
    genes.block <- genesBlockList # user-provided vector
  } else {
    genes.block <- NULL # No excluded genes
  }
  
  varfeat <- setdiff(varfeat, genes.block)
  #Also remove genes that are very poorly or always expressed
  means <- apply(GetAssayData(obj, assay=assay, slot="data")[varfeat,], 1, mean)
  removeGenes2 <- names(means[means<min.exp | means>max.exp])
  
  varfeat <- setdiff(varfeat, removeGenes2)
  n <- min(length(varfeat), nfeatures)
  
  VariableFeatures(obj) <- varfeat[1:n]
  
  return(obj)
}  

