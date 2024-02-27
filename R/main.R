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
#' @param do_centering Whether to center the data matrix
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
  nmf.res <- unlist(nmf.res, recursive = FALSE)
  
  return(nmf.res)
}  

#' Get list of genes for each NMF program
#'
#' Run it over a list of NMF models obtained using \code{multiNMF()}
#'
#' @param nmf.res A list of NMF models obtained using \code{multiNMF()}
#' @param method Parameter passed to \code{\link[NMF]{extractFeatures}} to
#'     obtain top genes for each program. When 'method' is a number between 0 
#'     and 1, it indicates
#'     the minimum relative basis contribution above which the feature is
#'     selected, i.e. how specific is a gene for a given program.
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
#' @importFrom NMF extractFeatures
#' @export  
getNMFgenes <- function(nmf.res, method=0.5, max.genes=200) {
  
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
#' @param method Parameter passed to \code{\link[NMF]{extractFeatures}} to 
#'     obtain top genes for each program
#' @param max.genes Max number of genes for each programs
#' @param hclust.method Method to build similarity tree between individual programs
#' @param nprograms Total number of meta-programs
#' @param min.confidence Percentage of programs in which a gene is seen 
#'      (out of programs in the corresponding program tree branch/cluster), to be 
#'      retained in the consensus metaprograms
#' @param remove.empty Whether to remove meta-programs with no genes above 
#' confidence threshold      
#' @return Returns a list with i) 'metaprograms.genes' top genes for each 
#'     meta-program; ii) 'metaprograms.metrics' dataframe with meta-programs 
#'     statistics: a) freq. of samples where the MP is present, b) average 
#'     silhouette width, c) mean Jaccard similarity, d) number of genes in MP, 
#'     e) number of gene programs in MP; iii) 'programs.jaccard': matrix of 
#'     Jaccard similarities between meta-programs; iv) 'programs.tree': 
#'     hierarchical clustering of meta-programs (hclust tree); v) 
#'     'programs.clusters': meta-program identity for each program
#'
#' @examples
#' library(Seurat)
#' data(sampleObj)
#' geneNMF_programs <- multiNMF(list(sampleObj), k=5)
#' geneNMF_metaprograms <- getMetaPrograms(geneNMF_programs, nprograms=3)
#' 
#' @importFrom NMF extractFeatures
#' @importFrom stats cutree dist
#' @importFrom cluster silhouette
#' @export  
getMetaPrograms <- function(nmf.res, method=0.5,
                            max.genes=200,
                            hclust.method="ward.D2",
                            nprograms=10,
                            min.confidence=0.2,
                            remove.empty=TRUE) {
  
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
  
  #Cluster programs by gene overlap (J index)
  tree <- hclust(Jdist, method=hclust.method)
  cl_members <- cutree(tree, k = nprograms)
  
  #Get consensus markers for MPs
  markers.consensus <- get_metaprogram_consensus(nmf.genes=nmf.genes,
                                                 nprograms=nprograms,
                                                 min.confidence=min.confidence,
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
  
  output.object <- list()
  output.object[["metaprograms.genes"]] <- markers.consensus
  output.object[["metaprograms.metrics"]] <- metaprograms.metrics
  output.object[["programs.jaccard"]] <- J
  output.object[["programs.tree"]] <- tree
  output.object[["programs.clusters"]] <- cl_members.new
  return(output.object)
}  

#' Visualizations for meta-programs
#'
#' Generates a clustered heatmap for meta-program similarities (by Jaccard 
#' index). This function is intended to be run on the object generated by 
#' \code{\link[GeneNMF]{getMetaPrograms}}, which contains a pre-calculated 
#' tree of pairwise similarities between clusters (as a 'hclust' object).
#'
#' @param mp.res The meta-programs object generated by \code{\link{getMetaPrograms}}
#' @param jaccard.cutoff Min and max values for plotting the Jaccard index
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
#' geneNMF_metaprograms <- getMetaPrograms(geneNMF_programs, nprograms=3)
#' plotMetaPrograms(geneNMF_metaprograms)
#' 
#' @importFrom pheatmap pheatmap
#' @importFrom viridis viridis
#' @importFrom stats as.dendrogram
#' @export  
plotMetaPrograms <- function(mp.res,
                            jaccard.cutoff=c(0,0.8),
                            scale = "none",
                            palette = viridis(100, option="A", direction=-1),
                            annotation_colors = NULL,
                            main = "Clustered Heatmap",
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            ...) {

  J <- mp.res[["programs.jaccard"]]
  tree <- mp.res[["programs.tree"]]
  cl_members <- mp.res[["programs.clusters"]]

  cl_names <- names(cl_members)
  cl_members[!is.na(cl_members)] <- paste0("MP",cl_members[!is.na(cl_members)])
  names(cl_members) <- cl_names

  #Recover order of MP clusters
  labs.order <- labels(as.dendrogram(tree))
  cluster.order <- unique(cl_members[labs.order])
  nprograms <- length(cluster.order)
  
  #Annotation column
  annotation_col <- as.data.frame(cl_members)
  colnames(annotation_col) <- "Metaprogram"
  annotation_col[["Metaprogram"]] <- factor(cl_members, levels=cluster.order)

  #Apply trimming to Jaccard for plotting  
  J[J<jaccard.cutoff[1]] <- jaccard.cutoff[1]
  J[J>jaccard.cutoff[2]] <- jaccard.cutoff[2]
  
  ph <- pheatmap(J,
                 scale = scale,
                 color = palette,
                 main = main,
                 cluster_rows = tree,
                 cluster_cols = tree,
                 cutree_rows = nprograms,
                 cutree_cols = nprograms,
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
#' Utility function to run Gene set enrichment analysis (GSEA) against gene 
#' sets from MSigDB.
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
#' gsea_res <- runGSEA(geneset, universe=rownames(sampleObj), category = "C8")
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
#' @param do_centering Whether to center the data matrix
#' @param L1 L1 regularization term for NMF
#' @param seed Random seed
#' @return Returns a Seurat object with a new dimensionality reduction (NMF)
#'
#' @examples
#' data(sampleObj)
#' sampleObj <- runNMF(sampleObj, k=8)
#' @importFrom RcppML nmf
#' @export  
runNMF <- function(obj, assay="RNA", slot="data", k=10,
                   new.reduction="NMF", seed=123,
                   L1=c(0,0), hvg=NULL,
                   do_centering=TRUE) {
  
  
  set.seed(seed)
  
  if (is.null(hvg)) {
    hvg <- VariableFeatures(obj, assay=assay)
  }
  if (is.null(hvg) | length(hvg)==0) {
    stop("No variable features found. Please run FindVariableFeatures() or specify genes with 'hvg' parameter")
  }
  
  mat <- getDataMatrix(obj=obj, assay=assay, slot=slot,
                    hvg=hvg, do_centering=do_centering)
  
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
#' @param do_centering Whether to center the data matrix
#'     
#' @return Returns a sparse data matrix (cells per genes), subset 
#' according to the given parameters
#'
#' @examples
#' data(sampleObj)
#' matrix <- getDataMatrix(sampleObj)
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

