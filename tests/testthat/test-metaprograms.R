mp.slots <- c("metaprograms.genes", "metaprograms.genes.weights",
              "metaprograms.metrics", "metaprograms.composition",
              "programs.similarity", "programs.tree", "programs.clusters")
mp.names <- paste0("MP", seq_len(test.nMP))
prog.names <- names(getNMFgenes(test.nmf))

test_that("getMetaPrograms returns the expected object structure", {
  expect_type(test.mp, "list")
  expect_named(test.mp, mp.slots)

  expect_named(test.mp$metaprograms.genes, mp.names)
  expect_named(test.mp$metaprograms.genes.weights, mp.names)
  expect_s3_class(test.mp$metaprograms.metrics, "data.frame")
  expect_equal(rownames(test.mp$metaprograms.metrics), mp.names)
  expect_named(test.mp$metaprograms.metrics,
               c("sampleCoverage", "silhouette", "meanSimilarity",
                 "numberGenes", "numberPrograms"))
  expect_equal(rownames(test.mp$metaprograms.composition), mp.names)
  expect_equal(colnames(test.mp$metaprograms.composition), names(test.objs))
  expect_s3_class(test.mp$programs.tree, "hclust")
})

test_that("getMetaPrograms genes and weights are consistent", {
  for (mp in mp.names) {
    genes <- test.mp$metaprograms.genes[[mp]]
    weights <- test.mp$metaprograms.genes.weights[[mp]]

    expect_gt(length(genes), 0)
    expect_true(all(genes %in% test.hvg))
    expect_equal(names(weights), genes)
    expect_equal(sum(weights), 1)
    expect_false(is.unsorted(rev(weights)))
  }
})

test_that("getMetaPrograms metrics are consistent with program clusters", {
  metrics <- test.mp$metaprograms.metrics
  clusters <- test.mp$programs.clusters
  comp <- test.mp$metaprograms.composition

  expect_setequal(names(clusters), prog.names)
  expect_setequal(unique(clusters), seq_len(test.nMP))
  expect_equal(sum(metrics$numberPrograms), length(prog.names))
  expect_equal(metrics$numberPrograms, as.numeric(table(clusters)))
  expect_equal(unname(rowSums(comp)), metrics$numberPrograms)
  expect_equal(metrics$numberGenes, unname(lengths(test.mp$metaprograms.genes)))

  expect_true(all(metrics$sampleCoverage >= 0 & metrics$sampleCoverage <= 1))
  expect_true(all(metrics$silhouette >= -1 & metrics$silhouette <= 1))
  #MPs are sorted by decreasing sample coverage
  expect_false(is.unsorted(rev(metrics$sampleCoverage)))
})

test_that("getMetaPrograms returns a valid similarity matrix", {
  J <- test.mp$programs.similarity

  expect_equal(dim(J), rep(length(prog.names), 2))
  expect_setequal(rownames(J), prog.names)
  expect_equal(unname(diag(J)), rep(1, nrow(J)))
  expect_true(isSymmetric(J))
  expect_true(all(J >= 0 & J <= 1 + 1e-8))
  expect_setequal(test.mp$programs.tree$labels, prog.names)
})

test_that("getMetaPrograms supports the jaccard metric", {
  mp <- getMetaPrograms(test.nmf, nMP = test.nMP, metric = "jaccard")

  expect_named(mp, mp.slots)
  expect_named(mp$metaprograms.genes, mp.names)
  expect_equal(unname(diag(mp$programs.similarity)),
               rep(1, length(prog.names)))
  expect_true(all(mp$programs.similarity >= 0 & mp$programs.similarity <= 1))

  expect_error(getMetaPrograms(test.nmf, nMP = test.nMP, metric = "euclidean"),
               "Unknown distance metric")
})

test_that("getMetaPrograms respects max.genes", {
  mp <- getMetaPrograms(test.nmf, nMP = test.nMP, max.genes = 3)

  expect_true(all(lengths(mp$metaprograms.genes) <= 3))
})

test_that("dropMetaPrograms removes the requested MPs", {
  mp <- dropMetaPrograms(test.mp, dropMP = "MP2")
  kept <- setdiff(mp.names, "MP2")
  dropped.progs <- names(test.mp$programs.clusters)[test.mp$programs.clusters == 2]

  expect_named(mp$metaprograms.genes, kept)
  expect_named(mp$metaprograms.genes.weights, kept)
  expect_equal(rownames(mp$metaprograms.metrics), kept)
  expect_equal(rownames(mp$metaprograms.composition), kept)
  expect_equal(mp$metaprograms.composition,
               test.mp$metaprograms.composition[kept, , drop = FALSE])
  expect_equal(mp$metaprograms.genes[kept], test.mp$metaprograms.genes[kept])

  #'programs.clusters' is an index to the remaining MPs
  expect_setequal(unique(mp$programs.clusters), seq_along(kept))
  expect_equal(clusterNames(mp),
               clusterNames(test.mp)[names(mp$programs.clusters)])
  expect_length(mp$programs.clusters,
                length(test.mp$programs.clusters) - length(dropped.progs))
  expect_false(any(dropped.progs %in% names(mp$programs.clusters)))
  expect_equal(rownames(mp$programs.similarity), names(mp$programs.clusters))
  expect_s3_class(mp$programs.tree, "hclust")
  expect_setequal(mp$programs.tree$labels, names(mp$programs.clusters))
})

test_that("dropMetaPrograms can drop down to a single MP", {
  mp <- dropMetaPrograms(test.mp, dropMP = c("MP2", "MP3"))

  expect_named(mp$metaprograms.genes, "MP1")
  expect_s3_class(mp$metaprograms.metrics, "data.frame")
  expect_equal(rownames(mp$metaprograms.metrics), "MP1")
  expect_equal(dim(mp$metaprograms.composition), c(1, length(test.objs)))
  expect_equal(rownames(mp$metaprograms.composition), "MP1")
})

test_that("dropMetaPrograms cannot drop all MPs", {
  expect_error(dropMetaPrograms(test.mp, dropMP = mp.names), "Cannot drop all")
})

test_that("dropMetaPrograms warns on unknown MPs", {
  expect_warning(mp <- dropMetaPrograms(test.mp, dropMP = "MP99"),
                 "not found in object: MP99")
  expect_named(mp$metaprograms.genes, mp.names)
  expect_equal(mp$programs.clusters, test.mp$programs.clusters)
})

test_that("renameMetaPrograms assigns new names", {
  new.names <- c("Bcell", "NK", "Other")
  mp <- renameMetaPrograms(test.mp, old.names = mp.names, new.names = new.names)

  expect_named(mp, mp.slots)
  expect_named(mp$metaprograms.genes, new.names)
  expect_named(mp$metaprograms.genes.weights, new.names)
  expect_equal(rownames(mp$metaprograms.metrics), new.names)
  expect_equal(rownames(mp$metaprograms.composition), new.names)
  expect_equal(unname(mp$metaprograms.genes), unname(test.mp$metaprograms.genes))
  expect_equal(mp$programs.clusters, test.mp$programs.clusters)
})

test_that("renameMetaPrograms reorders MPs following old.names", {
  ord <- rev(mp.names)
  mp <- renameMetaPrograms(test.mp, old.names = ord, new.names = mp.names)

  expect_named(mp$metaprograms.genes, mp.names)
  expect_equal(mp$metaprograms.genes[["MP1"]],
               test.mp$metaprograms.genes[[ord[1]]])
  expect_equal(mp$metaprograms.metrics$numberGenes,
               rev(test.mp$metaprograms.metrics$numberGenes))
  expect_equal(unname(mp$programs.clusters),
               unname(test.nMP + 1 - test.mp$programs.clusters))
})

test_that("renameMetaPrograms renumbers MPs after dropping", {
  dropped <- dropMetaPrograms(test.mp, dropMP = "MP2")
  mp <- renameMetaPrograms(dropped)

  expect_named(mp$metaprograms.genes, c("MP1", "MP2"))
  expect_equal(mp$metaprograms.genes[["MP2"]],
               test.mp$metaprograms.genes[["MP3"]])
  expect_equal(mp$metaprograms.metrics$numberPrograms,
               test.mp$metaprograms.metrics[c("MP1", "MP3"), "numberPrograms"])
  expect_equal(unname(mp$metaprograms.composition),
               unname(test.mp$metaprograms.composition[c("MP1", "MP3"), , drop = FALSE]))
  expect_equal(unname(rowSums(mp$metaprograms.composition)),
               mp$metaprograms.metrics$numberPrograms)
  expect_setequal(unique(mp$programs.clusters), 1:2)
  expect_equal(names(mp$programs.clusters), names(dropped$programs.clusters))
})

test_that("meta-programs with custom names can be renamed, dropped and plotted", {
  custom <- c("Bcell", "NK", "Other")
  mp <- renameMetaPrograms(test.mp, old.names = mp.names, new.names = custom)

  #rename again
  mp2 <- renameMetaPrograms(mp, old.names = rev(custom), new.names = c("X", "Y", "Z"))
  expect_false(anyNA(mp2$programs.clusters))
  expect_equal(unname(clusterNames(mp2)),
               unname(c(Bcell = "Z", NK = "Y", Other = "X")[clusterNames(mp)]))
  expect_equal(mp2$metaprograms.genes[["Z"]], test.mp$metaprograms.genes[["MP1"]])

  #drop
  dropped <- dropMetaPrograms(mp, dropMP = "NK")
  expect_named(dropped$metaprograms.genes, c("Bcell", "Other"))
  expect_setequal(clusterNames(dropped), c("Bcell", "Other"))
  expect_equal(sum(test.mp$programs.clusters != 2), length(dropped$programs.clusters))
  expect_equal(rownames(dropped$programs.similarity), names(dropped$programs.clusters))
  expect_setequal(dropped$programs.tree$labels, names(dropped$programs.clusters))
  expect_equal(unname(table(clusterNames(dropped))[c("Bcell", "Other")]),
               unname(table(test.mp$programs.clusters)[c("1", "3")]))

  #plot: annotations show the custom names
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  ph <- plotMetaPrograms(dropped, silent = TRUE, show_rownames = TRUE)
  expect_s3_class(ph, "pheatmap")
  legend <- ph$gtable$grobs[[which(ph$gtable$layout$name == "annotation_legend")]]
  labels <- unlist(lapply(legend$children, function(x) x$label))
  expect_true(all(c("Bcell", "Other") %in% labels))
})

test_that("renameMetaPrograms validates its input", {
  expect_error(renameMetaPrograms(test.mp, old.names = mp.names,
                                  new.names = c("A", "A", "B")),
               "new.names must be unique")
  expect_error(renameMetaPrograms(test.mp, old.names = c("MP1", "MP2")),
               "old.names must match")
  expect_error(renameMetaPrograms(test.mp, old.names = mp.names,
                                  new.names = c("A", "B")),
               "must have the same length")
})

test_that("plotMetaPrograms returns a heatmap", {
  #render to a null device, to avoid creating Rplots.pdf
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  ph <- plotMetaPrograms(test.mp, silent = TRUE)
  expect_s3_class(ph, "pheatmap")

  #downsampling the number of programs
  ph <- plotMetaPrograms(test.mp, downsample = 5, silent = TRUE)
  expect_s3_class(ph, "pheatmap")

  #downsampling is reproducible
  rows <- function(ph) {
    ph$gtable$grobs[[which(ph$gtable$layout$name == "row_names")]]$label
  }
  ph1 <- plotMetaPrograms(test.mp, downsample = 8, silent = TRUE, show_rownames = TRUE)
  ph2 <- plotMetaPrograms(test.mp, downsample = 8, silent = TRUE, show_rownames = TRUE)
  expect_length(rows(ph1), 8)
  expect_equal(rows(ph1), rows(ph2))

  #custom annotation colors
  cols <- setNames(c("red", "green", "blue"), mp.names)
  ph <- plotMetaPrograms(test.mp, annotation_colors = cols, silent = TRUE)
  expect_s3_class(ph, "pheatmap")
  ph <- plotMetaPrograms(test.mp, annotation_colors = list(Metaprogram = cols),
                         silent = TRUE)
  expect_s3_class(ph, "pheatmap")
})

test_that("plotMetaPrograms sorts programs by clustering when the tree is hidden", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  heatmap.rows <- function(ph) {
    ph$gtable$grobs[[which(ph$gtable$layout$name == "row_names")]]$label
  }
  tree.order <- labels(as.dendrogram(test.mp$programs.tree))

  ph <- plotMetaPrograms(test.mp, showtree = FALSE, silent = TRUE,
                         show_rownames = TRUE)
  expect_equal(heatmap.rows(ph), tree.order)
  #programs of the same MP are contiguous
  expect_length(rle(unname(test.mp$programs.clusters[heatmap.rows(ph)]))$values,
                test.nMP)

  #same order with or without the tree
  ph.tree <- plotMetaPrograms(test.mp, silent = TRUE, show_rownames = TRUE)
  expect_equal(heatmap.rows(ph.tree), tree.order)

  #also after dropping MPs
  dropped <- dropMetaPrograms(test.mp, dropMP = "MP2")
  ph <- plotMetaPrograms(dropped, silent = TRUE, show_rownames = TRUE)
  expect_equal(heatmap.rows(ph), labels(as.dendrogram(dropped$programs.tree)))
})

test_that("getMetaPrograms validates nMP", {
  nprog <- length(prog.names)
  expect_error(getMetaPrograms(test.nmf, nMP = nprog + 1), "cannot be larger")
  expect_error(getMetaPrograms(test.nmf, nMP = 0), "positive integer")
  expect_error(getMetaPrograms(test.nmf, nMP = 2.5), "positive integer")

  #a single MP, or one MP per program, are valid (silhouette is undefined)
  mp <- getMetaPrograms(test.nmf, nMP = 1)
  expect_named(mp$metaprograms.genes, "MP1")
  expect_true(is.na(mp$metaprograms.metrics$silhouette))
  expect_equal(mp$metaprograms.metrics$numberPrograms, nprog)

  mp <- getMetaPrograms(test.nmf, nMP = nprog, remove.empty = FALSE)
  expect_length(mp$metaprograms.genes, nprog)
})

test_that("getMetaPrograms applies an inclusive confidence threshold", {
  #with min.confidence=1, genes must be found in all programs of the MP
  strict <- getMetaPrograms(test.nmf, nMP = test.nMP, min.confidence = 1,
                            remove.empty = FALSE)
  loose <- getMetaPrograms(test.nmf, nMP = test.nMP, min.confidence = 0,
                           remove.empty = FALSE)
  expect_gt(sum(lengths(strict$metaprograms.genes)), 0)
  expect_gt(sum(lengths(loose$metaprograms.genes)),
            sum(lengths(strict$metaprograms.genes)))
})

test_that("plotMetaPrograms validates annotation colors", {
  expect_error(plotMetaPrograms(test.mp, annotation_colors = c("red", "blue")),
               "at least nMP elements")
  expect_error(plotMetaPrograms(test.mp, annotation_colors = list(other = "red")),
               "annotation_colors")
})
