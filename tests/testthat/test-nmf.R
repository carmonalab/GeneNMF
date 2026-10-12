test_that("multiNMF returns one model per sample and k", {
  expect_type(test.nmf, "list")
  expect_named(test.nmf, c("s1.k3", "s1.k4", "s2.k3", "s2.k4"))

  for (n in names(test.nmf)) {
    model <- test.nmf[[n]]
    k <- as.numeric(gsub(".*\\.k", "", n))
    obj <- test.objs[[gsub("\\.k\\d+$", "", n)]]

    expect_true(all(c("w", "d", "h") %in% names(model)))
    expect_equal(dim(model$w), c(length(test.hvg), k))
    expect_equal(dim(model$h), c(k, ncol(obj)))
    expect_equal(rownames(model$w), test.hvg)
    expect_equal(colnames(model$w), paste0("pattern", 1:k))
    expect_equal(rownames(model$h), paste0("pattern", 1:k))
    expect_equal(colnames(model$h), colnames(obj))
    expect_length(model$d, k)
    expect_true(all(model$w >= 0))
    expect_true(all(model$h >= 0))
  }
})

test_that("nmfInput returns a dense matrix unless Matrix is attached", {
  sparse <- getDataMatrix(test.objs$s1, hvg = test.hvg)
  dense <- as.matrix(sparse)
  expect_identical(nmfInput(dense), dense)

  if ("package:Matrix" %in% search()) {
    expect_s4_class(nmfInput(sparse), "sparseMatrix")
  } else {
    expect_identical(nmfInput(sparse), dense)
  }
})

test_that("multiNMF gives the same results on sparse and dense matrices", {
  skip_if("package:Matrix" %in% search(), "Matrix package is already attached")

  #sparse matrices are only used if the Matrix package is attached
  suppressPackageStartupMessages(library(Matrix))
  on.exit(detach("package:Matrix"))
  res <- quiet(multiNMF(test.objs, k = test.k, hvg = test.hvg))

  expect_equal(res, test.nmf, tolerance = 1e-8)
})

test_that("multiNMF calculates HVG when not provided", {
  res <- quiet(multiNMF(test.objs, k = 3, nfeatures = 100))

  expect_named(res, c("s1.k3", "s2.k3"))
  expect_equal(nrow(res[[1]]$w), 100)
  expect_equal(rownames(res[[1]]$w), rownames(res[[2]]$w))
})

test_that("multiNMF is reproducible for a fixed seed", {
  res1 <- quiet(multiNMF(test.objs["s1"], k = 3, hvg = test.hvg, seed = 42))
  res2 <- quiet(multiNMF(test.objs["s1"], k = 3, hvg = test.hvg, seed = 42))

  expect_equal(res1, res2)
})

test_that("multiNMF validates k", {
  expect_error(multiNMF(test.objs, k = "three", hvg = test.hvg),
               "k must be a numeric vector")

  expect_warning(res <- multiNMF(test.objs["s1"], k = 1:2, hvg = test.hvg),
                 "larger than 1")
  expect_named(res, "s1.k2")

  #zero, negative, non-integer and duplicated values
  expect_warning(res <- multiNMF(test.objs["s1"], k = c(0, -1, 2.5, 3, 3), hvg = test.hvg),
                 "larger than 1")
  expect_named(res, "s1.k3")

  expect_error(suppressWarnings(multiNMF(test.objs, k = c(0, 1), hvg = test.hvg)),
               "No valid values of k")
})

test_that("multiNMF skips k larger than the sample size", {
  objs <- list(s1 = test.objs$s1, small = sampleObj[, 1:12])
  expect_warning(res <- multiNMF(objs, k = c(3, 15), hvg = test.hvg),
                 "Sample small: k must be smaller.*k=15")
  expect_named(res, c("s1.k3", "s1.k15", "small.k3"))

  #no valid k for a sample
  expect_warning(res <- multiNMF(objs, k = 15, hvg = test.hvg), "Sample small")
  expect_named(res, "s1.k15")
  expect_error(suppressWarnings(multiNMF(objs["small"], k = 15, hvg = test.hvg)),
               "No NMF models could be calculated")
})

test_that("multiNMF ignores small samples", {
  objs <- c(test.objs, list(tiny = sampleObj[, 1:5]))
  res <- quiet(multiNMF(objs, k = 3, hvg = test.hvg, min.cells.per.sample = 10))

  expect_named(res, c("s1.k3", "s2.k3"))

  #samples with exactly the minimum number of cells are retained
  res <- quiet(multiNMF(list(a = sampleObj[, 1:10]), k = 2, hvg = test.hvg,
                        min.cells.per.sample = 10))
  expect_named(res, "a.k2")

  expect_error(multiNMF(test.objs, k = 3, hvg = test.hvg, min.cells.per.sample = 100),
               "No samples with at least 100 cells")
})

test_that("multiNMF assigns names to unnamed samples", {
  res <- quiet(multiNMF(unname(test.objs), k = 3, hvg = test.hvg))
  expect_named(res, c("sample1.k3", "sample2.k3"))

  #sample metrics are correct for unnamed lists
  res <- quiet(multiNMF(unname(test.objs), k = test.k, hvg = test.hvg))
  mp <- getMetaPrograms(res, nMP = test.nMP)
  expect_equal(colnames(mp$metaprograms.composition), c("sample1", "sample2"))
  expect_equal(mp$metaprograms.metrics, test.mp$metaprograms.metrics)

  objs <- setNames(test.objs, c("a", "a"))
  expect_error(multiNMF(objs, k = 3, hvg = test.hvg), "must be unique")
})

test_that("multiNMF ignores variable genes missing from the samples", {
  expect_warning(res <- multiNMF(test.objs["s1"], k = 3,
                                 hvg = c(test.hvg[1:50], "NOT_A_GENE")),
                 "1 genes in 'hvg' were not found")
  expect_equal(rownames(res[[1]]$w), test.hvg[1:50])

  expect_error(suppressWarnings(multiNMF(test.objs, k = 3, hvg = c("NOT_A_GENE", "NOR_THIS"))),
               "Not enough variable genes")
})

test_that("multiNMF calculates HVG on the specified assay", {
  genes <- rownames(sampleObj)[1:3000]
  objs <- lapply(test.objs, function(x) {
    x[["ALT"]] <- CreateAssayObject(counts = GetAssayData(x, layer = "counts")[genes, ])
    NormalizeData(x, assay = "ALT", verbose = FALSE)
  })
  res <- quiet(multiNMF(objs, assay = "ALT", k = 3, nfeatures = 100))

  expect_named(res, c("s1.k3", "s2.k3"))
  expect_equal(nrow(res[[1]]$w), 100)
  expect_true(all(rownames(res[[1]]$w) %in% genes))
})

test_that("getNMFgenes handles programs with the same number of genes", {
  w <- cbind(pattern1 = c(.4, .3, .2, .1, 0, 0),
             pattern2 = c(0, 0, .1, .2, .3, .4))
  rownames(w) <- paste0("g", 1:6)
  genes <- getNMFgenes(list(s.k2 = list(w = w)), specificity.weight = NULL,
                       weight.explained = 0.7)

  expect_named(genes, c("s.k2.1", "s.k2.2"))
  expect_equal(genes[["s.k2.1"]], c(g1 = 0.4, g2 = 0.3))
  expect_equal(genes[["s.k2.2"]], c(g6 = 0.4, g5 = 0.3))
})

test_that("getNMFgenes handles dominant genes and empty programs", {
  w <- cbind(pattern1 = c(.7, .1, .1, .1, 0, 0),
             pattern2 = c(0, 0, 0, 0, 0, 0),
             pattern3 = c(.1, .3, .25, .2, .1, .05))
  rownames(w) <- paste0("g", 1:6)
  model <- list(s.k3 = list(w = w))

  #the top gene explains more than 'weight.explained': keep it;
  #empty programs are dropped, others keep the index of their component
  genes <- getNMFgenes(model, specificity.weight = NULL, weight.explained = 0.5)
  expect_named(genes, c("s.k3.1", "s.k3.3"))
  expect_equal(genes[["s.k3.1"]], c(g1 = 0.7))
  expect_named(genes[["s.k3.3"]], c("g2", "g3"))

  expect_message(genes <- getNMFgenes(model, weight.explained = 0.5),
                 "Dropped 1 empty NMF programs")
  expect_named(genes, c("s.k3.1", "s.k3.3"))
})

test_that("getNMFgenes returns named gene weights for each program", {
  genes <- getNMFgenes(test.nmf)

  expect_type(genes, "list")
  expect_length(genes, 2 * sum(test.k))
  expect_true(all(grepl("^s[12]\\.k[34]\\.\\d+$", names(genes))))

  for (g in genes) {
    expect_type(g, "double")
    expect_gt(length(g), 0)
    expect_true(all(names(g) %in% test.hvg))
    #weights are sorted; all genes but the last explain less than 'weight.explained'
    expect_false(is.unsorted(rev(g)))
    expect_lt(sum(head(g, -1)), 0.5)
  }
})

test_that("getNMFgenes respects max.genes and weight.explained", {
  genes <- getNMFgenes(test.nmf, max.genes = 5)
  expect_true(all(lengths(genes) <= 5))

  few <- getNMFgenes(test.nmf, weight.explained = 0.2, max.genes = 1000)
  many <- getNMFgenes(test.nmf, weight.explained = 0.8, max.genes = 1000)
  expect_true(all(lengths(few[names(many)]) <= lengths(many)))
  expect_gt(sum(lengths(many)), sum(lengths(few)))
})

test_that("runNMF stores a new dimensionality reduction", {
  k <- 4
  obj <- runNMF(sampleObj, k = k)

  expect_s4_class(obj, "Seurat")
  expect_true("NMF" %in% Reductions(obj))
  expect_equal(dim(Embeddings(obj, "NMF")), c(ncol(sampleObj), k))
  expect_equal(dim(Loadings(obj, "NMF")), c(length(test.hvg), k))
  expect_equal(rownames(Embeddings(obj, "NMF")), colnames(sampleObj))
  expect_equal(colnames(Embeddings(obj, "NMF")), paste0("NMF_", 1:k))
  expect_true(all(Embeddings(obj, "NMF") >= 0))

  #a valid DimReduc object, with details of the NMF model
  dr <- obj[["NMF"]]
  expect_s4_class(dr, "DimReduc")
  expect_true(validObject(dr))
  expect_equal(DefaultAssay(dr), "RNA")
  expect_equal(Key(dr), "NMF_")
  expect_length(Stdev(dr), k)
  expect_named(Misc(dr), c("d", "tol", "iter"))
  expect_equal(Misc(dr)$d, Stdev(dr))
})

test_that("runNMF accepts custom genes and reduction name", {
  genes <- test.hvg[1:50]
  obj <- runNMF(sampleObj, k = 3, hvg = genes, new.reduction = "myNMF")

  expect_true("myNMF" %in% Reductions(obj))
  expect_equal(rownames(Loadings(obj, "myNMF")), genes)
  expect_equal(colnames(Embeddings(obj, "myNMF")), paste0("myNMF_", 1:3))
})

test_that("runNMF ignores variable genes missing from the object", {
  expect_warning(obj <- runNMF(sampleObj, k = 3, hvg = c(test.hvg[1:50], "NOT_A_GENE")),
                 "1 genes in 'hvg' were not found")
  expect_equal(rownames(Loadings(obj, "NMF")), test.hvg[1:50])
})

test_that("runNMF fails without variable features", {
  obj <- sampleObj
  VariableFeatures(obj) <- character(0)

  expect_error(runNMF(obj, k = 3), "No variable features found")
})
