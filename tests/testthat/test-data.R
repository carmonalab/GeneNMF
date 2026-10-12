test_that("getDataMatrix returns a genes-by-cells matrix", {
  mat <- getDataMatrix(sampleObj)

  expect_equal(dim(mat), dim(sampleObj))
  expect_equal(rownames(mat), rownames(sampleObj))
  expect_equal(colnames(mat), colnames(sampleObj))
  expect_true(all(mat >= 0))
})

test_that("getDataMatrix keeps sparse matrices sparse", {
  mat <- getDataMatrix(sampleObj, hvg = test.hvg)
  expect_s4_class(mat, "sparseMatrix")
  expect_equal(as.matrix(mat),
               as.matrix(GetAssayData(sampleObj, layer = "data")[test.hvg, ]))

  #centering or scaling return a dense matrix
  expect_true(is.matrix(getDataMatrix(sampleObj, hvg = test.hvg, center = TRUE)))
  expect_true(is.matrix(getDataMatrix(sampleObj, hvg = test.hvg, scale = TRUE)))
})

test_that("getDataMatrix enforces non-negative values on sparse matrices", {
  obj <- sampleObj
  mat <- GetAssayData(obj, layer = "data")
  mat[test.hvg[1], 1:5] <- -1
  obj <- SetAssayData(obj, layer = "data", new.data = mat)

  res <- getDataMatrix(obj, hvg = test.hvg[1:10])
  expect_s4_class(res, "sparseMatrix")
  expect_true(all(res >= 0))
  expect_equal(unname(res[1, 1:5]), rep(0, 5))

  res <- getDataMatrix(obj, hvg = test.hvg[1:10], non_negative = FALSE)
  expect_equal(unname(res[1, 1:5]), rep(-1, 5))
})

test_that("getDataMatrix subsets on hvg", {
  genes <- test.hvg[1:20]
  mat <- getDataMatrix(sampleObj, hvg = genes)

  expect_equal(rownames(mat), genes)
  expect_equal(ncol(mat), ncol(sampleObj))
})

test_that("getDataMatrix centering and scaling respect non-negativity", {
  genes <- test.hvg[1:20]

  mat <- getDataMatrix(sampleObj, hvg = genes, center = TRUE, scale = TRUE)
  expect_false(anyNA(mat))
  expect_true(all(mat >= 0))

  mat.neg <- getDataMatrix(sampleObj, hvg = genes, center = TRUE,
                           non_negative = FALSE)
  expect_true(any(mat.neg < 0))
  expect_equal(unname(rowMeans(mat.neg)), rep(0, length(genes)))
})

test_that("findVariableFeatures_wfilters returns the requested number of HVG", {
  obj <- quiet(findVariableFeatures_wfilters(sampleObj, nfeatures = 100))

  expect_s4_class(obj, "Seurat")
  hvg <- VariableFeatures(obj)
  expect_length(hvg, 100)
  expect_true(all(hvg %in% rownames(sampleObj)))
})

test_that("findVariableFeatures_wfilters applies expression filters", {
  min.exp <- 0.05
  max.exp <- 2
  obj <- quiet(findVariableFeatures_wfilters(sampleObj, nfeatures = 100,
                                             min.exp = min.exp, max.exp = max.exp))
  hvg <- VariableFeatures(obj)
  means <- rowMeans(as.matrix(GetAssayData(obj, layer = "data")[hvg, ]))

  expect_true(all(means >= min.exp & means <= max.exp))
})

test_that("findVariableFeatures_wfilters handles filters that exclude all genes", {
  #loess warnings from FindVariableFeatures are expected on tiny datasets
  warns <- character(0)
  obj <- withCallingHandlers(
    findVariableFeatures_wfilters(sampleObj, min.exp = 50, max.exp = 60),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_true(any(grepl("No variable features passed the filters", warns)))
  expect_length(VariableFeatures(obj), 0)
})

test_that("findVariableFeatures_wfilters excludes genes in blocklist", {
  ref <- VariableFeatures(quiet(findVariableFeatures_wfilters(sampleObj, nfeatures = 100)))
  block <- ref[1:10]

  #blocklist as a vector
  obj <- quiet(findVariableFeatures_wfilters(sampleObj, nfeatures = 100,
                                             genesBlockList = block))
  expect_length(intersect(VariableFeatures(obj), block), 0)

  #blocklist as a list of vectors
  obj <- quiet(findVariableFeatures_wfilters(sampleObj, nfeatures = 100,
                                             genesBlockList = list(a = block[1:5], b = block[6:10])))
  expect_length(intersect(VariableFeatures(obj), block), 0)
})
