test_that("jaccardIndex computes set overlap", {
  expect_equal(jaccardIndex(c("a", "b", "c"), c("b", "c", "d")), 0.5)
  expect_equal(jaccardIndex(c("a", "b"), c("a", "b")), 1)
  expect_equal(jaccardIndex(c("a", "b"), c("c", "d")), 0)
})

test_that("jaccardSimilarity agrees with pairwise jaccardIndex", {
  gl <- list(p1 = c(a = 1, b = 1, c = 1),
             p2 = c(b = 1, c = 1, d = 1),
             p3 = c(x = 1, y = 1))
  J <- jaccardSimilarity(gl)

  expect_equal(dim(J), c(3, 3))
  expect_equal(rownames(J), names(gl))
  expect_equal(colnames(J), names(gl))
  expect_true(isSymmetric(J))
  expect_equal(unname(diag(J)), rep(1, 3))
  expect_equal(J["p1", "p2"], jaccardIndex(names(gl$p1), names(gl$p2)))
  expect_equal(J["p1", "p3"], 0)
})

test_that("normVector rescales to unit sum", {
  expect_equal(normVector(c(1, 1, 2)), c(0.25, 0.25, 0.5))
  #all-zero vectors are returned unchanged (no division by zero)
  expect_equal(normVector(c(0, 0, 0)), c(0, 0, 0))
})

test_that("weightCumul keeps top elements up to the cumulative threshold", {
  x <- c(a = 1, b = 4, c = 3, d = 2)

  #sorted weights are 0.4, 0.3, 0.2, 0.1; cumulative 0.4, 0.7, 0.9, 1.0
  #the element that crosses the threshold is included
  expect_equal(weightCumul(x, weight.explained = 0.8), c(b = 0.4, c = 0.3, d = 0.2))
  expect_equal(weightCumul(x, weight.explained = 0.7), c(b = 0.4, c = 0.3))
  expect_equal(weightCumul(x, weight.explained = 1), c(b = 0.4, c = 0.3, d = 0.2, a = 0.1))
})

test_that("weightCumul always returns the top element, unless all are zero", {
  #a dominant element alone explains more than the threshold
  expect_equal(weightCumul(c(a = 8, b = 1, c = 1), weight.explained = 0.5), c(a = 0.8))
  #elements with zero weight are never returned
  expect_equal(weightCumul(c(a = 1, b = 0, c = 0), weight.explained = 1), c(a = 1))
  expect_length(weightCumul(c(a = 0, b = 0, c = 0)), 0)
})

test_that("nameSamples assigns unique sample names", {
  expect_named(nameSamples(list(1, 2)), c("sample1", "sample2"))
  expect_named(nameSamples(list(a = 1, 2)), c("a", "sample2"))
  expect_named(nameSamples(list(a = 1, b = 2)), c("a", "b"))
  expect_error(nameSamples(list(a = 1, a = 2)), "must be unique")
})

test_that("clusterNames maps programs to meta-program names", {
  mp <- list(metaprograms.genes = list(A = "g1", B = "g2"),
             programs.clusters = c(p1 = 2, p2 = 1, p3 = NA, p4 = 2))
  expect_equal(clusterNames(mp), c(p1 = "B", p2 = "A", p3 = NA, p4 = "B"))

  #objects from older versions, after dropping MP2
  mp <- list(metaprograms.genes = list(MP1 = "g1", MP3 = "g2"),
             programs.clusters = c(p1 = 3, p2 = 1, p3 = NA))
  expect_equal(clusterNames(mp), c(p1 = "MP3", p2 = "MP1", p3 = NA))
})

test_that("wgtLoad returns column-normalized loadings", {
  m <- matrix(c(5, 1, 1,
                1, 5, 1,
                2, 2, 2), nrow = 3, byrow = TRUE,
              dimnames = list(c("g1", "g2", "g3"), c("p1", "p2", "p3")))

  res0 <- wgtLoad(m, w = 0)
  expect_equal(dimnames(res0), dimnames(m))
  expect_equal(unname(colSums(res0)), rep(1, 3))
  #without specificity weighting, this is a simple column normalization
  expect_equal(res0, apply(m, 2, normVector))

  #with specificity weighting, the non-specific gene (g3) is penalized
  res5 <- wgtLoad(m, w = 5)
  expect_equal(unname(colSums(res5)), rep(1, 3))
  expect_lt(res5["g3", "p1"], res0["g3", "p1"])
  expect_gt(res5["g1", "p1"], res0["g1", "p1"])
})

test_that("geneList2table binds programs into a single table", {
  gl <- list(s1.k2 = matrix(1:6, nrow = 3, dimnames = list(c("a", "b", "c"), NULL)),
             s2.k2 = matrix(7:12, nrow = 3, dimnames = list(c("a", "b", "c"), NULL)))
  tab <- geneList2table(gl)

  expect_equal(dim(tab), c(3, 4))
  expect_equal(colnames(tab), c("s1.k2.1", "s1.k2.2", "s2.k2.1", "s2.k2.2"))
  expect_equal(rownames(tab), c("a", "b", "c"))
})

test_that("downsampleMin keeps at least one element per class", {
  v <- setNames(rep(1:4, each = 25), paste0("prog", 1:100))

  #no downsampling needed
  expect_identical(downsampleMin(v, size = 200), v)

  ds <- downsampleMin(v, size = 20)
  expect_length(ds, 20)
  expect_false(any(duplicated(ds)))
  expect_true(all(ds %in% names(v)))
  expect_setequal(unique(v[ds]), 1:4)

  #size smaller than the number of classes: one representative per class
  expect_length(downsampleMin(v, size = 2), 4)
})

test_that("downsampleMin is reproducible with a seed", {
  v <- setNames(rep(1:4, each = 25), paste0("prog", 1:100))

  expect_identical(downsampleMin(v, size = 20, seed = 1),
                   downsampleMin(v, size = 20, seed = 1))
  expect_false(identical(downsampleMin(v, size = 20, seed = 1),
                         downsampleMin(v, size = 20, seed = 2)))

  #the global random number stream is not altered
  set.seed(42)
  expected <- runif(3)
  set.seed(42)
  ds <- downsampleMin(v, size = 20, seed = 1)
  expect_identical(runif(3), expected)
})

test_that("programSample extracts sample names from program names", {
  expect_equal(programSample(c("s1.k3.1", "s2.k10.12")), c("s1", "s2"))
  #sample names that contain dots, or resemble the program suffix
  expect_equal(programSample("patient.A.k5.2"), "patient.A")
  expect_equal(programSample("run.k2.1_rep.k4.3"), "run.k2.1_rep")
})

test_that("getEntropy is maximal for uniform distributions", {
  expect_equal(getEntropy(c(5, 5, 5, 5), pseudo = 0), 2)
  expect_equal(getEntropy(c(10, 0), pseudo = 0), 0)
  expect_gt(getEntropy(c(5, 5, 5, 5)), getEntropy(c(17, 1, 1, 1)))
})
