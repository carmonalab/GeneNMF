#runGSEA depends on suggested packages; tests use a custom database of
#signatures, so that no internet access (or the msigdbr package) is required

test_that("runGSEA finds enriched signatures in a custom database", {
  skip_if_not_installed("fgsea")

  universe <- paste0("gene", 1:100)
  db <- list(setA = universe[1:10],
             setB = universe[51:60])
  res <- runGSEA(universe[1:8], universe = universe, custom.db = db)

  expect_s3_class(res, "data.frame")
  expect_equal(res$pathway, "setA")
  expect_equal(res$overlap, 8)
  expect_true(all(res$pval <= 0.05))

  #duplicated genes are ignored
  res.dup <- runGSEA(rep(universe[1:8], 2), universe = universe, custom.db = db)
  expect_equal(res.dup$overlap, res$overlap)

  #without a universe, all genes in the database are used
  expect_message(res.nouniv <- runGSEA(universe[1:8], custom.db = db, pval.thr = 1),
                 "No 'universe' specified")
  expect_setequal(res.nouniv$pathway, names(db))

  #relaxing the p-value threshold returns all signatures
  res.all <- runGSEA(universe[1:8], universe = universe, custom.db = db,
                     pval.thr = 1)
  expect_setequal(res.all$pathway, names(db))
})

test_that("runGSEA validates the custom database", {
  skip_if_not_installed("fgsea")

  expect_error(runGSEA(c("a", "b"), universe = letters,
                       custom.db = list(c("a", "b"))),
               "custom.db should be a named list")
  expect_error(runGSEA(c("a", "b"), universe = letters,
                       custom.db = c(a = "a", b = "b")),
               "custom.db should be a named list")
})
