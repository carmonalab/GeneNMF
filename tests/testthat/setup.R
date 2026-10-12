# Shared fixtures, computed once for the whole test suite.
# Everything runs on the 50-cell 'sampleObj' shipped with the package, split
# into two pseudo-samples, so the full suite takes a few seconds and needs no
# internet access.

# CRAN policy: do not use more than 2 cores
options(RcppML.threads = 1)

library(Seurat)
data(sampleObj, package = "GeneNMF")

# FindVariableFeatures emits harmless loess warnings on tiny datasets
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test.objs <- list(s1 = sampleObj[, seq(1, ncol(sampleObj), 2)],
                  s2 = sampleObj[, seq(2, ncol(sampleObj), 2)])

test.hvg <- VariableFeatures(sampleObj)

test.k <- 3:4
test.nMP <- 3

test.nmf <- quiet(multiNMF(test.objs, k = test.k, hvg = test.hvg))
test.mp <- getMetaPrograms(test.nmf, nMP = test.nMP)
