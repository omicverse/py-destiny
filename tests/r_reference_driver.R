#!/usr/bin/env Rscript
suppressMessages({library(destiny); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
fixture_path <- args[1]; output_path <- args[2]

data(guo, package = "destiny")
# guo is an ExpressionSet; pull the expression matrix
expr <- t(as.matrix(Biobase::exprs(guo)))  # cells × genes
cat("[ref] expr dims:", dim(expr), "\n")

# Save sidecar
write.csv(expr, sub("\\.rds$", "_expression.csv", fixture_path))

set.seed(42)
dm <- destiny::DiffusionMap(expr, sigma = "local", n_eigs = 5, verbose = FALSE)
cat("[ref] DiffusionMap:", length(dm@eigenvalues), "eigvals\n")

dpt <- destiny::DPT(dm, tips = 1L)
# DPT() returns a DPT object; convert to numeric pseudotime via dpt[, 1]
dpt_vec <- as.numeric(dpt[1, ])   # pseudotime from root cell 1 to all others
cat("[ref] DPT length:", length(dpt_vec), "\n")

out <- list(
  eigenvalues = as.numeric(dm@eigenvalues),
  eigenvectors = as.matrix(dm@eigenvectors),
  dpt = dpt_vec
)
write_json(out, output_path, auto_unbox = TRUE, digits = NA,
           matrix = "rowmajor", na = "null", pretty = FALSE)
cat("[ref] wrote", output_path, "\n")
