#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(destiny)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

expr <- as.matrix(read.csv(file.path(out_dir, "../../data/fixture_guo_expression.csv"), row.names=1))
rownames(expr) <- paste0("cell_", seq_len(nrow(expr)))

set.seed(42)
dm <- destiny::DiffusionMap(expr, sigma="local", n_eigs=5, k=20)
dpt <- destiny::DPT(dm, tips=1L)

# Plain DM
png(file.path(out_dir, "R_diffmap.png"), width=600, height=400, res=100)
plot(dm, 1:2)
dev.off()

# DPT
png(file.path(out_dir, "R_dpt.png"), width=600, height=400, res=100)
plot(dpt, root=1)
dev.off()

# Save eigenvectors and DPT vec so Py uses same numbers
write.csv(dm@eigenvectors, file.path(out_dir, "eigenvectors.csv"), row.names=FALSE)
write.csv(data.frame(dpt = dpt[1, ]), file.path(out_dir, "dpt.csv"), row.names=FALSE)
cat("R plots done\n")
