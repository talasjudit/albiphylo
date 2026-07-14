#!/usr/bin/env Rscript
# Plot pairwise SNP distance / count matrices from step 8 as heatmaps.
#
# Usage:
#   Rscript plot_pairwise_heatmap.R <distance_matrix.tsv> <snp_counts.tsv> <output_prefix>
#
# Output: <output_prefix>_distance_heatmap.pdf and <output_prefix>_snp_counts_heatmap.pdf
#
# Requires: pheatmap (install.packages("pheatmap"))

suppressPackageStartupMessages(library(pheatmap))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: plot_pairwise_heatmap.R <distance_matrix.tsv> <snp_counts.tsv> <output_prefix>")
}

dist_file   <- args[1]
counts_file <- args[2]
prefix      <- args[3]

read_matrix <- function(path) {
    m <- as.matrix(read.table(path, header = TRUE, sep = "\t",
                              row.names = 1, check.names = FALSE))
    storage.mode(m) <- "numeric"
    m
}

dist_mat   <- read_matrix(dist_file)
counts_mat <- read_matrix(counts_file)

# Dynamic plot dimensions based on sample count
n_samples <- nrow(dist_mat)
plot_dim <- max(6, 4 + n_samples * 0.25)

# Normalised distance: clustering by hierarchical average linkage
pheatmap(dist_mat,
         main          = "Pairwise SNP distance (Anderson 2023 dosage)",
         display_numbers = TRUE,
         number_format = "%.4f",
         clustering_distance_rows = as.dist(dist_mat),
         clustering_distance_cols = as.dist(dist_mat),
         clustering_method        = "average",
         filename      = paste0(prefix, "_distance_heatmap.pdf"),
         width = plot_dim, height = plot_dim)

# Raw differing-SNP counts: same clustering order as distance matrix for visual consistency
pheatmap(counts_mat,
         main          = "Pairwise differing SNPs (raw count)",
         display_numbers = TRUE,
         number_format = "%d",
         clustering_distance_rows = as.dist(dist_mat),
         clustering_distance_cols = as.dist(dist_mat),
         clustering_method        = "average",
         filename      = paste0(prefix, "_snp_counts_heatmap.pdf"),
         width = plot_dim, height = plot_dim)

cat("Wrote:",
    paste0(prefix, "_distance_heatmap.pdf"),
    paste0(prefix, "_snp_counts_heatmap.pdf"),
    sep = "\n")
