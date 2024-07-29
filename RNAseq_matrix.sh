#!/bin/bash 

# Load necessary libraries 
Rscript --vanilla RNAseq_counts.sh
library(org.Hs.eg.db)
library(AnnotationDbi)

# Set working directory to counts directory
setwd("path/to/file/counts")

# Function to read featureCounts output
read_featureCounts <- function(filename) {
  # Read the first few lines to determine the number of lines to skip
  lines <- readLines(filename, warn = FALSE)
  skip_lines <- sum(grepl("^#", lines))  # Count comment lines to skip
  
  # Read the file, skipping comment lines
  dat <- read.delim(filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = skip_lines)
  
  # Verify that 'Geneid' is among the column names
  if (!"Geneid" %in% colnames(dat)) {
    stop("Geneid column not found in ", filename)
  }
  
  # Extract 'Geneid' column and the counts column (assumed to be the last column)
  gene_id_col <- which(colnames(dat) == "Geneid")
  counts_col <- ncol(dat)
  dat_subset <- dat[, c(gene_id_col, counts_col)]
  
  # Rename the columns: Geneid and the sample name based on the file name
  colnames(dat_subset) <- c("Geneid", gsub("_counts\\.txt$", "", basename(filename)))
  return(dat_subset)
}

# List all count files
count_files <- list.files(pattern = "_counts.txt$")

# Read and process each count file
count_list <- lapply(count_files, read_featureCounts)

# Merge all count data into a single data frame by gene ID
gene_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)

# Map ENSEMBL gene IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_counts$Geneid,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Replace ENSEMBL IDs with gene symbols in the data frame
gene_counts$Geneid <- gene_symbols[match(gene_counts$Geneid, names(gene_symbols))]

# Aggregate counts by gene symbol, summing counts for duplicated symbols
aggregate_cols <- setdiff(names(gene_counts), "Geneid")
gene_counts_aggregated <- aggregate(. ~ Geneid, data = gene_counts, FUN = sum)

# Set Geneid column as row names
rownames(gene_counts_aggregated) <- gene_counts_aggregated$Geneid
gene_counts_aggregated$Geneid <- NULL  # Remove the now redundant 'Geneid' column

# Write the results to a CSV file 
write.csv(gene_counts_aggregated, file.path("path/to/file/Mapping", "gene_symbol_counts_matrix.csv"), row.names = TRUE)