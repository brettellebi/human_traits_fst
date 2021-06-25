#!/usr/bin/env Rscript

# Libraries

library(tidyverse)

# Load data

in_list = readRDS(snakemake@input[[1]])

# Extract SNP IDs

in_list@risk_alleles %>%
    dplyr::pull(variant_id) %>%
    writeLines(snakemake@output[[1]])

# Note: Some SNP IDs are written as e.g. chr5:17098189