#!/usr/bin/env Rscript

# Get variables

args = commandArgs(trailingOnly = TRUE)

efo_id = args[1]
out_file = args[2]

# Load libraries

library(gwasrapidd)

# Get traits

gcat_assoc = gwasrapidd::get_variants(efo_id = efo_id)

# Get list of associations

out = saveRDS(gcat_assoc, out_file)
