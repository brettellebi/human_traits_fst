#!/usr/bin/env Rscript

# Libraries

library(tidyverse)

# Read in data

assocs = readRDS(snakemake@input[[1]])