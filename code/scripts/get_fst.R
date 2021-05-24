#!/usr/bin/env Rscript

# Get variables

args = commandArgs(trailingOnly = TRUE)

in_vcf = args[1]
pop_file = args[2]
out_file = args[3]

# Load libraries

library(pegas)
library(tidyverse)

# Get VCF INFO

loci = pegas::VCFloci(in_vcf)
n_loci = nrow(loci) # Get number of loci

# Read VCF to get genotypes

vcf = pegas::read.vcf(in_vcf, to = n_loci)

# Read in population file and filter for those in `vcf`

pops = readr::read_csv(pop_file)
pops = pops[pops$Sample %in% rownames(vcf), ]

# Calculate Fst

fst = as.data.frame(pegas::Fst(vcf, pop = pops))

# Bind with loci info

out = cbind(dplyr::select(loci, CHROM, POS, REF, ALT),
            fst) %>%
        dplyr::mutate(LOCUS = paste(CHROM, POS, sep = ":")) %>%
        dplyr::select(LOCUS, everything())

# Write to file

readr::write_csv(out, out_file)



