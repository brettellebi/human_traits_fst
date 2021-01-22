#!/usr/bin/env Rscript

# Load libaries

library(here)
library(tidyverse)
library(pegas)
library(knitr)
library(plotly)
library(ggridges)
library(qqman)

# Functions

## Read `.afreq` files from `Plink`
read_afreq <- function(file){
  out = read.table(file, header = T, comment.char = "") %>% 
    dplyr::rename(CHR = X.CHROM,
                  SNP = ID)
  
  return(out)
}

## Get manhattan plot
get_man <- function(df, trait, chr, bp, snp, p){
  # Set title with number of SNPs
  title <- paste(trait, "\n", "SNP count:", nrow(df))
  # Plot
  manhattan(df, chr=chr, bp=bp, snp=snp, p=p,
            cex = 0.6,
            c(pal_primary[names(pal_primary) == trait],
              pal_secondary[names(pal_secondary) == trait]),
            main = title)
}

# Plotting parameters

# Create factor levels for `trait`
trait_levels = c("hei", "bmi", "edu", "int", "ibd", "pig")

# Create vector for recoding traits with full names
recode_vec = c("hei" = "Height",
               "bmi" = "BMI",
               "edu" = "Educational attainment",
               "int" = "Intelligence",
               "ibd" = "IBD",
               "pig" = "Pigmentation")

# Colour palettes
pal_primary = c("Height" = "#FC4E07",
                "BMI" = "#FFBF00",
                "Educational attainment" = "#0BC166",
                "Intelligence" = "#00AFBB",
                "IBD" = "#D84797",
                "Pigmentation" = "#360568")
pal_secondary = c("Height" = "#D11F1F",
                  "BMI" = "#cc7e08",
                  "Educational attainment" = "#39BFA2",
                  "Intelligence" = "#02395c",
                  "IBD" = "#82043B",
                  "Pigmentation" = "#960592")