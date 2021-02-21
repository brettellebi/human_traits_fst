#!/usr/bin/env Rscript

#############################
# Libaries
#############################

library(here)
library(tidyverse)
library(pegas)
library(knitr)
library(plotly)
library(ggridges)
library(qqman)

# Might be some problems loading `pegas` on MacOS Big Sur due to `gdal` and `sf` dependencies.
# Following instructions here to resolve: https://github.com/r-spatial/sf/issues/1536#issuecomment-727342736

#############################
# Paths
#############################

## Latest plot path
plot_path = here::here("plots", "20210221_batch")
dir.create(plot_path)


#############################
# Functions
#############################

## Read `.afreq` files from `Plink`
read_afreq <- function(file){
  out = read.table(file, header = T, comment.char = "") %>%
    dplyr::rename(CHR = X.CHROM,
                  SNP = ID) %>% 
    dplyr::mutate(CHR = as.integer(CHR))
  
  return(out)
}

## Get manhattan plot
get_man <- function(df, trait, title, chr, bp, snp, p){
  # Plot
  manhattan(df, chr=chr, bp=bp, snp=snp, p=p,
            cex = 0.6,
            c(pal_primary[names(pal_primary) == trait],
              pal_secondary[names(pal_secondary) == trait]),
            main = title)
}


#############################
# Parameters
#############################

## Clumping target parameter
clump_param = "r2-0.1_kb-1000"

# Factor levels for `trait` (or `PHENO`)
trait_levels = c("hei",
                 "bmi", 
                 "edu", 
                 "int", 
                 "ibd", 
                 "pig", 
                 "scz", 
                 "dep", 
                 "fgl", 
                 "myc", 
                 "ldl", 
                 "plt")
names(trait_levels) = trait_levels

trait_levels_verb =  c("Height",
                       "BMI",
                       "Educational attainment",
                       "Intelligence",
                       "IBD",
                       "Pigmentation",
                       "Schizophrenia",
                       "Depression",
                       "Fasting glucose",
                       "Myocardial infarction",
                       "LDL level",
                       "Platelet counts")
names(trait_levels_verb) = trait_levels_verb

# Factor levels for `HIT_CONTROL`
hit_control_levels = c("hit", "control")

# Create vectors for recoding traits with full names and vice versa
recode_vec = c("hei" = "Height",
               "bmi" = "BMI",
               "edu" = "Educational attainment",
               "int" = "Intelligence",
               "ibd" = "IBD",
               "pig" = "Pigmentation",
               "scz" = "Schizophrenia",
               "dep" = "Depression",
               "fgl" = "Fasting glucose",
               "myc" = "Myocardial infarction",
               "ldl" = "LDL level",
               "plt" = "Platelet counts")

rev_recode_vec = names(recode_vec)
names(rev_recode_vec) = recode_vec

# Colour palettes
pal_primary = c("Height" = "#FC4E07",
                "BMI" = "#FFBF00",
                "Educational attainment" = "#0BC166",
                "Intelligence" = "#00AFBB",
                "IBD" = "#D84797",
                "Pigmentation" = "#360568",
                "Schizophrenia" = "#4E30C2",
                "Depression" = "#C72084",
                "Fasting glucose" = "#0B63E6",
                "Myocardial infarction" = "#15B096",
                "LDL level" = "#F5990F",
                "Platelet counts" = "#E84141")
pal_secondary = c("Height" = "#B63502",
                  "BMI" = "#B88A00",
                  "Educational attainment" = "#07743D",
                  "Intelligence" = "#00727A",
                  "IBD" = "#972062",
                  "Pigmentation" = "#1F033A",
                  "Schizophrenia" = "#342183",
                  "Depression" = "#8C175D",
                  "Fasting glucose" = "#07439C",
                  "Myocardial infarction" = "#0D6E5D",
                  "LDL level" = "#B06D07",
                  "Platelet counts" = "#A31414")
