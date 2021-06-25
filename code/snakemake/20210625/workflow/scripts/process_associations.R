#!/usr/bin/env Rscript

# Libraries

library(gwasrapidd)
library(tidyverse)

# Read in data

assocs = readRDS(snakemake@input[[1]])

# Extract key data
## Pull `associations` and `risk_alleles` DFs
assocs_df = assocs@associations %>%
    dplyr::left_join(assocs@risk_alleles, by = "association_id")
## Add EFO_ID
assocs_df$efo_id = snakemake@params[["efo_id"]]

# Find study IDs

studies_key = gwasrapidd::association_to_study(unique(assocs_df$association_id))
print("studies_key success")

# Get studies data

studies = gwasrapidd::get_studies(study_id = unique(studies_key$study_id))
## Create list of slots
studies_list = list(studies@studies,
                    studies@genotyping_techs,
                    studies@platforms,
                    studies@ancestries,
                    studies@ancestral_groups,
                    studies@countries_of_origin,
                    studies@countries_of_recruitment,
                    studies@publications)
## Recursively join into single df
studies_comb = plyr::join_all(studies_list)
## Join with `studies_key`
studies_df = dplyr::right_join(studies_key,
                               studies_comb,
                               by = "study_id")
print("studies_df success")

# Get trait data

traits = gwasrapidd::get_traits(efo_id = snakemake@params[["efo_id"]])
traits = traits@traits
print("traits success")

# Join with associations df
## Studies
out = dplyr::left_join(assocs_df,
                       studies_df,
                       by = "association_id")
## Traits
out = dplyr::left_join(out,
                       traits,
                       by = "efo_id")
print("join success")

# Write to file

readr::write_csv(out, snakemake@output[[1]])