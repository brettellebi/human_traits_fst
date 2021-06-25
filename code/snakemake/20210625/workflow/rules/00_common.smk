#####################
# Libraries
#####################

import os.path
import pandas as pd

#####################
# Variables
#####################

# Config file
configfile: "code/snakemake/20210625/config/config.yaml"

CHRS  = [str(chrom) for chrom in list(range(1,23))]

CHRS_TEST = ['22']

# Date of GWAS data collection
DATE_OF_COLLECTION = config["date_of_gwas_collection"]

# Date of dbSNP data collection
DATE_OF_DBSNP_COLLECTION = config["date_of_dbsnp_collection"]

# Extensions
GZ_EXTS = config["gz_exts"]

# Traits
EFO_IDS = ["EFO_0004339",
           "EFO_0004340",
           "EFO_0004784",
           "EFO_0004337",
           "EFO_0003767",
           "EFO_0003784",
           "EFO_0007009",
           "EFO_0003949",
           "EFO_0009764",
           "EFO_0003924",
           "EFO_0007822",
           "EFO_0000692",
           "EFO_0003761",
           "EFO_0004465",
           "EFO_0000612",
           "EFO_0004611",
           "EFO_0004309"]

EFO_IDS_TEST = ["EFO_0000612"]

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()


#####################
# Rules
#####################

## Only to run once, before anything else, to create file with traits and IDs

rule trait_ids_file:
    output:
        "config/trait_ids.tsv"
    container:
        ""
    