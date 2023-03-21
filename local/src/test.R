#! /usr/bin/env Rscript

# 0. libraries ----
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))

# 1. create parser object ----
parser <- ArgumentParser(description = "Create lineplot from RNA counts data")

# add options
parser$add_argument("-c", "--count_matrix", action="store", required = T,
                    help = "feature count matrix path")
parser$add_argument("-m", "--metadata", action="store", required = T,
                    help = "metadata file path")
parser$add_argument("-f", "--feature", action="store", required = T,
                    help = "list of feature file path for count matrix subsetting")
parser$add_argument("-g", "--group_col", action="store", required = T,
                    help = "name of metadata column to use for averaging counts")
parser$add_argument("-d", "--condition_col", action="store", required = T,
                    help = "name of metadata column to use as line groups")
parser$add_argument("-t", "--time_point", action="store", required = T,
                    help = "name of metadata column to use as x axis in lineplot")
parser$add_argument("--time_order_col", action="store", default = NULL,
                    help = "name of metadata column to use for x axis ordering [Default: alphanumeric order]")
parser$add_argument("--pal", action="store", default=NULL,
                    help = "custom made palette [Default: RColorBrewer]")

# create args object
args <- parser$parse_args()
print(args)

# 2. read input files ----
counts <- read.delim(args$count_matrix, header = T, row.names = 1)
features <- scan(arguments$enst_file, character(), quote = "")
meta <- read.delim(arguments$metadata, header = T)

# filter TPM matrix with enst of interest
counts <- counts[row.names(counts) %in% features,]

print(counts)
print(features)
print(meta)
