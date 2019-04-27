#!/usr/bin/env R
library(scImpute)
args <- commandArgs(TRUE)
count.file <- as.character(args[1]) 
output <- as.character(args[2])
impute.dir <- as.character(args[3])
threads <- as.character(args[4])
scimpute(count_path = count.file, infile = "csv", outfile = "txt", type = "count",
        labeled = FALSE, drop_thre = 0.5, ncores = as.integer(threads)/2, Kcluster = 3, out_dir = impute.dir)
# data <- read.table(count.file, sep = "\t", header = T, stringsAsFactors = F)
q()
