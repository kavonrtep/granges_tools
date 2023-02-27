#!/usr/bin/env Rscript
library(optparse)
option_list <-list(
  make_option(c("-o", "--output_bigwig"), action = "store", type = "character",
              help = "output bigwig", default = NULL),
  make_option(c("-i", "--input_bedgraph"), action = "store", type = "character",
              help = "input bedgraph", default = NULL),
  make_option(c("-s", "--chrom_sizes"), action = "store", type = "character",
              help = "chrom_sizes", default = NULL)
)
description <- "Converts bedgraph to bigwig"
parser <- OptionParser(option_list = option_list, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))
if (any(is.null(opt$output_bigwig),
        is.null(opt$input_bedgraph),
        is.null(opt$chrom_sizes))){
  print_help(parser)
  stop('all argument are required!')
}
suppressPackageStartupMessages(library(rtracklayer))
bedgraph <- import(opt$input_bedgraph)
chrom_sizes <- read.table(opt$chrom_sizes, header = FALSE, stringsAsFactors = FALSE)
chrom_sizes <- chrom_sizes[,1:2]
chrom_sizes <- setNames(chrom_sizes[,2], chrom_sizes[,1])
# some ranges in bedgraph could be bigger that chrom_sizes, remove these ranges
bedgraph <- bedgraph[end(bedgraph) < chrom_sizes[as.character(seqnames(bedgraph))]]

seqlengths(bedgraph) <- chrom_sizes[seqlevels(bedgraph)]
export(bedgraph, format = "BigWig", con = opt$output_bigwig)
