#!/usr/bin/env Rscript
library(optparse)
parser <- OptionParser()
option_list <- list(
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output bigwig", default = NULL),
  make_option(c("-g", "--gff"), action = "store", type = "character",
              help = "input gff", default = NULL),
  make_option(c("-r", "--reference"), action = "store", type = "character",
              help = "reference sequence", default = NULL),
  make_option(c("-w", "--window"), action = "store", type = "numeric",
              help = "window size for binning, default 100000 ", default = 100000)

)
description <- "Creates bigwig from gff, density is calculated in provided window "

get_density <- function(x, tw=1000000){
  cvg <- coverage(x)
  bins <- tileGenome(seqlengths(x), tilewidth = tw)
  d <- binnedAverage(unlist(bins), cvg, "coverage")
  d
}



epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))

if (any(is.null(opt$output),
        is.null(opt$gff),
        is.null(opt$gff),
        is.null(opt$reference))){
  print_help(parser)
  stop('all argument are required!')
}
suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
})

r <- readDNAStringSet(opt$reference)
g <- import(opt$gff)
sl <- seqlengths(r)
rm(r)
seqlengths(g) <- sl
d <- get_density(g, tw=opt$window)
d$score=d$coverage

export(d, format = "BigWig", con = opt$output)