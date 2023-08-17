#!/usr/bin/env Rscript

library(optparse)

parser <- OptionParser()
option_list <- list(
  make_option(c("-g", "--gff3"), action = "store", type = "character",
              help = "gff3  with DANTE filtered output", default = NULL),
  make_option(c("-s", "--reference_sequence"), action = "store", type = "character",
              help = "reference sequence as fasta",
              default = NULL),
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "Output fasta file", default = NULL)
)

parser <- OptionParser(option_list = option_list)

opt <- parse_args(parser, args = commandArgs(TRUE))
suppressPackageStartupMessages({ library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
  library(parallel)
})




g <- import(opt$gff3)
s <- readDNAStringSet(opt$reference_sequence)
names(s) <- gsub(" .+$", "", names(s))
table(g$Final_Classification)
LINE <- g[g$Final_Classification == "Class_I|LINE"]

LINEplus <- LINE
offset <- round(runif(n=length(LINEplus), min = 0 , max= 200))
start(LINEplus) <- start(LINEplus) - 4000 - offset
end(LINEplus) <- end(LINEplus) + 4000 + offset
exc1 <- start(LINEplus) < 1
exc2 <- end(LINEplus) > seqlengths(s)[as.character(seqnames(LINEplus))]

LINEplus <- LINEplus[!(exc1 | exc2)]



LINEplus3 <- reduce(LINEplus)
s_line_plus3 <- getSeq(s, LINEplus3)

s_line_plus_parts5 <- unlist(strsplit(as.character(s_line_plus3), paste0("(?<=\\G.{", 100, "})"), perl=TRUE))

sout <- DNAStringSet((s_line_plus_parts5))
names(sout) <- seq_along(sout)
sout <- sout[nchar(sout)==100]

writeXStringSet(sout, filepath = opt$output)










