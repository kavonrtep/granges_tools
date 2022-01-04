#!/usr/bin/env Rscript
library(optparse)
parser <- OptionParser()
option_list <- list(
  make_option(c("-e", "--eggnog"), action = "store", type = "character",
              help = "oggnog output - out.emapper.annotations", default = NULL),
  make_option(c("-g", "--gff"), action = "store", type = "character",
              help = "input gff", default = NULL),
  make_option(c("-o", "--gff_out"), action = "store", type = "character",
              help = "output gff with eggnog annotation attached", default = NULL)

)
description <- ""

epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))
if (any(is.null(opt$eggnog),
        is.null(opt$gff),
        is.null(opt$gff_out))){
  print_help(parser)
  stop('all argument are required!')
}



suppressPackageStartupMessages(library(rtracklayer))
read_annot <- function(x) {
  l <- readLines(x)
  N <- max(grep("##", l)) - 8
  tbl <- read.table(x, sep = "\t", skip = 4, nrows = N, header = TRUE, comment.char = "", quote = "")
  colnames(tbl)[1] <- 'query'
  return(tbl)
}


if (file.exists(opt$gff_out)) {
  stop("output gff already exists")
}
g <- import(opt$gff)
egg <- read_annot(opt$eggnog)

index <- match(g$ID, egg$query)
g$seed_ortholog <- egg$seed_ortholog[index]
g$description <- egg$Description[index]
g$cog_category <- egg$COG_category[index]
g$pfam <- egg$PFAMs[index]
g$preferred_name <- egg$Preferred_name[index]
g$eggnog_ogs <- egg$eggNOG_OGs[index]
g[g$type == "mRNA"]
export(g, opt$gff_out, format = "gff3")
