#!/usr/bin/env Rscript
library(optparse)
parser <- OptionParser()
option_list <- list(
  make_option(c("-a", "--attribute_name"), action = "store", type = "character",
              help = "gff will be split by specified attribute", default = NULL),
  make_option(c("-g", "--input_gff"), action = "store", type = "character",
              help = "input gff", default = NULL),
  make_option(c("-w", "--min_width"), action = "store", type = "numeric",
              help = "minimal size of feature to be kept in output", default = 0),
  make_option(c("-o", "--output_prefix"), action = "store", type = "character",
              help = "output file base name", default = NULL)

)
description <- "Split gff file to multiple files based on given attribute of gff "
epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))


fn_sanitize <- function(filename, replacement = "_") {
  illegal <- "[/\\?<>\\:*|\":]"
  control <- "[[:cntrl:]]"
  reserved <- "^[.]+$"

  filename <- gsub(illegal, replacement, filename)
  filename <- gsub(control, replacement, filename)
  filename <- gsub(reserved, replacement, filename)
  return(filename)
}


suppressPackageStartupMessages(library(rtracklayer))
gff <- import(opt$input_gff)
BN <- opt$output_prefix
min_width <- as.numeric(opt$min_width)


gff_min_width <- sort(gff[width(gff) >= min_width], by = ~ seqnames * start)
head(gff_min_width)

if (is.null(opt$attribute_name)){
  export(gff_min_width, con = paste(BN, "gff3", sep = "."), format = "gff3")
}else{
  if (opt$attribute_name == "seqnames"){
    f <- seqnames(gff_min_width)
  }else{
    f <- mcols(gff_min_width)[,opt$attribute_name]
  }
  gff_min_width_parts <- split(gff_min_width, f = f)
  out_fname <- fn_sanitize(names(gff_min_width_parts))
  for (i in seq_along(out_fname)){
    export(gff_min_width_parts[[i]],
           format="gff3",
           con=paste0(opt$output_prefix, "_", out_fname[i], ".gff"
           )
    )
  }
}
