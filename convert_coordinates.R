#!/usr/bin/env Rscript
library(optparse)
parser <- OptionParser()
option_list <- list(
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output genomic track", default = NULL),
  make_option(c("-i", "--input"), action = "store", type = "character",
              help = "input genomic tracks", default = NULL),
  make_option(c("-f", "--format"), action = "store", type = "character",
              help = "format of genomic track: GFF3 (default), BED, WIG or BigWig) ", default = "GFF3"),
  make_option(c("-c", "--conversion_table"), action = "store", type = "character",
              help = "table with coordinates for conversion", default = NULL)
)
description <- ""
epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")

opt <- parse_args(parser, args = commandArgs(TRUE))
# Required options
if (any(is.null(opt$output),
        is.null(opt$input),
        is.null(opt$conversion_table))){
  print_help(parser)
  stop('all argument are required!')
}

suppressPackageStartupMessages({
  library(rtracklayer)
})

# for testing
if(FALSE){
  gin <-  import("/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/211209_CEN6_ANALYSIS/satellites/annotation/tandem_repeats_curated_track_220218.gff",format = "GFF")
  conversion_table <- read.table("/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/220329_CEN6_assembly_correction/coordinates_220406", as.is=TRUE, sep="\t", col.names = c("seqname", "start", "end", "strand"))
}


gin <-  import(opt$input, format = opt$format)
conversion_table <- read.table(opt$conversion_table, as.is=TRUE, sep="\t", col.names = c("seqname", "start", "end", "strand"))
# assume that new sequence is simple concatenation

conversion_table$width <- conversion_table$end - conversion_table$start + 1
N <-  nrow(conversion_table)
conversion_table$new_start <-  cumsum(c(1, conversion_table$width))[1:N]
conversion_table$new_end <-  cumsum(conversion_table$width)
ct <- makeGRangesFromDataFrame(conversion_table)

p <-  findOverlaps(gin, ct)
gin_new <- GRangesList()
for (i in seq_along(ct)){
  f <-  conversion_table$new_start[i] - conversion_table$start[i]
  gin_part <-  gin[from(p)[to(p) == i]]
  gin_part_trimmed <- restrict(gin_part, conversion_table$start[i], conversion_table$end[i])
  gin_new[[i]] <- shift(gin_part_trimmed, f)
}

gall <- unlist(gin_new)

gall_merged <- reduce(gall, with.revmap = TRUE)

if (length(gall) == length(gall_merged)){
  # no close ranges - gall is ok
  gout <- gall
}else{
  ok_ranges <- unlist(gall_merged$revmap[sapply(gall_merged$revmap, length) == 1])
  ranges_to_merge <- gall_merged$revmap[sapply(gall_merged$revmap, length) > 1]
  merged_ranges <- GRangesList()
  j <- 0
  for (i in seq_along(ranges_to_merge)){
    NR <- length(ranges_to_merge[[i]])
    if (NR != 2){
      message("something is wrong")
      stop()
    }
    r_attrs <- elementMetadata(gall[ranges_to_merge[[i]]])
    if (nrow(unique(r_attrs))==2){
      # ranges have different annotation - they will not be merged
      ok_ranges <- append(ok_ranges, unlist(ranges_to_merge[[i]]))
    }else{
      # ranges will be merged:
      j <- j + 1
      m <- reduce(gall[ranges_to_merge[[i]]])
      elementMetadata(m) <- elementMetadata(gall[ranges_to_merge[[i]]][1])
      merged_ranges[[j]] <- m
    }
  }
  gout <- c(gall[ok_ranges], unlist(merged_ranges))
}

export(gout, con = opt$output, format = opt$format)
