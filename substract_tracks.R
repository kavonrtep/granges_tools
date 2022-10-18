#!/usr/bin/env Rscript
library(optparse)
option_list <- list(
  make_option(c("-A", "--A_track"), action = "store", type = "character",
              help = "The A track in GFF3 format", default = NULL),
  make_option(c("-B", "--B_track"), action = "store", type = "character",
              help = "The B track in GFF3 format", default = NULL),
  make_option(c("-o", "--output_track"), action = "store", type = "character",
              help = "The output track name", default = NULL)
)

description <- "This script substracts track B from track A and outputs the result to GFF3 format, annotation from track A is kept in the output"

parser <- OptionParser(option_list = option_list, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))


suppressPackageStartupMessages(library(rtracklayer))
message("Reading A track")
t1 <- import(opt$A_track, format = "gff3")
message("Reading B track")
t2 <- import(opt$B_track, format = "gff3")

# for testing:
if (FALSE){
  t1 <- import("/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/220329_CEN6_assembly_correction/tracks_recalculated/DANTE.gff")
  t2 <- import("/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/220329_CEN6_assembly_correction/tracks_recalculated/All_LTR_TE.gff3")
  t1=GRanges(rep("chr1",5), IRanges(c(1, 10, 20, 30, 40), c(5, 15, 25, 35, 45)))
  t2=GRanges(rep("chr1",2), IRanges(c(42,55), c(43, 55)))


}


# it assume no overlap within the tracks
# finde overlapping regions, t1 s higher priority that t2

t1$priority <- 2
t2$priority <- 1
suppressWarnings({t12 <- c(t1, t2)})  # suppress warnings about unique seqlevels
message("Comparing tracks")
t12_disjoint <- disjoin(t12,with.revmap = TRUE, ignore.strand = TRUE)

# parts of t2 which ovelapts with t1 are removed

priority <- t12$priority
message("Removing overlapping regions")
p <- lapply(as.list(t12_disjoint$revmap), FUN = function(x)priority[x])
ok <- sapply(p, function(x)all(x==2))


t2ok = t12_disjoint[ok]
t2ok$revmapc = unlist(t2ok$revmap)

# add all atrributes from t2 to t2ok, use revmapc to find the right attributes
mcols(t2ok) <- mcols(t12)[t2ok$revmapc,, drop = FALSE]
t2ok$priority <- NULL
message("Sorting and writing output")
# sort by chr and start
t2ok_sorted <- sort(t2ok, by = ~ seqnames + start)
export(t2ok, opt$output_track, format = "gff3")
