#!/usr/bin/env Rscript
library(optparse)
parser <- OptionParser()
option_list <- list(
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output genomic track", default = NULL),
  make_option(c("-i", "--input"), action = "store", type = "character",
              help = "input genomic tracks", default = NULL),
  make_option(c("-f", "--format"), action = "store", type = "character",
              help = "format of genomic track: GFF3 (default), BED, WIG, BigWig or
              BedGraph) ", default = "GFF3"),
  make_option(c("-n", "--new_seqid"), action = "store", type = "character",
              help = "new seqid "),
  make_option(c("-c", "--conversion_table"), action = "store", type = "character",
              help = "table with coordinates for conversion", default = NULL)
)
description <- ""
epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description =
  description,
                       usage = "usage: %prog COMMAND [OPTIONS]")

opt <- parse_args(parser, args = commandArgs(TRUE))
# Required options
if (any(is.null(opt$output),
        is.null(opt$input),
        is.null(opt$new_seqid),
        is.null(opt$conversion_table))) {
  print_help(parser)
  stop('all argument are required!')
}

suppressPackageStartupMessages({
  library(rtracklayer)
})

# for testing
if (FALSE) {
  gin <- import("/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor
  /211209_CEN6_ANALYSIS/satellites/annotation/tandem_repeats_curated_track_220218.gff",
                format = "GFF")
  conversion_table <- read.table ("/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/220329_CEN6_assembly_correction/coordinates_220406",
                                  as.is = TRUE, sep = "\t", col.names = c("seqname", "start", "end", "strand"))

  gin <- import("tmp/C1P23_C1K_bs200.bigwig", format = "BigWig")
  gin <- gin[seqnames(gin) == "CEN6_ver_211209"]

  gin <- import("tmp/SAT_all_TideHunter.gff", format = "GFF3")
  conversion_table <- read.table("tmp/coordinates_220406", as.is = TRUE, sep = "\t",
col.names = c("seqname", "start", "end", "strand"))

  gin = makeGRangesFromDataFrame(data.frame(seqname = "CEN6_ver_211209", start=seq(1, 172063895, by=1000), end=seq(1, 172063895, by=1000)+1))
  gin$ID=start(gin)
}


gin <- import(opt$input, format = opt$format)
overlaping_tracks <- if (length(gin) == length(reduce(gin))) {
  FALSE
}else {
  TRUE
}

conversion_table <- read.table(opt$conversion_table, as.is = TRUE, sep = "\t",
                               col.names = c("seqname", "start", "end", "strand"))
# assume that new sequence is simple concatenation

conversion_table$width <- conversion_table$end - conversion_table$start + 1
N <- nrow(conversion_table)
conversion_table$new_start <- cumsum(c(1, conversion_table$width))[1:N]
conversion_table$new_end <- cumsum(conversion_table$width)
ct <- makeGRangesFromDataFrame(conversion_table)

p <- findOverlaps(gin, ct, ignore.strand = TRUE)
gin_new <- GRangesList()
for (i in seq_along(ct)) {
  f <- conversion_table$new_start[i] - conversion_table$start[i]
  gin_part <- gin[from(p)[to(p) == i]]
  gin_part_trimmed <- restrict(gin_part, conversion_table$start[i],
                               conversion_table$end[i])
  gin_new[[i]] <- shift(gin_part_trimmed, f)
  gin_new[[i]]$.BLOCK=i
}

gall <- unlist(gin_new)

# reduce only if format was gff or bed and original tracks does not overlap
gall_merged <- reduce(gall, with.revmap = TRUE)
if ((opt$format == "GFF3" | opt$format == "BED") & !overlaping_tracks) {
  if (length(gall) == length(gall_merged)) {
    # no close ranges - gall is ok
    gout <- gall
  }else {
    ok_ranges <- unlist(gall_merged$revmap[sapply(gall_merged$revmap, length) == 1])
    ranges_to_merge <- gall_merged$revmap[sapply(gall_merged$revmap, length) > 1]
    merged_ranges <- GRangesList()
    j <- 0
    for (i in seq_along(ranges_to_merge)) {
      NR <- length(ranges_to_merge[[i]])
      if (NR != 2) {
        message("something is wrong")
        stop()
      }
      r_attrs <- elementMetadata(gall[ranges_to_merge[[i]]])
      if (nrow(unique(r_attrs)) == 2) {
        # ranges have different annotation - they will not be merged
        ok_ranges <- append(ok_ranges, unlist(ranges_to_merge[[i]]))
      }else {
        # ranges will be merged:
        j <- j + 1
        m <- reduce(gall[ranges_to_merge[[i]]])
        elementMetadata(m) <- elementMetadata(gall[ranges_to_merge[[i]]][1])
        merged_ranges[[j]] <- m
      }
    }
    gout <- c(gall[ok_ranges], unlist(merged_ranges))
  }

}else {
  gout <- gall
}

seqlevels(gout)[seqlevels(gout) %in% seqlevels(ct)] <- opt$new_seqid
gout_sorted <- sort(gout, ignore.strand = TRUE)
if (!is.null(gout_sorted$ID)){
  if (any(duplicated(gout_sorted$ID))){
    message("some ID are duplicated, adding numerical suffix")
    gout_sorted$ID=paste0(gout_sorted$ID, "__", gout_sorted$.BLOCK)

  }
}
gout_sorted$.BLOCK=NULL
export(gout_sorted, con = opt$output, format = opt$format)
