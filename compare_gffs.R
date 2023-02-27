#!/usr/bin/env Rscript

library(optparse)


# parse command line arguments
parser <- OptionParser()
option_list <- list(
  make_option(c("-a", "--gff1"), action = "store", type = "character",
              help = "first gff", default = NULL),
  make_option(c("-b", "--gff2"), action = "store", type = "character",
              help = "second gff", default = NULL),
  make_option(c("-A", "--attr1"), action = "store", type = "character",
              help = "attribute to summarize for first gff", default = NULL),
  make_option(c("-B", "--attr2"), action = "store", type = "character",
              help = "attribute to summarize for second gff", default = NULL),
  make_option(c("-o", "--output_dir"), action = "store", type = "character",
              help = "output_directory", default = NULL)
)

description <- "Compare two gffs and create a table with overlaps"
epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))

suppressPackageStartupMessages({
  library(rtracklayer)
})





if (any(is.null(opt$gff1),
        is.null(opt$gff2),
        is.null(opt$attr1),
        is.null(opt$attr2),
        is.null(opt$output))){
  print_help(parser)
  stop('all argument are required!')
}

# just for testing
if (FALSE){
  f1 <- "/mnt/raid/454_data/MinION/analysis/centromere_assembly_Cameor/220329_CEN6_assembly_correction/tracks_recalculated/tandem_repeats_curated_track_220419.gff"
  f2 <- "/mnt/raid/users/petr/workspace/TideCluster/tmp2/CEN6_ver_220406_p40_P3000_c3_e25.tidecluster_02_mmseq_blast_cc18_75pid_w9_cc.gff"
  attr1 <- "Family"
  attr2 <- "Name"
  g1 <- import(f1, format = 'GFF3')
  g2 <- import(f2, format = "GFF3")
}

g1 <- import(opt$gff1, format = 'GFF3')
g2 <- import(opt$gff2, format = "GFF3")
attr1 <- opt$attr1
attr2 <- opt$attr2


# break it to non-overlaping regions
g1g2 <- disjoin(c(g1, g2))

g1_ovl <- findOverlaps(g1, g1g2)
g2_ovl <- findOverlaps(g2, g1g2)

g1g2$N1 <- "No_annotation"
g1g2$N1[to(g1_ovl)] <- mcols(g1)[, attr1][from(g1_ovl)]

g1g2$N2 <- "No_annotation"
g1g2$N2[to(g2_ovl)] <- mcols(g2)[, attr2][from(g2_ovl)]


match_tab <- sort(by(data = width(g1g2), INDICES = paste(g1g2$N1, g1g2$N2), FUN = sum), decreasing = FALSE)
match_df <- data.frame(gff1 = gsub(" .+", "", names(match_tab)),
                       gff2 = gsub("^.+ ","", names(match_tab)),
                       overlap = match_tab, stringsAsFactors = FALSE, row.names = NULL)
# sort by overlap
match_df <- match_df[order(match_df$overlap, decreasing = TRUE),]


# split into groups and sort by overlap
gff1_match <- split(match_df, match_df$gff1)
total1 <- sapply(gff1_match, function(x) sum(x$overlap))
order1 <- order(total1, decreasing = FALSE)
gff1_match <- gff1_match[order1]
total1 <- total1[order1]
stat1_out <-  do.call(rbind, lapply(gff1_match, function(x) rbind(data.frame(gff1="", gff2="", overlap=""), x, data.frame(gff1="", gff2="Total:", overlap=sum(x$overlap)))))





gff2_match <- split(match_df, match_df$gff2)
total2 <- sapply(gff2_match, function(x) sum(x$overlap))
order2 <- order(total2, decreasing = FALSE)
gff2_match <- gff2_match[order2]
total2 <- total2[order2]
stat2_out = do.call(rbind, lapply(gff2_match, function(x) rbind(data.frame(gff1="", gff2="", overlap=""), x, data.frame(gff1="", gff2="Total:", overlap=sum(x$overlap)))))
stat2_out <- stat2_out[, c("gff2", "gff1", "overlap")]

# export
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
write.table(match_df, file = file.path(opt$output_dir, "match_gff1_gff2_table.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(stat1_out, file = file.path(opt$output_dir, "match_by_gff1.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(stat2_out, file = file.path(opt$output_dir, "match_by_gff2.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

