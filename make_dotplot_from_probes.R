#!/usr/bin/env Rscript

library(optparse)
options_list <- list(
  make_option(c("-q","--query"), action = "store", type = "character",
              help = "query is blast output in tabular format", default = NULL),
  make_option(c("-s", "--subject"), action = "store", type = "character",
              help = "subject blast output in tabular format", default = NULL),
  make_option(c("-Q", "--query_chromosome_sizes"), action = "store", type = "character",
              help = "file with chromosome sizes for query(.fai format)", default = NULL),
  make_option(c("-S", "--subject_chromosome_sizes"), action = "store", type = "character",
              help = "file with chromosome sizes for subject(.fai format)", default = NULL),
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output png file", default = NULL)
)

opt_parser <- OptionParser(option_list = options_list, description = "Creates dotplot from two blast outputs")
opt <- parse_args(opt_parser, args = commandArgs(TRUE))

#library(tidyverse)

# blast tables are in format 6 with columns:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# for testing
if (FALSE){
  opt <- list(
    query = "/mnt/ceph/454_data/Pisum_pangenome/assemblies/JI2822__2023-06-22/analysis/painting_probes_CAM/all_oligos/JI2822_230622_x_oligos_CAMv2r2.blast_out",
    subject = "/mnt/ceph/454_data/Pisum_assembly_ver_2/assembly/230509_release_3/Analysis/painting_probes_CAM/all_oligos/CAMv2r3_x_oligos_CAMv2r2.blast_out",
    query_chromosome_sizes = "/mnt/ceph/454_data/Pisum_pangenome/assemblies/JI2822__2023-06-22/Pisum_sativum-JI2822-JIC_v1.0.fasta.fai",
    subject_chromosome_sizes = "/mnt/ceph/454_data/Pisum_assembly_ver_2/assembly/230509_release_3/Cameor_v2_release_3.fasta.fai",
    output = "test.png"
  )
}


query <- read.table(file = opt$query, sep="\t", header = FALSE,
                    col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                  "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
subject <- read.table(file = opt$subject, sep="\t", header = FALSE,
                      col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                    "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# chromosome sizes are in fai index format
chromosome_sizes_query <- read.table(file = opt$query_chromosome_sizes, sep="\t", header = FALSE,
                                     col.names = c('qseqid', 'qlen', 'V3','V4', 'V5'))
chromosome_sizes_subject <- read.table(file = opt$subject_chromosome_sizes, sep="\t", header = FALSE,
                                       col.names = c('sseqid', 'slen', 'V3','V4', 'V5'))


query$qstrand <- ifelse(query$sstart < query$send, "plus", "minus")
subject$sstrand <- ifelse(subject$sstart < subject$send, "plus", "minus")

# merge query and subject tables by qseqid - this is the probe name

QS <- merge(query, subject, by = "qseqid", suffixes = c("_query", "_subject"))

table(QS$qstrand, QS$sstrand)

QS$strand <- ifelse(QS$qstrand == QS$sstrand, "plus", "minus")

new_colnames <- c(
  sseqid_query = "qseqid",
  sseqid_subject = "sseqid",
  sstart_query = "qstart",
  send_query = "qend",
  sstart_subject = "sstart",
  send_subject = "send",
  strand = "sstrand"
)

blast_df <- QS[, names(new_colnames)]
colnames(blast_df) <- new_colnames

names(new_colnames) %in% names(QS)

# Create a color scheme for strand
color_scheme <- ifelse(blast_df$sstrand == "plus", "blue", "red")

# Create a factor with levels ordered by chromosome size
blast_df$qseqid <- factor(blast_df$qseqid, levels=chromosome_sizes_query$qseqid)
blast_df$sseqid <- factor(blast_df$sseqid, levels=chromosome_sizes_subject$sseqid)

# Calculate cumulative size for chromosome labels
chromosome_sizes_query$cumulative_qlen <- cumsum(as.numeric(chromosome_sizes_query$qlen))
chromosome_sizes_subject$cumulative_slen <- cumsum(as.numeric(chromosome_sizes_subject$slen))
chromosome_sizes_query$cumulative_start <- c(1, head(chromosome_sizes_query$cumulative_qlen, -1))
chromosome_sizes_subject$cumulative_start <- c(1, head(chromosome_sizes_subject$cumulative_slen, -1))



# Create a blank plot with appropriate limits and labels
query_label <- basename(opt$query)
subject_label <- basename(opt$subject)
png(opt$output, width=4000, height=4000, pointsize = 40)
plot(1, 1, xlim =c(1,  max(chromosome_sizes_query$cumulative_qlen)), ylim =c(1, max(chromosome_sizes_subject$cumulative_slen)),
     xlab = query_label, ylab =  subject_label, type = "n", axes = FALSE)

# Draw lines for each chromosome
par(lwd=5)
abline(v = c(1, chromosome_sizes_query$cumulative_qlen), lty = 2, col = "grey")


abline(h = c(1, chromosome_sizes_subject$cumulative_slen), lty = 2, col = "grey")


# Add text labels for chromosomes

axis(side = 1, at = chromosome_sizes_query$cumulative_qlen - chromosome_sizes_query$qlen / 2, labels = chromosome_sizes_query$qseqid, tick = FALSE)
axis(side = 2, at =chromosome_sizes_subject$cumulative_slen - chromosome_sizes_subject$slen / 2, labels = chromosome_sizes_subject$sseqid, tick = FALSE)

#
# Calculate the start and end of each hit relative to the start of the chromosome
blast_df$qstart_relative <- blast_df$qstart + chromosome_sizes_query$cumulative_start[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]
blast_df$qend_relative <- blast_df$qend + chromosome_sizes_query$cumulative_start[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]
blast_df$sstart_relative <- blast_df$sstart + chromosome_sizes_subject$cumulative_start[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]
blast_df$send_relative <- blast_df$send + chromosome_sizes_subject$cumulative_start[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]

# Draw lines for each BLAST hit
segments(blast_df$qstart_relative, blast_df$sstart_relative, blast_df$qend_relative, blast_df$send_relative, col = color_scheme)

# Add a legend
#legend("topright", legend = c("plus", "minus"), fill = c("red", "blue"), title = "Strand")
dev.off()


# export the data as table in PAF format
# columns in PAF format
# 1. Query sequence name - qseqid
# 2. Query sequence length - qlen
# 3. Query start coordinate (0-based; BED-like; closed) - qstart
# 4. Query end coordinate (0-based; BED-like; open) - qend
# 5. Relative strand: "+" or "-" convert from sstrand -> strand
# 6. Target sequence name - sseqid
# 7. Target sequence length - slen
# 8. Target start coordinate on original strand (0-based) - sstart
# 9. Target end coordinate on original strand (0-based) - send
# 10. Number of matching bases in the mapping - use length
# 11. Number bases that map as "M" - use length
# 12. Mapping quality (0-255 with 255 for missing) - use 255

blast_df$qlen <- chromosome_sizes_query$qlen[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]
blast_df$slen <- chromosome_sizes_subject$slen[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]
blast_df$length <- abs(blast_df$qend - blast_df$qstart)
blast_df$strand <- ifelse(blast_df$sstrand == "plus", "+", "-")
blast_df$mapq <- 255


if (is.null(opt$output_paf)) {
  opt$output_paf <- paste0(opt$output, ".paf")
}

PAF <- blast_df[, c("qseqid", "qlen", "qstart", "qend", "strand", "sseqid", "slen", "sstart", "send", "length", "length", "mapq")]
colnames(PAF) <- c("qseqid", "qlen", "qstart", "qend", "strand", "sseqid", "slen", "sstart", "send", "length", "m", "mapq")
write.table(PAF, file = opt$output_paf, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
