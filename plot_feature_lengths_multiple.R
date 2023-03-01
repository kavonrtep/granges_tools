#!/usr/bin/env Rscript
library(optparse)

parser <- OptionParser()

option_list = list(
  make_option(c("-g", "--gff"), action = "store", type = "character",
              help = "path to gff file", default = NULL),
  make_option(c("-o", "--output_pdf"), action = "store", type = "character",
              help = "path to output pdf file", default = NULL),
  make_option(c("-a", "--attribute"), action = "store", type = "character",
              help = "attribute name", default = NULL),
  make_option(c("-n", "--number_of_features"), action = "store", type = "character",
              help = "Number of the most abundand features to be ploted", default = 10),
  make_option(c("-N", "--nbins"), action = "store", type = "numeric",
              help = "number of bins for histogram ", default = 100)
)

description <- "Plot feature lengths histogram for the most abundand features"
epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
args <- parse_args(parser, args = commandArgs(TRUE))


if (is.null(args$gff) | is.null(args$output_pdf) | is.null(args$attribute)){
  message("ERROR: missing arguments")
  message("gff, output_png, attribute are required")
  quit()
}

suppressPackageStartupMessages({library(rtracklayer)})

g <- import(args$gff)
feature_totat_width <- sort(by(width(g), INDICES = mcols(g)[,args$attribute], FUN = sum), decreasing = TRUE)

pdf(args$output_pdf, width=15, height=5, pointsize=12)
for (i in 1:args$number_of_features){
  label <- names(feature_totat_width)[i]
  w <- width(g)[mcols(g)[,args$attribute] == label]
  W <- sum(w)
  if (length(w) == 0){
    break
  }
  N <- length(w)
  par(mfrow = c(1,4))
  # linear scale
  hist(w, breaks = args$nbins, col = "blue",
       main = paste0(args$attribute," = ",  label, ", N = ", N, ", ", "total width = ", W),
       xlab = "Feature length")
  # log scale
  h <- hist(log10(w), breaks = args$nbins, plot = FALSE)
  plot(h$mids, h$counts, type = "p", col = "blue",
       xlab = "Feature length - log10 scale", log='y', ylab="Frequency", cex=1.5, pch=16)
  stripchart(w, method = "jitter", vertical = TRUE, pch = 16, col = "#0000FF55",
             ylab = "Feature length", cex=1.5, jitter=0.5)
  stripchart(w, method = "jitter", vertical = TRUE, pch = 16, col = "#0000FF55",
             ylab = "Feature length", log='y', cex=1.5, jitter=0.5)
}
suppressMessages(dev.off())

