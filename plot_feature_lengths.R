#!/usr/bin/env Rscript

library(optparse)

parser <- OptionParser()

option_list = list(
  make_option(c("-g", "--gff"), action = "store", type = "character",
              help = "path to gff file", default = NULL),
  make_option(c("-o", "--output_png"), action = "store", type = "character",
              help = "path to output png file", default = NULL),
  make_option(c("-a", "--attribute"), action = "store", type = "character",
              help = "attribute name", default = NULL),
  make_option(c("-v", "--value"), action = "store", type = "character",
              help = "attribute value", default = NULL),
  make_option(c("-N", "--nbins"), action = "store", type = "numeric",
              help = "minimal size of feature to be kept in output", default = 100)
)

description <- "Plot feature lengths histogram"
epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
args <- parse_args(parser, args = commandArgs(TRUE))


if (is.null(args$gff) | is.null(args$output_png) | is.null(args$attribute) | is.null(args$value)){
  message("ERROR: missing arguments")
  message("gff, output_png, attribute, value are required")
  quit()
}

suppressPackageStartupMessages({library(rtracklayer)})


g <- import(args$gff)
w <- width(g)[mcols(g)[,args$attribute] == args$value]
if (length(w) == 0){
  message("ERROR: no features found")
  quit()
}

N = length(w)
png(args$output_png, width = 800, height = 600)
par(mfrow = c(1,2))
# linear scale
hist(w, breaks = args$nbins, col = "blue",
     main = paste0(args$attribute," = ",  args$value, " (N = ", N, ")"),
     xlab = "Feature length")
# log scale
hist(log10(w), breaks = args$nbins, col = "blue",
     main = paste0(args$attribute," = ",  args$value, " (N = ", N, ")"),
     xlab = "Feature length - log10 scale")
suppressMessages(dev.off())

