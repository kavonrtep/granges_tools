#!/usr/bin/env Rscript
library(Biostrings)
library(parallel)

# test data
if (TRUE){
  s <- readDNAStringSet("/mnt/raid/454_data/original_sequencing_data/0_submissions_to_databases/221003_Psativum_CEN6_assembly/CEN6_ver_220406.fasta.gz")
  kmer_size <- 31
  window_size <- 20000
}

kmers <- Views(s[[1]], start = 1:(nchar(s) - kmer_size  - window_size + 2 ), width = kmer_size)
windows <- Views(s[[1]], start = kmer_size:(nchar(s) - window_size +1 ), width = window_size)
interval=1:50000
#x = mcmapply(matchPattern, kmers[interval], windows[interval], max.mismatch = 3, SIMPLIFY = FALSE,
#             mc.cores = 10, mc.preschedule = FALSE)
xrc = mapply(matchPattern, reverseComplement(kmers[interval]), windows[interval], max.mismatch = 5, SIMPLIFY = FALSE)
x = mapply(matchPattern, kmers[interval], windows[interval], max.mismatch = 5, SIMPLIFY = FALSE)
# get pairs of matches with correct coordinates

S = which(sapply(x, length)> 0)
Xpos = x[sapply(x, length)> 0]
Epos = sapply(Xpos, start)
coords_list = mapply(cbind, S, Epos)
coordsF = do.call(rbind, coords_list)

S = which(sapply(xrc, length)> 0)
Xpos = x[sapply(xrc, length)> 0]
Epos = sapply(Xpos, start)
coords_list = mapply(cbind, S, Epos)
coordsRC = do.call(rbind, coords_list)

coords = rbind(coordsF, coordsRC)




png("test.png", width = 3000, height = 1000)
plot(coords[,1], coords[,2], pch = 19, cex = 1)
dev.off()


writeXStringSet(subseq(s,min(interval), max(interval)), "test.fasta")