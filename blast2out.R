#!/usr/bin/env Rscript

f = commandArgs(T)[1]
bl = read.table(f, col.names = c(
                     strsplit("qseqid qlen qstart qend sstrand sseqid slen sstart send nident length evalue bitscore"," ")[[1]]
                   ))
bl_out = bl[,1:9]
bl_out$sstrand = ifelse(bl$sstrand=="plus", "+", "-")
bl_out$sstart = ifelse(bl$sstart > bl$send, bl$send, bl$sstart)
bl_out$send = ifelse(bl$sstart < bl$send, bl$send, bl$sstart)

bl_out$PID = round(bl$nident/bl$length * 100,1)

write.table(bl_out, paste0(f,".out"), sep=" ", col.names = FALSE, row.names = FALSE, quote = FALSE)
L = bl$length
write.table(bl_out[L>100,], paste0(f, "_100.out"), sep=" ", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(bl_out[L>5000,], paste0(f, "_5k.out"), sep=" ", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(bl_out[L>20000,], paste0(f, "_20k.out"), sep=" ", col.names = FALSE, row.names = FALSE, quote = FALSE)

