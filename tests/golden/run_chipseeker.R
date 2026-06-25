#!/usr/bin/env Rscript
# Generate the value-exact golden: ChIPseeker annotatePeak on anno.gtf + peaks.bed.
suppressMessages({library(ChIPseeker); library(txdbmaker); library(GenomicRanges)})
args <- commandArgs(trailingOnly=TRUE)
gtf <- ifelse(length(args)>=1, args[1], "anno.gtf")
bed <- ifelse(length(args)>=2, args[2], "peaks.bed")
out <- ifelse(length(args)>=3, args[3], "golden.tsv")
txdb <- suppressWarnings(makeTxDbFromGFF(gtf, format="gtf"))
pk <- readPeakFile(bed)
ann <- annotatePeak(pk, TxDb=txdb, tssRegion=c(-3000,3000), level="transcript", verbose=FALSE)
df <- as.data.frame(ann)
cols <- intersect(c("seqnames","start","end","V4","annotation","geneChr","geneStart",
  "geneEnd","geneLength","geneStrand","geneId","transcriptId","distanceToTSS"), colnames(df))
write.table(df[,cols], out, sep="\t", quote=FALSE, row.names=FALSE)
cat("wrote", out, "with", nrow(df), "rows;", ncol(df), "cols\n")
cat("annotations:\n"); print(df$annotation)
cat("distanceToTSS:\n"); print(df$distanceToTSS)
