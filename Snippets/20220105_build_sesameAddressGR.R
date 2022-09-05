#!/usr/bin/env Rscript-4.2.0

suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
args <- commandArgs(trailingOnly = TRUE)
tsv_fn <- args[1]
out_gr_fn <- args[2]
decoy <- args[3]

cat(sprintf("=== build GRanges %s ===\n", out_gr_fn))
guess_chrmorder <- function(chrms) {
    chrms1 <- chrms[!(chrms %in% c("chrX","chrY","chrM"))]
    paste0("chr",c(as.character(seq_len(max(as.integer(str_replace(
        sort(unique(chrms1)), "chr", "")), na.rm=TRUE))), c("X","Y","M")))
}

suppressMessages(df <- read_tsv(tsv_fn))
for (col in c("CpG_beg","CpG_end","address_A","address_B")) {
    if (col %in% colnames(df)) {
        df[[col]] <- as.integer(df[[col]])
    }
}

idx <- is.na(df$CpG_chrm)
df$CpG_chrm[idx] <- "*"
df$CpG_beg[idx] <- -1
df$CpG_end[idx] <- 0

chrms <- df$CpG_chrm
chrms <- chrms[!is.na(chrms)]
if (decoy == "decoy") {
    suppressWarnings(chrms <- c(guess_chrmorder(chrms[!grepl("_", chrms)]),
        sort(unique(chrms[grepl("_", chrms)]))))
} else if (decoy=="nodecoy") {
    chrms <- guess_chrmorder(chrms[!grepl("_", chrms)])
} else { # general
    chrms <- sort(unique(chrms))
}
chrms <- c(chrms, "*")
## df <- df[!is.na(df$CpG_chrm) & !is.na(df$CpG_beg) & !is.na(df$CpG_end),]
## df <- df[df$CpG_chrm %in% chrms,]
gr <- GRanges(df$CpG_chrm,
    IRanges(df$CpG_beg+1, df$CpG_end),
    strand = ifelse(is.na(df$mapFlag_A), "*", ifelse(df$mapFlag_A=="0", "+", "-")),
    seqinfo = Seqinfo(chrms))
## mcols(gr) <- df[,columns]
names(gr) <- df$Probe_ID
gr <- sort(gr, ignore.strand = TRUE)
saveRDS(gr, file=out_gr_fn)

cat(sprintf("%d probes in GRanges.\n", length(gr)))
cat(sprintf("%d probes belong to chromosome *.\n", sum(seqnames(gr)=="*")))
cat(sprintf("%d probes belong to decoy chromosomes.\n", sum(grepl("_", seqnames(gr)))))

cat(sprintf("==============================\n"))
