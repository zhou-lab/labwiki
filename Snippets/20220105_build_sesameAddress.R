#!/usr/bin/env Rscript-4.2.0

suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
args <- commandArgs(trailingOnly = TRUE)
tsv_fn <- args[1]
out_order_rds_fn <- args[2]

create_mask <- function(df) {
    unmapped <- ifelse(is.na(df$mapAS_A) | df$mapAS_A < 35 | (!is.na(df$mapAS_B) & df$mapAS_B < 35), 1, 0)
    masks <- with(df, data.frame(
        Probe_ID = Probe_ID,
        nonunique = ifelse((!unmapped) & (mapQ_A == 0 | (!is.na(mapQ_B) & mapQ_B == 0)), 1, 0),
        missing_target = ifelse((!unmapped) & (is.na(target) | (target != "CG")) & grepl("^cg", Probe_ID), 1, 0)))
    masks$control = ifelse(grepl("^ctl", df$Probe_ID), 1, 0)
    masks$design_issue = ifelse(grepl("^uk", df$Probe_ID), 1, 0)
    masks$unmapped <- ifelse(unmapped == 1 & masks$control != 1 & masks$design_issue != 1, 1, 0)
    masks$low_mapq <- with(df, ifelse((!is.na(mapQ_A)) & (mapQ_A < 30 | (!is.na(mapQ_B) & mapQ_B < 30)), 1, 0))
    masks$ref_issue <- with(masks, ifelse(unmapped == 1 | missing_target == 1, 1, 0))
    masks[c("Probe_ID","unmapped","missing_target","ref_issue","nonunique","low_mapq","control","design_issue")]
}

build_refManifest <- function(in_fn) {
    dfref <- read_tsv(in_fn, col_types=cols(CpG_chrm = col_character(), CpG_beg = col_integer(), CpG_end = col_integer(),
        address_A = col_integer(), address_B = col_integer(),
        target = col_character(), nextBase = col_character(), channel = col_character(),
        Probe_ID = col_character(), mapFlag_A = col_integer(), mapChrm_A = col_character(),
        mapPos_A = col_integer(), mapQ_A = col_integer(), mapCigar_A = col_character(),
        AlleleA_ProbeSeq = col_character(),
        mapNM_A = col_character(), mapAS_A = col_integer(), mapYD_A = col_character(),
        mapFlag_B = col_integer(), mapChrm_B = col_character(), mapPos_B = col_integer(),
        mapQ_B = col_integer(), mapCigar_B = col_character(), AlleleB_ProbeSeq = col_character(),
        mapNM_B = col_character(), mapAS_B = col_integer(), mapYD_B = col_character(), type = col_character()))
    ordering <- data.frame(Probe_ID = dfref$Probe_ID, M=dfref$address_B, U=dfref$address_A, col=factor(dfref$channel, level=c("G","R")), mask=FALSE)
    ordering$mask <- (create_mask(dfref)$ref_issue == 1)
    ordering
}

cat(sprintf("=== build address order %s ===\n", out_order_rds_fn))

## source("~/repo/labwiki/Snippets/20220117_buildSpeciesManifest.R")
## source("https://raw.githubusercontent.com/zhou-lab/labwiki/master/Snippets/20220302_sesameDataExtra_utils.R")
a <- list()
a$ordering <- build_refManifest(tsv_fn)
a$controls <- NULL
cat(sprintf("%d probes masked\n", sum(a$ordering$mask)))
saveRDS(a, file=out_order_rds_fn)

cat(sprintf("%d probes/rows in ordering\n", nrow(a$ordering)))
cat(sprintf("%d probes masked\n", sum(a$ordering$mask)))
cat(sprintf("%d red probes\n", sum(na.omit(a$ordering$col=="R"))))
cat(sprintf("%d grn probes\n", sum(na.omit(a$ordering$col=="G"))))
cat(sprintf("========================================\n"))
