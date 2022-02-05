build_GENCODE_gtf <- function(x) {

    ## https://www.gencodegenes.org/pages/data_format.html
    download.file(x, sprintf("%s/gtf.gz", tempdir()), mode="wb")
    gtf = read_tsv(sprintf("%s/gtf.gz", tempdir()), comment="#",
        col_names = c("chrm", "source", "feature_type", "start", "end",
            "score_not_used", "strand", "cds_phase", "additional"),
        col_types=cols(
            start=col_integer(),
            end=col_integer(),
            .default=col_character()))

    ## transcript
    g1 <- gtf[gtf$feature_type == "transcript", ]
    g1$transcript_id = str_match(
        g1$additional, 'transcript_id "([^"]*)"')[,2]
    g1$transcript_name = str_match(
        g1$additional, 'transcript_name "([^"]*)"')[,2]
    g1$transcript_type = str_match(
        g1$additional, 'transcript_type "([^"]*)"')[,2]
    g1$gene_id = str_match(g1$additional, 'gene_id "([^"]*)"')[,2]
    g1$gene_name = str_match(g1$additional, 'gene_name "([^"]*)"')[,2]
    g1$gene_type = str_match(
        g1$additional, 'gene_type "([^"]*)"')[,2]
    ## there is also transcript_name, is it useful? not included
    g1$level = str_match(g1$additional, 'level ([^;]*)')[,2]
    ## gene_status and transcript_status are obsolete after 25 and M11
    stopifnot(length(g1$transcript_id) == length(unique(g1$transcript_id)))

    ## exon
    g2 <- gtf[gtf$feature_type == "exon", ]
    g2$transcript_id <- str_match(
        g2$additional, 'transcript_id "([^"]*)"')[,2]
    g2$exon_id <- str_match(
        g2$additional, 'exon_id "([^"]*)"')[,2]
    g2$exon_number <- str_match(g2$additional, 'exon_number ([^;]*)')[,2]
    chrms <- paste0("chr",c(as.character(seq_len(max(as.integer(str_replace(
        sort(unique(g2$chrm)), "chr", "")), na.rm=TRUE))), c("X","Y","M")))
    gr <- GRanges(
        seqnames = g2$chrm, ranges=IRanges(g2$start, g2$end),
        strand = g2$strand, seqinfo=Seqinfo(chrms))
    mcols(gr)$exon_number <- as.integer(g2$exon_number)
    names(gr) <- g2$exon_id
    grl <- GRangesList(split(gr, g2$transcript_id)) # slow
    stopifnot(length(grl) == length(g1$transcript_id))
    stopifnot(all(sort(names(grl)) == sort(g1$transcript_id)))

    ## CDS
    g3 <- gtf[gtf$feature_type == "CDS", ]
    g3$transcript_id = str_match(
        g3$additional, 'transcript_id "([^"]*)"')[,2]
    g1$cdsStart <- g3$start[match(g1$transcript_id, g3$transcript_id)]
    g1$cdsEnd <- g3$end[match(g1$transcript_id, g3$transcript_id)]

    ## put together
    g1 <- g1[order(factor(g1$chrm, level=chrms), g1$start),]
    grl <- grl[g1$transcript_id]
    
    mcl <- g1[match(names(grl), g1$transcript_id), c(
        "chrm", "start", "end", "strand",
        "transcript_id", "transcript_type", "transcript_name",
        "gene_name", "gene_id", "gene_type", "source", "level",
        "cdsStart", "cdsEnd")]
    colnames(mcl)[2] <- "transcript_start"
    colnames(mcl)[3] <- "transcript_end"
    colnames(mcl)[4] <- "transcript_strand"
    mcols(grl) <- mcl
    grl
}

getTranscriptGRanges <- function(grl) {
    mcl <- mcols(grl)
    gr <- GRanges(
        seqnames = mcl$chrm, ranges = IRanges(
            mcl$transcript_start, mcl$transcript_end),
        strand = mcl$transcript_strand, seqinfo = seqinfo(grl))
    names(gr) <- mcl$transcript_id
    mcols(gr) <- mcl[,colnames(mcl)[!(colnames(mcl) %in% c(
        "chrm","transcript_id", "transcript_start",
        "transcript_end","transcript_strand"))]]
    gr <- sort(gr, ignore.strand = TRUE)
    gr
}



