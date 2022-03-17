## check https://github.com/zwdzwd/sesame/blob/e5fa29de3f895f5e37675400dc2ac09797a5f7d1/R/cnv.R#L252
## see also https://github.com/anastasia-lucas/hudson

#' Make a Manhattan plot to summarize EWAS
#'
#' @param pvals named vector of p-values, vector name is Probe ID.
#' @param
sesamePlot_Manhattan <- function(pvals) {

    pvals

    stopifnot(is(pvals, "numeric"))
    GenomicRanges::values(bin.coords)$bin.mids <-
        (start(bin.coords) + end(bin.coords)) / 2
    GenomicRanges::values(bin.coords)$bin.x <-
        seqstart[as.character(seqnames(bin.coords))] + bin.coords$bin.mids

    ## plot bin
    p <- ggplot2::qplot(bin.coords$bin.x / totlen,
        bin.signals, color=bin.signals, alpha=I(0.8))

    ## plot segment
    seg.beg <- (seqstart[sigs$chrom] + sigs$loc.start) / totlen
    seg.end <- (seqstart[sigs$chrom] + sigs$loc.end) / totlen
    p <- p + ggplot2::geom_segment(ggplot2::aes(x = seg.beg, xend = seg.end,
        y = sigs$seg.mean, yend=sigs$seg.mean), size=1.5, color='blue')

    ## chromosome boundary
    p <- p + ggplot2::geom_vline(xintercept=seqstart[-1]/totlen, alpha=I(0.5))

    ## chromosome label
    p <- p + ggplot2::scale_x_continuous(
        labels=seq.names, breaks=(seqstart+seqlen/2)/totlen) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=0.5))

    p <- p + ggplot2::scale_colour_gradient2(
        limits=c(-0.3,0.3), low='red', mid='grey', high='green',
        oob=scales::squish) + ggplot2::xlab('') + ggplot2::ylab('') +
        ggplot2::theme(legend.position="none")
    p
}
