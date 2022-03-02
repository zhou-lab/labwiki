## a convenience function, only works on Mac
sesameDataClearHub <- function() {
    unlink("~/Library/Caches/org.R-project.R/R/ExperimentHub/", recursive=TRUE)
}

manifest_base <- "https://github.com/zhou-lab/InfiniumManifestsV"
manifest_base_default_version <- 1

#' download Infinium manifest from Github repositories
#'
#' @param platform Mammal40, MM285, EPIC, and HM450
#' @param genome hg38, mm10 etc.
#' @param version release version, default is the latest
#' @return tibble
#' @importFrom readr read_tsv
#' @importFrom readr cols
#' @importFrom readr col_integer
#' @importFrom readr col_character
#' @examples
#' mft <- sesameData_getManifestDF("Mammal40")
#' @export
sesameData_getManifestDF <- function(platform, genome=NULL,
    version = manifest_base_default_version) {
    
    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)
    
    title <- sprintf("InfiniumManifestV%d_%s_%s", version, platform, genome)
    data <- sesameDataGet_checkEnv(title)
    if (is.null(data)) {
        u1 <- sprintf(
            "%s%d/raw/main/%s/%s.tsv.gz",
            manifest_base, version, platform, genome)
        stopifnot(valid_url(u1))
        data <- read_tsv(u1,
            col_types=cols(CpG_beg=col_integer(), CpG_end=col_integer(),
                address_A=col_integer(), address_B=col_integer(),
                .default=col_character()))
        sesameDataGet_assignEnv(title, data)
    }
    data
}

#' download Infinium manifest from repositories and return GRanges
#'
#' @param platform Mammal40, MM285, EPIC, and HM450
#' @param genome hg38, mm10 etc.
#' @param decoy whether to include probes mapped to decoy sequence
#' @param version release version, default is the latest
#' @param check_EH return the ExperimentHub version if present
#' @param columns additional columns to add from the manifest to mcols
#' @return GRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @examples
#' gr <- sesameData_getManifestGRanges("Mammal40")
#' @export
sesameData_getManifestGRanges <- function(
    platform, genome = NULL, version = manifest_base_default_version,
    check_EH = TRUE, decoy = FALSE, columns = NULL) {

    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)

    ## check EH first
    addr <- sesameDataGet(sprintf("%s.address", platform))
    if (genome %in% names(addr)) { return(addr[[genome]]) }
    
    df <- sesameData_getManifestDF(platform, genome=genome, version=version)

    chrms <- df$CpG_chrm
    chrms <- chrms[!is.na(chrms)]
    if (genome %in% c("mm10","mm39","hg19","hg38")) {
        if (decoy) {
            chrms <- c(
                guess_chrmorder(chrms[!grepl("_", chrms)]),
                sort(unique(chrms[grepl("_", chrms)])))
        } else {
            chrms <- guess_chrmorder(chrms[!grepl("_", chrms)])
        }
    } else {
        chrms <- sort(unique(chrms))
    }
    df <- df[!is.na(df$CpG_chrm) & !is.na(df$CpG_beg) & !is.na(df$CpG_end),]
    df <- df[df$CpG_chrm %in% chrms,]
    gr <- GRanges(df$CpG_chrm,
        IRanges(df$CpG_beg+1, df$CpG_end),
        strand = ifelse(df$mapFlag_A=="0", "+", "-"),
        seqinfo = Seqinfo(chrms))
    if (length(columns) > 0) {
        mcols(gr) <- df[,columns] }
    names(gr) <- df$Probe_ID
    sort(gr, ignore.strand = TRUE)
}
