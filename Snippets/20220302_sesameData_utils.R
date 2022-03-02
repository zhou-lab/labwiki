## a convenience function, only works on Mac
sesameDataClearHub <- function() {
    unlink("~/Library/Caches/org.R-project.R/R/ExperimentHub/", recursive=TRUE)
}

manifest_base <- "https://github.com/zhou-lab/InfiniumManifestsV"
manifest_base_default_version <- 1
anno_base <- "https://github.com/zhou-lab/InfiniumAnnotationV"
anno_base_default_version <- 1

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


download_file <- function(title, version, dest_dir) {
    dest_file <- sprintf("%s/%s", dest_dir, title)
    dir.create(sprintf("%s/%s", dest_dir, dirname(title)),
        recursive = TRUE, showWarnings = FALSE)
    url <- sprintf("%s%d/raw/main/%s", anno_base, version, title)
    stopifnot(valid_url(url))
    download.file(url, dest_file, mode="wb")
    stopifnot(file.exists(dest_file) && file.info(dest_file)$size > 0)
    
    invisible(list(
        url = url,
        dest_dir = dest_dir,
        dest_file = dest_file,
        file_name = title))
}

#' Retrieve additional annotation files
#'
#' From the Infinium annotation website associated github repo
#' e.g., 
#' https://github.com/zhou-lab/InfiniumAnnotationV1
#'
#' The default version number should always work. One need to
#' refer to the actual repo to see which one of the other versions
#' also work.
#' 
#' See also
#' http://zwdzwd.github.io/InfiniumAnnotation
#'
#' @param title title of the annotation file
#' @param version version number
#' @param dest_dir if not NULL, download to this directory
#' @return annotation file
#' @examples
#'
#' ## avoided testing these function as they use external resources
#' if (FALSE) {
#' mft <- sesameData_getAnno("HM27/HM27.hg19.manifest.tsv.gz")
#' annoI <- sesameData_getAnno("EPIC/EPIC.hg19.typeI_overlap_b151.rds")
#' annoS <- sesameData_getAnno("EPIC/EPIC.hg19.snp_overlap_b151.rds")
#' sesameData_getAnno("test/3999492009_R01C01_Grn.idat", dest_dir = tempdir())
#' }
#' 
#' @export
sesameData_getAnno <- function(
    title, version = anno_base_default_version, dest_dir = NULL) {

    if (!is.null(dest_dir)) {
        return(download_file(title, version, dest_dir))
    }
    
    download_path <- sprintf("%s%d/raw/main/%s", anno_base, version, title)
    if (!valid_url(download_path)) {
        message(sprintf("File not available %s.", download_path))
        return(NULL)
    }

    if (endsWith(title, ".tsv.gz")) {
        z <- gzcon(url(download_path))
        raw <- textConnection(readLines(z))
        close(z)
        message("Retrieving annotation from ",download_path,
            "... ", appendLF = FALSE)
        anno <- read.table(raw, header=TRUE, sep="\t")
        close(raw)
        message("Done.")
    } else if (endsWith(title, ".rds")) {
        message("Retrieving annotation from ",download_path,
            "... ", appendLF = FALSE)
        anno <- readRDS(url(download_path))
        message("Done.")
    }
    anno
}


