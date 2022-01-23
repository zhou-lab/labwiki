create_mask <- function(df) {
    unmapped <- ifelse(is.na(df$mapAS_A) | df$mapAS_A < 35 | (!is.na(df$mapAS_B) & df$mapAS_B < 35), 1, 0)
    masks <- with(df, data.frame(
         Probe_ID = Probe_ID,
         nonunique = ifelse((!unmapped) & (mapQ_A == 0 | (!is.na(mapQ_B) & mapQ_B == 0)), 1, 0),
         missing_target = ifelse((!unmapped) & (target != "CG") & grepl("^cg", Probe_ID), 1, 0)))
    masks$control = ifelse(grepl("^ctl", df$Probe_ID), 1, 0)
    masks$design_issue = ifelse(grepl("^uk", df$Probe_ID), 1, 0)
    masks$unmapped <- ifelse(unmapped == 1 & masks$control != 1 & masks$design_issue != 1, 1, 0)
    masks$low_mapq <- with(df, ifelse((!is.na(mapQ_A)) & (mapQ_A < 30 | (!is.na(mapQ_B) & mapQ_B < 30)), 1, 0))
    masks$ref_issue <- with(masks, as.integer(unmapped * missing_target))
    masks
}

build_speciesManifest <- function(in_dir, reference) {
  dfsp <- read_excel('~/samplesheets/2020/20201202_310_species_EnsemblVertebrates3.xlsx')

  dfref <- read_tsv(sprintf("%s/%s.tsv.gz", in_dir, reference), col_types=cols(CpG_chrm = col_character(), CpG_beg = col_integer(), CpG_end = col_integer(), address_A = col_integer(), address_B = col_integer(), target = col_character(), nextBase = col_character(), channel = col_character(), Probe_ID = col_character(), mapFlag_A = col_integer(), mapChrm_A = col_character(), mapPos_A = col_integer(), mapQ_A = col_integer(), mapCigar_A = col_character(), AlleleA_ProbeSeq = col_character(), mapNM_A = col_character(), mapAS_A = col_integer(), mapYD_A = col_character(), mapFlag_B = col_integer(), mapChrm_B = col_character(), mapPos_B = col_integer(), mapQ_B = col_integer(), mapCigar_B = col_character(), AlleleB_ProbeSeq = col_character(), mapNM_B = col_character(), mapAS_B = col_integer(), mapYD_B = col_character(), type = col_character()))
  ordering <- data.frame(Probe_ID = dfref$Probe_ID, M=dfref$address_B, U=dfref$address_A, col=factor(dfref$channel, level=c("G","R")), mask=FALSE)
  ordering$mask <- create_mask(dfref)$ref_issue
    
  species <- mclapply(seq_along(dfsp$scientificName), function(ii) {
    nm <- dfsp$scientificName[ii]
    df <- read_tsv(sprintf("%s/%s.tsv.gz", in_dir, nm), col_types=cols(CpG_chrm = col_character(), CpG_beg = col_integer(), CpG_end = col_integer(), address_A = col_integer(), address_B = col_integer(), target = col_character(), nextBase = col_character(), channel = col_character(), Probe_ID = col_character(), mapFlag_A = col_integer(), mapChrm_A = col_character(), mapPos_A = col_integer(), mapQ_A = col_integer(), mapCigar_A = col_character(), AlleleA_ProbeSeq = col_character(), mapNM_A = col_character(), mapAS_A = col_integer(), mapYD_A = col_character(), mapFlag_B = col_integer(), mapChrm_B = col_character(), mapPos_B = col_integer(), mapQ_B = col_integer(), mapCigar_B = col_character(), AlleleB_ProbeSeq = col_character(), mapNM_B = col_character(), mapAS_B = col_integer(), mapYD_B = col_character(), type = col_character()))
    df <- df[match(ordering$Probe_ID, df$Probe_ID),]
    o1 <- list(col = factor(df$channel, level=c("G","R")), mask = rep(FALSE, nrow(df)), AS = with(df, pmax(mapAS_A, mapAS_B, na.rm=TRUE)))
    o1$mask <- create_mask(df)$ref_issue
    o1$scientificName <- nm
    o1$taxonID <- dfsp$taxonID[ii]
    o1$commonName <- dfsp$commonName[ii]
    o1$assembly <- dfsp$assembly[ii]
    o1
  }, mc.cores=20)
  names(species) <- dfsp$scientificName

  list(ordering = ordering, controls = data.frame(), species = species)
}
