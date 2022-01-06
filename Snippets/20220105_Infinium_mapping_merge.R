
mergeI_and_II_1 <- function(x) {
  dir.create("mapping", showWarnings=FALSE)
  
  p2addr = read_tsv("tmp/probe2address.txt", col_names = c("Probe_ID","addrA","addrB"))
  p2addrA = with(p2addr, setNames(addrA, Probe_ID))
  p2addrB = with(p2addr, setNames(addrB, Probe_ID))
  p2seq = read_tsv("tmp/probe2originalseq.txt", col_names = c("Probe_ID","seqA","seqB"))
  p2seqA = with(p2seq, setNames(seqA, Probe_ID))
  p2seqB = with(p2seq, setNames(seqB, Probe_ID))
  
  df1 <- read_tsv(sprintf('tmp/%s_AB_final', x), col_names=c(
    'chrmA','begA','endA','lastA','probeID','flagA','samChrmA','samPosA','mapqA','cigarA','samSeqA','nmA','asA','ydA',
    'chrmB','begB','endB','lastB','probeID.B','flagB','samChrmB','samPosB','mapqB','cigarB','samSeqB','nmB','asB','ydB',
    'ext','tgt','col'), guess_max=300000);
  df1$type='I'; df1$species=x;
  df2 <- read_tsv(sprintf('tmp/%s_II_final', x), col_names=c(
    'chrm','beg','end','last','probeID','flag','samChrm','samPos','mapq','cigar','samSeq','nm','as','yd','tgt'), guess_max=300000);
  df22 <- tibble(chrmA=df2$chrm, begA=df2$beg, endA=df2$end, 
                 lastA=df2$last, probeID=df2$probeID, flagA=df2$flag, samChrmA=df2$samChrm, samPosA=df2$samPos, 
                 mapqA=df2$mapq, cigarA=df2$cigar, samSeqA=df2$samSeq, nmA=df2$nm, asA=df2$as, ydA=df2$yd, 
                 chrmB=NA, begB=NA, endB=NA, lastB=NA, probeID.B=NA, flagB=NA, samChrmB=NA, samPosB=NA, mapqB=NA, 
                 cigarB=NA, samSeqB=NA, nmB=NA, asB=NA, ydB=NA, ext=NA, tgt=df2$tgt, col=NA, type='II', species=x)
  df3 <- rbind(df1,df22)
  
  df4 <- with(df3, tibble(
    CpG_chrm = chrmA, CpG_beg = begA, CpG_end = endA, 
    address_A = p2addrA[probeID],
    address_B = p2addrB[probeID],
    target = tgt,
    nextBase = ext,
    channel = col,
    Probe_ID = probeID,
    lastBase_A = lastA, mapFlag_A = flagA, mapChrm_A = samChrmA, mapPos_A = samPosA, mapQ_A = mapqA, 
    mapCigar_A = cigarA, AlleleA_ProbeSeq = p2seqA[probeID],
    mapNM_A = nmA, mapAS_A = asA, mapYD_A = ydA, 
    lastBase_B = lastB, mapFlag_B = flagB, mapChrm_B = samChrmB, mapPos_B = samPosB, mapQ_B = mapqB, 
    mapCigar_B = cigarB, AlleleB_ProbeSeq = p2seqB[probeID],
    mapNM_B = nmB, mapAS_B = asB, mapYD_B = ydB,
    design = NA,
    mask = NA, type = type))
  saveRDS(df4, file=sprintf('mapping/%s.rds', x))
}

mergeI_and_II <- function(dir) {
  setwd(dir)
    
  dir.create("mapping", showWarnings=FALSE)
  df <- read_excel('~/samplesheets/2020/20201202_310_species_EnsemblVertebrates.xlsx')
  mfts <- mclapply(df$species, mergeI_and_II_1, mc.cores=20)

  ## count number of probes
  ns <- simplify2array(mclapply(df$species, function(x) nrow(readRDS(sprintf('mapping/%s.rds', x))), mc.cores=20))
  names(ns) <- df$species
  write_tsv(data.frame(ns, species=names(ns)), file=sprintf('%s/num_probes.tsv', dir))
  ns
}
