library(sesame)
pfx = searchIDATprefixes('IDATs')

sset = readIDATpair(pfx[1])
head(IG(sset))
head(oobR(sset))
head(IR(sset))
head(oobG(sset))


ssets = lapply(pfx[1:3], readIDATpair)

library(parallel)
ssets = mclapply(pfx[1:3], readIDATpair, mc.cores=4)

library(tidyverse)
betas = do.call(cbind, mclapply(pfx[1:3], function(x) {
  readIDATpair(x) %>% noob %>% dyeBiasCorrTypeINorm %>% detectionMask %>% getBetas
} , mc.cores=4))

?noob
?detectionMask
noob

# how pvals are calculated
pOOBAH

# beta == M / (M+U)
deIdentify("IDATs/204617710004_R01C01_Red.idat", out_path='test_Red.idat')
betas2 = getBetas(readIDATpair('test'))
head(betas2[grep('^rs', names(betas2))])


