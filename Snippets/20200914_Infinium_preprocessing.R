## version 1
library(sesame)
source('https://raw.githubusercontent.com/zhou-lab/tbmate/master/scripts/tbmate.R')

## change your input folder here
base_dir <- '/mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE122126/'
pfxs     <- searchIDATprefixes(file.path(base_dir, 'IDATs'))
idx_dir  <- '/mnt/isilon/zhou_lab/projects/20191221_references/InfiniumArray'

tmp <- mclapply(seq_along(pfxs), function(i) {
    sset <- readIDATpair(pfxs[i])
    cat(pfxs[i],sset@platform,'\n')
    betas <- getBetas(dyeBiasCorrTypeINorm(noob(sset)))
    pvals <- pval(sset)[names(betas)]
    tbk_pack(betas, data2 = pvals, 
             out_dir = sprintf('%s/tbk_%s/', base_dir, sset@platform), 
             out_fname = str_split(names(pfxs)[i],'_')[[1]][1], 
             idx_fname = sprintf('%s/%s/%s.idx.gz', idx_dir, sset@platform, sset@platform), 
             dtype='FLOAT_FLOAT')
} , mc.cores=15)
