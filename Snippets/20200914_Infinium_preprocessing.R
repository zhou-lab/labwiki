###################################
## version 1 with name cleaning
###################################
library(sesame)
source('https://raw.githubusercontent.com/zhou-lab/tbmate/master/scripts/tbmate.R')

## change your input folder here
##  ./IDATs
##  ./tbk_EPIC   # output
##  ./tbk_HM450  # output
base_dir <- '/mnt/isilon/zhoulab/labprojects/20200822_TCGA_HM450/'
pfxs     <- searchIDATprefixes(file.path(base_dir, '/IDATs'))
idx_dir  <- '/mnt/isilon/zhou_lab/projects/20191221_references/InfiniumArray'

tmp <- mclapply(seq_along(pfxs), function(i) {
    sset <- readIDATpair(pfxs[i])
    cat(pfxs[i],sset@platform,'\n')
    betas <- getBetas(dyeBiasCorrTypeINorm(noob(sset)))
    pvals <- pval(sset)[names(betas)]
    tbk_pack(betas, data2 = pvals, 
             out_dir = sprintf('%s/tbk_%s/', base_dir, sset@platform), 
             out_fname = names(pfxs)[i], # str_split(names(pfxs)[i],'_')[[1]][1], 
             idx_fname = sprintf('%s/%s/%s.idx.gz', idx_dir, sset@platform, sset@platform), 
             dtype='FLOAT_FLOAT')
} , mc.cores=20)

###################################
## version 2 with name cleaning
###################################
library(sesame)
source('https://raw.githubusercontent.com/zhou-lab/tbmate/master/scripts/tbmate.R')

base_dir <- '/mnt/isilon/zhoulab/labprojects/20200822_TCGA_HM450/'
pfxs     <- searchIDATprefixes(file.path(base_dir, '/IDATs'))
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
} , mc.cores=20)


###################################
## version 3 beta values only
###################################
library(parallel)
source('https://raw.githubusercontent.com/zhou-lab/tbmate/master/scripts/tbmate.R')
base_dir <- '/mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE41826'
idx_fname  <- '/mnt/isilon/zhou_lab/projects/20191221_references/InfiniumArray/HM450/HM450.idx.gz'

load('legacy/betas.rda')

tmp <- mclapply(seq_len(ncol(betas)), function(i) {
    tbk_pack(betas[,i],
             out_dir = sprintf('%s/tbk_HM450/', base_dir), 
             out_fname = colnames(betas)[i], 
             idx_fname = idx_fname, 
             dtype='FLOAT')
}, mc.cores=20)




