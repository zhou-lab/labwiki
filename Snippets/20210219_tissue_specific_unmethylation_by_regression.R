coefs <- mclapply(seq_len(ncol(betasT)), function(i) {
    a <- with(meta3, lm(betasT[,i] ~ Stage + Tissue + Sex + Strain))
    a$coefficients
}, mc.cores=20)

length(unique(lapply(coefs, names))) # names are identical expect 1
coefs <- do.call(cbind, coefs)
colnames(coefs) <- colnames(betasT)

ticfs <- grep('Tissue',rownames(coefs),value=T)
tissue_coefs <- apply(coefs, 2, function(x) {xx <- x[ticfs]; delta <- -sum(xx)/(length(xx)+1); xx <- xx+delta; xx['TissueBrain'] <- delta; xx})

saveRDS(tissue_coefs, file='/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/tissue_coefs.rds')
saveRDS(coefs, file='/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/allcoefs.rds')
## scp hpc2:/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/tissue_coefs.rds ~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/tissue_coefs.rds
## scp hpc2:/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/allcoefs.rds ~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/allcoefs.rds

tissue_coefs <- readRDS('~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/tissue_coefs.rds')
rownames(tissue_coefs) <- str_replace(rownames(tissue_coefs), 'Tissue','')
tissue_coefs <- tissue_coefs[,intersect(getAutosomeProbes('MM285','mm10'),colnames(tissue_coefs))]
coefs <- readRDS('~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210217_integrativeAnalysis/allcoefs.rds')

head(tissue_coefs[,1:4])
dim(tissue_coefs)

pdf('~/gallery/20210217_integrativeAnalysis_tissue_density.pdf', width=10, height=5, onefile=FALSE)
wzPlotDens1d.fromMatrixGG(t(tissue_coefs)) + xlab('Delta Methylation Level') + scale_color_discrete(name='Tissue')
dev.off()

most_u <- apply(tissue_coefs, 2, function(x) names(sort(x))[1])
disc_power <- apply(tissue_coefs, 2, function(x) {xx <- sort(x); xx[2]-xx[1];})
## ti_u <- lapply(split(disc_power, most_u), function(x) names(head(sort(x[x>0.2], decreasing=TRUE),100)))
ti_u <- lapply(split(disc_power, most_u), function(x) names(head(sort(x[x>0.2], decreasing=TRUE), 200)))
