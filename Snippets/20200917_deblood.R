deblood <- function(x) {

    blood_invar0 <- readRDS('/mnt/isilon/zhoulab/labprojects/20200916_germlayers/blood_invariants0.rds')
    blood_invar1 <- readRDS('/mnt/isilon/zhoulab/labprojects/20200916_germlayers/blood_invariants1.rds')
    bloodMean <- readRDS('/mnt/isilon/zhoulab/labprojects/20200916_germlayers/bloodMean.rds')

    xx0 <- na.omit(x[blood_invar0])
    isTop0 <- (xx0 >= quantile(xx0, 1-20/length(xx0), na.rm=T))

    if (sum(!is.na(xx0[isTop0])) < 5) {
        purity0 <- NA
    } else {
        ds <- density(na.omit(xx0[isTop0]))
        purity0 <- ds$x[which.max(ds$y)]
    }

    xx1 <- na.omit(x[blood_invar1])
    isTop1 <- (xx1 <= quantile(xx1, 20/length(xx1), na.rm=T))

    if (sum(!is.na(xx1[isTop1])) < 5) {
        purity1 <- NA
    } else {
        ds <- density(na.omit(xx1[isTop1]))
        purity1 <- 1-ds$x[which.max(ds$y)]
    }

    nonblood_purity = max(purity0, purity1, na.rm=T)
    list(
        nonblood_purity = nonblood_purity,
        x = x - bloodMean * (1-nonblood_purity) / nonblood_purity
    )
}
