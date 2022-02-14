
## source: DescTools https://github.com/AndriSignorell/DescTools/blob/f46b633486104c963137f574cbb826f0bf5e182d/R/StatsAndCIs.r
UncertCoef <- function(x, y = NULL, direction = c("symmetric", "row", "column"),
                       p.zero.correction = 1/sum(x)^2, ... ) {
  # Theil's UC (1970)
  # slightly nudge zero values so that their logarithm can be calculated (cf. Theil 1970: x->0 => xlogx->0)
  if(!is.null(y)) x <- table(x, y, ...)

  x[x == 0] <- p.zero.correction

  n <- sum(x)
  rsum <- rowSums(x)
  csum <- colSums(x)

  hx <- -sum((apply(x, 1, sum) * log(apply(x, 1, sum)/n))/n)
  hy <- -sum((apply(x, 2, sum) * log(apply(x, 2, sum)/n))/n)
  hxy <- -sum(apply(x, c(1, 2), sum) * log(apply(x, c(1, 2), sum)/n)/n)

  switch( match.arg( arg = direction, choices = c("symmetric", "row", "column") )
          , "symmetric" = { res <- 2 * (hx + hy - hxy)/(hx + hy) }
          , "row" = { res <- (hx + hy - hxy)/hx }
          , "column" = { res <- (hx + hy - hxy)/hy }
  )
  return(res)
}


Entropy <- function(x, y = NULL, base = 2, ...) {

  # x is either a table or a vector if y is defined

  if(!is.null(y)) { x <- table(x, y, ...) }
  x <- as.matrix(x)

  ptab <- x / sum(x)
  H <- - sum( ifelse(ptab > 0, ptab * log(ptab, base=base), 0) )
  return(H)

}


MutInf <- function(x, y = NULL, base = 2, ...){
  # ### Ref.:  http://en.wikipedia.org/wiki/Cluster_labeling

  if(!is.null(y)) { x <- table(x, y, ...) }
  x <- as.matrix(x)

  return(
    Entropy(rowSums(x), base=base) +
      Entropy(colSums(x), base=base) - Entropy(x, base=base)
  )

}
