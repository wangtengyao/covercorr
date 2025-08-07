kmac.test <- function(x, y, B){
  y <- as.matrix(y)
  stats <- rep(0, B+1)
  stats[1] <- KMAc(y, x)
  for (b in 1:B) stats[b+1] <- KMAc(y[sample(nrow(y)), ], x)
  return(list(stat=stats[1], pval=mean(stats[1] <= stats)))
}

dcor.test <- function(x, y, B){
  y <- as.matrix(y)
  stats <- rep(0, B+1)
  stats[1] <- Rfast::dcor(x, y)$dcor
  for (b in 1:B) stats[b+1] <- Rfast::dcor(x, y[sample(nrow(y)), ])$dcor
  return(list(stat=stats[1], pval=mean(stats[1] <= stats)))
}
