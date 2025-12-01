rdcor.mat <- function(x, B = 1) {

  dm <- dim(x)
  n <- dm[1]  ;  p <- dm[2]
  mat <- matrix(1, p, p)
  Rx <- Rfast::colRanks(x) / n
  for ( i in 1:c(p-1) ) {
    mat[i, -c(1:i)] <- mat[-c(1:i), i] <- dcov::mdcor(Rx[, i], Rx[, -c(1:i)], type = "V")
  }

  pval <- NULL
  if ( B > 1 ) {
    pval <- mat
    for ( i in 1:c(p-1) ) {
      pr <- matrix( nrow = B, ncol = p - i)
      py <- replicate( B, Rfast2::Sample(x[, i], n) )
      Rpy <- Rfast::colRanks(py) / n
      for ( j in 1:B )  pr[j, ] <- dcov::mdcor(Rpy[, j], Rx[, -c(1:i)], type = "V")
      pval[i, -c(1:i)] <- pval[-c(1:i), i] <- ( Rfast::rowsums( t(pr) > mat[i, -c(1:i)] ) + 1 ) / (B + 1)
    }
  }

  list(r = mat, pvalue = pval)
}


