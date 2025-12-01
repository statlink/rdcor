rdcor.test <- function(y, x, B = 499) {

  if ( !is.matrix(x) ) {
    n <- length(y)
    Ry <- Rfast::Rank(y) / n
    Rx <- Rfast::Rank(x) / n
    r <- dcov::dcor(Ry, Rx, type = "V")
    px <- replicate( B, Rfast2::Sample(x, n) )
    Rpx <- Rfast::colRanks(px) / n
    pr <- dcov::mdcor(y, Rpx, type = "V")
    res <- c( r, ( sum( pr > r ) + 1 ) / (B + 1) )
    names(res) <- c("rdcor", "permutation p-value")

  } else  {
    n <- dim(x)[1]
    Ry <- Rfast::Rank(y) / n
    Rx <- Rfast::colRanks(x) / n
    r <- as.vector( dcov::mdcor(Ry, Rx, type = "V") )
    py <- replicate( B, Rfast2::Sample(y, n) )
    Rpy <- Rfast::colRanks(py) / n
    pr <- matrix(nrow = B, ncol = dim(x)[2])
    for (i in 1:B)  pr[i, ] <- dcov::mdcor(Rpy[, i], Rx, type = "V")
    res <- cbind( r, ( Rfast::rowsums( t(pr) > r ) + 1 ) / (B + 1) )
    colnames(res) <- c("rdcor", "permutation p-value")
  }

  res
}






