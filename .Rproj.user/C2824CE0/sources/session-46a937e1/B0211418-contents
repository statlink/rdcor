rdcor <- function(y, x) {

  if ( !is.matrix(x) ) {
    n <- length(y)
    Ry <- Rfast::Rank(y) / n
    Rx <- Rfast::Rank(x) / n
    r <- dcov::dcor2d(Ry, Rx, type = "V")

  } else {
    n <- dim(x)[1]
    Ry <- Rfast::Rank(y) / n
    Rx <- Rfast::colRanks(x) / n
    r <- dcov::mdcor(Ry, Rx, type = "V")
  }

  r
}
