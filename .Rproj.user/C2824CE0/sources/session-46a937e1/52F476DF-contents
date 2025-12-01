mpdcor <- function(y, x, z) {
  a1 <- dcov::mdcor(y, x, type = "U")
  a2 <- dcov::mdcor(z, x, type = "U")
  a3 <- dcov::dcor(y, z, type = "U")
  r <-  up <- a1 - a2 * a3
  down <- sqrt(1 - a2^2) * sqrt(1 - a3^2)
  as.vector( up/down )
}
