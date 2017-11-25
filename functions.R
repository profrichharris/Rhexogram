sppoly2hex <- function(sppoly, bins = 45) {
  
  pts <- coordinates(sppoly)
  hex <- hexBinning(pts, bins = bins)
  return(hex)
  
}

hexconflicts <- function(sppoly, bins = 45) {
  
  pts <- coordinates(sppoly)
  hex <- sppoly2hex(sppoly, bins)
  kn <- nn2(query = pts, data = cbind(hex$x, hex$y), k=1)[[1]]
  subset <- as.logical(duplicated(kn) + duplicated(kn, fromLast = TRUE))
  N <- length(subset[subset])
  return(list(N, bins, subset, hex, kn))
  
}

binmin <- function(sppoly, N) {
  
  bins <- 2
  target <- N
  while(N > target) {
    bins <- bins + 1
    hex <- hexconflicts(sppoly, bins)
    N <- hex[[1]]
  }
  return(hex)
  
}