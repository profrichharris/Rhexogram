

sppoly2hex <- function(sppoly, bins) {
  
  pts <- sp::coordinates(sppoly)
  hex <- fMultivar::hexBinning(pts, bins = bins)
  return(hex)
  
}

hexconflicts <- function(sppoly, bins) {
  
  pts <- sp::coordinates(sppoly)
  hex <- fMultivar::hexBinning(pts, bins = bins)
  kn <- RANN::nn2(query = pts, data = cbind(hex$x, hex$y), k=1)[[1]]
  subset <- as.logical(duplicated(kn) + duplicated(kn, fromLast = TRUE))
  N <- length(subset[subset])
  return(list(N, bins, subset, hex, kn))
  
}

binmin <- function(sppoly, N) {
  
  bins <- 1
  N2 <- N + 1
  while(N2 > N) {
    bins <- bins + 1
    hex <- hexconflicts(sppoly, bins)
    N2 <- hex[[1]]
  }
  return(hex)
  
}

shift <- function(sppoly, hex, trials) {
  
  N <- hex[[1]]
  bins <- hex[[2]]
  subset <- hex[[3]]
  
  pts <- sp::coordinates(sppoly)
  
  hex.save <- hex
  sppoly.save <- sppoly
  
  for(i in 1:trials) {
  
    sppoly <- sppoly.save  
    submap <- sppoly[subset,]
    n <- length(submap)
    rpnts <- sapply(1:n, function(x) {
      sp::coordinates(sp::spsample(submap[x,], 1, type="random", iter=1000))
    })
    rpnts <- t(rpnts)
    pts[subset,] <- rpnts
    for(i in 1:length(sppoly)) {
      sppoly@polygons[[i]]@labpt <- pts[i,]
    }
    hex <- hexconflicts(sppoly, bins)
    N2 <- hex[[1]]
    if(N2 < N) {
      hex.save <- hex
      subset <- hex[[3]]
      sppoly.save <- sppoly
    }
    if(N2 == 0) break
    
  }
  return(list(hex.save, sppoly.save))
}



hex2sphex <- function(hex) {
  
  x <- hex
  
  X = x$x
  Y = x$y
  rx = min(diff(unique(sort(X))))
  ry = min(diff(unique(sort(Y))))
  rt = 2 * ry
  u = c(rx, 0, -rx, -rx, 0, rx)
  v = c(ry, rt, ry, -ry, -rt, -ry)/3
  
  polys <- lapply(1:length(X), function(i) {
    coords <- cbind(u + X[i], v + Y[i])
    poly <- sp::Polygon(coords)
    sp::Polygons(list(poly), ID=i)
  })
  sphex <- sp::SpatialPolygons(polys)
  return(sphex)
  
}


hexogram <- function(sppoly, bins = 45, maxSizeError = 1.25,
                     binsmax = bins + 10, trials = 500,
                     set.seed = 101) {
  
  if(!is.null(set.seed)) set.seed(set.seed)
  carto <- sppoly

  areas <- rgeos::gArea(carto, byid = TRUE)
  carto$scaleby <- areas
  
  bbx <- sp::bbox(carto)
  A <- (bbx[1,2] - bbx[1,1]) * (bbx[2,2] - bbx[2,1])
  su <- 2 * A / bins^2
  
  subset <- carto$scaleby <= su
  carto$scaleby[subset] <- su

  cat("\nStep 1, creating initial cartogram")
  carto <- cartogram::cartogram(carto, "scaleby", itermax=10,
                     maxSizeError = maxSizeError, prepare="none")
  
  # Find the number of conflicts
  hex <- hexconflicts(carto, bins)
  N <- hex[[1]]
  
  cat("\nStep 2, resolving remaining conflicts")
  
  while(N > 0) {
    
    carto.bak <- carto
    N.bak <- N
    
    bins <- hex[[2]]
    if(bins > binsmax) stop("Number of bins exceeds maximum")
    
    cat("\nNumber of conflicts = ", N)
    cat("\nNumber of bins =", bins)

    # See if the same number can be achieved with fewer bins
    hex <- binmin(carto, N)
    N <- hex[[1]]
    
    # Try moving the centroids
    if(N > 0) {
      newhex <- shift(carto, hex, trials)
      hex <- newhex[[1]]
      carto <- newhex[[2]]
      N <- hex[[1]]
    }
    
    # Try enlarging the conflicted areas
    if(N > 0) {
      subset <- hex[[3]]
      su <- min(areas[!subset & areas > max(areas[subset])])
      carto$scaleby[subset] <- su
    
      carto <- cartogram::cartogram(carto, "scaleby", itermax=10,
                       maxSizeError = maxSizeError, prepare="none")
    
      # Find the number of conflicts
      hex <- hexconflicts(carto, bins)
      N <- hex[[1]]
    }
    
    # If all else fails
    if(N > N.bak) {
      carto <- carto.bak
      N <- N.bak
      hex <- binmin(carto, N-1)
      N <- hex[[1]]
    }
    
  }
  
  hbins <- hex[[4]]
  kn <- hex[[5]][,1]

  sphex <- hex2sphex(hbins)
  sphex <- sphex[kn,]
  return(list(carto, sphex))
  
}


binN <- function(sppoly, bin.max = 65) {
  
  x <- 1:bin.max
  y <- sapply(x, function(bins) {
    hexconflicts(sppoly, bins)[[1]]
  })
  N <- length(sppoly)
  py <- y/N
  plot(x, py, xlab = "bins", ylab = "% areas conflicted", type="l", las=1)
  lo <- loess(y ~ x)
  y2 <- predict(lo)
  dy <- diff(y2) 
  dy <- (dy - min(dy)) / (max(dy) - min(dy)) * (max(py) - min(py)) + min(py)
  points(1:length(dy), dy, type="l", lty="dashed")
  rug(0: bin.max)
  rug(seq(0, bin.max, by = 5), lwd=2)
  
}

reallocate <- function(map, hexoutput) {
  
  pts.hex <- sp::coordinates(hexoutput[[2]])
  pts.map <- sp::coordinates(map)
  
  X = pts.hex[,1]
  Y = pts.hex[,2]
  rx = diff(unique(sort(X)))
  rx <- min(rx[rx > 0.001])
  ry = diff(unique(sort(Y)))
  ry <- min(ry[ry > 0.001])  
  
  D <- sqrt(rx^2 + ry^2)
  
  xx <- pts.map[,1] - pts.hex[,1]
  yy <- pts.map[,2] - pts.hex[,2]
  xsign <- sign(xx)
  ysign <- sign(yy)
  xx <- abs(xx)
  yy <- abs(yy)
  dd <- sqrt(xx^2 + yy^2)
  
  n <- nrow(pts.map)
  saved <- rep(FALSE, n)
  pts.bak <- pts.hex
  pts.save <- pts.hex
  
  for(i in seq(1,0.1,by=-0.05)) {
    
    theta <- atan(xx/yy)
    xx <- sin(theta) * i *dd
    yy <- cos(theta) * i *dd
    
    pts.try <- pts.bak
    pts.try[,1] <- pts.try[,1] + xsign * xx
    pts.try[,2] <- pts.try[,2] + ysign * yy
    
    subset <- !complete.cases(pts.try)
    pts.try[subset,] <- pts.bak[subset,]
    
    for(i in 1:n) {
      pts.hex <- pts.save
      pts.hex[i,] <- pts.try[i,] 
      qry <- cbind(pts.hex[i,1], pts.hex[i,2])
      if(RANN::nn2(pts.hex, qry, k=2)$nn.dists[,2] > D) {
        if(!saved[i]) {
          saved[i] <- TRUE
          pts.save[i,] <- pts.try[i,]
        }
      }
    }
    
  }
  
  
  X <- pts.save[,1]
  Y <- pts.save[,2]
  rt = 2 * ry
  u = c(rx, 0, -rx, -rx, 0, rx)
  v = c(ry, rt, ry, -ry, -rt, -ry)/3
  
  polys <- lapply(1:length(X), function(i) {
    coords <- cbind(u + X[i], v + Y[i])
    poly <- sp::Polygon(coords)
    sp::Polygons(list(poly), ID=i)
  })
  sppolys <- sp::SpatialPolygons(polys)
  return(sppolys)
  
}