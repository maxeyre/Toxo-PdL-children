# loc.primary are the locations for which you want to perform the minimum distance calculation for every loc.secondary locations.

dist.3D.fn <- function(loc.primary, loc.secondary){
  dist.fn <- function(loc1, loc2){
    dist <- sqrt((loc1$X-loc2$X)^2+(loc1$Y-loc2$Y)^2+
                   (loc1$Z-loc2$Z)^2)
  }
  dist.in.fn2 <- function(i){
    d <- min(dist.fn(loc.primary[i,], loc.secondary))
  }
  n <- nrow(loc.primary)
  d.out <- sapply(1:n, dist.in.fn2)
  return(d.out)
}