# Code to sample land cover 
# Raster: lc.raster with three types: 1 = veg, 2 = soil, 3 = imperv.
# Locations: a set of unique locations unique.locs with X and Y columns for coordinates
# Raster and locations must be in the same UTM coordinate system
# Buffer: a vector of the different buffers sizes to be extracted (radius in m)

lc.sampler.fn <- function(unique.locs, lc.raster, buffer, name=FALSE){
  library(raster)
  
  ifelse(name==FALSE, name <- "lc_", name <- name)
  
  coords <- unique.locs[,c("X","Y")]
  n <- nrow(coords)
  bn <- length(buffer)
  
  output <- matrix(NA, nrow=n, ncol=0)
  
  for(i in 1:bn){
    y <- raster::extract(lc.raster, coords, buffer=buffer[i])
    count <- matrix(NA, nrow=n, ncol = 4)
    count[,1] <- sapply(y, function(i) sum(i==1)) # vegetation
    count[,2] <- sapply(y, function(i) sum(i==2)) # soil
    count[,3] <- sapply(y, function(i) sum(i==3)) # impervious
    count[,4] <- rowSums(count[,1:3])
    
    out <- matrix(NA, nrow=n, ncol=3)
    out[,1] <- count[,1]/count[,4] # vegetation
    out[,2] <- count[,2]/count[,4] # soil
    out[,3] <- count[,3]/count[,4] # impervious
    out <- as.data.frame(out)
    
    lc.types <- c("veg", "soil", "imperv")
    names.col <- paste0(name, buffer[i], "m_", lc.types)
    colnames(out) <- names.col
    
    output <- cbind(output, out)
  }
  
  output <- cbind(coords, output)
  return(output)
}

