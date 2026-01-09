# Code to sample Topographic Wetness Index (TWI)
# Raster: twi.raster with continuous values
# Locations: a set of unique locations unique.locs with X and Y columns for coordinates
# Raster and locations must be in the same UTM coordinate system
# Buffer: a vector of the different buffers sizes to be extracted (radius in m)

twi.sampler.fn <- function(unique.locs, twi.raster, buffer, name=FALSE, mean=TRUE, median=TRUE, range_low=FALSE, range_upp=FALSE){
  library(raster)
  
  ifelse(name==FALSE, name <- "twi_", name <- name)
  
  coords <- unique.locs[,c("X","Y")]
  n <- nrow(coords)
  bn <- length(buffer)
  
  vars <- c(mean, median, range_low, range_upp)
  
  output <- matrix(NA, nrow=n, ncol=0)
  
  for(i in 1:bn){
    y <- raster::extract(twi.raster, coords, buffer=buffer[i])
    out <- matrix(NA, nrow=n, ncol = 4)
    out[,1] <- sapply(y, function(x) mean(x)) # mean
    out[,2] <- sapply(y, function(x) median(x)) # median
    out[,3] <- sapply(y, function(x) range(x))[1,] # lower range
    out[,4] <- sapply(y, function(x) range(x))[2,] # upper range
    
    out <- as.data.frame(out)
    
    twi.cols <- c("mean", "median", "range_low", "range_upp")
    names.col <- paste0(name, buffer[i], "m_", twi.cols)
    colnames(out) <- names.col
    out <- out[,c(names.col[vars])]
    
    output <- cbind(output, out)
  }
  
  output <- cbind(coords, output)
  return(output)
}

