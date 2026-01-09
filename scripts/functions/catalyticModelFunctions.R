
mcmcModelRandomSampler <- function(numberOfSamples, mcmcMatrix, ageVector){
  for (i in 1:numberOfSamples){
    # Specify random number to sample from mcmc chain for parameter estimates
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    lambdaSample <- mcmcMatrix[randomNumber,2]
    deltaSample <- mcmcMatrix[randomNumber,1]
    # Run the model with new esitmates
    newRow <- (lambdaSample / (lambdaSample+deltaSample)) *(1 - exp(-ageVector*(lambdaSample+deltaSample)))
    # Add this to a new matrix, where each row is a sample, and each column is an age
    outDf[i,] <- newRow
  }
  outDf
}

mcmcRandomSamplerCat <- function(numberOfSamples, mcmcMatrix, ageVector, ageTotals){
  
  outDf <- matrix(NA,nrow=numberOfSamples, ncol = length(ageVector))
  
  
  for (i in 1:numberOfSamples){
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    
    lambdaSample <- mcmcMatrix[randomNumber,1]
    
    newRow <- 1 - exp(-ageVector*(lambdaSample))
    updateRow <- c()
    for(j in 1:length(ageTotals)){
      randomlySampleBinomialDis <- rbinom(1,size = ageTotals[j],prob = newRow[j])
      if(randomlySampleBinomialDis > 0){
        result <- randomlySampleBinomialDis / ageTotals[j]
      } else {
        result <- 0
      }
      updateRow[j] <- result
    }
    rbinom(1,20,0.03)
    outDf[i,] <- updateRow
  }
  outDf
}


mcmcRandomSampler <- function(numberOfSamples, mcmcMatrix, ageVector, ageTotals){
  
  outDf <- matrix(NA,nrow=numberOfSamples, ncol = length(ageVector))
  
  
  for (i in 1:numberOfSamples){
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    
    lambdaSample <- mcmcMatrix[randomNumber,2]
    deltaSample <- mcmcMatrix[randomNumber,1]
    
    newRow <- (lambdaSample / (lambdaSample+deltaSample)) *(1 - exp(-ageVector*(lambdaSample+deltaSample)))
    updateRow <- c()
    for(j in 1:length(ageTotals)){
      randomlySampleBinomialDis <- rbinom(1,size = ageTotals[j],prob = newRow[j])
      if(randomlySampleBinomialDis > 0){
        result <- randomlySampleBinomialDis / ageTotals[j]
      } else {
        result <- 0
      }
      updateRow[j] <- result
    }
    
    outDf[i,] <- updateRow
  }
  outDf
}

ageQuantiles <- function(mcmcDF){
  quantileMatrix <- matrix(NA,nrow=ncol(mcmcDF), ncol = 3)
  for(i in 1:ncol(mcmcDF)){
    quantiles <- mcmcDF[,i] %>% quantile(probs=c(.5,.025,.975))
    quantileMatrix[i,] <- quantiles
  }
  quantileMatrix
}

extend_sampling_to_bounds <- function(df_sampling,
                                      lower_age = 4,
                                      upper_age = 18) {
  stopifnot(all(c("midpoint", "prev", "ci_lower", "ci_upper") %in% names(df_sampling)))
  df_sampling <- df_sampling[order(df_sampling$midpoint), ]
  
  # If already includes bounds, return as-is
  if (lower_age %in% df_sampling$midpoint &&
      upper_age %in% df_sampling$midpoint) {
    return(df_sampling)
  }
  
  full_midpoint_values <- sort(unique(c(lower_age, df_sampling$midpoint, upper_age)))
  
  linear_extrapolate <- function(x, y, new_x) {
    fit <- lm(y ~ x)
    as.numeric(predict(fit, newdata = data.frame(x = new_x)))
  }
  
  interp_prev <- approx(df_sampling$midpoint, df_sampling$prev,
                        xout = full_midpoint_values, rule = 2)$y
  interp_lower <- approx(df_sampling$midpoint, df_sampling$ci_lower,
                         xout = full_midpoint_values, rule = 2)$y
  interp_upper <- approx(df_sampling$midpoint, df_sampling$ci_upper,
                         xout = full_midpoint_values, rule = 2)$y
  
  # Left boundary extrapolation
  if (lower_age < min(df_sampling$midpoint)) {
    interp_prev[1]  <- linear_extrapolate(df_sampling$midpoint[1:2],
                                          df_sampling$prev[1:2], lower_age)
    interp_lower[1] <- linear_extrapolate(df_sampling$midpoint[1:2],
                                          df_sampling$ci_lower[1:2], lower_age)
    interp_upper[1] <- linear_extrapolate(df_sampling$midpoint[1:2],
                                          df_sampling$ci_upper[1:2], lower_age)
  }
  
  # Right boundary extrapolation
  if (upper_age > max(df_sampling$midpoint)) {
    n <- nrow(df_sampling)
    interp_prev[length(interp_prev)]  <-
      linear_extrapolate(df_sampling$midpoint[(n-1):n],
                         df_sampling$prev[(n-1):n], upper_age)
    interp_lower[length(interp_lower)] <-
      linear_extrapolate(df_sampling$midpoint[(n-1):n],
                         df_sampling$ci_lower[(n-1):n], upper_age)
    interp_upper[length(interp_upper)] <-
      linear_extrapolate(df_sampling$midpoint[(n-1):n],
                         df_sampling$ci_upper[(n-1):n], upper_age)
  }
  
  data.frame(
    midpoint = full_midpoint_values,
    prev     = interp_prev,
    ci_lower = interp_lower,
    ci_upper = interp_upper
  )
}


