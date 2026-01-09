## Includes functions for standardising/unstandardising regression coefficients

# Please note that this works for a model.matrix set-up to ensure all covariates are correctly scaled.

## standardise.coeff.fn
# beta: set of regression coefficients
# intercept: logical for whether including intercept or not (TRUE = included/FALSE = not included). 
# (please note that intercept must be included if you have discrete covariates)
# x: model matrix of unscaled covariate values (unscaled) with columns in the same order as regression coefficients
standardise.coeff.fn <- function(beta, intercept=TRUE, x){
  if(is.logical(intercept)==FALSE){stop("intercept must be logical: TRUE (with intercept) or FALSE (without intercept) ")}
  
  n <- length(beta)
  
  if(intercept==TRUE){
    x.s <- scale(x[,2:n])
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta0.s <- beta[1] + sum(beta[2:n]*center)
    beta.s <- beta[2:n]*scaled
    beta.full.s <- c(beta0.s,beta.s)
  }
  
  if(intercept==FALSE){
    x.s <- scale(x)
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta.full.s <- beta[1:n]*scaled
  }
  
  return(beta.full.s)
}

## unstandardise.coeff.fn
# beta.s: set of regression coefficients estimated on standardised data
# intercept: logical for whether including intercept or not (TRUE = included/FALSE = not included). 
# (please note that intercept must be included if you have discrete covariates)
# x: model matrix of unscaled covariate values (unscaled) with columns in the same order as regression coefficients
unstandardise.coeff.fn <- function(beta.s, intercept = TRUE, x){
  if(is.logical(intercept)==FALSE){stop("intercept must be logical: TRUE (with intercept) or FALSE (without intercept) ")}
  
  n <- length(beta.s)
  
  if(intercept==TRUE){
    x.s <- scale(x[,2:n])
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta0.us <- beta.s[1] - sum(beta.s[2:n]*center/scaled)
    beta.us <- beta.s[2:n]/scaled
    beta.full.us <- c(beta0.us,beta.us)
  }
  
  if(intercept==FALSE){
    x.s <- scale(x)
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta.full.us <- beta.s[1:n]/scaled
  }
  
  return(beta.full.us)
}