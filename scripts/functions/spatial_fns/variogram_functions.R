ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords))/3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "distance", 
                    show_nbins = F, envelop=F, theo=F, my.l=1) {
  require(geoR)
  
  coords <- as.matrix(coords)
  min_dist <- min(dist(coords))
  if(dim(coords)[1] < 2){
    # if(is.null(uvec))  uvec <- seq(min_dist, maxdist, l = bins)
    
    fit.dat <- data.frame(xx= seq(0, maxdist, length.out = 101), 
                          yy = theo.var.fit(u= seq(0, maxdist, length.out = 101), my.l = my.l))
    maxofy <- max(fit.dat$yy)
    pp <- ggplot() + geom_line(data = fit.dat, aes(x = xx, y = yy), inherit.aes = FALSE) +
      scale_y_continuous(name = "semivariance", limits = c(0,  maxofy))
    return(pp)
  } 
  else{
    if(min_dist < maxdist){
      if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
    }else{
      if(is.null(uvec))  uvec <- seq(maxdist, min_dist, l = bins)
    }
  }
  
  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
  # plot(empvario)
  if(envelop==FALSE){
    dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                          nbins = empvario$n)
    p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
      geom_point(col = "black", fill = color, shape = 21, size = 3) +
      scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                         breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
      # scale_y_continuous(name = "semivariance",
      #                    limits = c(0,  max(dfvario$empirical))) +
      ggtitle("Empirical semivariogram") 
    # theme_classic()
    p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
    
    if(theo == T){
      fit.dat <- data.frame(xx= seq(0, maxdist, length.out = 101), 
                            yy = theo.var.fit(u= seq(0, maxdist, length.out = 101), my.l = my.l))
      maxofy <- max(max(fit.dat$yy), max(dfvario$empirical))
      if(show_nbins){
        p3 <- p2 + geom_line(data = fit.dat , aes(x = xx, y = yy), inherit.aes = FALSE) +
          scale_y_continuous(name = "semivariance", limits = c(0,  maxofy)) 
        p3
      } else{
        p3 <- p1 + geom_line(data = fit.dat, aes(x = xx, y = yy), inherit.aes = FALSE) +
          scale_y_continuous(name = "semivariance", limits = c(0,  maxofy))
        p3
      } 
    }else{
      if(show_nbins){
        p3 <- p2 +  scale_y_continuous(name = "semivariance", limits = c(0,  max(dfvario$empirical))) 
        p3
      } else{
        p3 <- p1 +  scale_y_continuous(name = "semivariance", limits = c(0,  max(dfvario$empirical))) 
        p3
      }
    }
  }else{
    envmc <- variog_envelope(coords = coords, data = data, 
                             obj.variog = empvario, nsim = nsim, messages = F)
    dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                          lowemp = envmc$v.lower, upemp = envmc$v.upper, 
                          nbins = empvario$n)
    p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
      geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = .3) +
      geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
      scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                         breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
      scale_y_continuous(name = "semivariance", 
                         #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1), 
                         limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
      ggtitle("Empirical semivariogram") 
    # theme_classic()
    p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
    if(show_nbins) p2 else p1
  }
}


variog_envelope <- function (geodata, coords = geodata$coords, data = geodata$data, 
                             obj.variog, nsim = 99, save.sim = FALSE, messages) 
{
  call.fc <- match.call()
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  obj.variog$v <- NULL
  if ((is.matrix(data) | is.data.frame(data))) 
    if (ncol(data) > 1) 
      stop("envelops can be computed for only one data set at once")
  if (!is.null(obj.variog$estimator.type)) 
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    if (abs(obj.variog$lambda) < 1e-04) 
      data <- log(data)
    else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
  }
  xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ xmat + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim, "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data) {
    return(data[sample(1:n.data)])
  }
  simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, 
                       n.data = n.data)
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      temp <- .C("binit", as.integer(obj.variog$n.data), 
                 as.double(as.vector(coords[, 1])), as.double(as.vector(coords[, 
                                                                               2])), as.double(as.vector(sim)), as.integer(nbins), 
                 as.double(as.vector(obj.variog$bins.lim)), as.integer(estimator.type == 
                                                                         "modulus"), as.double(max(obj.variog$u)), as.double(cbin), 
                 vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin), 
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  }
  else {
    variog.vbin <- function(x, ...) {
      variog(geodata = geodata, 
             data = x, uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type, 
             nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist, 
             pairs.min = obj.variog$pairs.min, direction = obj.variog$direction, 
             tolerance = obj.variog$tolerance, messages.screen = FALSE,...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (save.sim == FALSE) 
    simula$data <- NULL
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], 
                  v.upper = limits[2, ])
  if (save.sim) 
    res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}