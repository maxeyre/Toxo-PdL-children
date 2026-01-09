################################################################################
# TOXO-PdL-children
################################################################################

# -------------------------------------------------------------------------
# Script: 04_spatial_analysis.R
#
# Purpose:
#  - Variable selection for prediction model (GLMM dredge)
#  - Build prediction grid
#  - Fit intercept-only binomial geostatistical model (MCML)
#  - Predict prevalence + exceedance
#  - Fit full geostatistical model + predict S(x)
#  - Make maps and panel figure
# 
# Note: coordinates are not provided in de-identified data and 
# fitting geostatistical model without them is not be possible
# -------------------------------------------------------------------------

## ------------------------------------------------------
## 0. Packages & global options ####
## ------------------------------------------------------

required_pkgs <- c(
  "tidyverse",
  "lme4",
  "sf",
  "MuMIn",
  "PrevMap",
  "raster",
  "terra",
  "scales",
  "RColorBrewer",
  "cowplot",
  "magick"
)

invisible(lapply(required_pkgs, library, character.only = TRUE))

set.seed(01092025)

source("scripts/functions/spatial_fns/variogram_functions.R")

# ------------------------------------------------------
## 1. Load cleaned individual-level data ####
## ------------------------------------------------------

dat <- readRDS("data/derived/Toxo2003_full_cleaned_data.rds")
# dat <- readRDS("data/derived/Toxo2003_full_cleaned_data_deid.rds")


## ------------------------------------------------------
## 2. Variable selection ####
## ------------------------------------------------------

# Candidate variables (from multivariable results, p<0.05)
vars <- c(
  "agegroup", "sex", "scale(income_pcap)", "scale(elevation)", "scale(dist_road)",
  "cat", "contact_sewerwater"
)

f_pred <- as.formula(paste("toxo_igg ~", paste(vars, collapse = " + "), "+ (1|hh_id)"))

options(na.action = "na.fail")
m0 <- glmer(
  f_pred, data = dat, family = "binomial", nAGQ = 10,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

m0_d <- dredge(m0, m.lim = c(4, 7))
saveRDS(m0_d, "outputs/models/pred_model_glmm_dredge.rds")

# reload
# m0_d <- readRDS("outputs/models/pred_model_glmm_dredge.rds")

## ------------------------------------------------------
## 3. Create prediction grid  ####
## ------------------------------------------------------

# Study area outline (UTM/meters)
aoi <- st_read("data/spatial/PdL_cleaned_study_area.shp", quiet = TRUE) |>
  st_make_valid()

# Extend outline by 5 m buffer (so don't end up with empty cells near edges) ---
aoi_buf <- st_buffer(aoi, 5)

# Create grid points (4m by 4m)
grid_pts <- st_make_grid(aoi_buf, cellsize = 4, what = "centers") |>
  st_as_sf(crs = st_crs(aoi_buf)) |>
  st_filter(aoi_buf, .predicate = st_within) |>
  mutate(id = row_number())

xy <- st_coordinates(grid_pts)
grid_pts <- bind_cols(grid_pts, tibble(X = xy[,1], Y = xy[,2])) |>
  dplyr::select(id, X, Y)

write.csv(st_drop_geometry(grid_pts), "data/derived/prediction_grid_points_4m.csv", row.names = FALSE)

pred_grid <- read_csv("data/derived/prediction_grid_points_4m.csv", show_col_types = FALSE)

## ------------------------------------------------------
## 4. Intercept-only geostatistical model ####
## ------------------------------------------------------

# MCMC control for MCML iterations
cmcmc <- list(
  control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8),
  control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8),
  control.mcmc.MCML(n.sim = 65000, burnin = 5000, thin = 6)
)

dat_geo <- dat %>% 
  dplyr::select(
    toxo_igg, X, Y,
    agegroup, sex, income_pcap,
    elevation, dist_road, cat, contact_sewerwater) %>% 
  as.data.frame() %>%
  mutate(agegroup = relevel(factor(agegroup), ref = "10-12"),
         units.m = 1, # 1 individual per row
         ID = create.ID.coords(., coords = ~ X + Y)) # location IDs

# Starting values from simple GLM
f0 <- toxo_igg ~ 1
fit_binom0 <- glm(f0, family = "binomial", data = dat_geo)

par0 <- as.numeric(coef(fit_binom0))
p <- length(par0)

# par0: c(beta, sigma2, phi, tau2)
par0 <- c(par0, 1, 25, 1)

theta <- list()
theta[[1]] <- par0

# Fitting models with nugget effect (fixed.rel.nugget = NULL)
# ---- Iteration 1 ----
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])

geo_binomial_intercept <- binomial.logistic.MCML(formula = f0, 
                                                 units.m = ~ units.m,
                                                 coords = ~ X + Y,
                                                 data = dat_geo,
                                                 par0 = par0,
                                                 ID.coords = dat_geo$ID,
                                                 control.mcmc = cmcmc[[1]], 
                                                 kappa = 0.5,
                                                 start.cov.pars = init_cov_pars,
                                                 method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial_intercept))
theta[[2]] <- par0
saveRDS(geo_binomial_intercept, file = "outputs/models/geostat_intercept.rds")
print(summary(geo_binomial_intercept, l = F))

# ---- Iteration 2 ----
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])

geo_binomial_intercept <- binomial.logistic.MCML(
  formula = f0,
  units.m = ~ units.m,
  coords  = ~ X + Y,
  data    = dat_geo,
  par0    = par0,
  ID.coords = dat_geo$ID,
  control.mcmc = cmcmc[[2]],
  kappa = 0.5,
  start.cov.pars = init_cov_pars,
  method = "nlminb",
  messages = TRUE
)

par0 <- as.numeric(coef(geo_binomial_intercept))
theta[[3]] <- par0
saveRDS(geo_binomial_intercept, "outputs/models/geostat_intercept.rds")
print(summary(geo_binomial_intercept, l = FALSE))

# ---- Iteration 3 ----
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])

geo_binomial_intercept <- binomial.logistic.MCML(
  formula = f0,
  units.m = ~ units.m,
  coords  = ~ X + Y,
  data    = dat_geo,
  par0    = par0,
  ID.coords = dat_geo$ID,
  control.mcmc = cmcmc[[3]],
  kappa = 0.5,
  start.cov.pars = init_cov_pars,
  method = "nlminb",
  messages = TRUE
)

par0 <- as.numeric(coef(geo_binomial_intercept))
theta[[4]] <- par0
saveRDS(geo_binomial_intercept, "outputs/models/geostat_intercept.rds")
print(summary(geo_binomial_intercept, l = FALSE))

# Parameter table (95% CI)
se  <- sqrt(diag(geo_binomial_intercept$covariance))
uci <- geo_binomial_intercept$estimate + 1.96 * se
lci <- geo_binomial_intercept$estimate - 1.96 * se

# exponentiate spatial params stored on log scale
idx <- c("log(sigma^2)", "log(phi)", "log(nu^2)")
uci[idx] <- exp(uci[idx])
lci[idx] <- exp(lci[idx])

par_est <- tibble(
  par = names(coef(geo_binomial_intercept)),
  est = round(coef(geo_binomial_intercept), 3),
  lci = round(lci, 3),
  uci = round(uci, 3)
)

write_csv(par_est, "outputs/models/geostat_intercept_par.csv")
write_csv(par_est, "outputs/tables/tableS5_geostat_intercept_par.csv")

## Predictions (intercept-only)
cmcmc_pred <- control.mcmc.MCML(n.sim = 25000, burnin = 5000, thin = 10)

geo_pred_intercept <- spatial.pred.binomial.MCML(
  geo_binomial_intercept,
  grid.pred = pred_grid[, c("X", "Y")],
  control.mcmc = cmcmc_pred,
  type = "marginal",
  messages = TRUE,
  plot.correlogram = TRUE
)

saveRDS(geo_pred_intercept, "outputs/models/geostat_intercept_pred.rds")

geo_pred_intercept <- readRDS("outputs/models/geostat_intercept_pred.rds")

# make into a dataframe then a raster
prev_preds_intercept <- tibble(
  X    = geo_pred_intercept$grid$X,
  Y    = geo_pred_intercept$grid$Y,
  prev = geo_pred_intercept$prevalence$predictions,
  prev_u95ci = apply(geo_pred_intercept$samples, 1, function(draws) plogis(quantile(draws, 0.975))),
  prev_l95ci = apply(geo_pred_intercept$samples, 1, function(draws) plogis(quantile(draws, 0.025))),
  ex50 = rowMeans(geo_pred_intercept$samples > 0.5)
)

pred_rast_intercept <- rasterFromXYZ(prev_preds_intercept, crs = crs(aoi))

# quick plot check
plot(pred_rast_intercept[["prev"]])
plot(pred_rast_intercept[["prev_u95ci"]])
plot(pred_rast_intercept[["prev_l95ci"]])
plot(pred_rast_intercept[["ex50"]])

# Save raster (GeoTIFF)
terra::writeRaster(rast(pred_rast_intercept), "outputs/models/geostat_intercept_pred_raster.tif", overwrite = TRUE)

## ------------------------------------------------------
## 4. Full geostatistical model ####
## ------------------------------------------------------
f_full <- toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + I(dist_road/10) + cat + contact_sewerwater
fit_binom_full <- glm(f_full, family = "binomial", data = dat_geo)

par_covar <- as.numeric(coef(fit_binom_full))
p <- length(par_covar)

# Start from GLM betas + spatial params from intercept-only model
par0_full <- c(par_covar, theta[[4]][(length(theta[[4]]) - 2):length(theta[[4]])])  # last 3 are sigma2/phi/tau2-ish
theta_full <- list()
theta_full[[1]] <- par0_full

# ---- Iteration 1 ----
init_cov_pars <- c(par0_full[p + 2], par0_full[p + 3] / par0_full[p + 1])

geo_binomial_full <- binomial.logistic.MCML(
  formula = f_full,
  units.m  = ~ units.m,
  coords   = ~ X + Y,
  data     = dat_geo,
  par0     = par0_full,
  ID.coords = dat_geo$ID,
  control.mcmc = cmcmc[[1]],
  kappa = 0.5,
  start.cov.pars = init_cov_pars,
  method = "nlminb",
  messages = TRUE
)

par0_full <- as.numeric(coef(geo_binomial_full))
theta_full[[2]] <- par0_full
saveRDS(geo_binomial_full, "outputs/models/geostat_full.rds")
print(summary(geo_binomial_full, l = FALSE))

# ---- Iteration 2 ----
init_cov_pars <- c(par0_full[p + 2], par0_full[p + 3] / par0_full[p + 1])

geo_binomial_full <- binomial.logistic.MCML(
  formula = f_full,
  units.m  = ~ units.m,
  coords   = ~ X + Y,
  data     = dat_geo,
  par0     = par0_full,
  ID.coords = dat_geo$ID,
  control.mcmc = cmcmc[[2]],
  kappa = 0.5,
  start.cov.pars = init_cov_pars,
  method = "nlminb",
  messages = TRUE
)

par0_full <- as.numeric(coef(geo_binomial_full))
theta_full[[3]] <- par0_full
saveRDS(geo_binomial_full, "outputs/models/geostat_full.rds")
print(summary(geo_binomial_full, l = FALSE))

# ---- Iteration 3 ----
init_cov_pars <- c(par0_full[p + 2], par0_full[p + 3] / par0_full[p + 1])

geo_binomial_full <- binomial.logistic.MCML(
  formula = f_full,
  units.m  = ~ units.m,
  coords   = ~ X + Y,
  data     = dat_geo,
  par0     = par0_full,
  ID.coords = dat_geo$ID,
  control.mcmc = cmcmc[[3]],
  kappa = 0.5,
  start.cov.pars = init_cov_pars,
  method = "nlminb",
  messages = TRUE
)

par0_full <- as.numeric(coef(geo_binomial_full))
theta_full[[4]] <- par0_full
saveRDS(geo_binomial_full, "outputs/models/geostat_full.rds")
print(summary(geo_binomial_full, l = FALSE))

geo_binomial_full <- readRDS("outputs/models/geostat_full.rds")

# Parameter table (95% CI)
se <- sqrt(diag(geo_binomial_full$covariance))
uci <- geo_binomial_full$estimate + 1.96*se
lci<- geo_binomial_full$estimate - 1.96*se

par_est_full <- tibble(
  par = names(coef(geo_binomial_full)),
  est = round(exp(geo_binomial_full$estimate), 3),
  lci = round(exp(lci), 3),
  uci = round(exp(uci), 3)
)

write_csv(par_est_full, "outputs/tables/tableS6_geostat_full_par.csv")
write_csv(par_est_full, "outputs/models/geostat_full_par.csv")

# Compute the proportion of residual variance attributable to spatial process vs. nugget
sigma2 <- exp(geo_binomial_full$estimate["log(sigma^2)"])
tau2   <- exp(geo_binomial_full$estimate["log(nu^2)"])
prop_spat <- sigma2 / (sigma2 + tau2)
print(prop_spat) # i.e. About 31% of the unexplained variation is spatially structured over the range defined by φ (~30 m), and 69% occurs at very fine scales captured by the nugget.

## Predictions (full model): predict on logit scale using dummy covariates, then subtract the linear predictor for these covariates values from the prediction
# Create dummy predictor values
newdat <- pred_grid %>%
  mutate(
    agegroup = relevel(factor(levels(factor(dat_geo$agegroup))[1], levels = levels(factor(dat_geo$agegroup))), ref = "10-12"),
    sex      = factor(levels(dat_geo$sex)[1],      levels = levels(dat_geo$sex)),
    cat      = factor(levels(dat_geo$cat)[1],      levels = levels(dat_geo$cat)),
    contact_sewerwater = factor(levels(dat_geo$contact_sewerwater)[1], levels = levels(dat_geo$contact_sewerwater)),
    income_pcap = mean(dat_geo$income_pcap, na.rm = TRUE),
    elevation  = mean(dat_geo$elevation,  na.rm = TRUE),
    dist_road  = mean(dat_geo$dist_road,  na.rm = TRUE),
    `I(elevation/10)` = mean(dat_geo$elevation / 10, na.rm = TRUE),
    `I(dist_road/10)` = mean(dat_geo$dist_road / 10, na.rm = TRUE),
    units.m = 1
  ) %>%
  as.data.frame()

# Joint predictions on logit scale (eta = Xβ + S(x)), using dummy predictor values
cmcmc_pred_full <- control.mcmc.MCML(n.sim = 15000, burnin = 5000, thin = 40)

# run in chunks as high computational cost
chunk_size <- 5000
n_grid <- nrow(newdat)
chunk_id <- ceiling(seq_len(n_grid) / chunk_size)
idx_list <- split(seq_len(n_grid), chunk_id)

length(idx_list)  # how many chunks

# compute Xβ on the grid using the SAME formula/contrasts
rhs_terms  <- delete.response(terms(f_full))
Xgrid_full <- model.matrix(rhs_terms, data = newdat)

# extract fixed effects (drop spatial params), align to Xgrid_full
b_all <- coef(geo_binomial_full)
b_fix <- b_all[!names(b_all) %in% c("sigma^2","phi","tau^2")]

stopifnot(all(colnames(Xgrid_full) %in% names(b_fix)))

beta_hat  <- b_fix[colnames(Xgrid_full)]
fixed_hat <- as.numeric(Xgrid_full %*% beta_hat)

# run joint predictions by chunk
pred_list <- vector("list", length(idx_list))

for (k in seq_along(idx_list)) {
  ii <- idx_list[[k]]
  message("Chunk ", k, "/", length(idx_list), " (n=", length(ii), ")")
  
  geo_pred_k <- spatial.pred.binomial.MCML(
    geo_binomial_full,
    grid.pred         = newdat[ii, c("X","Y")],
    predictors        = newdat[ii, ],
    control.mcmc      = cmcmc_pred_full,
    type              = "joint",
    scale.predictions = "logit",
    messages          = FALSE,
    plot.correlogram  = FALSE
  )
  
  # logit predictions (Xβ + S(x)) for this chunk
  eta_k      <- as.numeric(geo_pred_k$logit$predictions)
  eta_l95_k  <- as.numeric(geo_pred_k$logit$quantiles[, 1])
  eta_u95_k  <- as.numeric(geo_pred_k$logit$quantiles[, 2])
  
  fixed_hat_k <- fixed_hat[ii]
  
  pred_list[[k]] <- tibble(
    X = newdat$X[ii],
    Y = newdat$Y[ii],
    
    eta_logit = eta_k,
    eta_l95   = eta_l95_k,
    eta_u95   = eta_u95_k,
    
    # (Only interpret these if newdat has meaningful covariates)
    prev      = plogis(eta_k),
    prev_l95  = plogis(eta_l95_k),
    prev_u95  = plogis(eta_u95_k),
    
    fixed_hat = fixed_hat_k,
    Sx_mean   = eta_k - fixed_hat_k
  )
  
  rm(geo_pred_k, eta_k, eta_l95_k, eta_u95_k)
  gc()
}

geo_pred_full <- bind_rows(pred_list) 

# quick checks
summary(geo_pred_full$Sx_mean)
hist(geo_pred_full$Sx_mean, breaks = 40)

saveRDS(geo_pred_full, "outputs/models/geostat_full_pred_chunked.rds")
geo_pred_full <- readRDS("outputs/models/geostat_full_pred_chunked.rds")

# Make a data.frame for rasterFromXYZ
# Include whichever layers you want to rasterise:
prev_preds_full <- geo_pred_full %>%
  dplyr::select(X, Y, prev, prev_l95, prev_u95, Sx_mean)

# Rasterise (creates a RasterBrick with multiple layers)
pred_rast_full <- raster::rasterFromXYZ(prev_preds_full, crs = raster::crs(aoi))

# Check plots
plot(pred_rast_full[["prev"]])
plot(pred_rast_full[["Sx_mean"]])

# Save raster (GeoTIFF)
terra::writeRaster(rast(pred_rast_full),
                   "outputs/models/geostat_full_pred_raster.tif",
                   overwrite = TRUE)

## ------------------------------------------------------
## 5. Map making ####
## ------------------------------------------------------

# helper: crop + mask (RETURN the object)
crop.fn <- function(r, outline) {
  r <- crop(r, extent(outline))
  r <- terra::mask(r, outline)
  r <- raster::raster(r) #  (change RasterLayer from Terra SpatRaster for better plotting)
  r
}

# param
pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(200)
lwd <- 1.5
res <- 300
xdim <- 8
ydim <- xdim * (7.94/10)
par(mar = c(3, 3, 5, 6)) 

legend_w   <- 1.3 # wider colourbar
legend_sh  <- 0.85 # longer colourbar
cex_leg    <- 1.5 # legend tick text size
cex_title  <- 1.8 # panel title size
title_line <- 1 # distance of title from plot
title_font <- 2 # bold

# if need to read in again
pred_rast_intercept <- rast("outputs/models/geostat_intercept_pred_raster.tif")
pred_rast_full <- rast("outputs/models/geostat_full_pred_raster.tif")

# crop rasters
p_intercept <- crop.fn(pred_rast_intercept[["prev"]], aoi)
p_intercept_lci <- crop.fn(pred_rast_intercept[["prev_l95ci"]], aoi)
p_intercept_uci <- crop.fn(pred_rast_intercept[["prev_u95ci"]], aoi)
p_intercept_ex50 <- crop.fn(pred_rast_intercept[["ex50"]], aoi)

Sx_full <- crop.fn(pred_rast_full[["Sx_mean"]], aoi)

# --- Predicted seroprevalence (intercept model)
p_intercept.min <- 0
p_intercept.max <- 0.8
at_A  <- seq(p_intercept.min, p_intercept.max, by = 0.2)
lab_A <- paste0(at_A * 100, "%")

tiff("outputs/figures/Fig4A_pred_p_intercept.tiff", units = "in", width = xdim, height = ydim, res = res, bg="white")

plot(st_geometry(aoi), lwd = lwd)
plot(p_intercept, col = rev(pal), xaxt = "n", yaxt = "n", axes = FALSE, box = FALSE,
     zlim = c(p_intercept.min, p_intercept.max), legend = FALSE, add = TRUE)

mtext("Predicted seroprevalence", side = 3, line = title_line, adj = 0.5, cex = cex_title, font = title_font)

plot(p_intercept, legend.only = TRUE, col = rev(pal),
     zlim = c(p_intercept.min, p_intercept.max),
     legend.width = legend_w, legend.shrink = legend_sh,
     axis.args = list(cex.axis = cex_leg, at = at_A, labels = lab_A))

plot(st_geometry(aoi), add = TRUE, lwd = lwd)
dev.off()

# --- Exceedance prob (>50%)
p_ex50_min <- floor(raster::cellStats(p_intercept_ex50, "min", na.rm = TRUE) * 100) / 100
p_ex50_max <- ceiling(raster::cellStats(p_intercept_ex50, "max", na.rm = TRUE) * 100) / 100

at_ex50  <- seq(p_ex50_min, p_ex50_max, by = 0.2)
lab_ex50 <- paste0(at_ex50 * 100, "%")

tiff("outputs/figures/Fig4B_pred_p_intercept_ex50.tiff", units = "in", width = xdim, height = ydim, res = res, bg="white")

plot(st_geometry(aoi), lwd = lwd)
plot(p_intercept_ex50, col = rev(pal), xaxt = "n", yaxt = "n", axes = FALSE, box = FALSE,
     zlim = c(p_ex50_min, p_ex50_max), legend = FALSE, add = TRUE)

mtext("Exceedance probability (>50%)", side = 3, line = title_line, adj = 0.5, cex = cex_title, font = title_font)

plot(p_intercept_ex50, legend.only = TRUE, col = rev(pal),
     zlim = c(p_ex50_min, p_ex50_max),
     legend.width = legend_w, legend.shrink = legend_sh,
     axis.args = list(cex.axis = cex_leg, at = at_ex50, labels = lab_ex50))

plot(st_geometry(aoi), add = TRUE, lwd = lwd)
dev.off()

# --- Residual spatial effect S(x)
Sx_min <- floor(raster::cellStats(Sx_full, "min", na.rm = TRUE))
Sx_max <- ceiling(raster::cellStats(Sx_full, "max", na.rm = TRUE))
lim <- max(abs(c(Sx_min, Sx_max)))
Sx_min <- -lim
Sx_max <-  lim

at_Sx  <- seq(Sx_min, Sx_max, by = 0.5)
lab_Sx <- format(at_Sx, trim = TRUE)

tiff("outputs/figures/Fig4C_pred_Sx_full.tiff", units = "in", width = xdim, height = ydim, res = res, bg="white")

plot(st_geometry(aoi), lwd = lwd)
plot(Sx_full, col = rev(pal), xaxt = "n", yaxt = "n", axes = FALSE, box = FALSE,
     zlim = c(Sx_min, Sx_max), legend = FALSE, add = TRUE)

mtext("Residual spatial effect S(x)", side = 3, line = title_line, adj = 0.5, cex = cex_title, font = title_font)

plot(Sx_full, legend.only = TRUE, col = rev(pal),
     zlim = c(Sx_min, Sx_max),
     legend.width = legend_w, legend.shrink = legend_sh,
     axis.args = list(cex.axis = cex_leg, at = at_Sx, labels = lab_Sx))

plot(st_geometry(aoi), add = TRUE, lwd = lwd)
dev.off()

# --- Panel figure (A/B/C)
read_panel <- function(path) {
  image_read(path) |>
    image_background("white", flatten = TRUE) |>
    image_border("white", "20x20")
}

imgA <- read_panel("outputs/figures/Fig4A_pred_p_intercept.tiff")
imgB <- read_panel("outputs/figures/Fig4B_pred_p_intercept_ex50.tiff")
imgC <- read_panel("outputs/figures/Fig4C_pred_Sx_full.tiff")

# Make a blank panel the same size as A
infoA <- image_info(imgA)
blank <- image_blank(width = infoA$width, height = infoA$height, color = "white")

# Compose 2x2
top <- image_append(c(imgA, imgB), stack = FALSE)
bot <- image_append(c(imgC, blank), stack = FALSE)
panel_img <- image_append(c(top, bot), stack = TRUE)


infoP <- image_info(panel_img)
W <- infoP$width
H <- infoP$height

# panel dimensions (2x2 layout)
halfW <- W / 2
halfH <- H / 2

# margins inside each panel (tune if needed)
pad_x <- round(0.045 * W)
pad_y <- round(0.01 * H)

label_size <- round(0.09 * halfH)  # scales nicely with figure size

panel_img <- panel_img |>
  image_annotate(
    "A",
    size = label_size, weight = 700, color = "black",
    location = paste0("+", pad_x, "+", pad_y)
  ) |>
  image_annotate(
    "B",
    size = label_size, weight = 700, color = "black",
    location = paste0("+", halfW + pad_x, "+", pad_y)
  ) |>
  image_annotate(
    "C",
    size = label_size, weight = 700, color = "black",
    location = paste0("+", pad_x, "+", halfH + pad_y)
  )
image_write(panel_img, "outputs/figures/Fig4_pred_panel.tiff", format = "tiff", density = "400x400")


