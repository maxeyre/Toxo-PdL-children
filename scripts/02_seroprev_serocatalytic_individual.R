################################################################################
# TOXO-PdL-children
################################################################################

# -------------------------------------------------------------------------
# Script: 02_seroprev_serocatalytic_individual.R
#
# Purpose:
# Estimate age- and sex-specific seroprevalence and fit individual-level Bayesian serocatalytic
# models to estimate the force of infection.
#
# Inputs:
# - data/derived/Toxo2003_full_cleaned_data.rds  (from 01_setup_import_clean.R)
# -------------------------------------------------------------------------

## ------------------------------------------------------
## 0. Packages & global options ####
## ------------------------------------------------------

required_pkgs <- c(
  "tidyverse",
  "binom",
  "rjags",
  "coda",
  "MCMCvis",
  "loo",
  "cowplot"
)
invisible(lapply(required_pkgs, library, character.only = TRUE))

source("scripts/functions/catalyticModelFunctions.R")
set.seed(01092025)

## ------------------------------------------------------
## 1. Load cleaned individual-level data ####
## ------------------------------------------------------

dat <- readRDS("data/derived/Toxo2003_full_cleaned_data.rds")
# dat <- readRDS("data/derived/Toxo2003_full_cleaned_data_deid.rds")

## ------------------------------------------------------
## 2. Seroprevalence by sex x agegroup (binomial CIs) ####
## ------------------------------------------------------

seroDatGrouped <- dat %>%
  group_by(sex, age) %>%
  summarise(
    n = n(),
    n_pos = sum(toxo_igg, na.rm = TRUE),
    prev = n_pos / n,
    ci_lower = binom.confint(n_pos, n, conf.level = 0.95, methods = "wilson")$lower,
    ci_upper = binom.confint(n_pos, n, conf.level = 0.95, methods = "wilson")$upper,
    .groups = "drop"
  ) %>%
  arrange(sex, age) %>%
  mutate(sex = relevel(sex, ref="Male"))


# Plot: seroprevalence bars with CIs
col_male   <- "#2c7bb6"
col_female <- "#d7191c"

seroprev_plot <- ggplot(seroDatGrouped, aes(x = agegroup, y = 100 * prev, fill = sex)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.5, colour = "black") +
  geom_errorbar(
    aes(ymin = 100 * ci_lower, ymax = 100 * ci_upper),
    width = 0.3,
    position = position_dodge(width = 0.9)
  ) +
  geom_text(aes(label = paste0(sprintf("%.1f", 100 * prev), "%","\n","n=",n_pos,"/", n),
                y = 100*ci_upper),
            position = position_dodge(0.9), vjust = -.5, size = 4) +
  scale_fill_manual(values = c("Male" = col_male, "Female" = col_female)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  labs(x = "Age (years)", y = "Seroprevalence (%)") +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

seroprev_plot

ggsave("outputs/figures/Fig2A_seroprev_by_age_sex.tiff", seroprev_plot, width = 10, height = 6, dpi = 500)


## ------------------------------------------------------
## 3. JAGS serocatalytic model definition ####
## ------------------------------------------------------

# Catalytic model (constant FOI):
# P(seropos | age) = 1 - exp(-lambda * age)

jcode <- "
model{
  for (i in 1:N){
    n_pos[i] ~ dbin(p[i], n[i])
    
    # Catalytic model
    p[i] <- 1 - exp(-lambda * age[i])
    
    # Calculate likelihood (used for WAIC)
    loglik[i] <- logdensity.bin(n_pos[i], p[i], n[i])
  }
  ## Priors
  lambda ~ dunif(0, 0.5)
}"

## ------------------------------------------------------
## 4. Helper: run model + diagnostics + uncertainty bands ####
## ------------------------------------------------------

run_serocatalytic <- function(df_sex,
                              sex_label,
                              mcmc_length = 10000,
                              n_chains = 4,
                              n_adapt = 2000,
                              n_burn = 2000,
                              n_samples_uncert = 1000,
                              age_range = 4:18,
                              out_prefix = NULL) {
  
  #stopifnot(all(c("midpoint", "n", "n_pos") %in% names(df_sex)))
  if (is.null(out_prefix)) out_prefix <- paste0("sero_", tolower(sex_label))
  
  # JAGS data
  jdat <- list(
    n_pos = df_sex$n_pos,
    age   = df_sex$age,
    n     = df_sex$n,
    N = nrow(df_sex)
  )
  
  # Compile + burn-in
  jmod <- rjags::jags.model(
    file = textConnection(jcode),
    data = jdat,
    n.chains = n_chains,
    n.adapt = n_adapt
  )
  update(jmod, n.iter = n_burn)
  
  jpos <- coda.samples(
    model = jmod,
    variable.names = c("lambda", "loglik"),
    n.iter = mcmc_length
  )
  
  ## ---- Diagnostics ----
  # Trace/density for lambda + logliks
  # (If many logliks, plot can be busy; still useful for trace mixing.)
  pdf(file = paste0("outputs/models/", out_prefix, "_mcmc_trace.pdf"), width = 10, height = 7)
  plot(jpos)
  dev.off()
  
  mcmc_summary <- MCMCvis::MCMCsummary(jpos, round = 3)
  write.csv(mcmc_summary, paste0("outputs/models/", out_prefix, "_MCMCsummary.csv"), row.names = FALSE)
  
  # Matrix of posterior draws
  mcmcMatrix <- as.matrix(jpos)
  
  # Posterior density plots (facetted)
  mcmcDF <- tibble::as_tibble(mcmcMatrix) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "param", values_to = "value")
  
  dens_plot <- ggplot(mcmcDF, aes(value)) +
    geom_density() +
    facet_wrap(~ param, scales = "free") +
    theme_bw(base_size = 14)
  
  ggsave(paste0("outputs/models/", out_prefix, "_posterior_densities.png"),
         dens_plot, width = 10, height = 6, dpi = 300)
  
  # Point estimates for lambda
  lambda_draws <- mcmcMatrix[, "lambda"]
  lambda_q <- quantile(lambda_draws, probs = c(0.5, 0.025, 0.975))
  names(lambda_q) <- c("median", "l95", "u95")
  
  # DIC
  dic <- dic.samples(jmod, n.iter = mcmc_length)
  
  # WAIC + LOO from pointwise loglik columns (robust extraction)
  logLik <- mcmcMatrix[, grepl("^loglik\\[", colnames(mcmcMatrix)), drop = FALSE]
  waic_out <- loo::waic(logLik)
  loo_out  <- loo::loo(logLik)
  
  saveRDS(list(dic = dic, waic = waic_out, loo = loo_out),
          file = paste0("outputs/models/", out_prefix, "_fit_metrics.rds"))
  
  ## ---- Uncertainty bands for catalytic curve ----
  # 1) Model uncertainty: use quantiles of lambda
  df_mod <- data.frame(
    midpoint = age_range,
    prev  = 1 - exp(-lambda_q["median"] *age_range),
    ci_lower = 1 - exp(-lambda_q["l95"]    *age_range),
    ci_upper = 1 - exp(-lambda_q["u95"]    *age_range)
  )
  
  # 2) Sampling uncertainty: requires helper functions in catalyticModelFunctions.R
  #    - mcmcRandomSamplerCat(n, mcmcMatrix, midpoints, totals)
  #    - ageQuantiles(samples)
  midpoints <- df_sex$age
  totals    <- df_sex$n
  
  # sample the chain + add binomial sampling variability
  sampled <- mcmcRandomSamplerCat(n_samples_uncert, mcmcMatrix, midpoints, totals)
  aq <- ageQuantiles(sampled)
  
  df_sampling <- data.frame(
    midpoint = df_sex$age,
    prev  = 1 - exp(-df_sex$age * lambda_q["median"]),
    ci_lower = aq[, 2],
    ci_upper = aq[, 3]
  )
  
  df_sampling <- extend_sampling_to_bounds(df_sampling, lower_age = 4, upper_age = 18)
  
  # Save outputs for downstream plotting/panels
  saveRDS(df_mod,      paste0("outputs/models/", out_prefix, "_modelUncertainty.rds"))
  saveRDS(df_sampling, paste0("outputs/models/", out_prefix, "_samplingUncertainty.rds"))
  saveRDS(lambda_q,    paste0("outputs/models/", out_prefix, "_lambda_quantiles.rds"))
  
  list(
    sex = sex_label,
    jmod = jmod,
    jpos = jpos,
    mcmcMatrix = mcmcMatrix,
    lambda_q = lambda_q,
    dic = dic,
    waic = waic_out,
    loo = loo_out,
    df_mod = df_mod,
    df_sampling = df_sampling
  )
}

## ------------------------------------------------------
## 5. Run serocatalytic model for males + females ####
## ------------------------------------------------------
sero_male   <- seroDatGrouped %>% filter(sex == "Male")
sero_female <- seroDatGrouped %>% filter(sex == "Female")

## MALE
fit_male <- run_serocatalytic(
  df_sex = sero_male,
  sex_label = "Male",
  mcmc_length = 10000,
  n_chains = 4,
  n_adapt = 2000,
  n_burn = 2000,
  n_samples_uncert = 1000,
  age_range = 4:18,
  out_prefix = "serocatalytic_male"
)

# diagnostic checks
plot(fit_male$jpos) # check convergence
summary(fit_male$jpos)
MCMCsummary(fit_male$jpos, round = 2) ## Check ESS and Rhat

# point estimate and CIs
lambdaEst_male <- fit_male$lambda_q

# Loo and WAIC
fit_male$loo
fit_male$waic

## FEMALE
fit_female <- run_serocatalytic(
  df_sex = sero_female,
  sex_label = "Female",
  mcmc_length = 10000,
  n_chains = 4,
  n_adapt = 2000,
  n_burn = 2000,
  n_samples_uncert = 1000,
  age_range = 4:18,
  out_prefix = "serocatalytic_female"
)

# diagnostic checks
plot(fit_female$jpos) # check convergence
summary(fit_female$jpos)
MCMCsummary(fit_female$jpos, round = 2) ## Check ESS and Rhat

# point estimate and CIs
lambdaEst_female <- fit_female$lambda_q
lambdaEst_female
# Loo and WAIC
fit_female$loo
fit_female$waic


## ------------------------------------------------------
## 6. Plot serocatalytic curves (male + female) + combine with seroprev bar plot #####
## ------------------------------------------------------
make_serocatalytic_plot <- function(df_mod, df_sampling, df_obs, lambda_q, sex_label, col_fill) {
  
  FOI <- lambda_q
  AIP <- (1 - exp(-FOI)) * 100
  
  label <- paste0(
    "FOI = ", sprintf("%.3f", FOI["median"]),
    " (95% CrI: ", sprintf("%.3f", FOI["l95"]), ", ", sprintf("%.3f", FOI["u95"]), ")",
    "\nAIP = ", sprintf("%.2f", AIP["median"]), "% (95% CrI: ",
    sprintf("%.2f", AIP["l95"]), "%, ", sprintf("%.2f", AIP["u95"]), "%)"
  )
  
  ggplot(df_mod, aes(x = midpoint, y = 100 * prev)) +
    geom_ribbon(aes(ymin = 100 * ci_lower, ymax = 100 * ci_upper), alpha = 0.3, fill = col_fill) +
    geom_line(colour = "black") +
    geom_ribbon(
      data = df_sampling,
      aes(x = midpoint, ymin = 100 * ci_lower, ymax = 100 * ci_upper),
      alpha = 0.30,
      fill = col_fill
    ) +
    geom_point(data = df_obs, aes(x = midpoint, y = 100 * prev), inherit.aes = FALSE) +
    geom_linerange(data = df_obs, aes(x = midpoint, ymin = 100 * ci_lower, ymax = 100 * ci_upper),
                   inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(4, 18, by = 2), limits = c(4, 18)) +
    scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
    labs(x = "Age (years)", y = "Seroprevalence (%)") +
    annotate("text", x = 4.3, y = 85, label = label, hjust = 0, vjust = 0, size = 4, color = "black") +
    theme_bw(base_size = 16)
}

catPlot_male <- make_serocatalytic_plot(
  df_mod      = fit_male$df_mod,
  df_sampling = fit_male$df_sampling,
  df_obs      = sero_male,
  lambda_q    = fit_male$lambda_q,
  sex_label   = "Male",
  col_fill    = col_male
)

catPlot_female <- make_serocatalytic_plot(
  df_mod      = fit_female$df_mod,
  df_sampling = fit_female$df_sampling,
  df_obs      = sero_female,
  lambda_q    = fit_female$lambda_q,
  sex_label   = "Female",
  col_fill    = col_female
)

# Combine into one figure (A: seroprev bars; B/C: serocatalytic)
top_row <- cowplot::plot_grid(
  seroprev_plot, labels = "A", ncol = 1,
  label_fontface = "bold", label_size = 24, label_x = -0.01
)

bottom_row <- cowplot::plot_grid(
  catPlot_male, catPlot_female,
  labels = c("B", "C"), ncol = 2,
  label_fontface = "bold", label_size = 24, label_x = -0.01
)

Fig_sero_panel <- cowplot::plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
Fig_sero_panel
ggsave("outputs/figures/Fig2_seroprev_serocatalytic_panel.tiff", Fig_sero_panel,
       width = 10, height = 10, dpi = 500)

## ------------------------------------------------------ 
## 7. Save key numeric outputs for manuscript text ####
## ------------------------------------------------------
key_out <- tibble(
  sex = c("Male", "Female"),
  FOI_median = c(fit_male$lambda_q["median"], fit_female$lambda_q["median"]),
  FOI_l95    = c(fit_male$lambda_q["l95"],    fit_female$lambda_q["l95"]),
  FOI_u95    = c(fit_male$lambda_q["u95"],    fit_female$lambda_q["u95"])
) %>%
  mutate(
    AIP_median = (1 - exp(-FOI_median)) * 100,
    AIP_l95    = (1 - exp(-FOI_l95)) * 100,
    AIP_u95    = (1 - exp(-FOI_u95)) * 100
  )

write_csv(key_out, "outputs/models/serocatalytic_key_estimates.csv")

