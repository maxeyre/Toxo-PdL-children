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

# Group by 1-year age
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

# Function to make age grouped data for neater plotting (e.g. for 3-year age groups)
make_obs_age_groups <- function(df, width = 3, min_age = 4, max_age = 18) {
  
  breaks <- seq(min_age, max_age + 1, by = width)
  
  if (tail(breaks, 1) < max_age + 1) {
    breaks <- c(breaks, max_age + 1)
  }
  
  age_levels <- tibble(
    age_min = head(breaks, -1),
    age_max = tail(breaks, -1) - 1,
    midpoint = (age_min + age_max) / 2,
    agegroup = paste0(age_min, "-", age_max)
  ) %>%
    filter(age_min <= max_age) %>%
    mutate(
      age_max = pmin(age_max, max_age),
      midpoint = (age_min + age_max) / 2,
      agegroup = paste0(age_min, "-", age_max)
    ) %>%
    arrange(midpoint)
  
  df %>%
    mutate(
      agegroup = cut(
        age,
        breaks = breaks,
        right = FALSE,
        include.lowest = TRUE,
        labels = age_levels$agegroup
      )
    ) %>%
    left_join(age_levels, by = "agegroup") %>%
    group_by(sex, agegroup, midpoint) %>%
    summarise(
      n = sum(n),
      n_pos = sum(n_pos),
      prev = n_pos / n,
      ci_lower = binom.confint(n_pos, n, conf.level = 0.95, methods = "wilson")$lower,
      ci_upper = binom.confint(n_pos, n, conf.level = 0.95, methods = "wilson")$upper,
      .groups = "drop"
    ) %>%
    arrange(sex, midpoint) %>%
    mutate(
      agegroup = factor(agegroup, levels = age_levels$agegroup),
      sex = relevel(sex, ref = "Male")
    )
}

sero_obs_3yr <- make_obs_age_groups(seroDatGrouped, width = 3)
# or:
# sero_obs_5yr <- make_obs_age_groups(seroDatGrouped, width = 5)

sero_male_obs_plot <- sero_obs_3yr %>% filter(sex == "Male")
sero_female_obs_plot <- sero_obs_3yr %>% filter(sex == "Female")

# Plot: seroprevalence bars with CIs
col_male   <- "#2c7bb6"
col_female <- "#d7191c"

seroprev_plot <- ggplot(sero_obs_3yr, aes(x = agegroup, y = 100 * prev, fill = sex)) +
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

# ggsave("outputs/figures/Fig2A_seroprev_by_age_sex.tiff", seroprev_plot, width = 10, height = 6, dpi = 500)


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
  
  # Save outputs for downstream plotting/panels
  saveRDS(df_mod,      paste0("outputs/models/", out_prefix, "_modelUncertainty.rds"))
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
    df_mod = df_mod
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

# Get binomial sampling uncertainty for light shaded ribbon in plot
make_sampling_uncertainty_for_plot <- function(df_obs_plot,
                                               mcmcMatrix,
                                               lambda_q,
                                               n_samples_uncert = 1000,
                                               lower_age = 4,
                                               upper_age = 18) {
  
  sampled <- mcmcRandomSamplerCat(
    n_samples_uncert,
    mcmcMatrix,
    df_obs_plot$midpoint,
    df_obs_plot$n
  )
  
  aq <- ageQuantiles(sampled)
  
  df_sampling <- data.frame(
    midpoint = df_obs_plot$midpoint,
    prev = 1 - exp(-lambda_q["median"] * df_obs_plot$midpoint),
    ci_lower = aq[, 2],
    ci_upper = aq[, 3]
  ) %>%
    mutate(
      lower_width = prev - ci_lower,
      upper_width = ci_upper - prev
    )
  
  boundary_ages <- c(lower_age, upper_age)
  boundary_prev <- 1 - exp(-lambda_q["median"] * boundary_ages)
  
  boundary_lower_width <- approx(
    x = df_sampling$midpoint,
    y = df_sampling$lower_width,
    xout = boundary_ages,
    rule = 2
  )$y
  
  boundary_upper_width <- approx(
    x = df_sampling$midpoint,
    y = df_sampling$upper_width,
    xout = boundary_ages,
    rule = 2
  )$y
  
  boundary_rows <- data.frame(
    midpoint = boundary_ages,
    prev = boundary_prev,
    ci_lower = pmax(0, boundary_prev - boundary_lower_width),
    ci_upper = pmin(1, boundary_prev + boundary_upper_width)
  )
  
  bind_rows(
    df_sampling %>% select(midpoint, prev, ci_lower, ci_upper),
    boundary_rows
  ) %>%
    arrange(midpoint) %>%
    distinct(midpoint, .keep_all = TRUE)
}

df_sampling_male_plot <- make_sampling_uncertainty_for_plot(
  df_obs_plot = sero_male_obs_plot,
  mcmcMatrix = fit_male$mcmcMatrix,
  lambda_q = fit_male$lambda_q,
  n_samples_uncert = 1000
)

df_sampling_female_plot <- make_sampling_uncertainty_for_plot(
  df_obs_plot = sero_female_obs_plot,
  mcmcMatrix = fit_female$mcmcMatrix,
  lambda_q = fit_female$lambda_q,
  n_samples_uncert = 1000
)

make_serocatalytic_plot <- function(df_mod,
                                    df_obs,
                                    lambda_q,
                                    sex_label,
                                    col_fill,
                                    df_sampling = NULL,
                                    show_sampling_uncertainty = FALSE) {
  
  FOI <- lambda_q
  AIP <- (1 - exp(-FOI)) * 100
  
  label <- paste0(
    "FOI = ", sprintf("%.3f", FOI["median"]),
    " (95% CrI: ", sprintf("%.3f", FOI["l95"]), ", ", sprintf("%.3f", FOI["u95"]), ")",
    "\nAIP = ", sprintf("%.2f", AIP["median"]), "% (95% CrI: ",
    sprintf("%.2f", AIP["l95"]), "%, ", sprintf("%.2f", AIP["u95"]), "%)"
  )
  
  p <- ggplot(df_mod, aes(x = midpoint, y = 100 * prev))
  
  if (show_sampling_uncertainty && !is.null(df_sampling)) {
    p <- p +
      geom_ribbon(
        data = df_sampling,
        aes(x = midpoint, ymin = 100 * ci_lower, ymax = 100 * ci_upper),
        inherit.aes = FALSE,
        alpha = 0.2,
        fill = col_fill
      )
  }
  
  p +
    geom_ribbon(
      aes(ymin = 100 * ci_lower, ymax = 100 * ci_upper),
      alpha = 0.35,
      fill = col_fill
    ) +
    geom_line(colour = "black") +
    geom_point(
      data = df_obs,
      aes(x = midpoint, y = 100 * prev),
      inherit.aes = FALSE
    ) +
    geom_linerange(
      data = df_obs,
      aes(x = midpoint, ymin = 100 * ci_lower, ymax = 100 * ci_upper),
      inherit.aes = FALSE
    ) +
    scale_x_continuous(breaks = seq(4, 18, by = 2), limits = c(4, 18)) +
    scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
    labs(x = "Age (years)", y = "Seroprevalence (%)") +
    annotate("text", x = 4.3, y = 85, label = label, hjust = 0, vjust = 0, size = 4, color = "black") +
    theme_bw(base_size = 16)
}

catPlot_male <- make_serocatalytic_plot(
  df_mod      = fit_male$df_mod,
  df_obs      = sero_male_obs_plot,
  lambda_q    = fit_male$lambda_q,
  sex_label   = "Male",
  col_fill    = col_male,
  df_sampling = df_sampling_male_plot,
  show_sampling_uncertainty = TRUE # change to FALSE to remove light shaded binomial uncertainty
)

catPlot_female <- make_serocatalytic_plot(
  df_mod      = fit_female$df_mod,
  df_obs      = sero_female_obs_plot,
  lambda_q    = fit_female$lambda_q,
  sex_label   = "Female",
  col_fill    = col_female,
  df_sampling = df_sampling_female_plot,
  show_sampling_uncertainty = TRUE
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
# ggsave("outputs/figures/Fig2_seroprev_serocatalytic_panel.tiff", Fig_sero_panel,
       # width = 10, height = 10, dpi = 500)

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

