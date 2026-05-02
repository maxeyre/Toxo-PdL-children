################################################################################
# TOXO-PdL-children
################################################################################

# -------------------------------------------------------------------------
# Script: 03_female_piecewise_foi_sensitivity.R
#
# Purpose:
# Sensitivity analysis for reviewer comment on age-dependent exposure.
# Fits individual-level serocatalytic models among female participants:
#   1. Constant FOI model
#   2. Piecewise FOI model allowing FOI to change from age 12 onwards
#
# Model:
#   Constant:
#     p_i = 1 - exp(-lambda * age_i)
#
#   Piecewise:
#     H_i = lambda_pre * min(age_i, 12) +
#           lambda_post * max(age_i - 12, 0)
#     p_i = 1 - exp(-H_i)
#
# Inputs:
# - data/derived/Toxo2003_full_cleaned_data.rds
#
# -------------------------------------------------------------------------

## ------------------------------------------------------
## 0. Packages and global options ####
## ------------------------------------------------------

required_pkgs <- c(
  "tidyverse",
  "rjags",
  "coda",
  "MCMCvis",
  "loo",
  "binom"
)

invisible(lapply(required_pkgs, library, character.only = TRUE))

set.seed(01092025)

## ------------------------------------------------------
## 1. Load cleaned individual-level data ####
## ------------------------------------------------------

dat <- readRDS("data/derived/Toxo2003_full_cleaned_data.rds")

## Check required variables exist
required_vars <- c("sex", "age", "toxo_igg")
missing_vars <- setdiff(required_vars, names(dat))

if (length(missing_vars) > 0) {
  stop(
    "The following required variables are missing from dat: ",
    paste(missing_vars, collapse = ", ")
  )
}

## Restrict to female participants with non-missing age and T. gondii IgG status
dat_female <- dat %>%
  filter(
    sex == "Female",
    !is.na(age),
    !is.na(toxo_igg)
  ) %>%
  mutate(
    toxo_igg = as.integer(toxo_igg),
    age = as.numeric(age)
  )

summary(dat_female$age)
table(dat_female$toxo_igg)

## ------------------------------------------------------
## 2. Observed seroprevalence by age group for plotting ####
## ------------------------------------------------------

## Use existing agegroup if present; otherwise create 3-year groups
if (!"agegroup" %in% names(dat_female)) {
  dat_female <- dat_female %>%
    mutate(
      agegroup = cut(
        age,
        breaks = c(3, 6, 9, 12, 15, 18),
        labels = c("4-6", "7-9", "10-12", "13-15", "16-18"),
        right = TRUE
      )
    )
}

age_midpoints <- tibble(
  agegroup = factor(c("4-6", "7-9", "10-12", "13-15", "16-18"),
                    levels = c("4-6", "7-9", "10-12", "13-15", "16-18")),
  midpoint = c(5, 8, 11, 14, 17)
)

sero_female_grouped <- dat_female %>%
  mutate(
    agegroup = factor(agegroup, levels = c("4-6", "7-9", "10-12", "13-15", "16-18"))
  ) %>%
  left_join(age_midpoints, by = "agegroup") %>%
  group_by(agegroup, midpoint) %>%
  summarise(
    n = n(),
    n_pos = sum(toxo_igg, na.rm = TRUE),
    prev = n_pos / n,
    ci_lower = binom::binom.confint(n_pos, n, conf.level = 0.95, methods = "wilson")$lower,
    ci_upper = binom::binom.confint(n_pos, n, conf.level = 0.95, methods = "wilson")$upper,
    .groups = "drop"
  ) %>%
  arrange(midpoint)

print(sero_female_grouped)

## ------------------------------------------------------
## 3. JAGS model definitions ####
## ------------------------------------------------------

## Individual-level constant FOI model
jcode_constant <- "
model {
  for (i in 1:N) {
    y[i] ~ dbern(p[i])

    p[i] <- 1 - exp(-lambda * age[i])

    loglik[i] <- logdensity.bern(y[i], p[i])
  }

  lambda ~ dunif(0, 0.5)
}
"

## Individual-level piecewise FOI model
## FOI allowed to change from cutoff_age onwards
jcode_piecewise <- "
model {
  for (i in 1:N) {
    y[i] ~ dbern(p[i])

    H[i] <- lambda_pre * pre_time[i] + lambda_post * post_time[i]
    p[i] <- 1 - exp(-H[i])

    loglik[i] <- logdensity.bern(y[i], p[i])
  }

  lambda_pre  ~ dunif(0, 0.5)
  lambda_post ~ dunif(0, 0.5)

  lambda_diff  <- lambda_post - lambda_pre
  lambda_ratio <- lambda_post / lambda_pre

  ## This is 1 in posterior samples where post-12 FOI is lower than pre-12 FOI
  post_lower_indicator <- step(lambda_pre - lambda_post)
}
"

## ------------------------------------------------------
## 4. Fit piecewise FOI model with change at age 12 ####
## ------------------------------------------------------

run_piecewise_foi_individual <- function(dat_individual,
                                         cutoff_age = 12,
                                         mcmc_length = 20000,
                                         n_chains = 4,
                                         n_adapt = 3000,
                                         n_burn = 5000,
                                         out_prefix = "female_individual_piecewise12_foi") {
  
  dat_model <- dat_individual %>%
    mutate(
      pre_time = pmin(age, cutoff_age),
      post_time = pmax(age - cutoff_age, 0)
    )
  
  jdat <- list(
    y = dat_model$toxo_igg,
    age = dat_model$age,
    pre_time = dat_model$pre_time,
    post_time = dat_model$post_time,
    N = nrow(dat_model)
  )
  
  jmod <- rjags::jags.model(
    file = textConnection(jcode_piecewise),
    data = jdat,
    n.chains = n_chains,
    n.adapt = n_adapt
  )
  
  update(jmod, n.iter = n_burn)
  
  jpos <- rjags::coda.samples(
    model = jmod,
    variable.names = c(
      "lambda_pre",
      "lambda_post",
      "lambda_diff",
      "lambda_ratio",
      "post_lower_indicator",
      "loglik"
    ),
    n.iter = mcmc_length
  )
  
  mcmc_mat <- as.matrix(jpos)
  
  ## Diagnostics
  pdf(paste0("outputs/models/", out_prefix, "_trace.pdf"), width = 10, height = 7)
  plot(jpos)
  dev.off()
  
  mcmc_summary <- MCMCvis::MCMCsummary(jpos, round = 4)
  write.csv(
    mcmc_summary,
    paste0("outputs/models/", out_prefix, "_MCMCsummary.csv"),
    row.names = FALSE
  )
  
  ## Fit metrics
  loglik <- mcmc_mat[, grepl("^loglik\\[", colnames(mcmc_mat)), drop = FALSE]
  waic_out <- loo::waic(loglik)
  loo_out <- loo::loo(loglik)
  dic_out <- rjags::dic.samples(jmod, n.iter = mcmc_length)
  
  ## Posterior summaries
  lambda_pre_draws <- mcmc_mat[, "lambda_pre"]
  lambda_post_draws <- mcmc_mat[, "lambda_post"]
  lambda_diff_draws <- mcmc_mat[, "lambda_diff"]
  lambda_ratio_draws <- mcmc_mat[, "lambda_ratio"]
  post_lower_draws <- mcmc_mat[, "post_lower_indicator"]
  
  key_estimates <- tibble(
    parameter = c(
      "lambda_pre_12",
      "lambda_post_12",
      "lambda_difference_post_minus_pre",
      "lambda_ratio_post_over_pre",
      "posterior_probability_post_foi_lower",
      "annual_incidence_probability_pre_12",
      "annual_incidence_probability_post_12"
    ),
    median = c(
      median(lambda_pre_draws),
      median(lambda_post_draws),
      median(lambda_diff_draws),
      median(lambda_ratio_draws),
      mean(post_lower_draws),
      median((1 - exp(-lambda_pre_draws)) * 100),
      median((1 - exp(-lambda_post_draws)) * 100)
    ),
    l95 = c(
      quantile(lambda_pre_draws, 0.025),
      quantile(lambda_post_draws, 0.025),
      quantile(lambda_diff_draws, 0.025),
      quantile(lambda_ratio_draws, 0.025),
      NA,
      quantile((1 - exp(-lambda_pre_draws)) * 100, 0.025),
      quantile((1 - exp(-lambda_post_draws)) * 100, 0.025)
    ),
    u95 = c(
      quantile(lambda_pre_draws, 0.975),
      quantile(lambda_post_draws, 0.975),
      quantile(lambda_diff_draws, 0.975),
      quantile(lambda_ratio_draws, 0.975),
      NA,
      quantile((1 - exp(-lambda_pre_draws)) * 100, 0.975),
      quantile((1 - exp(-lambda_post_draws)) * 100, 0.975)
    )
  )
  
  write_csv(
    key_estimates,
    paste0("outputs/models/", out_prefix, "_key_estimates.csv")
  )
  
  list(
    cutoff_age = cutoff_age,
    jmod = jmod,
    jpos = jpos,
    mcmc_mat = mcmc_mat,
    key_estimates = key_estimates,
    waic = waic_out,
    loo = loo_out,
    dic = dic_out
  )
}

fit_piecewise12_female_individual <- run_piecewise_foi_individual(
  dat_individual = dat_female,
  cutoff_age = 12
)

fit_piecewise12_female_individual$key_estimates
fit_piecewise12_female_individual$waic
fit_piecewise12_female_individual$loo
