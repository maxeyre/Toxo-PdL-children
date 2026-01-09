################################################################################
# TOXO-PdL-children
################################################################################

# -------------------------------------------------------------------------
# Script: 03_descriptive_gams_regression.R
#
# Purpose:
# (1) Fit and plot univariable GAM smooths for selected continuous predictors.
# (2) Fit univariable household-random-intercept logistic mixed models for a
#     broad set of candidate risk factors and export a summary table.
# (3) Fit DAG-informed multivariable mixed models, calculate E-values, and
#     export a multivariable table + forest plot.
#
# Inputs:
# - data/derived/Toxo2003_full_cleaned_data.rds  (from 01_setup_import_clean.R)
# -------------------------------------------------------------------------

## ------------------------------------------------------
## 0. Packages & global options ####
## ------------------------------------------------------

required_pkgs <- c(
  "tidyverse",
  "lme4",
  "performance",
  "broom.mixed",
  "mgcv",
  "mgcViz",
  "car",
  "binom",
  "EValue",
  "cowplot",
  "purrr", 
  "officer", 
  "flextable"
)

invisible(lapply(required_pkgs, library, character.only = TRUE))

set.seed(01092025)

## ------------------------------------------------------
## 1. Load cleaned individual-level data ####
## ------------------------------------------------------

dat <- readRDS("data/derived/Toxo2003_full_cleaned_data.rds")
# dat <- readRDS("data/derived/Toxo2003_full_cleaned_data_deid.rds")

## ------------------------------------------------------
## 2. Table 1 descriptives (Total + Seropositive) by domain ####
## ------------------------------------------------------

# Helper: median (IQR) formatter
fmt_median_iqr <- function(x, digits = 1) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  paste0(
    round(q[2], digits), " (",
    round(q[1], digits), ", ",
    round(q[3], digits), ")"
  )
}

# Helper: n (%) formatter
fmt_n_pct <- function(n, denom, digits = 1) {
  paste0(n, " (", round(100 * n / denom, digits), "%)")
}

# Helper: make categorical block (works for Yes/No, sex, race, etc.)
cat_block <- function(data, var, domain, label_map = NULL) {
  v <- rlang::ensym(var)
  
  out <- data %>%
    dplyr::filter(!is.na(!!v)) %>%
    dplyr::group_by(!!v) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_pos = sum(toxo_igg),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      domain = domain,
      variable = as.character(!!v),
      total = fmt_n_pct(n, nrow(data)),
      seropos = fmt_n_pct(n_pos, n)
    ) %>%
    dplyr::select(domain, variable, total, seropos)
  
  # relabel levels
  if (!is.null(label_map)) {
    out <- out %>%
      dplyr::mutate(variable = dplyr::recode(variable, !!!label_map))
  }
  
  out
}

# Helper: continuous block with optional scaling label (e.g., show per 50m later in ORs)
cont_block <- function(data, var, domain, label, digits = 1) {
  v <- rlang::ensym(var)
  tibble::tibble(
    domain = domain,
    variable = label,
    total = fmt_median_iqr(dplyr::pull(data, !!v), digits = digits),
    seropos = fmt_median_iqr(dplyr::pull(dplyr::filter(data, toxo_igg == 1), !!v), digits = digits)
  )
}

N_total <- nrow(dat)

## --- 1) Demographic & socioeconomic ---
tab_demo <- dplyr::bind_rows(
  
  cont_block(dat, age, "Demographic & socioeconomic", "Age (continuous)", digits = 0),
  
  dat %>%
    dplyr::group_by(agegroup) %>%
    dplyr::summarise(n = dplyr::n(), n_pos = sum(toxo_igg), .groups = "drop") %>%
    dplyr::mutate(
      domain = "Demographic & socioeconomic",
      variable = as.character(agegroup),
      total = fmt_n_pct(n, N_total),
      seropos = fmt_n_pct(n_pos, n)
    ) %>%
    dplyr::select(domain, variable, total, seropos),
  
  cat_block(dat, sex,  "Demographic & socioeconomic"),
  cat_block(dat, race, "Demographic & socioeconomic"),
  
  cont_block(dat, income_pcap, "Demographic & socioeconomic",
             "Per capita daily household income (US$)", digits = 2),
  
  cat_block(dat, house_rented, "Demographic & socioeconomic",
            label_map = c("0" = "House is rented: No", "1" = "House is rented: Yes",
                          "No" = "House is rented: No", "Yes" = "House is rented: Yes")),
  
  cat_block(dat, house_title, "Demographic & socioeconomic",
            label_map = c("0" = "Owns the title to the household: No", "1" = "Owns the title to the household: Yes",
                          "No" = "Owns the title to the household: No", "Yes" = "Owns the title to the household: Yes"))
)

## --- 2) Household animals ---
tab_animals <- dplyr::bind_rows(
  cat_block(dat, cat,     "Household animals",
            label_map = c("0"="Cat in household: No","1"="Cat in household: Yes","No"="Cat in household: No","Yes"="Cat in household: Yes")),
  cat_block(dat, dog,     "Household animals",
            label_map = c("0"="Dog in household: No","1"="Dog in household: Yes","No"="Dog in household: No","Yes"="Dog in household: Yes")),
  cat_block(dat, chicken, "Household animals",
            label_map = c("0"="Raise chickens: No","1"="Raise chickens: Yes","No"="Raise chickens: No","Yes"="Raise chickens: Yes")),
  cat_block(dat, rats_observed, "Household animals",
            label_map = c("0"="Observed rats in or near house: No","1"="Observed rats in or near house: Yes",
                          "No"="Observed rats in or near house: No","Yes"="Observed rats in or near house: Yes"))
)

## --- 3) Household & peridomestic environment ---
tab_env <- dplyr::bind_rows(
  cont_block(dat, elevation,  "Household & peridomestic environment", "Household elevation (m)", digits = 1),
  cont_block(dat, dist_road,  "Household & peridomestic environment", "Distance to the main road (m)", digits = 1),
  cont_block(dat, dist_trash, "Household & peridomestic environment", "Distance to nearest trash dump (m)", digits = 1),
  cont_block(dat, dist_sewer, "Household & peridomestic environment", "Distance to nearest open sewer (m)", digits = 1),
  
  cat_block(dat, hh_floods, "Household & peridomestic environment",
            label_map = c("0"="House flooded in last 6 months: No","1"="House flooded in last 6 months: Yes",
                          "No"="House flooded in last 6 months: No","Yes"="House flooded in last 6 months: Yes")),
  
  cat_block(dat, veg, "Household & peridomestic environment",
            label_map = c("0"="Vegetation within 10 m of the house: No","1"="Vegetation within 10 m of the house: Yes",
                          "No"="Vegetation within 10 m of the house: No","Yes"="Vegetation within 10 m of the house: Yes"))
)

## --- 4) Contact with environment ---
tab_contact <- dplyr::bind_rows(
  cat_block(dat, contact_sewerwater, "Contact with environment",
            label_map = c("0"="Contact with sewer water: No","1"="Contact with sewer water: Yes",
                          "No"="Contact with sewer water: No","Yes"="Contact with sewer water: Yes")),
  cat_block(dat, contact_trash, "Contact with environment",
            label_map = c("0"="Contact with trash: No","1"="Contact with trash: Yes",
                          "No"="Contact with trash: No","Yes"="Contact with trash: Yes")),
  cat_block(dat, contact_floodwater, "Contact with environment",
            label_map = c("0"="Contact with flood water: No","1"="Contact with flood water: Yes",
                          "No"="Contact with flood water: No","Yes"="Contact with flood water: Yes"))
)

## Combine
tab1_desc <- dplyr::bind_rows(tab_demo, tab_animals, tab_env, tab_contact)

# Order domains for printing/export
tab1_desc <- tab1_desc %>%
  dplyr::mutate(
    domain = factor(
      domain,
      levels = c(
        "Demographic & socioeconomic",
        "Household animals",
        "Household & peridomestic environment",
        "Contact with environment"
      )
    )
  ) %>%
  dplyr::arrange(domain)

write_csv(tab1_desc, "outputs/tables/table1_descriptives.csv")

# write as a word table
ft <- flextable(tab1_desc) %>%
  autofit()

save_as_docx(
  "Table 1. Study population characteristics and univariable associations" = ft,
  path = "outputs/tables/table1_descriptives.docx"
)

## ------------------------------------------------------
## 3. Univariable GAM smooths (visual exploration) ####
## ------------------------------------------------------

# Fit the models for each variable
g_age <- getViz(gam(data = dat, toxo_igg ~ s(age), family="binomial",seWithMean = TRUE))
g_road <- getViz(gam(data = dat, toxo_igg ~ s(dist_road), family="binomial",seWithMean = TRUE))
g_trash <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, dist_trash) %>% na.omit(), toxo_igg ~ s(dist_trash), family="binomial",seWithMean = TRUE))
g_sewer <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, dist_sewer) %>% na.omit(), toxo_igg ~ s(dist_sewer), family="binomial",seWithMean = TRUE))
g_elevation <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, elevation) %>% na.omit(), toxo_igg ~ s(elevation), family="binomial",seWithMean = TRUE))
g_income <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, income_pcap) %>% na.omit(), toxo_igg ~ s(income_pcap), family="binomial",seWithMean = TRUE))

# Extract the ggplot objects from the gamViz objects
p_age <- plot(sm(g_age, 1),seWithMean = TRUE) + labs(x="Age (years)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p_road <- plot(sm(g_road, 1),seWithMean = TRUE) + labs(x="Distance to the main road (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p_trash <- plot(sm(g_trash, 1),seWithMean = TRUE) + labs(x="Distance to nearest trash dump (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p_sewer <- plot(sm(g_sewer, 1),seWithMean = TRUE) + labs(x="Distance to nearest sewer (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p_elevation <- plot(sm(g_elevation,1),seWithMean = TRUE) + labs(x="Household elevation (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p_income <- plot(sm(g_income, 1),seWithMean = TRUE) + labs(x="Per capita daily household income (USD)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))

# save plots
save_mgcviz_tiff <- function(p, filename, width = 8, height = 4, dpi = 500) {
  grDevices::tiff(filename, units = "in", width = width, height = height,
                 res = dpi)
  print(p)
  grDevices::dev.off()
}

save_mgcviz_tiff(p_age,       "outputs/figures/gam_age.tiff")
save_mgcviz_tiff(p_road,      "outputs/figures/gam_road.tiff")
save_mgcviz_tiff(p_trash,     "outputs/figures/gam_trash.tiff")
save_mgcviz_tiff(p_sewer,     "outputs/figures/gam_sewer.tiff")
save_mgcviz_tiff(p_elevation, "outputs/figures/gam_elevation.tiff")
save_mgcviz_tiff(p_income,    "outputs/figures/gam_income.tiff")

p_age_i   <- ggdraw() + draw_image("outputs/figures/gam_age.tiff")
p_road_i  <- ggdraw() + draw_image("outputs/figures/gam_road.tiff")
p_trash_i <- ggdraw() + draw_image("outputs/figures/gam_trash.tiff")
p_sewer_i <- ggdraw() + draw_image("outputs/figures/gam_sewer.tiff")
p_elev_i  <- ggdraw() + draw_image("outputs/figures/gam_elevation.tiff")
p_inc_i   <- ggdraw() + draw_image("outputs/figures/gam_income.tiff")

top_row <- cowplot::plot_grid(p_age_i, p_road_i, ncol = 2, labels = c("A","B"),
                              label_size = 40, label_y = 1.01)
mid_row <- cowplot::plot_grid(p_trash_i, p_sewer_i, ncol = 2, labels = c("C","D"),
                              label_size = 40, label_y = 1.01)
bot_row <- cowplot::plot_grid(p_elev_i, p_inc_i, ncol = 2, labels = c("E","F"),
                              label_size = 40, label_y = 1.01)

fig_gam <- cowplot::plot_grid(top_row, mid_row, bot_row, ncol = 1, rel_heights = c(1,1,1))

ggsave("outputs/figures/FigS2_gam_all.tiff", fig_gam, units = "mm", width = 320, height = 320, dpi = 300)

## ------------------------------------------------------
## 4. Univariable logistic mixed models (random hh intercept) ####
## ------------------------------------------------------

# Robust tidy extractor for OR + Wald CI

extract_or_table <- function(model, model_label) {
  # Fixed effects summary
  s <- summary(model)$coefficients
  s <- as.data.frame(s)
  s$term <- rownames(s)
  
  # Drop intercept
  s <- s[s$term != "(Intercept)", , drop = FALSE]
  
  # Wald CI from confint
  ci <- suppressMessages(suppressWarnings(confint(model, method = "Wald")))
  ci <- as.data.frame(ci)
  ci$term <- rownames(ci)
  ci <- ci[ci$term != "(Intercept)", , drop = FALSE]
  
  out <- dplyr::left_join(
    s %>% dplyr::select(term,
                        estimate = Estimate,
                        std.error = `Std. Error`,
                        statistic = `z value`,
                        p.value = `Pr(>|z|)`),
    ci %>% dplyr::rename(conf.low = `2.5 %`, conf.high = `97.5 %`),
    by = "term"
  ) %>%
    dplyr::mutate(
      OR = exp(estimate),
      conf.low = exp(conf.low),
      conf.high = exp(conf.high),
      psig = dplyr::case_when(
        p.value < 0.001 ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      OR_CI = paste0(
        sprintf("%.2f", OR), " (",
        sprintf("%.2f", conf.low), ", ",
        sprintf("%.2f", conf.high), ")"
      ),
      model = model_label
    ) %>%
    dplyr::select(term, model, OR, conf.low, conf.high, p.value, psig, OR_CI,
                  estimate, std.error, statistic)
  
  tibble::as_tibble(out)
}

# Fit a list of univariable models

# 1 parameter models
one_param_terms <- c(
  "age", "sex", "income_pcap", "house_rented", "house_title", "I(elevation/10)",
  "I(dist_trash/10)", "I(dist_sewer/10)", "I(dist_road/50)", "hh_floods", "veg",
  "cat", "dog", "chicken", "rats_observed", "contact_trash", "contact_floodwater", "contact_sewerwater"
)

uni_1p <- purrr::map_dfr(one_param_terms, function(var) {
  f <- as.formula(paste0("toxo_igg ~ ", var, " + (1|hh_id)"))
  m <- lme4::glmer(
    f, data = dat, family = binomial, nAGQ = 10,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  )
  extract_or_table(m, model_label = var)
})

# multi-parameter models
m_race <- glmer(
  toxo_igg ~ as.factor(race) + (1|hh_id),
  data = dat, family = binomial, nAGQ = 10,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

m_agegroup <- glmer(
  toxo_igg ~ agegroup + (1|hh_id),
  data = dat, family = binomial, nAGQ = 10,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

uni_race <- extract_or_table(m_race, model_label = "race")
uni_ageg <- extract_or_table(m_agegroup, model_label = "agegroup")

univar_table <- dplyr::bind_rows(uni_1p, uni_race, uni_ageg)

write_csv(univar_table, "outputs/tables/table1_univar.csv")

# write as a word table
ft_uni <- flextable(univar_table) %>%
  autofit()

save_as_docx(
  "Table 1. Study population characteristics and univariable associations" = ft_uni,
  path = "outputs/tables/table1_univar.docx"
)


## ------------------------------------------------------
## 5. DAG-informed multivariable regression models ####
## ------------------------------------------------------

## Collinearity check (VIF) 
glm_model <- glm(
  toxo_igg ~ sex + agegroup + race +
    scale(income_pcap) + scale(elevation) +
    I(dist_trash/10) + I(dist_road/10) + I(dist_sewer/10) +
    hh_floods + veg + cat + dog + chicken + rats_observed +
    contact_trash + contact_floodwater + contact_sewerwater,
  data = dat,
  family = binomial
)

vif_values <- car::vif(glm_model)
print(vif_values)
print(vif_values[vif_values > 5]) # Check if any VIF values exceed the common threshold of 5


## Helper functions

# E-value for a specific fixed-effect term (Wald CI)
get_evalue <- function(model, term, rare = FALSE) {
  cf <- tryCatch(
    suppressMessages(confint(model, parm = term, method = "Wald")),
    error = function(e) NA
  )
  if (any(is.na(cf))) return(NA_real_)
  
  or_est <- exp(lme4::fixef(model)[term])
  or_lo  <- exp(cf[1])
  or_hi  <- exp(cf[2])
  
  # returns a matrix; take point-estimate E-value
  EValue::evalues.OR(or_est, lo = or_lo, hi = or_hi, rare = rare)[2, 1]
}

# Fit glmer with consistent controls
fit_glmer <- function(formula, data) {
  lme4::glmer(
    formula, data = data, family = binomial, nAGQ = 10,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  )
}

# Extract OR table (drops intercept, returns consistent schema)
extract_or_table <- function(model, model_id, domain, pretty_map = NULL, keep_terms = NULL) {
  
  tt <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      OR = exp(estimate),
      conf.low = exp(conf.low),
      conf.high = exp(conf.high),
      OR_CI = paste0(sprintf("%.2f", OR), " (", sprintf("%.2f", conf.low), ", ", sprintf("%.2f", conf.high), ")"),
      psig = dplyr::case_when(
        p.value < 0.001 ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      evalue = purrr::map_dbl(term, ~ get_evalue(model, .x, rare = FALSE)),
      model_id = model_id,
      domain = domain
    )
  
  # map to nice labels
  tt <- if (!is.null(pretty_map)) {
    tt %>% dplyr::mutate(variable = dplyr::recode(term, !!!pretty_map, .default = term))
  } else {
    tt %>% dplyr::mutate(variable = term)
  }
  
  # keep only terms of interest (if supplied)
  if (!is.null(keep_terms)) {
    tt <- tt %>% dplyr::filter(term %in% keep_terms)
  }
  
  tt %>%
    dplyr::select(domain, model_id, term, variable,
                  estimate, std.error, statistic, p.value,
                  OR, conf.low, conf.high, OR_CI, psig, evalue)
}


## Model specifications (one per “exposure of interest”)
## Each entry: list(domain, model_id, formula, data_filter, pretty_map)

pretty_common <- c(
  "agegroup7-9"   = "Age: 7-9 vs. 4-6",
  "agegroup10-12" = "Age: 10-12 vs. 4-6",
  "agegroup13-15" = "Age: 13-15 vs. 4-6",
  "agegroup16-18" = "Age: 16-18 vs. 4-6",
  "sexMale"       = "Sex: Male vs. Female",
  "raceBlack" = "Race: Black vs. Pardo",
  "raceWhite" = "Race: White vs. Pardo",
  "raceOther" = "Race: Other vs. Pardo",
  "income_pcap"   = "Per-capita household income (per US$)",
  "I(elevation/10)"   = "Household elevation (per 10 m)",
  "I(dist_road/50)"   = "Distance to main road (per 50 m)",
  "I(dist_trash/10)"  = "Distance to trash dump (per 10 m)",
  "I(dist_sewer/10)"  = "Distance to open sewer (per 10 m)",
  "hh_floodsYes"      = "House flooded in last 6 months",
  "vegYes"            = "Vegetation within 10 m",
  "catYes"            = "Cat in household",
  "dogYes"            = "Dog in household",
  "chickenYes"        = "Raise chickens",
  "rats_observedYes"  = "Rats observed in or near house",
  "contact_trashYes"        = "Contact with trash",
  "contact_floodwaterYes"   = "Contact with flood water",
  "contact_sewerwaterYes"   = "Contact with sewer water"
)

models_spec <- list(
  
  # Demographic & socioeconomic:
  list(
    domain = "Demographic & socioeconomic",
    model_id = "m0_age_sex_race",
    formula = toxo_igg ~ sex + agegroup + race + (1|hh_id),
    data = dat,
    pretty_map = pretty_common,
    keep_terms = c(
      "agegroup7-9","agegroup10-12","agegroup13-15","agegroup16-18",
      "sexMale",
      "raceBlack","raceWhite","raceOther"
    )
  ),
  
  list(
    domain = "Demographic & socioeconomic",
    model_id = "m1_income",
    formula = toxo_igg ~ income_pcap + agegroup + race + (1|hh_id),
    data = dat,
    pretty_map = pretty_common,
    keep_terms = c("income_pcap")
  ),
  
  # Household & peridomestic environment
  list(
    domain = "Household & peridomestic environment",
    model_id = "m2_elevation",
    formula = toxo_igg ~ I(elevation/10) + agegroup + race + income_pcap + (1|hh_id),
    data = dat,
    pretty_map = pretty_common,
    keep_terms = c("I(elevation/10)")
  ),
  list(
    domain = "Household & peridomestic environment",
    model_id = "m3_road",
    formula = toxo_igg ~ I(dist_road/50) + agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat,
    pretty_map = pretty_common,
    keep_terms = c("I(dist_road/50)")
  ),
  list(
    domain = "Household & peridomestic environment",
    model_id = "m4_trash",
    formula = toxo_igg ~ I(dist_trash/10) + scale(dist_road) + agegroup + race +
      scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(dist_trash)),
    pretty_map = pretty_common,
    keep_terms = c("I(dist_trash/10)")
  ),
  list(
    domain = "Household & peridomestic environment",
    model_id = "m5_sewer",
    formula = toxo_igg ~ I(dist_sewer/10) + agegroup + race +
      scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(dist_sewer)),
    pretty_map = pretty_common,
    keep_terms = c("I(dist_sewer/10)")
  ),
  list(
    domain = "Household & peridomestic environment",
    model_id = "m6_flooding",
    formula = toxo_igg ~ hh_floods + veg + scale(dist_sewer) + scale(dist_road) +
      agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(hh_floods)),
    pretty_map = pretty_common,
    keep_terms = c("hh_floodsYes")
  ),
  list(
    domain = "Household & peridomestic environment",
    model_id = "m7_veg",
    formula = toxo_igg ~ veg + scale(dist_road) + agegroup + race +
      scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(veg)),
    pretty_map = pretty_common,
    keep_terms = c("vegYes")
  ),
  
  # Household animals
  list(
    domain = "Household animals",
    model_id = "m8_cat",
    formula = toxo_igg ~ cat + agegroup + race + scale(income_pcap) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(cat)),
    pretty_map = pretty_common,
    keep_terms = c("catYes")
  ),
  list(
    domain = "Household animals",
    model_id = "m9_dog",
    formula = toxo_igg ~ dog + agegroup + race + scale(income_pcap) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(dog)),
    pretty_map = pretty_common,
    keep_terms = c("dogYes")
  ),
  list(
    domain = "Household animals",
    model_id = "m10_chicken",
    formula = toxo_igg ~ chicken + scale(dist_road) + agegroup + race +
      scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(chicken)),
    pretty_map = pretty_common,
    keep_terms = c("chickenYes")
  ),
  list(
    domain = "Household animals",
    model_id = "m11_rats",
    formula = toxo_igg ~ rats_observed + hh_floods + veg +
      scale(dist_sewer) + scale(dist_road) + scale(dist_trash) +
      agegroup + race + scale(income_pcap) + scale(elevation) +
      dog + chicken + cat + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(rats_observed)),
    pretty_map = pretty_common,
    keep_terms = c("rats_observedYes")
  ),
  
  # Contact with environment
  list(
    domain = "Contact with environment",
    model_id = "m12_contact_trash",
    formula = toxo_igg ~ contact_trash + scale(dist_trash) + scale(dist_road) +
      agegroup + sex + race + scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(contact_trash)),
    pretty_map = pretty_common,
    keep_terms = c("contact_trashYes")
  ),
  list(
    domain = "Contact with environment",
    model_id = "m13_contact_flood",
    formula = toxo_igg ~ contact_floodwater + hh_floods + veg + scale(dist_sewer) + scale(dist_road) +
      agegroup + sex + race + scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(contact_floodwater)),
    pretty_map = pretty_common,
    keep_terms = c("contact_floodwaterYes")
  ),
  list(
    domain = "Contact with environment",
    model_id = "m14_contact_sewer",
    formula = toxo_igg ~ contact_sewerwater + hh_floods + veg + scale(dist_sewer) + scale(dist_road) +
      agegroup + sex + race + scale(income_pcap) + scale(elevation) + (1|hh_id),
    data = dat %>% dplyr::filter(!is.na(contact_sewerwater)),
    pretty_map = pretty_common,
    keep_terms = c("contact_sewerwaterYes")
  )
)

## Fit all models + extract results
multivar_table <- purrr::map_dfr(models_spec, function(s) {
  m <- fit_glmer(s$formula, s$data)
  extract_or_table(
    m,
    model_id = s$model_id,
    domain = s$domain,
    pretty_map = s$pretty_map,
    keep_terms = s$keep_terms
  )
})

readr::write_csv(multivar_table, "outputs/tables/tableS2_S3_multivar_table.csv")


## G-computation (marginal predicted seroprevalence) for selected binary exposures
gcomp_binary <- function(model, data, var, yes_level = "Yes", no_level = "No") {
  # Prepare datasets
  dat_yes <- data %>% dplyr::mutate(!!rlang::ensym(var) := yes_level)
  dat_no  <- data %>% dplyr::mutate(!!rlang::ensym(var) := no_level)
  
  # Predict marginal probabilities
  p_yes <- predict(model, newdata = dat_yes, type = "response", re.form = NA)
  p_no  <- predict(model, newdata = dat_no,  type = "response", re.form = NA)
  
  tibble::tibble(
    variable = var,
    risk_yes = mean(p_yes, na.rm = TRUE),
    risk_no  = mean(p_no,  na.rm = TRUE),
    risk_diff = mean(p_yes, na.rm = TRUE) - mean(p_no, na.rm = TRUE)
  )
}

# G-computation for the cat variable
m_cat <- fit_glmer(toxo_igg ~ cat + agegroup + race + scale(income_pcap) + (1|hh_id),
                   dat %>% filter(!is.na(cat)))
gcomp_cat <- gcomp_binary(m_cat, dat %>% filter(!is.na(cat)), var = "cat")
print(gcomp_cat)


## Forest plot
domain_cols <- c(
  "Demographic & socioeconomic" = "#9e2a2b",
  "Household animals" = "#457b9d",
  "Household & peridomestic environment" = "#2a9d8f",
  "Contact with environment" = "#e76f51"
)

# Choose the terms you actually want in the forest plot (drop nuisance adjustment terms)
# Here we keep one “headline” exposure per model + key demographic terms
keep_terms <- c(
  # agegroup contrasts
  "agegroup7-9","agegroup10-12","agegroup13-15","agegroup16-18",
  # sex + race contrasts
  "sexMale","raceBlack","raceWhite",
  # income
  "income_pcap",
  # environment exposures
  "I(elevation/10)","I(dist_road/50)","I(dist_trash/10)","I(dist_sewer/10)","hh_floodsYes","vegYes",
  # animals
  "catYes","dogYes","chickenYes","rats_observedYes",
  # contact
  "contact_trashYes","contact_floodwaterYes","contact_sewerwaterYes"
)

forest_df <- multivar_table %>%
  dplyr::filter(term %in% keep_terms) %>%
  dplyr::mutate(
    domain = factor(domain, levels = c(
      "Demographic & socioeconomic",
      "Household animals",
      "Household & peridomestic environment",
      "Contact with environment"
    ))
  )

# Order variables (using the pretty label you already created)
variable_levels <- c(
  "Age: 7-9 vs. 4-6",
  "Age: 10-12 vs. 4-6",
  "Age: 13-15 vs. 4-6",
  "Age: 16-18 vs. 4-6",
  "Sex: Male vs. Female",
  "Race: Black vs. Pardo",
  "Race: White vs. Pardo",
  "Per-capita household income (per US$)",
  "Cat in household",
  "Dog in household",
  "Raise chickens",
  "Rats observed in or near house",
  "Household elevation (per 10 m)",
  "Distance to main road (per 50 m)",
  "House flooded in last 6 months",
  "Distance to trash dump (per 10 m)",
  "Distance to open sewer (per 10 m)",
  "Vegetation within 10 m",
  "Contact with sewer water",
  "Contact with trash",
  "Contact with flood water"
)

forest_df <- forest_df %>%
  dplyr::mutate(variable = factor(variable, levels = rev(variable_levels)))

fplot <- ggplot(forest_df, aes(x = OR, y = variable, colour = domain)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25, linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~domain, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = domain_cols, guide = "none") +
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20),
    limits = c(0.1, 40),
    labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10", "20")
  ) +
  labs(x = "Adjusted odds ratio (log scale)", y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.y = unit(0.6, "lines"),
    strip.background = element_rect(fill = "grey90", colour = "grey40", linewidth = 0.8),
    strip.text.x     = element_text(face = "bold", size = 15, margin = margin(t = 4, b = 4)),
    strip.placement  = "outside",
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.text.x = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 15, colour = "black")
  )

ggsave("outputs/figures/Fig4_multivar_forest_plot.tiff", fplot, width = 12, height = 13, dpi = 400)

