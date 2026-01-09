################################################################################
# TOXO-PdL-children
################################################################################

# -------------------------------------------------------------------------
# Script: 01_setup_import_clean.R
#
# Purpose:
# Load required packages, define project settings, import raw toxoplasmosis
# serosurvey data, clean and harmonise variables, create household IDs, and
# generate a single analysis-ready dataset used by all downstream scripts.
#
# -------------------------------------------------------------------------

## ------------------------------------------------------
## 0. Packages & global options ####
## ------------------------------------------------------

required_pkgs <- c(
  "tidyverse"
)

missing <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

set.seed(01092025)
options(stringsAsFactors = FALSE)

## ------------------------------------------------------
## 1. Import raw data ####
## ------------------------------------------------------

dat_raw <- readr::read_csv(
  "data/raw/toxoplasma database 728 - binomial RES16.csv"
)

## ------------------------------------------------------
## 2. Variable cleaning & harmonisation ####
## ------------------------------------------------------

dat_clean <- dat_raw %>%
  rename(
    toxo_igg   = RES16_bin,
    income_pcap = renduss,
    elevation  = height,
    dist_road  = NEAR_Main,
    dist_trash = trash,
    dist_sewer = waste,
    X          = este,
    Y          = norte,
    house_title = tit, 
    house_rented = alug
  ) %>%
  mutate(
    # Sex
    sex = factor(
      if_else(sexfim == 1, "Male", "Female"),
      levels = c("Female", "Male")
    ),
    
    # Age groups
    agegroup = factor(
      case_when(
        agegroup == 1 ~ "4-6",
        agegroup == 2 ~ "7-9",
        agegroup == 3 ~ "10-12",
        agegroup == 4 ~ "13-15",
        agegroup == 5 ~ "16-18"
      ),
      levels = c("4-6","7-9","10-12","13-15","16-18")
    ),
    
    # Race/ethnicity
    race = factor(
      case_when(
        racare == 0 ~ "White",
        racare == 1 ~ "Pardo",
        racare == 2 ~ "Black",
        TRUE        ~ "Other"
      ),
      levels = c("Pardo","White","Black","Other")
    ),
    
    # Household animals
    cat      = factor(if_else(gat  == 1, "Yes", "No")),
    dog      = factor(if_else(cach == 1, "Yes", "No")),
    chicken = factor(if_else(galc == 1, "Yes", "No")),
    
    # Environmental exposure
    contact_trash      = factor(if_else(clxd   == 1, "Yes", "No")),
    contact_floodwater = factor(if_else(caalad == 1, "Yes", "No")),
    contact_sewerwater = factor(if_else(caesgd == 1, "Yes", "No")),
    
    # Household flooding
    hh_floods = factor(if_else(alaca == 1, "Yes", "No")),
    
    # Vegetation
    veg  = factor(if_else(veg  == 1, "Yes", "No")),
    
    # Rodents
    rats_observed = factor(if_else(ratd == 1, "Yes", "No"))
  ) %>%
  dplyr::select(
    toxo_igg, age, agegroup, sex, race,
    income_pcap, elevation,
    dist_road, dist_trash, dist_sewer,
    hh_floods, veg, cat, dog, chicken,
    rats_observed,
    contact_trash, contact_floodwater, contact_sewerwater, house_title, 
    house_rented,
    X, Y, casano
  )

## ------------------------------------------------------
## 3. Household ID creation ####
## ------------------------------------------------------

hh_ref <- dat_clean %>%
  distinct(casano) %>%
  arrange(casano) %>%
  mutate(hh_id = row_number())

dat <- dat_clean %>%
  left_join(hh_ref, by = "casano") %>%
  dplyr::select(-casano)

## ------------------------------------------------------
## 4. Save analysis-ready dataset ####
## ------------------------------------------------------

write_csv(dat, "data/derived/Toxo2003_full_cleaned_data.csv")
saveRDS(dat, "data/derived/Toxo2003_full_cleaned_data.rds")

## ------------------------------------------------------
## 5. Create deidentified dataset ####
## ------------------------------------------------------

data_deid <- dat %>%
  dplyr::select(-age, -X, -Y) %>%
  mutate(race = factor(case_when(
    race == "Other" ~ "White/Other",
    race == "White" ~ "White/Other",
    TRUE ~ race
  ), levels = c("Pardo","White/Other","Black"))) # collapsed other and race categories due to small number of observations

write_csv(data_deid, "data/derived/Toxo2003_full_cleaned_data_deid.csv")
saveRDS(data_deid, "data/derived/Toxo2003_full_cleaned_data_deid.rds")
