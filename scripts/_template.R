################################################################################
# PROJECT ACRONYM: PROJECT TITLE
################################################################################

# -------------------------------------------------------------------------
# Script: 01_setup_import_clean.R
#
# Purpose:
# Load required packages, define project settings, import raw data, 
# clean and harmonise variables, and generate a single analysis-ready 
# dataset used by all downstream scripts.
#
# -------------------------------------------------------------------------

## ------------------------------------------------------
## 0. Packages & global options ####
## ------------------------------------------------------

proj_title <- "X"

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

dat_raw <- read_csv(
  "data/raw/...csv"
)

## ------------------------------------------------------
## 2. Variable cleaning & harmonisation ####
## ------------------------------------------------------

dat <- dat_raw %>%
  select()


## ------------------------------------------------------
## 3. Save analysis-ready dataset ####
## ------------------------------------------------------

write_csv(dat, paste0("data/derived/",proj_title,"cleaned_data.csv"))
saveRDS(dat, paste0("data/derived/",proj_title,"cleaned_data.rds"))

## ------------------------------------------------------
## 4. Create deidentified dataset ####
## ------------------------------------------------------

data_deid <- dat %>%
  select()

write_csv(data_deid, paste0("data/derived/",proj_title,"cleaned_data_deid.csv"))
saveRDS(data_deid, paste0("data/derived/",proj_title,"cleaned_data_deid.rds"))  
