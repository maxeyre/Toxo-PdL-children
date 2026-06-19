# Toxo-PdL-children

This repository contains the analysis code supporting the manuscript:

**“Social marginalisation, environmental degradation and Toxoplasma gondii exposure in urban informal settlements in Brazil”**

Citation: Eyre MT, Wang JY, Carneiro IdO, Reis RB, Wunder Jr EA, Júnior NN, et al. (2026) Social marginalisation, environmental degradation and Toxoplasma gondii exposure in urban informal settlements in Brazil. PLoS Negl Trop Dis 20(6): e0014453. https://doi.org/10.1371/journal.pntd.0014453

The study investigates demographic, socioeconomic, environmental, and behavioural risk factors for *T. gondii* IgG seropositivity among children and adolescents living in a marginalised urban community, using serological data and a combination of regression, serocatalytic, and spatial statistical methods.

---

## Repository structure

```text
├── data/
│   ├── raw/                # Original input data (not public)
│   ├── spatial/            # Study area shapefiles (not public)
│   └── derived/            # All derived analysis datasets (not public except deid)
│       └── deid/           # De-identified analysis datasets (public)
├── scripts/
│   ├── 01_setup_import_clean.R
│   ├── 02_seroprev_serocatalytic.R
│   ├── 03_female_piecewise_foi_sensitivity.R
│   ├── 04_descriptive_gams_regression.R
│   └── 05_spatial_analysis.R
├── outputs/                # Model outputs, figures, tables (not public)
│   ├── figures/
│   ├── tables/
│   └── models/
├── documents/              # Manuscripts / internal docs (not public)
└── README.md
```

---

## Script overview

### `01_setup_import_clean.R`
**Purpose**
- Import raw serological, household, and environmental data  
- Harmonise variable coding  
- Construct derived variables (e.g. age groups, distances, binary indicators)  
- Save cleaned individual-level datasets for downstream analyses  

**Outputs**
- `data/derived/Toxo2003_full_cleaned_data.rds`  
- `data/derived/Toxo2003_full_cleaned_data_deid.rds` (public-use version)

---

### `02_seroprev_serocatalytic.R`
**Purpose**
- Estimate age-specific seroprevalence  
- Fit catalytic (force-of-infection) models stratified by sex  
- Quantify uncertainty in seroprevalence and FOI estimates  

**Key methods**
- Bayesian serocatalytic modelling (JAGS)  
- Model-based and sampling uncertainty bands  

**Outputs**
- FOI estimates and uncertainty intervals  
- Seroprevalence curves and supporting tables/figures  

---

### `02_seroprev_serocatalytic_individual.R`
Same as 02_seroprev_serocatalytic.R but for individual-level age data.

---

### `03_female_piecewise_foi_sensitivity.R`
**Purpose**
- Conduct a sensitivity analysis for age-dependent exposure among female participants
- Fit an individual-level piecewise serocatalytic model
- Allow the force of infection (FOI) to differ before and after age 12
- Estimate pre- and post-12 FOI, annual infection probabilities, and posterior evidence for reduced exposure after age 12
- Generate model diagnostics and fit metrics 

**Key methods**
- Individual-level Bayesian serocatalytic model fitted using JAGS
- Piecewise cumulative hazard model with separate FOI parameters before and after age 12
- Posterior estimation of the difference and ratio between post-12 and pre-12 FOI
- Posterior probability that post-12 FOI is lower than pre-12 FOI
- Model assessment using WAIC, LOO, DIC, and MCMC diagnostics

**Outputs**
- Posterior summaries for pre- and post-12 FOI
- Annual infection probabilities before and after age 12
- Posterior probability that FOI declines after age 12
- MCMC trace plots and summary tables
- WAIC, LOO, and DIC model fit metrics
  
---

### `04_descriptive_gams_regression.R`
**Purpose**
- Generate descriptive Table 1 (overall and seropositive populations)
- Explore non-linear associations using univariable GAM smooths
- Fit univariable mixed-effects logistic regression models
- Fit DAG-informed multivariable mixed-effects models
- Compute E-values and generate forest plots  

**Key methods**
- Logistic mixed-effects models with household random intercepts  
- Generalised additive models (GAMs) for exploratory analysis  
- DAG-informed covariate adjustment  
- E-values for unmeasured confounding  

**Outputs**
- Descriptive tables (CSV / Word)  
- Univariable and multivariable regression tables  
- Forest plots and GAM figures  

---

### `05_spatial_analysis.R`
**Purpose**
- Select predictors for spatial prediction  
- Build a fine-resolution prediction grid  
- Fit binomial geostatistical models using MCML  
- Generate spatial predictions of seroprevalence  
- Estimate and map residual spatial effects \(S(x)\)  

**Key methods**
- Model-based geostatistics (PrevMap)  
- Intercept-only and covariate-adjusted spatial models  
- Marginal and joint spatial prediction  
- Exceedance probability mapping  

⚠️ **Note**  
Geostatistical analyses require household coordinates, which are **not included** in the public de-identified dataset. This script will therefore not run end-to-end using only the shared data.

---

## Variable codebook

### Outcome

| Variable | Description |
|--------|-------------|
| `toxo_igg` | *T. gondii* IgG serostatus (0 = negative, 1 = positive) |

---

### Demographic & socioeconomic

| Variable | Description |
|--------|-------------|
| `age` | Age in years (removed in de-identified data) |
| `agegroup` | Age categories: 4–6, 7–9, 10–12, 13–15, 16–18 |
| `sex` | Sex (Male / Female) |
| `race` | Self-reported race (Pardo, Black, White, Other\*) |
| `income_pcap` | Per-capita daily household income (USD) |
| `house_rented` | Household is rented (Yes/No) |
| `house_title` | Household has legal title (Yes/No) |

\* In the de-identified dataset, **“Other” race was combined with “White”** due to small cell counts.

---

### Household animals

| Variable | Description |
|--------|-------------|
| `cat` | Cat present in household (Yes/No) |
| `dog` | Dog present in household (Yes/No) |
| `chicken` | Chickens raised at household (Yes/No) |
| `rats_observed` | Rats observed in or near household in last 6 months(Yes/No) |

---

### Household & peridomestic environment

| Variable | Description |
|--------|-------------|
| `elevation` | Household elevation (m) |
| `dist_road` | Distance to nearest main road (m) |
| `dist_trash` | Distance to nearest trash dump (m) |
| `dist_sewer` | Distance to nearest open sewer (m) |
| `hh_floods` | Household flooded in last 6 months (Yes/No) |
| `veg` | Vegetation within 10 m of household (Yes/No) |

---

### Contact with environment

| Variable | Description |
|--------|-------------|
| `contact_trash` | Contact with trash near house in last 6 months (Yes/No) |
| `contact_floodwater` | Contact with flood water near house in last 6 months (Yes/No) |
| `contact_sewerwater` | Contact with sewer water near house in last 6 months (Yes/No) |

---

## De-identified public-use data

A de-identified dataset is provided for transparency and reuse:

**`data/derived/Toxo2003_full_cleaned_data_deid.rds`**

To protect participant confidentiality:
- Household **coordinates have been removed**
- Exact **age has been removed** (age groups retained)
- Race category **“Other” has been merged with “White”**

As a result:
- Descriptive, regression, GAM, and serocatalytic analyses can be reproduced
- Geostatistical analyses **cannot** be fully reproduced without restricted-access data

---

## Software and packages

Analyses were conducted in **R (≥ 4.2)** using key packages including:

- `tidyverse`
- `lme4`, `broom.mixed`
- `mgcv`, `mgcViz`
- `PrevMap`
- `sf`, `raster`, `terra`
- `cowplot`, `flextable`, `EValue`

See individual scripts for full dependency lists.

---

## Citation

If you use this code or data, please cite the associated manuscript (details to be updated upon publication).

---

## Contact

For questions about the analysis or data:

**Max Eyre**  
London School of Hygiene & Tropical Medicine  
📧 max.eyre@lshtm.ac.uk
