# -------------------------------------------------------------------------
# Script accompanying the manuscript:
# "Looking beyond household cats: environmental deficiencies in public spaces 
#  are associated with Toxoplasma gondii exposure in urban Brazilian children"
#
# Short title:
# "Environmental deficiencies and Toxoplasma gondii exposure in children"
#
# Authors:
# Max T. Eyre, Joyce Y. Wang, Ianei de O. Carneiro, Renato B. Reis, Elsio A. Wunder Jr, 
# Nivison N. Júnior, Guilherme S. Ribeiro, Fabio N. Souza, Juliet O. Santana, 
# Ellie Delight, Ridalva D. M. Felzemburgh, Francisco S. Santana, Sharif Mohr, 
# Astrid X. T. O. Melendez, Adriano Queiroz, Andreia C. Santos, Meghan Owens, 
# Claudia Muñoz-Zanzi, Mitermayer G. Reis, Bruno Martorelli Di Genova, 
# Peter J. Diggle, Federico Costa, Albert I. Ko
#
# Description:
# This script contains the statistical analyses and geospatial modelling conducted 
# for the above manuscript. It includes data preparation, model fitting, and 
# production of figures and tables reported in the paper. All analyses were 
# performed in R.
#
# Notes:
# - De-identified data only are provided here to protect participant confidentiality. 
# - The code is provided to enable transparency and reproducibility of analyses. 
# - For questions regarding the code, please contact the corresponding author.
# -------------------------------------------------------------------------

# List of required packages
packages <- c(
  "terra",
  "plyr",
  "tidyverse",
  "PrevMap",
  "lme4",
  "MuMIn",
  "performance",
  "binom",
  "mgcv",
  "mgcViz",
  "cowplot",
  "car",
  "sf",
  "RColorBrewer",
  "rjags",
  "varhandle",
  "MCMCvis",
  "loo"
)

# Install missing packages
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

# Load functions
source("scripts/functions/variogram_functions.R")
source("scripts/functions/catalyticModelFunctions.R")

set.seed(01092025)

#### 1. Load data ####
dat <- read_csv("data/old/toxoplasma database 728 - binomial RES16.csv") %>%
  rename(house_title = tit, house_rented = alug, dist_road = NEAR_Main, dist_trash = trash, dist_sewer = waste, elevation = height, income_pcap = renduss, dog = cach,
         X = este, Y = norte) %>%
  mutate(
    race = relevel(as.factor(case_when(
      racare == 0 ~ "white",
      racare == 1 ~ "mixed",
      racare == 2 ~ "black",
      racare == 3 ~ "other",
    )),"mixed"),
    agegroup = factor(case_when(
      agegroup == 1 ~ "4-6",
      agegroup == 2 ~ "7-9",
      agegroup == 3 ~ "10-12",
      agegroup == 4 ~ "13-15",
      agegroup == 5 ~ "16-18",
    ),c("4-6","7-9","10-12","13-15","16-18")),
    sex = case_when(
      sexfim == 0 ~ "Female",
      sexfim == 1 ~ "Male"
    ),
    hh_floods = case_when(
      alaca == 0 ~ "No",
      alaca == 1 ~ "Yes"
    ),
    cat = case_when(
      gat == 0 ~ "No",
      gat == 1 ~ "Yes"
    ),
    chicken = case_when(
      galc == 0 ~ "No",
      galc == 1 ~ "Yes"
    ),
    contact_floodwater = case_when(
      caalad == 0 ~ "No",
      caalad == 1 ~ "Yes"
    ),
    contact_sewerwater = case_when(
      caesgd == 0 ~ "No",
      caesgd == 1 ~ "Yes"
    ),
    ratd = case_when(
      ratd == 0 ~ "No",
      ratd == 1 ~ "Yes"
    ),
    qratd = case_when(
      qratd == 0 ~ "No",
      qratd == 1 ~ "Yes"
    ),
    contact_trash = case_when(
      clxd == 0 ~ "No",
      clxd == 1 ~ "Yes"
    )
  )

dat <- dat %>% 
  select(X, Y, casano, toxo_igg = RES16_bin, age, agegroup, sex, race,
         income_pcap, elevation, house_title, house_rented, dist_road, 
         dist_trash, dist_sewer, hh_floods, veg, dog, cat, chicken, 
         contact_floodwater, contact_sewerwater, contact_trash, ratd)

# create new household ids
ref_table  <- dat %>% 
  distinct(casano) %>%             
  arrange(casano) %>%
  mutate(hh_id = row_number())

dat <- dat %>%
  left_join(ref_table, by = "casano") %>%
  select(-casano)

write.csv(dat, "data/cleaned/Toxo2003_full_cleaned_data.csv")


dat <- read.csv("data/cleaned/Toxo2003_full_cleaned_data.csv") %>%
  mutate(
    race = factor(race),
    race = relevel(race, ref = "mixed"),
    agegroup = factor(
      agegroup,
      levels = c("4-6", "7-9", "10-12", "13-15", "16-18")
    )
  )

# COME BACK AND DE-IDENTIFY !! (remove: age, X, Y)

#### 2. Check household clustering of cases #### 
# descriptive
cluster_check <- dat %>% group_by(hh_id) %>%
  summarise(n= n(), n_pos = sum(toxo_igg), case_per_hh = n_pos/n)
mean(cluster_check$n) # mean number of participants per hh
mean(cluster_check$n_pos) # mean number of positive participants per hh
mean(cluster_check$case_per_hh) # mean % positive participants per hh
mean(cluster_check$case_per_hh[cluster_check$n_pos>0]) # mean % positive participants per hh (positive hh only)
100*prop.table(table(cluster_check$n_pos)) # percentage of hhs with n cases

# get ICC
mod_icc <- glmer(toxo_igg ~ 1 + (1|hh_id), dat, family="binomial", nAGQ =10, 
              control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
icc(mod_icc)

#### 3. Seroprevalence ####
# 3a. Figure 3 - seroprevalence plot ----

results <- dat %>%
  group_by(sex, agegroup) %>%
  dplyr::summarise(
    n = n(),
    count = sum(toxo_igg),
    proportion = count / n,
    ci_lower = binom.confint(count, n, conf.level = 0.95, methods = "wilson")$lower,
    ci_upper = binom.confint(count, n, conf.level = 0.95, methods = "wilson")$upper,
    .groups = "drop")

# plotting
results$sex <- factor(results$sex, levels = rev(levels(factor(results$sex))))
results$agegroup<- factor(results$agegroup, levels = c("4-6","7-9","10-12","13-15","16-18"))

# DELETE
# seroprev_plot <- ggplot(data = results, aes(x = agegroup, y = 100*proportion, fill = sex)) +
#   geom_bar(stat = "identity", position = position_dodge(), alpha = 0.5, colour = "black") +
#   geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), width = 0.3, position = position_dodge(0.9)) +
#   #facet_wrap(~sex) +
#   labs(x = "Age (years)",
#        y = "Seroprevalence (%)") +
#   geom_text(aes(label = paste0(sprintf("%.1f", 100 * proportion), "%","\n","n=",count,"/", n), y = 100*ci_upper),
#             position = position_dodge(0.9), vjust = -.5, size = 4) +  # Adjust vjust and size as needed
#   ylim(c(0,100)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 0),
#         legend.title = element_blank(),
#         strip.text = element_blank(),
#         text = element_text(size = 16))
# ggsave(seroprev_plot, "figures/seroprev_plot.tiff", dpi=300, width = 10, height = 6)

# 3b. serocatalytic models ----
# Reading in 5-year age grouped data
seroDatGrouped <- read_csv("Data/toxoSeroDatGrouped.csv")

# rename to new naming convention
seroDatGrouped <- seroDatGrouped %>%
  dplyr::select(age_name, 
         midpoint = age,
         sex,
         total = n,
         seropositive = n_pos,
         mean = prev)

cis <- sapply(1:nrow(seroDatGrouped), function(i) 
{binom.test(seroDatGrouped$seropositive[i], seroDatGrouped$total[i])$conf.int})

seroDatGrouped$lower <- cis[1,]
seroDatGrouped$upper <- cis[2,]

seroDatGrouped_all <- seroDatGrouped

### define model catalytic model

#define model
jcode <- "model{ 
for (i in 1:length(n)){

n.pos[i] ~ dbin(seropos_est[i],n[i]) 

#catalytic model
seropos_est[i] = 1-exp(-lambda*age[i]) 

#calculate likelihood (used for WAIC)
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],n[i]) 

}
## Priors
lambda ~ dunif(0,0.5) 
}"

##### MALE 
seroDatGrouped_male <- seroDatGrouped_all %>% filter(sex=="Male")

fiveYearAgeTotals <- seroDatGrouped_male$total

# Running the model

## Number of model iterations
mcmc.length=10000

## Specify my data
jdat <- list(n.pos = seroDatGrouped_male$seropositive,
             age = seroDatGrouped_male$midpoint,
             n=seroDatGrouped_male$total)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=10, n.adapt=2000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","loglik"), n.iter=mcmc.length)

plot(jpos) # check convergence

summary(jpos)
MCMCsummary(jpos, round = 2) ## Check ESS and Rhat

#convert mcmc.list to a matrix
mcmcMatrix <- as.matrix(jpos)

# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

# Get point estimates
lambdaPointEst <- mcmcMatrix[,"lambda"] %>% quantile(probs=c(.5,.025,.975))

# calculate DIC 
dic.samples(jmod, n.iter = mcmc.length)

## extract just the log liklihood
logLik <- mcmcMatrix[,2:6]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)

# waic = 30.3, looic = 30.3 

## Get data in correct format for plotting
## Sample from mcmc chains for credible intervals
## Add binomial sampling uncertainty

midpoints <- seroDatGrouped_male$midpoint
ager=4:18

## 1. Sample from mcmc chain to get 95% credible intervals (model uncertainty)
numSamples = 1000
outDf <- matrix(NA, nrow=numSamples, ncol=1)

for (i in 1:numSamples ) {
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  outDf[i,] <- mcmcMatrix[randomNumber,][[1]]
}

lambdaEst <- outDf %>% quantile(probs=c(.5,.025,.975))
lambdaEst %>% round(3) %>% print

## Create a df with model uncertainty
df_mod=data.frame(midpoint=ager, 
                  mean=1-exp(-lambdaEst[1]*ager),
                  lower=1-exp(-lambdaEst[2]*ager),
                  upper=1-exp(-lambdaEst[3]*ager))


## 2. Binomial sample uncertainty - accounts for the sample size of the underlying data
randomlySampleMcmcChain <- mcmcRandomSamplerCat(1000,mcmcMatrix,midpoints,fiveYearAgeTotals)
ageQuantilesSamplingUncertainty <- ageQuantiles(randomlySampleMcmcChain)

## Create a df with sample uncertainty
df_sampling = data.frame(
  midpoint = seroDatGrouped$midpoint,
  mean = 1 - exp(-seroDatGrouped$midpoint*(lambdaEst[1])),
  upper = ageQuantilesSamplingUncertainty[,3],
  lower = ageQuantilesSamplingUncertainty[,2]
)

# get lambda and calculate annual attack rate
lambdaEst_male <- lambdaEst
AIP_male <- (1-exp(-lambdaEst_male))*100

## Save the dataframes for plotting (along with reverse catalytic model)
saveRDS(df_sampling,"outputs/samplingUncertaintyCat_5_male.rds")
saveRDS(df_mod, "outputs/modelUncertaintyCat_5_male.rds" )
saveRDS(lambdaEst_male, "outputs/lambdaEst_male.rds" )

##  Plot figures ####

## Read in dataframes for catalytic and reverse catalytic models (sampling and model uncertainty)
samplingUncertaintyCat_male <- readRDS("outputs/samplingUncertaintyCat_5_male.rds")
modelUncertaintyCat_male <- readRDS("outputs/modelUncertaintyCat_5_male.rds" )
lambdaEst_male <- as.numeric(readRDS("outputs/lambdaEst_male.rds" ))
AIP_male <- (1-exp(-lambdaEst_male))*100

## Read in 3-year age grouped data
seroDat_male <- seroDatGrouped_male ## seroprevalence binned into 3 year age groups (used for plotting)

## Catalytic plot
# Extend the data with interpolated values for x = 4 and x = 18

data <- samplingUncertaintyCat_male

# Define the range of x values including the points for extrapolation
full_midpoint_values <- c(4, data$midpoint, 18)

# Function to perform linear extrapolation
linear_extrapolate <- function(x, y, new_x) {
  fit <- lm(y ~ x)
  predict(fit, newdata = data.frame(x = new_x))
}

# Interpolated and extrapolated values
interpolated_mean <- approx(data$midpoint, data$mean, xout = full_midpoint_values, rule = 2)$y
interpolated_lower <- approx(data$midpoint, data$lower, xout = full_midpoint_values, rule = 2)$y
interpolated_upper <- approx(data$midpoint, data$upper, xout = full_midpoint_values, rule = 2)$y

# Extrapolate for the points outside the original range
extrapolated_mean_left <- linear_extrapolate(data$midpoint[1:2], data$mean[1:2], 4)
extrapolated_mean_right <- linear_extrapolate(data$midpoint[(nrow(data)-1):nrow(data)], data$mean[(nrow(data)-1):nrow(data)], 18)

extrapolated_lower_left <- linear_extrapolate(data$midpoint[1:2], data$lower[1:2], 4)
extrapolated_lower_right <- linear_extrapolate(data$midpoint[(nrow(data)-1):nrow(data)], data$lower[(nrow(data)-1):nrow(data)], 18)

extrapolated_upper_left <- linear_extrapolate(data$midpoint[1:2], data$upper[1:2], 4)
extrapolated_upper_right <- linear_extrapolate(data$midpoint[(nrow(data)-1):nrow(data)], data$upper[(nrow(data)-1):nrow(data)], 18)

# Replace the boundary values with extrapolated values
interpolated_mean[1] <- extrapolated_mean_left
interpolated_mean[length(interpolated_mean)] <- extrapolated_mean_right

interpolated_lower[1] <- extrapolated_lower_left
interpolated_lower[length(interpolated_lower)] <- extrapolated_lower_right

interpolated_upper[1] <- extrapolated_upper_left
interpolated_upper[length(interpolated_upper)] <- extrapolated_upper_right

# Create interpolated data frame
samplingUncertaintyCat_male <- data.frame(
  midpoint = full_midpoint_values,
  mean = interpolated_mean,
  lower = interpolated_lower,
  upper = interpolated_upper
)

lambdaEst_male <- round(lambdaEst_male, 3)
lambdaEst_male<- format(lambdaEst_male, nsmall = 2) # to keep trailing zeroes

AIP_male <- round(AIP_male, 2)
AIP_male<- format(AIP_male, nsmall = 2) # to keep trailing zeroes

label_male <- paste0("FOI = ", lambdaEst_male[1], " (95% CrI: ", lambdaEst_male[2], ", ",lambdaEst_male[3],
                     ") \nAIP = ", AIP_male[1], "% (95% CrI: ", AIP_male[2], "%, ",AIP_male[3],"%)")

catPlot_male <- ggplot(modelUncertaintyCat_male, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.3, fill = "#457b9d")+
  geom_line()+
  geom_point(data=seroDat_male)+
  geom_linerange(data=seroDat_male) +
  geom_ribbon(data=samplingUncertaintyCat_male, alpha=0.3,fill = "#457b9d")+
  scale_y_continuous(breaks=seq(0,1,by=0.2), lim=c(0,.9))+
  xlab("Age (years)") + ylab("Proportion seropositive") +
  scale_x_continuous(breaks=seq(4,18,by=2), lim= c(4,18)) +
  theme_bw() +
  annotate("text", x = 5, y = 0.815, label = label_male,
           size = 3.5, color = "black", hjust= 0, vjust=0)
catPlot_male


##### FEMALE #####
seroDatGrouped_female <- seroDatGrouped_all %>% filter(sex=="Female")

fiveYearAgeTotals <- seroDatGrouped_female$total


# Running the model ####

## Number of model iterations
mcmc.length=10000

## Specify my data
jdat <- list(n.pos = seroDatGrouped_female$seropositive,
             age = seroDatGrouped_female$midpoint,
             n=seroDatGrouped_female$total)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=10, n.adapt=2000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","loglik"), n.iter=mcmc.length)

plot(jpos) # check convergence

summary(jpos)
MCMCsummary(jpos, round = 2) ## Check ESS and Rhat

#convert mcmc.list to a matrix
mcmcMatrix <- as.matrix(jpos)

# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

# Get point estimates
lambdaPointEst <- mcmcMatrix[,"lambda"] %>% quantile(probs=c(.5,.025,.975))

# calculate DIC 
dic.samples(jmod, n.iter = mcmc.length)

## extract just the log liklihood
logLik <- mcmcMatrix[,2:6]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)


################################################################################
## Get data in correct format for plotting
## Sample from mcmc chains for credible intervals
## Add binomial sampling uncertainty
################################################################################

midpoints <- seroDatGrouped_female$midpoint
ager=4:18

## 1. Sample from mcmc chain to get 95% credible intervals (model uncertainty)
numSamples = 1000
outDf <- matrix(NA, nrow=numSamples, ncol=1)

for (i in 1:numSamples ) {
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  outDf[i,] <- mcmcMatrix[randomNumber,][[1]]
}

lambdaEst <- outDf %>% quantile(probs=c(.5,.025,.975))
lambdaEst %>% round(3) %>% print

## Create a df with model uncertainty
df_mod=data.frame(midpoint=ager, 
                  mean=1-exp(-lambdaEst[1]*ager),
                  lower=1-exp(-lambdaEst[2]*ager),
                  upper=1-exp(-lambdaEst[3]*ager))


## 2. Binomial sample uncertainty - accounts for the sample size of the underlying data
randomlySampleMcmcChain <- mcmcRandomSamplerCat(1000,mcmcMatrix,midpoints,fiveYearAgeTotals)
ageQuantilesSamplingUncertainty <- ageQuantiles(randomlySampleMcmcChain)

## Create a df with sample uncertainty
df_sampling = data.frame(
  midpoint = seroDatGrouped_female$midpoint,
  mean = 1 - exp(-seroDatGrouped_female$midpoint*(lambdaEst[1])),
  upper = ageQuantilesSamplingUncertainty[,3],
  lower = ageQuantilesSamplingUncertainty[,2]
)

# get lambda and calculate annual attack rate
lambdaEst_female <- lambdaEst
AIP_female <- (1-exp(-lambdaEst))*100

## Save the dataframes for plotting (along with reverse catalytic model)
saveRDS(df_sampling,"outputs/samplingUncertaintyCat_5_female.rds")
saveRDS(df_mod, "outputs/modelUncertaintyCat_5_female.rds" )
saveRDS(lambdaEst_female, "outputs/lambdaEst_female.rds" )

### Plot Figures ####

## Read in dataframes for catalytic and reverse catalytic models (sampling and model uncertainty)
samplingUncertaintyCat_female <- readRDS("outputs/samplingUncertaintyCat_5_female.rds")
modelUncertaintyCat_female <- readRDS("outputs/modelUncertaintyCat_5_female.rds" )
lambdaEst_female <- readRDS("outputs/lambdaEst_female.rds" )
AIP_female <- (1-exp(-lambdaEst_female))*100

## Read in 5 year age grouped data
seroDat_female <- seroDatGrouped_female ## seroprevalence binned into 5 year age groups (used for plotting)

## Catalytic plot
# Extend the data with interpolated values for x = 4 and x = 18

data <- samplingUncertaintyCat_female

# Define the range of x values including the points for extrapolation
full_midpoint_values <- c(4, data$midpoint, 18)

# Function to perform linear extrapolation
linear_extrapolate <- function(x, y, new_x) {
  fit <- lm(y ~ x)
  predict(fit, newdata = data.frame(x = new_x))
}

# Interpolated and extrapolated values
interpolated_mean <- approx(data$midpoint, data$mean, xout = full_midpoint_values, rule = 2)$y
interpolated_lower <- approx(data$midpoint, data$lower, xout = full_midpoint_values, rule = 2)$y
interpolated_upper <- approx(data$midpoint, data$upper, xout = full_midpoint_values, rule = 2)$y

# Extrapolate for the points outside the original range
extrapolated_mean_left <- linear_extrapolate(data$midpoint[1:2], data$mean[1:2], 4)
extrapolated_mean_right <- linear_extrapolate(data$midpoint[(nrow(data)-1):nrow(data)], data$mean[(nrow(data)-1):nrow(data)], 18)

extrapolated_lower_left <- linear_extrapolate(data$midpoint[1:2], data$lower[1:2], 4)
extrapolated_lower_right <- linear_extrapolate(data$midpoint[(nrow(data)-1):nrow(data)], data$lower[(nrow(data)-1):nrow(data)], 18)

extrapolated_upper_left <- linear_extrapolate(data$midpoint[1:2], data$upper[1:2], 4)
extrapolated_upper_right <- linear_extrapolate(data$midpoint[(nrow(data)-1):nrow(data)], data$upper[(nrow(data)-1):nrow(data)], 18)

# Replace the boundary values with extrapolated values
interpolated_mean[1] <- extrapolated_mean_left
interpolated_mean[length(interpolated_mean)] <- extrapolated_mean_right

interpolated_lower[1] <- extrapolated_lower_left
interpolated_lower[length(interpolated_lower)] <- extrapolated_lower_right

interpolated_upper[1] <- extrapolated_upper_left
interpolated_upper[length(interpolated_upper)] <- extrapolated_upper_right

# Create interpolated data frame
samplingUncertaintyCat_female <- data.frame(
  midpoint = full_midpoint_values,
  mean = interpolated_mean,
  lower = interpolated_lower,
  upper = interpolated_upper
)

modelUncertaintyCat_female <- modelUncertaintyCat

lambdaEst_female <- round(lambdaEst_female, 3)
lambdaEst_female<- format(lambdaEst_female, nsmall = 2) # to keep trailing zeroes

AIP_female <- round(AIP_female, 2)
AIP_female<- format(AIP_female, nsmall = 2) # to keep trailing zeroes

label_female <- paste0("FOI = ", lambdaEst_female[1], " (95% CrI: ", lambdaEst_female[2], ", ",lambdaEst_female[3],
                       ") \nAIP = ", AIP_female[1], "% (95% CrI: ", AIP_female[2], "%, ",AIP_female[3],"%)")

# catPlot_female <- ggplot(modelUncertaintyCat_female, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
#   geom_ribbon(alpha=0.3, fill = "#9e2a2b")+
#   geom_line()+
#   geom_point(data=seroDat_female)+
#   geom_linerange(data=seroDat_female) +
#   geom_ribbon(data=samplingUncertaintyCat_female, alpha=0.3,fill = "#9e2a2b")+
#   scale_y_continuous(breaks=seq(0,1,by=0.2), lim=c(0,.9))+
#   xlab("Age (years)") + ylab("Proportion seropositive") +
#   scale_x_continuous(breaks=seq(4,18,by=2), lim= c(4,18)) +
#   theme_bw() + 
#   annotate("text", x = 5, y = 0.815, label = label_female,
#            size = 3.5, color = "black", hjust= 0, vjust=0)
# catPlot_female
# 
# catPlot_male <- ggplot(modelUncertaintyCat_male, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
#   geom_ribbon(alpha=0.3, fill = "#457b9d")+
#   geom_line()+
#   geom_point(data=seroDat_male)+
#   geom_linerange(data=seroDat_male) +
#   geom_ribbon(data=samplingUncertaintyCat_male, alpha=0.3,fill = "#457b9d")+
#   scale_y_continuous(breaks=seq(0,1,by=0.2), lim=c(0,.9))+
#   xlab("Age (years)") + ylab("Proportion seropositive") +
#   scale_x_continuous(breaks=seq(4,18,by=2), lim= c(4,18)) +
#   theme_bw() +
#   annotate("text", x = 5, y = 0.815, label = label_male,
#            size = 3.5, color = "black", hjust= 0, vjust=0)
# catPlot_male

# seroprev_plot <- ggplot(data = results, aes(x = agegroup, y = 100*proportion, fill = sex)) +
#   geom_bar(stat = "identity", position = position_dodge(), alpha = 0.5, colour = "black") +
#   geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), width = 0.3, position = position_dodge(0.9)) +
#   #facet_wrap(~sex) +
#   labs(x = "Age (years)",
#        y = "Seroprevalence (%)") +
#   geom_text(aes(label = paste0(sprintf("%.1f", 100 * proportion), "%","\n","n=",count,"/", n), y = 100*ci_upper),
#             position = position_dodge(0.9), vjust = -.5, size = 4) +  # Adjust vjust and size as needed
#   ylim(c(0,100)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 0),
#         legend.title = element_blank(),
#         strip.text = element_blank(),
#         text = element_text(size = 16))
# ggsave(seroprev_plot, "figures/seroprev_plot.tiff", dpi=300, width = 10, height = 6)
# 
# ## Create panel plot
# # top row: just the seroprev_plot, full width
# top_row <- plot_grid(seroprev_plot, labels = "A", ncol = 1)
# 
# # bottom row: male and female plots side by side
# bottom_row <- plot_grid(catPlot_male, catPlot_female,
#                         labels = c("B", "C"),
#                         ncol = 2)
# 
# # combine: top row stacked above bottom row
# catPlot <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
# catPlot
# 
# catPlot <- plot_grid(seroprev_plot, NULL, catPlot_male, catPlot_female, labels = c('A',NULL, 'B','C', rows = 2))
# catPlot
# # save plot
# ggsave(plot = catPlot, "outputs/serocatalytic_plot_label.png", width = 10, height = 5, dpi = 300)

####
# Define consistent colours
col_male   <-  "#2c7bb6" # blue
col_female <-  "#d7191c" # red

# Female plot
catPlot_female <- ggplot(modelUncertaintyCat_female, 
                         aes(x=midpoint, y=100*mean, ymin=100*lower, ymax=100*upper)) +
  geom_ribbon(alpha=0.3, fill = col_female) +
  geom_line(color = "black") +
  geom_point(data=seroDat_female, color = "black") +
  geom_linerange(data=seroDat_female, color = "black") +
  geom_ribbon(data=samplingUncertaintyCat_female, alpha=0.3, fill = col_female) +
  scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
  scale_x_continuous(breaks=seq(4,18,by=2), limits=c(4,18)) +
  labs(x = "Age (years)", y = "") +
  annotate("text", x = 4.3, y = 100*0.815, label = label_female,
           size = 4, color = "black", hjust = 0, vjust = 0) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.title.x = element_text(size = 16))

# Male plot
catPlot_male <- ggplot(modelUncertaintyCat_male, 
                       aes(x=midpoint, y=100*mean, ymin=100*lower, ymax=100*upper)) +
  geom_ribbon(alpha=0.3, fill = col_male) +
  geom_line(color = "black") +
  geom_point(data=seroDat_male, color = "black") +
  geom_linerange(data=seroDat_male, color = "black") +
  geom_ribbon(data=samplingUncertaintyCat_male, alpha=0.3, fill = col_male) +
  scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
  scale_x_continuous(breaks=seq(4,18,by=2), limits=c(4,18)) +
  labs(x = "Age (years)", y = "Seroprevalence (%)") +
  annotate("text", x = 4.3, y = 100*0.815, label = label_male,
           size = 4, color = "black", hjust = 0, vjust = 0) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.title.x = element_text(size = 16))

# Seroprevalence plot
seroprev_plot <- ggplot(data = results, aes(x = agegroup, y = 100*proportion, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.5, colour = "black") +
  geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), 
                width = 0.3, position = position_dodge(0.9)) +
  geom_text(aes(label = paste0(sprintf("%.1f", 100 * proportion), "%","\n","n=",count,"/", n),
                y = 100*ci_upper),
            position = position_dodge(0.9), vjust = -.5, size = 4) +
  scale_fill_manual(values = c("Male" = col_male, "Female" = col_female)) +
  labs(x = "Age (years)", y = "Seroprevalence (%)") +
  scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        legend.title = element_blank(),
        strip.text = element_blank(),
        text = element_text(size = 16))

# Combine into final figure
top_row <- plot_grid(seroprev_plot, labels = "A", ncol = 1, label_fontface = "bold", label_size = 25,label_x = -0.01)
bottom_row <- plot_grid(catPlot_male, catPlot_female, labels = c("B", "C"), ncol = 2, label_fontface = "bold", label_size = 25, label_x = -0.01)
catPlot <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
catPlot
# Save
ggsave(plot = catPlot, "outputs/serocatalytic_plot_label.tiff", width = 10, height = 10, dpi = 500)





#### 4. GAMs ####
# Fit the models for each variable
g0 <- getViz(gam(data = dat, toxo_igg ~ s(age), family="binomial",seWithMean = TRUE))
g1 <- getViz(gam(data = dat, toxo_igg ~ s(dist_road), family="binomial",seWithMean = TRUE))
g2 <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, dist_trash) %>% na.omit(), toxo_igg ~ s(dist_trash), family="binomial",seWithMean = TRUE))
g3 <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, dist_sewer) %>% na.omit(), toxo_igg ~ s(dist_sewer), family="binomial",seWithMean = TRUE))
g4 <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, elevation) %>% na.omit(), toxo_igg ~ s(elevation), family="binomial",seWithMean = TRUE))
g5 <- getViz(gam(data = dat %>% dplyr::select(toxo_igg, income_pcap) %>% na.omit(), toxo_igg ~ s(income_pcap), family="binomial",seWithMean = TRUE))

# Generate plots manually using ggplot for each variable

# Extract the ggplot objects from the gamViz objects
p0 <- plot(sm(g0, 1),seWithMean = TRUE) + labs(x="Age (years)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p1 <- plot(sm(g1, 1),seWithMean = TRUE) + labs(x="Distance to the main road (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p2 <- plot(sm(g2, 1),seWithMean = TRUE) + labs(x="Distance to nearest trash dump (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p3 <- plot(sm(g3, 1),seWithMean = TRUE) + labs(x="Distance to nearest sewer (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p4 <- plot(sm(g4, 1),seWithMean = TRUE) + labs(x="Household elevation (m)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))
p5 <- plot(sm(g5, 1),seWithMean = TRUE) + labs(x="Per capita daily household income (USD)", y="Log-odds of seropositivity risk") + theme(text = element_text(size = 16))

# Make joint figure
jpeg("figures/gam_age.tiff", units="in", width=8, height=4, res=300)
p0
dev.off()

jpeg("figures/gam_road.tiff", units="in", width=8, height=4, res=300)
p1
dev.off()

jpeg("figures/gam_trash.tiff", units="in", width=8, height=4, res=300)
p2
dev.off()

jpeg("figures/gam_sewer.tiff", units="in", width=8, height=4, res=300)
p3
dev.off()

jpeg("figures/gam_elev.tiff", units="in", width=8, height=4, res=300)
p4
dev.off()

jpeg("figures/gam_income.tiff", units="in", width=8, height=4, res=300)
p5
dev.off()

p0 <- ggdraw() + draw_image("figures/gam_age.tiff")
p1 <- ggdraw() + draw_image("figures/gam_road.tiff")
p2 <- ggdraw() + draw_image("figures/gam_trash.tiff")
p3 <- ggdraw() + draw_image("figures/gam_sewer.tiff")
p4 <- ggdraw() + draw_image("figures/gam_elev.tiff")
p5 <- ggdraw() + draw_image("figures/gam_income.tiff")

top_row <- cowplot::plot_grid(p0, p1, ncol=2, labels = c("A","B"), label_size = 40, label_y = 1.01)
middle_row <- cowplot::plot_grid(p2, p3, ncol=2, labels = c("C","D"), label_size = 40, label_y = 1.01)
bottom_row <- cowplot::plot_grid(p4, p5,ncol=2, labels = c("E","F"), label_size = 40, label_y = 1.01)

fig2 <- cowplot::plot_grid(top_row, middle_row, bottom_row, ncol=1, rel_heights=c(1,1,1))

jpeg("figures/gam_all.tiff", units="mm", width=320, height=320, res=300)
print(fig2)
dev.off()

#### 5. Univariable regression models ####

# 1 parameter variables
options(na.action = "na.omit")
OR.1p <- sapply(c(
  "age", "sex","income_pcap", "house_rented", "house_title",  "I(elevation/10)",
  "I(dist_trash/10)", "I(dist_sewer/10)", "I(dist_road/50)", "hh_floods","veg", "cat", "dog", "chicken",
  "ratd", "contact_trash", "contact_floodwater","contact_sewerwater"),
  
  function(var) {
    
    formula    <- as.formula(paste("toxo_igg ~ ", var, "+ (1|hh_id) "))
    res.logist <- glmer(formula, data = dat, family = binomial, nAGQ = 10, 
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    n <- nrow(coef(summary(res.logist)))
    x <- confint(res.logist,method = c("Wald"))
    
    out <- c(coef(summary(res.logist))[2:n,],x[(2+1):(n+1),1:2])
    return(out)
  })

# race (4 parameter values)
race <- glmer(toxo_igg ~ as.factor(race) + (1|hh_id), dat, family="binomial", nAGQ =10, 
              control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
x <- confint(race,method = c("Wald"))
out1 <- cbind(coef(summary(race))[2:4,],x[3:5,1:2])

# age group (5 parameter values)
age <- glmer(toxo_igg ~ agegroup + (1|hh_id), dat, family="binomial", nAGQ =10, 
              control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
x <- confint(age,method = c("Wald"))
out2 <- cbind(coef(summary(age))[2:5,],x[3:6,1:2])

# combine together
OR.full <- bind_rows(as.data.frame(t(OR.1p)),as.data.frame(out1),as.data.frame(out2)) %>%
  mutate(Estimate = round(exp(Estimate),2), 
         `2.5 %` = round(exp(`2.5 %`),2), 
         `97.5 %` = round(exp(`97.5 %`),2),
         variable = rownames(.)) %>%
  dplyr::select(variable, Estimate, `2.5 %`, `97.5 %`, pval = `Pr(>|z|)`) %>%
  mutate(psig = case_when(
    pval < 0.05 & pval >= 0.001~ "*",
    pval < 0.001 ~ "**",
    TRUE ~ ""
  ),
  OR = paste0(Estimate, " (",`2.5 %`,", ",`97.5 %`,")")) %>%
  write_csv("outputs/univar_table.csv")

#### 6. DAG-informed multivariable regression models ####
# check for collinearity:
# Fit a glm (without random effects) for VIF calculation
glm_model <- glm(toxo_igg ~ sex + agegroup + race + scale(income_pcap) + scale(elevation) + 
                 I(dist_trash/10) + I(dist_road/10) + I(dist_sewer/10) + hh_floods + veg + 
                 cat + dog + chicken + ratd + contact_trash + contact_floodwater + contact_sewerwater, 
               data = dat, family="binomial")

# Calculate VIF for the linear model
vif_values <- vif(glm_model)

# Print the VIF values
print(vif_values)

# Check if any VIF values exceed the common threshold of 5
high_vif <- vif_values[vif_values > 5]
print(high_vif)

# Function for getting e-values
get_evalue <- function(model, term, rare = FALSE) {
  cf <- tryCatch(
    confint(model, parm = term, method = "Wald"),
    error = function(e) NA
  )
  
  if (any(is.na(cf))) return(NA_real_)
  
  or_est <- exp(fixef(model)[term])
  or_lo  <- exp(cf[1])
  or_hi  <- exp(cf[2])
  
  # point-estimate E-value
  EValue::evalues.OR(or_est, lo = or_lo, hi = or_hi, rare = rare)[2, 1]
}


# DEMOGRAPHIC & SOCIOECONOMIC VARIABLES
# Age, sex and race
multivar_results <- tibble()

m0 <- glmer(toxo_igg ~ sex + agegroup + race + (1|hh_id), dat, family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m0
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)

multivar_results <- rbind(multivar_results, x[3:10,])

# SES (income)
m1 <- glmer(toxo_igg ~ income_pcap + agegroup + race + (1|hh_id), dat, family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m1
summary(mod)
exp(fixef(mod))

x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)

multivar_results <- rbind(multivar_results, x[3,])



# PUBLIC ENVIRONMENT
# Household elevation
m2 <- glmer(toxo_igg ~ I(elevation/10) + agegroup + race + income_pcap + (1|hh_id), dat, family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m2
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# Distance to road 
m3 <- glmer(toxo_igg ~ I(dist_road/50) + agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id), dat, family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m3
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# Distance to trash 
m4 <- glmer(toxo_igg ~ I(dist_trash/10) + scale(dist_road) + agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m4
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# Distance to open sewer 
m5 <- glmer(toxo_igg ~ I(dist_sewer/10) + agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m5
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# Household flooding
m6 <- glmer(toxo_igg ~ hh_floods + veg + scale(dist_sewer) + scale(dist_road) + agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m6
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# Vegetation
m7 <- glmer(toxo_igg ~ veg + scale(dist_road) + agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m7
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# HOUSEHOLD ENVIRONMENT
# cat ownership
m8 <- glmer(toxo_igg ~ cat + agegroup + race + scale(income_pcap) + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m8
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

## G-computation for marginal predicted seroprevalence

# Prepare datasets with cat = Yes and cat = No
dat_yes <- dat %>% mutate(cat = "Yes")
dat_no  <- dat %>% mutate(cat = "No")

# Predict *marginal* probabilities (re.form = NA)
p_yes <- predict(m8, newdata = dat_yes, type = "response", re.form = NA)
p_no  <- predict(m8, newdata = dat_no,  type = "response", re.form = NA)

risk_yes <- mean(p_yes)
risk_no  <- mean(p_no)

cat("Marginal risk (cat=Yes):", risk_yes, "\n")
cat("Marginal risk (cat=No): ", risk_no,  "\n")


# dog ownership
m9 <- glmer(toxo_igg ~ dog + agegroup + race + scale(income_pcap) + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m9
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# chicken ownership
m10 <- glmer(toxo_igg ~ chicken + scale(dist_road) + agegroup + race + scale(income_pcap) + scale(elevation) + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m10
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# rats observed
m11 <- glmer(toxo_igg ~ ratd + hh_floods + veg + scale(dist_sewer) + scale(dist_road) + scale(dist_trash) + agegroup + race + scale(income_pcap) + scale(elevation) + dog + chicken + cat + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m11
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# CONTACT WITH ENVIRONMENT
# contact w trash
m12 <- glmer(toxo_igg ~ contact_trash + scale(dist_trash) + scale(dist_road) + agegroup + sex +  race + scale(income_pcap) + scale(elevation)  + (1|hh_id), dat %>% filter(!is.na(dist_trash)), family="binomial", nAGQ =10, 
             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m12
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# contact w flood water
m13 <- glmer(toxo_igg ~ contact_floodwater + hh_floods + veg + scale(dist_sewer) + scale(dist_road) + agegroup + sex +  race + scale(income_pcap) + scale(elevation)  + (1|hh_id), dat , family="binomial", nAGQ =10, 
             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m13
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])

# contact w sewer water
m14 <- glmer(toxo_igg ~ contact_sewerwater + hh_floods + veg + scale(dist_sewer) + scale(dist_road) + agegroup + sex +  race + scale(income_pcap) + scale(elevation)  + (1|hh_id), dat, family="binomial", nAGQ =10, 
             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
mod <- m14
summary(mod)
exp(fixef(mod))
x <- tibble(
  parameter = names(exp(confint(mod, method="Wald")[,1])),
  est = round(c(NA, exp(fixef(mod))), 2),
  lci = round(exp(confint(mod, method="Wald")[,1]), 2),
  uci = round(exp(confint(mod, method="Wald")[,2]), 2),
  evalue = sapply(
    names(exp(confint(mod, method="Wald")[,1])),
    function(t) get_evalue(mod, t, rare = FALSE)
  )
)
multivar_results <- rbind(multivar_results, x[3,])


## G-computation for marginal predicted seroprevalence

# Prepare datasets with cat = Yes and cat = No
dat_yes <- dat %>% mutate(contact_sewerwater = "Yes")
dat_no  <- dat %>% mutate(contact_sewerwater = "No")

# Predict *marginal* probabilities (re.form = NA)
p_yes <- predict(m14, newdata = dat_yes, type = "response", re.form = NA)
p_no  <- predict(m14, newdata = dat_no,  type = "response", re.form = NA)

risk_yes <- mean(p_yes)
risk_no  <- mean(p_no)

cat("Marginal risk (contact_sewerwater=Yes):", risk_yes, "\n")
cat("Marginal risk (contact_sewerwater=No): ", risk_no,  "\n")

check <- multivar_results %>%
  mutate(
    OR = paste0(
      sprintf("%.2f", est),
      " (",
      sprintf("%.2f", lci),
      ", ",
      sprintf("%.2f", uci),
      ")"
    ),
    evalue = round(evalue, 2)
  ) %>%
  write_csv("outputs/multivar_table.csv")

# forest plot
domain_cols <- c(
  "Demographic & socioeconomic" = "#9e2a2b",
  "Household animals" = "#457b9d",
  "Household & peridomestic environment" = "#2a9d8f",
  "Contact with environment" = "#e76f51"
)


df <- data.frame(
  domain = c(
    rep("Demographic & socioeconomic", 8),
    rep("Household animals", 4),
    rep("Household & peridomestic environment", 6),
    rep("Contact with environment", 3)
  ),
  
  variable = c(
    # Age
    "Age: 7–9 vs. 4–6",
    "Age: 10–12 vs. 4–6",
    "Age: 13–15 vs. 4–6",
    "Age: 16–18 vs. 4–6",
    
    # Sex
    "Sex: Male vs. Female",
    
    # Race
    "Race: Black vs. Pardo",
    "Race: White vs. Pardo",

    # Income
    "Per-capita household income (per US$)",
    
    # Household
    "Cat in household",
    "Dog in household",
    "Raise chickens",
    "Rats observed in or near house",
    
    # Public environment
    "House flooded in last 6 months",
    "Household elevation (per 10 m)",
    "Distance to trash dump (per 10 m)",
    "Distance to open sewer (per 10 m)",
    "Distance to main road (per 50 m)",
    "Vegetation within 10 m",
    
    # Contact
    "Contact with trash",
    "Contact with flood water",
    "Contact with sewer water"
  ),
  
  aOR = c(
    3.26, 5.05, 12.96, 12.30,
    2.46,
    1.52, 0.33,
    0.54,
    1.93, 1.06, 1.28, 1.31,
    1.04, 0.66, 1.00, 0.91, 1.16, 0.85,
    1.28, 1.01, 2.54
  ),
  
  lower = c(
    1.58, 2.32, 5.58, 5.43,
    1.59,
    0.95, 0.11,
    0.38,
    1.08, 0.67, 0.82, 0.79,
    0.52, 0.55, 0.94, 0.79, 1.04, 0.49,
    0.76, 0.65, 1.50
  ),
  
  upper = c(
    6.74, 11.00, 30.08, 27.85,
    3.81,
    2.45, 1.02,
    0.78,
    3.44, 1.67, 1.99, 2.18,
    2.09, 0.80, 1.06, 1.05, 1.30, 1.47,
    2.15, 1.58, 4.33
  )
)


domain_levels <- c(
  "Demographic & socioeconomic",
  "Household animals",
  "Household & peridomestic environment",
  "Contact with environment"
)

variable_levels <- c(
  # Demographic & socioeconomic
  "Age: 7–9 vs. 4–6",
  "Age: 10–12 vs. 4–6",
  "Age: 13–15 vs. 4–6",
  "Age: 16–18 vs. 4–6",
  "Sex: Male vs. Female",
  "Race: Black vs. Pardo",
  "Race: White vs. Pardo",
  "Per-capita household income (per US$)",
  
  # Household animals
  "Cat in household",
  "Dog in household",
  "Raise chickens",
  "Rats observed in or near house",
  
  # Household & peridomestic environment
  "Household elevation (per 10 m)",
  "Distance to main road (per 50 m)",
  "House flooded in last 6 months",
  "Distance to trash dump (per 10 m)",
  "Distance to open sewer (per 10 m)",
  "Vegetation within 10 m",
  
  # Environmental contact
  "Contact with sewer water",
  "Contact with trash",
  "Contact with flood water"
)

df <- df %>%
  mutate(
    domain   = factor(domain, levels = domain_levels),
    variable = factor(variable, levels = rev(variable_levels))
  )

fplot <- ggplot(df, aes(x = aOR, y = variable, colour = domain)) +
  
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.25,
    linewidth = 0.9
  ) +
  
  geom_point(size = 3) +
  
  facet_wrap(
   ~domain, ncol = 1,
    scales = "free_y",
  ) +

  scale_colour_manual(values = domain_cols, guide = "none") +
  
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20),
    limits = c(0.1, 40),
    labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10", "20")
  ) +
  
  labs(
    x = "Adjusted odds ratio (log scale)",
    y = NULL
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    
    panel.spacing.y = unit(0.6, "lines"),
    
      strip.background = element_rect(fill = "grey90", colour = "grey40", linewidth = 0.8),
    strip.text.x     = element_text(
      face = "bold",
      size = 15,
      margin = margin(t= 4, b = 4)
    ),
    strip.placement  = "outside",
    
    axis.text.y = element_text(size = 14, colour="black"),
    axis.text.x = element_text(size = 15,colour="black"),
    axis.title.x = element_text(size = 15, colour="black")
  )



ggsave("figures/multivar_forest_plot.jpg", fplot, width = 12, height = 13, dpi = 400)

#### 7. Prediction model: variable selection ####

# First identify the best predictive model (all variables in multivariable analysis with p<0.05)
var <- c("agegroup", "sex", "scale(income_pcap)", "scale(elevation)", "scale(dist_road)",
         "cat", "contact_sewerwater") 

# glmer dredge
f <- as.formula(paste("toxo_igg ~", paste(var, collapse = " + "), "+ (1|hh_id)"))

options(na.action = "na.fail")
m0 <- glmer(f, dat, family="binomial", nAGQ =10, 
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
m0_d <- dredge(m0, m.lim = c(4,7))
saveRDS(m0_d, "outputs/pred_model_glmm_dredge.rds")


m0_d <- readRDS("outputs/pred_model_glmm_dredge.rds")

# Output the model parameter estimates in a table
pred_model <- glmer(toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + I(dist_road/10) + cat + contact_sewerwater + (1|hh_id), dat, family="binomial", nAGQ =10, 
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
icc(pred_model)
summary(pred_model)
pred_model_out <- tibble(parameter = names(exp(confint(pred_model,method="Wald")[,1])), 
                         est = round(c(NA,exp(fixef(pred_model))),2), 
                         lci = round(exp(confint(pred_model,method="Wald")[,1]),2), 
                         uci = round(exp(confint(pred_model,method="Wald")[,2]),2)) %>%
  mutate(OR = paste0(paste0(sprintf("%.2f", est)), 
                     " (",
                     paste0(sprintf("%.2f", lci)),
                     ", ",
                     paste0(sprintf("%.2f", uci)),
                     ")")) %>%
  write_csv("outputs/pred_table.csv")

#### 8. Spatial analysis ####

# 8a. Create prediction grid [DELETE FOR REPO] ---- 

# Study area outline (already UTM/meters)
aoi <- st_read("data/PdL_study_area_2003/PdL_cleaned_study_area.shp", quiet = TRUE) |>
  st_make_valid()

# Extend outline by 5 m buffer (to make sure don't end up with empty cells near edges) ---
aoi_buf <- st_buffer(aoi, 5)

# Create 5m by 5m grid points
grid_pts <- st_make_grid(aoi_buf, cellsize = 4, what = "centers") |>
  st_as_sf(crs = st_crs(aoi_buf)) |>
  st_filter(aoi_buf, .predicate = st_within) |>
  mutate(id = row_number())

# Coordinates as columns (if you need X/Y for functions)
xy <- st_coordinates(grid_pts)
grid_pts <- bind_cols(grid_pts, tibble(X = xy[,1], Y = xy[,2])) |>
  dplyr::select(id, X, Y)

write.csv(st_drop_geometry(grid_pts), "data/cleaned/prediction_grid_points_4m.csv", row.names = FALSE)

## 8b. Intercept-only geostatistical model ----

# create IDs for individuals at same locations
dat$ID <- create.ID.coords(dat, coords= ~X+Y)

## Fitting binomial geostatistical model ----
cmcmc <- list()
cmcmc[[1]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
cmcmc[[2]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
cmcmc[[3]] <- control.mcmc.MCML(n.sim = 65000, burnin = 5000, thin = 6)

# set up number of plug-ins and mcmc samples
control.mcmc <- control.mcmc.MCML(n.sim=1000,burnin=200,thin=8)
dat$units.m <- 1 # 1 individual per row

dat_geo <- dat %>% dplyr::select(toxo_igg, X, Y, units.m, ID, 
                                 agegroup, sex, income_pcap,
                                 elevation, dist_road, cat, contact_sewerwater) %>%
  as.data.frame()

# Fit a model with nugget effect (fixed.rel.nugget = NULL)
f <- toxo_igg~  1
fit_binom <- glm(f, family = "binomial", data = dat_geo)
par0 <- as.numeric(coef(fit_binom))
p <- length(par0)
par0 <- c(par0, 1, 25, 1) # c(beta,sigma2,phi,tau2) tau2 is variance of nugget effect; 

theta <- list()
theta[[1]] <- par0

# iteration 1
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial_intercept <- binomial.logistic.MCML(formula = f, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       ID.coords = dat_geo$ID,
                                       control.mcmc = cmcmc[[1]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial_intercept))
print(summary(geo_binomial_intercept, l = F))
theta[[2]] <- par0
saveRDS(geo_binomial_intercept, file = "outputs/geostat_intercept.rds")

# iteration 2
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial_intercept <- binomial.logistic.MCML(formula = f, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       ID.coords = dat_geo$ID,
                                       control.mcmc = cmcmc[[2]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial_intercept))
print(summary(geo_binomial_intercept, l = F))
theta[[3]] <- par0
saveRDS(geo_binomial_intercept, file = "outputs/geostat_intercept.rds")

# iteration 3
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial_intercept <- binomial.logistic.MCML(formula = f, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       ID.coords = dat_geo$ID,
                                       control.mcmc = cmcmc[[3]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial_intercept))
print(summary(geo_binomial_intercept, l = F))
theta[[4]] <- par0
saveRDS(geo_binomial_intercept, file = "outputs/geostat_intercept.rds")

print(summary(geo_binomial_intercept, l = F))

# get parameter estimates out with 95%CIs
se <- sqrt(diag(geo_binomial_intercept$covariance))
uci <- geo_binomial_intercept$estimate + 1.96*se
uci[c("log(sigma^2)","log(phi)","log(nu^2)")] <- exp(uci[c("log(sigma^2)","log(phi)","log(nu^2)")])
lci<- geo_binomial_intercept$estimate - 1.96*se
lci[c("log(sigma^2)","log(phi)","log(nu^2)")] <- exp(lci[c("log(sigma^2)","log(phi)","log(nu^2)")])

par_est <- tibble(par = names(coef(geo_binomial_intercept)),
                  est = round(coef(geo_binomial_intercept),3),
                  lci = round(lci,3),
                  uci = round(uci,3))
write.csv(par_est, "outputs/geostat_intercept_par.csv")

# geo_binomial_intercept <- readRDS("outputs/geostat_intercept.rds")
# par0 <- as.numeric(coef(geo_binomial_intercept))

## Generate predictions ----
# load grid 
pred_grid <- read.csv("data/cleaned/prediction_grid_points_4m.csv")

cmcmc_pred <- control.mcmc.MCML(n.sim = 25000, burnin = 5000, thin = 10)
geo_pred_intercept <- spatial.pred.binomial.MCML(geo_binomial_intercept,
                                          grid.pred = pred_grid[,c("X","Y")], 
                                          # predictors = predictors,  
                                          control.mcmc = cmcmc_pred,
                                          type = "marginal",
                                          messages = T,
                                          plot.correlogram = T)
saveRDS(geo_pred_intercept, "outputs/geostat_intercept_pred.rds")
geo_pred_intercept <- readRDS("outputs/geostat_intercept_pred.rds")

# make into a dataframe then a raster
prev_preds_intercept <- tibble(
  X    = geo_pred_intercept$grid$X,
  Y    = geo_pred_intercept$grid$Y,
  prev = geo_pred_intercept$prevalence$predictions,
  prev_u95ci = apply(geo_pred_intercept$samples, 1, function(draws) {
    plogis(quantile(draws, 0.975))
  }),
  prev_l95ci = apply(geo_pred_intercept$samples, 1, function(draws) {
    plogis(quantile(draws, 0.025))
  }),
  ex50 = rowMeans(geo_pred_intercept$samples > 0.5)
)

pred_rast_intercept <- rasterFromXYZ(prev_preds_intercept, crs = crs(aoi))

plot(pred_rast_intercept[["prev"]])
plot(pred_rast_intercept[["prev_u95ci"]])
plot(pred_rast_intercept[["prev_l95ci"]])
plot(pred_rast_intercept[["ex50"]])


## 8c. Full geostatistical model ----

## Fitting binomial geostatistical model ----
# Fit a model with nugget effect (fixed.rel.nugget = NULL)
f_full <- toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + I(dist_road/10) + cat + contact_sewerwater
fit_binom <- glm(f_full, family = "binomial", data = dat_geo)
par_covar <- as.numeric(coef(fit_binom))
p <- length(par_covar)
par0 <- c(fixef(pred_model),par0[2:4]) # c(beta,sigma2,phi,tau2) tau2 is variance of nugget effect; using estimates from intercept-only model here

theta <- list()
theta[[1]] <- par0

# iteration 1
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial_full <- binomial.logistic.MCML(formula = f_full, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       ID.coords = dat_geo$ID,
                                       control.mcmc = cmcmc[[1]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial_full))
print(summary(geo_binomial_full, l = F))
theta[[2]] <- par0
saveRDS(geo_binomial_full, file = "outputs/geostat_full.rds")

# iteration 2
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial_full <- binomial.logistic.MCML(formula = f_full, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       ID.coords = dat_geo$ID,
                                       control.mcmc = cmcmc[[2]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial_full))
print(summary(geo_binomial_full, l = F))
theta[[3]] <- par0
saveRDS(geo_binomial_full, file = "outputs/geostat_full.rds")

# iteration 3
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial_full <- binomial.logistic.MCML(formula = f_full, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       ID.coords = dat_geo$ID,
                                       control.mcmc = cmcmc[[3]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial_full))
print(summary(geo_binomial_full, l = F))
theta[[4]] <- par0
saveRDS(geo_binomial_full, file = "outputs/geostat_full.rds")

print(summary(geo_binomial_full, l = F))
geo_binomial_full <- readRDS(file = "outputs/geostat_full.rds")

# get parameter estimates out with 95%CIs
se <- sqrt(diag(geo_binomial_full$covariance))
uci <- geo_binomial_full$estimate + 1.96*se
lci<- geo_binomial_full$estimate - 1.96*se

par_est <- tibble(par = names(coef(geo_binomial_full)),
                  est = round(exp(geo_binomial_full$estimate),3),
                  lci = round(exp(lci),3),
                  uci = round(exp(uci),3))
write.csv(par_est, "outputs/geostat_full_par.csv")

# compute spatial ICC analogue for geostat model (the proportion of residual variance attributable to spatial vs. household-level processes)
sigma2 <- par_est$est[par_est$par=="sigma^2"]
tau2 <- par_est$est[par_est$par=="tau^2"]

icc_spat <- sigma2/(sigma2 + tau2)
icc_spat

# i.e. About 31% of the unexplained variation is spatially structured over the range defined by φ (~30 m), and 69% occurs at very fine scales captured by the nugget.


## Generate predictions ----
# To get S(x) we first predict on logit scale based on a given set of predictor values
# then we subtract the linear predictor for these predictor values from the prediction

# Create dummy predictor values
newdat <- pred_grid %>%
  mutate(
    agegroup = relevel(factor(levels(factor(dat_geo$agegroup))[1], levels = levels(factor(dat_geo$agegroup))), ref = "10-12"),
    sex      = factor(levels(factor(dat_geo$sex))[1],      levels = levels(factor(dat_geo$sex))),
    cat      = factor(levels(factor(dat_geo$cat))[1],      levels = levels(factor(dat_geo$cat))),
    contact_sewerwater = factor(levels(factor(dat_geo$contact_sewerwater))[1],
                                levels = levels(factor(dat_geo$contact_sewerwater))),
    income_pcap = mean(dat_geo$income_pcap, na.rm = TRUE),
    elevation = mean(dat_geo$elevation,   na.rm = TRUE),
    dist_road = mean(dat_geo$dist_road,   na.rm = TRUE),
    `I(elevation/10)`   = mean(dat_geo$elevation/10,   na.rm = TRUE),
    `I(dist_road/10)`   = mean(dat_geo$dist_road/10,   na.rm = TRUE),
    units.m     = 1   # if your fit included units.m
  )

# Marginal predictions on the LOGIT scale (eta = Xβ + S(x)), using dummy predictor values
cmcmc_pred <- control.mcmc.MCML(n.sim = 25000, burnin = 5000, thin = 20)

geo_pred_full <- spatial.pred.binomial.MCML(
  geo_binomial_full,
  grid.pred         = newdat[, c("X","Y")],
  predictors        = newdat ,
  control.mcmc      = cmcmc_pred,
  type              = "marginal",
  scale.predictions = "logit",
  messages          = TRUE,
  plot.correlogram  = FALSE
)

saveRDS(geo_pred_full, "outputs/geostat_full_pred.rds")

geo_pred_full <- readRDS("outputs/geostat_full_pred.rds")

# make into a dataframe then a raster
prev_preds_full <- tibble(
  X    = geo_pred_full$grid$X,
  Y    = geo_pred_full$grid$Y,
  eta_logit = geo_pred_full$logit$predictions,
  prev_u95ci_logit = geo_pred_full$logit$quantiles[,2],
  prev_l95ci_logit = geo_pred_full$logit$quantiles[,1],
  prev = plogis(geo_pred_full$logit$predictions),
  prev_u95ci = plogis(geo_pred_full$logit$quantiles[,2]),
  prev_l95ci = plogis(geo_pred_full$logit$quantiles[,1])
)

# now need to remove the linear predictor part:
eta_hat <- as.numeric(prev_preds_full$eta_logit)  

# Compute Xβ on the grid using the SAME formula/contrasts, then subtract
rhs_terms  <- delete.response(terms(f_full))
Xgrid_full <- model.matrix(rhs_terms, data = newdat)  # includes intercept column

# Extract fixed effects (drop spatial params), align to Xgrid_full
b_all    <- coef(geo_binomial_full)
sp_names <- intersect(names(b_all), c("sigma^2","phi","tau^2"))
b_fix    <- b_all[setdiff(names(b_all), sp_names)]

beta0    <- if ("(Intercept)" %in% names(b_fix)) unname(b_fix["(Intercept)"]) else 0
b_slopes <- b_fix[setdiff(names(b_fix), "(Intercept)")]
b_slopes <- b_slopes[colnames(Xgrid_full)[-1]]   # align to columns (no intercept)

fixed_hat <- as.numeric(beta0 + Xgrid_full[, -1, drop = FALSE] %*% b_slopes)

# Isolate S(x)
prev_preds_full$Sx_mean <- eta_hat - fixed_hat -1.88528501 

# make into a raster
pred_rast_full <- rasterFromXYZ(prev_preds_full, crs = crs(aoi))

plot(pred_rast_full[["prev"]])
plot(pred_rast_full[["prev_u95ci"]])
plot(pred_rast_full[["prev_l95ci"]])
plot(pred_rast_full[["Sx_mean"]])



## 8d. Full geostatistical model minus distance to road ----

## Fitting binomial geostatistical model ----
# Fit a model with nugget effect (fixed.rel.nugget = NULL)

# get spatial parameter estimates from full mode
# par0spat <- par0[(p+1):(p+3)]
# 
# f_noroad <- toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + cat + contact_sewerwater
# fit_binom <- glm(f_noroad, family = "binomial", data = dat_geo)
# par_covar <- as.numeric(coef(fit_binom))
# p <- length(par_covar)
# par0 <- c(par_covar,par0spat) # c(beta,sigma2,phi,tau2) tau2 is variance of nugget effect; using estimates from full model here
# 
# theta <- list()
# theta[[1]] <- par0
# 
# # iteration 1
# init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
# geo_binomial_no_road <- binomial.logistic.MCML(formula = f_noroad, units.m = ~ units.m,
#                                        coords = ~ X + Y,
#                                        data = dat_geo,
#                                        par0 = par0,
#                                        ID.coords = dat_geo$ID,
#                                        control.mcmc = cmcmc[[1]], 
#                                        kappa = 0.5,
#                                        start.cov.pars = init_cov_pars,
#                                        method = "nlminb", messages = T)
# par0 <- as.numeric(coef(geo_binomial_no_road))
# print(summary(geo_binomial_no_road, l = F))
# theta[[2]] <- par0
# saveRDS(geo_binomial_no_road, file = "outputs/geostat_no_road.rds")
# 
# # iteration 2
# init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
# geo_binomial_no_road <- binomial.logistic.MCML(formula = f_noroad, units.m = ~ units.m,
#                                        coords = ~ X + Y,
#                                        data = dat_geo,
#                                        par0 = par0,
#                                        ID.coords = dat_geo$ID,
#                                        control.mcmc = cmcmc[[2]], 
#                                        kappa = 0.5,
#                                        start.cov.pars = init_cov_pars,
#                                        method = "nlminb", messages = T)
# par0 <- as.numeric(coef(geo_binomial_no_road))
# print(summary(geo_binomial_no_road, l = F))
# theta[[3]] <- par0
# saveRDS(geo_binomial_no_road, file = "outputs/geostat_no_road.rds")
# 
# # iteration 3
# init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
# geo_binomial_no_road <- binomial.logistic.MCML(formula = f_noroad, units.m = ~ units.m,
#                                        coords = ~ X + Y,
#                                        data = dat_geo,
#                                        par0 = par0,
#                                        ID.coords = dat_geo$ID,
#                                        control.mcmc = cmcmc[[3]], 
#                                        kappa = 0.5,
#                                        start.cov.pars = init_cov_pars,
#                                        method = "nlminb", messages = T)
# par0 <- as.numeric(coef(geo_binomial_no_road))
# print(summary(geo_binomial_no_road, l = F))
# theta[[4]] <- par0
# saveRDS(geo_binomial_no_road, file = "outputs/geostat_no_road.rds")
# geo_binomial_no_road <- readRDS("outputs/geostat_no_road.rds")
# 
# print(summary(geo_binomial_no_road, l = F))
# 
# # Marginal predictions on the LOGIT scale (eta = Xβ + S(x)), using dummy predictor values
# cmcmc_pred <- control.mcmc.MCML(n.sim = 25000, burnin = 5000, thin = 20)
# 
# geo_pred_no_road <- spatial.pred.binomial.MCML(
#   geo_binomial_no_road,
#   grid.pred         = newdat[, c("X","Y")],
#   predictors        = newdat,
#   control.mcmc      = cmcmc_pred,
#   type              = "marginal",
#   scale.predictions = "logit",
#   messages          = TRUE,
#   plot.correlogram  = FALSE
# )
# 
# saveRDS(geo_pred_no_road, "outputs/geostat_no_road_pred.rds")
# 
# # make into a dataframe then a raster
# prev_preds_no_road <- tibble(
#   X    = geo_pred$grid$X,
#   Y    = geo_pred$grid$Y,
#   eta_logit = geo_pred$logit$predictions,
#   prev_u95ci_logit = geo_pred$logit$quantiles[,2],
#   prev_l95ci_logit = geo_pred$logit$quantiles[,1],
#   prev = plogis(geo_pred$logit$predictions),
#   prev_u95ci = plogis(geo_pred$logit$quantiles[,2]),
#   prev_l95ci = plogis(geo_pred$logit$quantiles[,1])
# )
# 
# # now need to remove the linear predictor part:
# eta_hat <- as.numeric(prev_preds_no_road$eta_logit)  
# 
# # Compute Xβ on the grid using the SAME formula/contrasts, then subtract
# rhs_terms  <- delete.response(terms(f_noroad))
# Xgrid_full <- model.matrix(rhs_terms, data = newdat)  # includes intercept column
# 
# # Extract fixed effects (drop spatial params), align to Xgrid_full
# b_all    <- coef(geo_binomial_no_road)
# sp_names <- intersect(names(b_all), c("sigma^2","phi","tau^2"))
# b_fix    <- b_all[setdiff(names(b_all), sp_names)]
# 
# beta0    <- if ("(Intercept)" %in% names(b_fix)) unname(b_fix["(Intercept)"]) else 0
# b_slopes <- b_fix[setdiff(names(b_fix), "(Intercept)")]
# b_slopes <- b_slopes[colnames(Xgrid_full)[-1]]   # align to columns (no intercept)
# 
# fixed_hat <- as.numeric(beta0 + Xgrid_full[, -1, drop = FALSE] %*% b_slopes)
# 
# # Isolate S(x)
# prev_preds_no_road$Sx_mean <- eta_hat - fixed_hat
# 
# # make into a raster
# pred_rast_no_road <- rasterFromXYZ(prev_preds_no_road, crs = crs(aoi))
# 
# plot(pred_rast_no_road[["prev"]])
# plot(pred_rast_no_road[["prev_u95ci"]])
# plot(pred_rast_no_road[["prev_l95ci"]])
# plot(pred_rast_no_road[["Sx_mean"]])
# 
# pred_rast_no_road <- rast(pred_rast_no_road)
# terra::writeRaster(pred_rast_no_road, 
#                    paste0("outputs/geostat_no_road_pred_raster.tif"), 
#                    overwrite = TRUE)

## 8e. map making ----
# Colours schemes
pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(200)

# function for cropping to fit cleanly within study area
crop.fn <- function(raster, outline){
  raster <- crop(raster, extent(outline))
  raster <- terra::mask(raster, outline)
}

# crop rasters
p_intercept <- crop.fn(pred_rast_intercept[["prev"]], aoi)
plot(p_intercept)

p_intercept_lci <- crop.fn(pred_rast_intercept[["prev_l95ci"]], aoi)
plot(p_intercept_lci)

p_intercept_uci <- crop.fn(pred_rast_intercept[["prev_u95ci"]], aoi)
plot(p_intercept_uci)

p_intercept_ex50 <- crop.fn(pred_rast_intercept[["ex50"]], aoi)
plot(p_intercept_ex50)

Sx_full <- crop.fn(pred_rast_full[["Sx_mean"]], aoi)
plot(Sx_full)

# choose line width
lwd <- 1.5
# choose resolution
res <- 300

# choose x dimension of output
xdim <- 8
ydim <- xdim*(7.94/10)

# Figure 4. Spatial predictions
# function for identifying potential break sizes
breaks.fn <- function(min,max,level=0.05){
  message(paste0("Difference = ", max-min))
  x <- seq(level,20,by=level)
  out <- x[which(as.integer((max-min)/x) == (max-min)/x)]
  return(out)
}

legend_w   <- 1.3   # wider colourbar
legend_sh  <- 0.85  # longer colourbar
cex_leg    <- 1.5   # legend tick text size
cex_title  <- 1.8   # panel title size
title_line <- 1   # distance of title from plot
title_font <- 2     # bold

par(mar = c(3, 3, 5, 6)) 

# seroprevalence - mean
p_intercept.min <- round_any(raster::cellStats(p_intercept, "min", na.rm = TRUE), 0.01, floor)
p_intercept.max <- round_any(raster::cellStats(p_intercept, "max", na.rm = TRUE), 0.01, ceiling)

p_intercept.min <- 0
p_intercept.max <- 0.8

# Legend ticks + labels (0–80%)
at_A  <- seq(p_intercept.min, p_intercept.max, by = 0.2)
lab_A <- paste0(at_A * 100, "%")

jpeg("outputs/pred_p_intercept.jpeg", units = "in", width = xdim, height = ydim, res = res)

plot(st_geometry(aoi), lwd = lwd)

plot(
  p_intercept,
  col = rev(pal),
  xaxt = "n", yaxt = "n",
  axes = FALSE, box = FALSE,
  zlim = c(p_intercept.min, p_intercept.max),
  legend = FALSE,
  add = TRUE
)

mtext(
  "Predicted seroprevalence",
  side = 3, line = title_line, adj = 0.5,
  cex = cex_title, font = 2
)

plot(
  p_intercept,
  legend.only = TRUE,
  col = rev(pal),
  zlim = c(p_intercept.min, p_intercept.max),
  legend.width = legend_w,
  legend.shrink = legend_sh,
  axis.args = list(cex.axis = cex_leg, at = at_A, labels = lab_A)
)

plot(st_geometry(aoi), add = TRUE, lwd = lwd)

dev.off()

# seroprevalence - ex50
p_intercept_ex50.min <- round_any(raster::cellStats(p_intercept_ex50, 'min', na.rm = TRUE), 0.01, floor) # limits for legend
p_intercept_ex50.max <- round_any(raster::cellStats(p_intercept_ex50, 'max', na.rm = TRUE), 0.01, ceiling)

at_ex50  <- seq(p_intercept_ex50.min, p_intercept_ex50.max, by = 0.2)
lab_ex50 <- paste0(at_ex50 * 100, "%")  # 0–80%

jpeg("outputs/pred_p_intercept_ex50.jpeg", units = "in", width = xdim, height = ydim, res = res)

plot(st_geometry(aoi), lwd = lwd)

plot(p_intercept_ex50, col = rev(pal),
     xaxt = "n", yaxt = "n", axes = FALSE, box = FALSE,
     zlim = c(p_intercept_ex50.min, p_intercept_ex50.max),
     legend = FALSE, add = TRUE)

mtext("Exceedance probability (>50%)",
      side = 3, line = title_line, adj = 0.5,
      cex = cex_title, font = title_font)

plot(p_intercept_ex50, legend.only = TRUE, col = rev(pal),
     zlim = c(p_intercept_ex50.min, p_intercept_ex50.max),
     legend.width = legend_w,
     legend.shrink = legend_sh,
     axis.args = list(cex.axis = cex_leg, at = at_ex50, labels = lab_ex50))

plot(st_geometry(aoi), add = TRUE, lwd = lwd)

dev.off()



# S(x) full model - mean
Sx_full.min <- floor(raster::cellStats(Sx_full, "min", na.rm = TRUE))
Sx_full.max <- ceiling(raster::cellStats(Sx_full, "max", na.rm = TRUE))

lim <- max(abs(c(Sx_full.min, Sx_full.max)))
Sx_full.min <- -lim
Sx_full.max <-  lim

at_Sx  <- seq(Sx_full.min, Sx_full.max, by = 0.5)
lab_Sx <- format(at_Sx, trim = TRUE)

jpeg("outputs/pred_Sx_full.jpeg", units = "in", width = xdim, height = ydim, res = res)

plot(st_geometry(aoi), lwd = lwd)

plot(Sx_full, col = rev(pal),
     xaxt = "n", yaxt = "n", axes = FALSE, box = FALSE,
     zlim = c(Sx_full.min, Sx_full.max),
     legend = FALSE, add = TRUE)

mtext("Residual spatial effect S(x)",
      side = 3, line = title_line, adj = 0.5,
      cex = cex_title, font = title_font)

plot(Sx_full, legend.only = TRUE, col = rev(pal),
     zlim = c(Sx_full.min, Sx_full.max),
     legend.width = legend_w,
     legend.shrink = legend_sh,
     axis.args = list(cex.axis = cex_leg, at = at_Sx, labels = lab_Sx))

plot(st_geometry(aoi), add = TRUE, lwd = lwd)

dev.off()

# Bring plots together as panels
library(cowplot)
library(magick)

read_panel <- function(path) image_read(path) |> image_border("white", "20x20")  # no trim

imgA <- read_panel("outputs/pred_p_intercept.jpeg")
imgB <- read_panel("outputs/pred_p_intercept_ex50.jpeg")
imgC <- read_panel("outputs/pred_Sx_full.jpeg")

pA <- ggdraw() + draw_image(imgA)
pB <- ggdraw() + draw_image(imgB)
pC <- ggdraw() + draw_image(imgC)

# Put C under A (and keep an empty placeholder under B)
blank <- ggdraw()

panel <- plot_grid(
  pA, pB,
  pC, blank,
  ncol = 2,
  rel_widths  = c(1, 1),
  rel_heights = c(1, 1)
)

# Add panel letters in consistent positions
panel <- ggdraw(panel) +
  draw_label("A", x = 0.04, y = 0.966,
             hjust = 0, vjust = 1,
             fontface = "bold", size = 30) +
  draw_label("B", x = 0.54, y = 0.966,
             hjust = 0, vjust = 1,
             fontface = "bold", size = 30) +
  draw_label("C", x = 0.04, y = 0.466,
             hjust = 0, vjust = 1,
             fontface = "bold", size = 30)

ggsave("outputs/Fig4_pred_panel.jpg", panel, width = 12, height = 10, dpi = 400)


# save rasters
pred_rast <- rast(pred_rast)

terra::writeRaster(pred_rast, 
                   paste0("outputs/geostat_intercept_pred_raster.tif"), 
                   overwrite = TRUE)

pred_rast_full <- rast(pred_rast_full)
terra::writeRaster(pred_rast_full, 
                   paste0("outputs/geostat_full_pred_raster.tif"), 
                   overwrite = TRUE)


# Create forest plot of spatial parameter estimates

# -----------------------------
# Spatial parameter data
# -----------------------------

spatial_data <- tribble(
  ~model, ~parameter, ~estimate, ~lower, ~upper,
  "Intercept-only", "φ (range, m)",               58.72, 24.68, 139.71,
  "Intercept-only", "σ² (spatial variance)",       0.76,  0.46,   1.25,
  "Intercept-only", "τ² (nugget variance)",        0.75,  0.44,   2.17,
  "Full model",     "φ (range, m)",               29.95, 15.48,  57.95,
  "Full model",     "σ² (spatial variance)",       0.66,  0.37,   1.18,
  "Full model",     "τ² (nugget variance)",        1.48,  0.62,   3.55
)

# -----------------------------
# Compute ICC proxy (proportion spatially structured)
# -----------------------------
icc <- tribble(
  ~model, ~icc_value,
  "Intercept-only", 0.76 / (0.76 + 0.75),
  "Full model",     0.66 / (0.66 + 1.48)
) %>%
  mutate(icc_label = sprintf("%s = %.2f", model, icc_value))

icc_text <- paste("Proportion spatially structured:\n",
                  paste0("  ", icc$icc_label, collapse = "\n"))

# -----------------------------
# Forest plot with facets
# -----------------------------

pD <- ggplot(spatial_data,
             aes(x = estimate,
                 y = model,
                 colour = model)) +
  
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.2, linewidth = 0.9) +
  
  facet_wrap(~parameter, scales = "free_x", ncol = 1, strip.position = "top") +
  
  scale_colour_manual(values = c("Intercept-only" = "#1b9e77",
                                 "Full model" = "#d95f02")) +
  
  labs(x = "Estimate (95% CI)", y = NULL,
       title = "Comparison of Spatial Parameters Between Models") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  
  # annotation: ICC proxy displayed once
  annotate("text", x = -Inf, y = -Inf,
           label = icc_text,
           hjust = -0.05, vjust = -1.2,
           size = 4)

pD



#### 8. GAM plot of residual spatial variability with all variables ####

# GAM with dist_road
gam_model <- gam(toxo_igg ~ agegroup + sex + income_pcap + elevation + contact_sewerwater + dist_road,
                 data = dat, 
                 family = binomial(link = "logit"))

# Get predictions and residuals
dat$predicted <- predict(gam_model, type = "response")
dat$residuals <- residuals(gam_model, type = "response")

# Fit a New GAM to Model the Spatial Residual Pattern
residual_gam <- gam(residuals ~ s(X, Y), data = dat)

# Step 1: Create a regular grid over the study area
x_range <- seq(min(dat$X), max(dat$X), length.out = 100)  # Adjust grid density
y_range <- seq(min(dat$Y), max(dat$Y), length.out = 100)

grid <- expand.grid(X = x_range, Y = y_range)

# Step 2: Predict the residuals at the grid points
grid$predicted_residuals <- predict(residual_gam, newdata = grid, type = "response")

# Reshape the predicted values into a matrix for contour plotting
z_matrix <- matrix(grid$predicted_residuals, nrow = length(x_range), ncol = length(y_range))

# Convert z_matrix into a format suitable for ggplot2
# We need to convert the grid and the predicted values into a data frame
grid_df <- expand.grid(X = x_range, Y = y_range)
grid_df$z <- as.vector(z_matrix)  # Add the predicted values to the grid

# Create the contour plot using ggplot2
ggplot() +
  geom_contour(data = grid_df, aes(x = X, y = Y, z = z), color = "black", bins = 10) +  # Contour lines
  geom_point(data = dat, aes(x = X, y = Y, color = factor(toxo_igg)), alpha = 0.5, size = 2) +  # Observed data points
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Negative", "Positive")) +  # Custom legend labels
  labs(x = "Longitude", y = "Latitude", color = "Serostatus") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border around plot
        panel.grid = element_blank())
ggsave("figures/spatial_resid_gam_all_var.jpeg", dpi=300, width = 10, height = 6)

#### 9. GAM plot of residual spatial variability without distance to road ####
# GAM without dist_road
gam_model <- gam(toxo_igg ~ agegroup + sex + income_pcap + elevation +  contact_sewerwater,
                 data = dat, 
                 family = binomial(link = "logit"))

# Get predictions and residuals
dat2$predicted <- predict(gam_model, type = "response")
dat2$residuals <- residuals(gam_model, type = "response")

# Fit a New GAM to Model the Spatial Residual Pattern
residual_gam <- gam(residuals ~ s(X, Y), data = dat2)

# Assuming the gam_model has already been fitted and the 'dat2' dataset is loaded

# Step 1: Create a regular grid over the study area
x_range <- seq(min(dat2$X), max(dat2$X), length.out = 100)  # Adjust grid density
y_range <- seq(min(dat2$Y), max(dat2$Y), length.out = 100)

grid <- expand.grid(X = x_range, Y = y_range)

# Step 2: Predict the residuals at the grid points
grid$predicted_residuals <- predict(residual_gam, newdata = grid, type = "response")

# Reshape the predicted values into a matrix for contour plotting
z_matrix <- matrix(grid$predicted_residuals, nrow = length(x_range), ncol = length(y_range))

# Convert z_matrix into a format suitable for ggplot2
# We need to convert the grid and the predicted values into a data frame
grid_df <- expand.grid(X = x_range, Y = y_range)
grid_df$z <- as.vector(z_matrix)  # Add the predicted values to the grid

# Create the contour plot using ggplot2
ggplot() +
  geom_contour(data = grid_df, aes(x = X, y = Y, z = z), color = "black", bins = 10) +  # Contour lines
  geom_point(data = dat2, aes(x = X, y = Y, color = factor(toxo_igg)), alpha = 0.5, size = 2) +  # Observed data points
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Negative", "Positive")) +  # Custom legend labels
  labs(x = "Longitude", y = "Latitude", color = "Serostatus") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        legend.background = element_rect(color = "black", size = 0.5),  # Add a box around the legend
        legend.box = "vertical",
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border around plot
        panel.grid = element_blank())
ggsave("figures/spatial_resid_gam_no_road.jpeg", dpi=300, width = 10, height = 6)

# make a joint figure
gamp0 <- ggdraw() + draw_image("figures/spatial_resid_gam_all_var.jpeg")
gamp1 <- ggdraw() + draw_image("figures/spatial_resid_gam_no_road.jpeg")

fig_gams <- cowplot::plot_grid(gamp0, gamp1, labels = c("A","B"), label_size = 40, label_y = 1.01, ncol=2)

jpeg("figures/spatial_resid_gam_figure.jpeg", units="in", width=20, height=6, res=300)
print(fig_gams)
dev.off()

## get variograms ----
data_resid <- dat %>% select(X, Y)

# full model
pred_model <- glmer(toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + I(dist_road/10) + cat + contact_sewerwater + (1|hh_id), dat, family="binomial", nAGQ =10, 
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# fit model without road
pred_model_no_road <- glmer(toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + cat + contact_sewerwater + (1|hh_id), dat, family="binomial", nAGQ =10, 
                            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
summary(pred_model_no_road)

# residuals for full model and model without road out
data_resid$res_full <- residuals(pred_model)
data_resid$res_no_road <- residuals(pred_model_no_road)


vario_full <- ggvario(coords = data_resid[,c("X","Y")], 
        data = data_resid$res_full, xlab = "distance (m)", 
        show_nbins = F, nsim = 1000, envelop = T)
vario_full

vario_no_road <- ggvario(coords = data_resid[,c("X","Y")], 
                      data = data_resid$res_no_road, xlab = "distance (m)", 
                      show_nbins = F, nsim = 1000, envelop = T)
vario_no_road


#
p_full    <- predict(pred_model, type="response", re.form = NA)  # no random effects
p_no_road <- predict(pred_model_no_road, type="response", re.form = NA)
r_full    <- dat$toxo_igg - p_full
r_no_road <- dat$toxo_igg - p_no_road


data_resid$res_full <- r_full
data_resid$res_no_road <- r_no_road

vario_full <- ggvario(coords = data_resid[,c("X","Y")], 
                      data = data_resid$res_full, xlab = "distance (m)", 
                      show_nbins = F, nsim = 1000, envelop = T)
vario_full

vario_no_road <- ggvario(coords = data_resid[,c("X","Y")], 
                         data = data_resid$res_no_road, xlab = "distance (m)", 
                         show_nbins = F, nsim = 1000, envelop = T)
vario_no_road

# intercept-only geostatistical model

# create IDs for individuals at same locations (important for geostat fitting models)
dat <- as.data.frame(dat)
dat$ID <- create.ID.coords(dat, coords= ~X+Y)

# # Fitting binomial geostatistical model
cmcmc <- list()
cmcmc[[1]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
cmcmc[[2]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
cmcmc[[3]] <- control.mcmc.MCML(n.sim = 65000, burnin = 5000, thin = 6)

# set up number of plug-ins and mcmc samples
control.mcmc <- control.mcmc.MCML(n.sim=1000,burnin=200,thin=8)
dat$units.m <- 1 # 1 individual per row

dat_geo <- dat %>% dplyr::select(toxo_igg, X, Y, units.m, ID, 
                          agegroup, sex, income_pcap,
                          elevation, dist_road, contact_sewerwater)

# 4. Fit a model with nugget effect (fixed.rel.nugget = NULL)
f <- toxo_igg~  1
fit_binom <- glm(f, family = "binomial", data = dat_geo)
par0 <- as.numeric(coef(fit_binom))
p <- length(par0)
par0 <- c(par0, 1, 25, 1) # c(beta,sigma2,phi,tau2) tau2 is variance of nugget effect; 

theta <- list()
theta[[1]] <- par0

# iteration 1
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial <- binomial.logistic.MCML(formula = f, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       control.mcmc = cmcmc[[1]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial))
print(summary(geo_binomial, l = F))
theta[[2]] <- par0

# iteration 2
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial <- binomial.logistic.MCML(formula = f, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       control.mcmc = cmcmc[[2]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial))
print(summary(geo_binomial, l = F))
theta[[3]] <- par0

# iteration 3
init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
geo_binomial <- binomial.logistic.MCML(formula = f, units.m = ~ units.m,
                                       coords = ~ X + Y,
                                       data = dat_geo,
                                       par0 = par0,
                                       control.mcmc = cmcmc[[3]], 
                                       kappa = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
par0 <- as.numeric(coef(geo_binomial))
print(summary(geo_binomial, l = F))
theta[[4]] <- par0


# Get phi, sigma2, tau2 out for all iterations into tibble:
theta.hat <- tibble(phi=unlist(lapply(theta, '[[', length(x)-1)),
                    sigma2=unlist(lapply(theta, '[[', length(x)-2)),
                    tau2=unlist(lapply(theta, '[[', length(x))))





## 7a. Get variograms of residuals ----

# full model
pred_model <- glmer(toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + I(dist_road/10) + cat + contact_sewerwater + (1|hh_id), dat, family="binomial", nAGQ =10, 
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# fit model without road
pred_model_no_road <- glmer(toxo_igg ~ agegroup + sex + income_pcap + I(elevation/10) + cat + contact_sewerwater + (1|hh_id), dat, family="binomial", nAGQ =10, 
                            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


# get residuals for full model and model without road out
# 
data_resid <- dat %>% select(X, Y, hh_id) %>%
  distinct()

data_resid$res_full <- residuals(pred_model, type = "pearson")
data_resid$res_no_road <- residuals(pred_model_no_road, type = "pearson")

data_resid$res_full <- plogis(ranef(pred_model)$hh_id)
data_resid$res_no_road <- residuals(pred_model_no_road, type = "pearson")


vario_full <- ggvario(coords = data_resid[,c("X","Y")], 
                      data = data_resid$res_full, xlab = "distance (m)", 
                      show_nbins = F, nsim = 1000, envelop = T)
vario_full

vario_no_road <- ggvario(coords = data_resid[,c("X","Y")], 
                         data = data_resid$res_no_road, xlab = "distance (m)", 
                         show_nbins = F, nsim = 1000, envelop = T)
vario_no_road


#
p_full    <- predict(pred_model, type="response", re.form = NA)  # no random effects
p_no_road <- predict(pred_model_no_road, type="response", re.form = NA)
r_full    <- dat$toxo_igg - p_full
r_no_road <- dat$toxo_igg - p_no_road


data_resid$res_full <- r_full
data_resid$res_no_road <- r_no_road

vario_full <- ggvario(coords = data_resid[,c("X","Y")], 
                      data = data_resid$res_full, xlab = "distance (m)", 
                      show_nbins = F, nsim = 1000, envelop = T)
vario_full

vario_no_road <- ggvario(coords = data_resid[,c("X","Y")], 
                         data = data_resid$res_no_road, xlab = "distance (m)", 
                         show_nbins = F, nsim = 1000, envelop = T)
vario_no_road



## 8. Additional analyses  ----
library(lme4)
library(EValue)
library(dplyr)
library(purrr)

## 8a. E-values for ORs  ----
# Sensitivity analysis of "mechanistic" (conditional) effects for unmeasured confounding

# 

# Fit logistic mixed model
m8 <- glmer(
  toxo_igg ~ cat + agegroup + race + scale(income_pcap) + (1|hh_id),
  dat %>% filter(!is.na(dist_trash)),
  family = "binomial",
  nAGQ = 10,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 100000))
)

# Extract conditional OR and CI
or_est <- exp(fixef(m8)["catYes"])   # assuming "No" is reference
ci_or  <- confint(m8, parm = "catYes", method = "Wald")
or_lo  <- exp(ci_or[1])
or_hi  <- exp(ci_or[2])

cat("Conditional OR:", or_est, "CI:", or_lo, "-", or_hi, "\n")

# E-values for conditional OR (non-rare outcome)
e_or <- EValue::evalues.OR(or_est, lo = or_lo, hi = or_hi, rare = FALSE)
e_or

## 8b. G-computation for marginal RRs & E-values  ----

# Prepare datasets with cat = Yes and cat = No
dat_yes <- dat %>% mutate(cat = "Yes")
dat_no  <- dat %>% mutate(cat = "No")

# Predict *marginal* probabilities (re.form = NA)
p_yes <- predict(m8, newdata = dat_yes, type = "response", re.form = NA)
p_no  <- predict(m8, newdata = dat_no,  type = "response", re.form = NA)

risk_yes <- mean(p_yes)
risk_no  <- mean(p_no)

cat("Marginal risk (cat=Yes):", risk_yes, "\n")
cat("Marginal risk (cat=No): ", risk_no,  "\n")

# Marginal RR
RR_marg <- risk_yes / risk_no
RR_marg

## Bootstrap CI for marginal RR ----
set.seed(123)
B <- 10

boot_RR <- replicate(B, {
  # resample data
  idx <- sample(seq_len(nrow(dat)), replace = TRUE)
  dat_b <- dat[idx, ] %>% filter(!is.na(dist_trash))
  
  # refit model
  m_b <- try(
    glmer(
      toxo_igg ~ cat + agegroup + race + scale(income_pcap) + (1|hh_id),
      dat_b, family = "binomial", nAGQ = 10,
      control = glmerControl(optimizer = "bobyqa",
                             optCtrl = list(maxfun = 100000))
    ),
    silent = TRUE
  )
  
  if (inherits(m_b, "try-error")) return(NA)
  
  # prepare counterfactual datasets
  dat_yes_b <- dat_b %>% mutate(cat = "Yes")
  dat_no_b  <- dat_b %>% mutate(cat = "No")
  
  # marginal predictions
  p_yes_b <- predict(m_b, newdata = dat_yes_b, type = "response", re.form = NA)
  p_no_b  <- predict(m_b, newdata = dat_no_b,  type = "response", re.form = NA)
  
  mean(p_yes_b) / mean(p_no_b)
})

boot_RR <- boot_RR[!is.na(boot_RR)]
RR_CI <- quantile(boot_RR, c(0.025, 0.975))
RR_CI

## E-values for marginal RR ----
EValue::evalues.RR(RR_marg, lo = RR_CI[1], hi = RR_CI[2])

## 8c. Estimate PAF ----

# Step 1: predicted risk under actual exposure distribution
p_obs <- predict(m8, newdata = dat,
                 type = "response", re.form = NA)
risk_obs <- mean(p_obs)

# Step 2: predicted risk if no one had cats (cat = "No")
p_no <- predict(m8, newdata = dat_no, type = "response", re.form = NA)
risk_no <- mean(p_no)

# PAF
PAF <- (risk_obs - risk_no) / risk_obs
PAF

## Bootstrap CI for PAF ----
boot_PAF <- replicate(B, {
  idx <- sample(seq_len(nrow(dat)), replace = TRUE)
  dat_b <- dat[idx, ] %>% filter(!is.na(dist_trash))
  
  m_b <- try(
    glmer(
      toxo_igg ~ cat + agegroup + race + scale(income_pcap) + (1|hh_id),
      dat_b, family = "binomial", nAGQ = 10,
      control = glmerControl(optimizer = "bobyqa",
                             optCtrl = list(maxfun = 100000))
    ),
    silent = TRUE
  )
  if (inherits(m_b, "try-error")) return(NA)
  
  # Counterfactual datasets
  dat_no_b <- dat_b %>% mutate(cat = "No")
  
  # risks
  p_obs_b <- predict(m_b, newdata = dat_b, type = "response", re.form = NA)
  p_no_b  <- predict(m_b, newdata = dat_no_b, type = "response", re.form = NA)
  
  risk_obs_b <- mean(p_obs_b)
  risk_no_b  <- mean(p_no_b)
  
  (risk_obs_b - risk_no_b) / risk_obs_b
})

boot_PAF <- boot_PAF[!is.na(boot_PAF)]
PAF_CI <- quantile(boot_PAF, c(0.025, 0.975))
PAF_CI

# The adjusted association between cat ownership and Toxoplasma gondii IgG seropositivity demonstrated moderate robustness to unmeasured confounding.
# The E-value for the point estimate (OR = 1.93) was 2.12, indicating that an unmeasured confounder associated with both cat ownership and seropositivity by more than two-fold each would be required to fully explain away the observed effect.
# he E-value for the lower confidence bound was 1.25, suggesting that a relatively modest unmeasured confounder could shift the confidence interval to include the null.
# Thus, while the point estimate is unlikely to be fully attributable to weak confounding, moderate unmeasured confounding remains a plausible explanation for the lower-bound effect.

# The association between cat ownership and T. gondii seropositivity showed moderate robustness to unmeasured confounding.
# The point-estimate E-value was 2.12, indicating that a confounder associated with both cat ownership and seropositivity by more than two-fold would be required to fully explain away the observed effect.
# The E-value for the lower confidence bound was 1.25, suggesting that modest unmeasured confounding could attenuate the statistical significance of the association but would not be sufficient to explain away the point estimate.
# Overall, the analysis indicates a plausible association that is not highly fragile but remains sensitive to moderate unmeasured confounding, which is typical for environmental and behavioural exposures.