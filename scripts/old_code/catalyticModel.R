############################################################################
##  Catalytic model 
## M Eyre adapted from E Rees
############################################################################

## Read in packages
require(tidyverse)
require(rjags)
require(binom)
require(varhandle)
require(MCMCvis)
require(loo)
require(cowplot)


################################################################################
### Read in data and source functions
################################################################################

# read in model functions
source("Catalytic_model/catalyticModelfunctions.R")

# Reading in 5-year age grouped data
seroDatGrouped <- read_csv("Data/toxoSeroDatGrouped.csv")

# rename to new naming convention
seroDatGrouped <- seroDatGrouped %>%
  select(age_name, 
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

######################################
### define model catalytic model
######################################

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

##### MALE #####
seroDatGrouped <- seroDatGrouped_all %>% filter(sex=="Male")

fiveYearAgeTotals <- seroDatGrouped$total


################################################################################
# Running the model ####
################################################################################

## Number of model iterations
mcmc.length=10000

## Specify my data
jdat <- list(n.pos = seroDatGrouped$seropositive,
             age = seroDatGrouped$midpoint,
             n=seroDatGrouped$total)

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

# waic = 30.3, looic = 30.3 ---> this is lower than the reverse catalytic model


################################################################################
## Get data in correct format for plotting
## Sample from mcmc chains for credible intervals
## Add binomial sampling uncertainty
################################################################################

midpoints <- seroDatGrouped$midpoint
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
AR_male <- (1-exp(-lambdaEst))*100

## Save the dataframes for plotting (along with reverse catalytic model)
saveRDS(df_sampling,"Catalytic_model/samplingUncertaintyCat_5_male.rds")
saveRDS(df_mod, "Catalytic_model/modelUncertaintyCat_5_male.rds" )

##  Plot figures ####

## Read in dataframes for catalytic and reverse catalytic models (sampling and model uncertainty)
samplingUncertaintyCat <- readRDS("Catalytic_model/samplingUncertaintyCat_5_male.rds")
modelUncertaintyCat <- readRDS("Catalytic_model/modelUncertaintyCat_5_male.rds" )

## Read in 5 year age grouped data
seroDat <- seroDatGrouped ## seroprevalence binned into 5 year age groups (used for plotting)

## Catalytic plot
# Extend the data with interpolated values for x = 4 and x = 18

data <- samplingUncertaintyCat

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

modelUncertaintyCat_male <- modelUncertaintyCat
seroDat_male <- seroDat

lambdaEst_male <- round(lambdaEst_male, 3)
lambdaEst_male<- format(lambdaEst_male, nsmall = 2) # to keep trailing zeroes

AR_male <- round(AR_male, 2)
AR_male<- format(AR_male, nsmall = 2) # to keep trailing zeroes

label_male <- paste0("FOI = ", lambdaEst_male[1], " (95% CrI: ", lambdaEst_male[2], ", ",lambdaEst_male[3],
                       ") \nAR = ", AR_male[1], "% (95% CrI: ", AR_male[2], "%, ",AR_male[3],"%)")

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
seroDatGrouped <- seroDatGrouped_all %>% filter(sex=="Female")

fiveYearAgeTotals <- seroDatGrouped$total


################################################################################
# Running the model ####
################################################################################

## Number of model iterations
mcmc.length=10000

## Specify my data
jdat <- list(n.pos = seroDatGrouped$seropositive,
             age = seroDatGrouped$midpoint,
             n=seroDatGrouped$total)

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

midpoints <- seroDatGrouped$midpoint
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
lambdaEst_female <- lambdaEst
AR_female <- (1-exp(-lambdaEst))*100

## Save the dataframes for plotting (along with reverse catalytic model)
saveRDS(df_sampling,"Catalytic_model/samplingUncertaintyCat_5_female.rds")
saveRDS(df_mod, "Catalytic_model/modelUncertaintyCat_5_female.rds" )

### Plot Figures ####

## Read in dataframes for catalytic and reverse catalytic models (sampling and model uncertainty)
samplingUncertaintyCat <- readRDS("Catalytic_model/samplingUncertaintyCat_5_female.rds")
modelUncertaintyCat <- readRDS("Catalytic_model/modelUncertaintyCat_5_female.rds" )

## Read in 5 year age grouped data
seroDat <- seroDatGrouped ## seroprevalence binned into 5 year age groups (used for plotting)

## Catalytic plot
# Extend the data with interpolated values for x = 4 and x = 18

data <- samplingUncertaintyCat

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
seroDat_female <- seroDat

lambdaEst_female <- round(lambdaEst_female, 3)
lambdaEst_female<- format(lambdaEst_female, nsmall = 2) # to keep trailing zeroes

AR_female <- round(AR_female, 2)
AR_female<- format(AR_female, nsmall = 2) # to keep trailing zeroes

label_female <- paste0("FOI = ", lambdaEst_female[1], " (95% CrI: ", lambdaEst_female[2], ", ",lambdaEst_female[3],
                       ") \nAR = ", AR_female[1], "% (95% CrI: ", AR_female[2], "%, ",AR_female[3],"%)")

catPlot_female <- ggplot(modelUncertaintyCat_female, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.3, fill = "#9e2a2b")+
  geom_line()+
  geom_point(data=seroDat_female)+
  geom_linerange(data=seroDat_female) +
  geom_ribbon(data=samplingUncertaintyCat_female, alpha=0.3,fill = "#9e2a2b")+
  scale_y_continuous(breaks=seq(0,1,by=0.2), lim=c(0,.9))+
  xlab("Age (years)") + ylab("Proportion seropositive") +
  scale_x_continuous(breaks=seq(4,18,by=2), lim= c(4,18)) +
  theme_bw() + 
  annotate("text", x = 5, y = 0.815, label = label_female,
                        size = 3.5, color = "black", hjust= 0, vjust=0)
catPlot_female

## Create panel plot
catPlot <- plot_grid(catPlot_male, catPlot_female, labels = c('A', 'B'))

# save plot
ggsave(plot = catPlot, "Catalytic_model/serocatalytic_plot_label.png", width = 10, height = 5, dpi = 300)
