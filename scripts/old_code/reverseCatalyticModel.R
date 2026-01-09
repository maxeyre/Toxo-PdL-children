############################################################################
## Reverse catalytic model 
## E Rees
## Due to data sharing constraints, this code used data binned into 5 year age 
## categories rather than the individual data
## Therefore, the results are slightly different to those that appear in the paper
############################################################################

## Load in packages
require(tidyverse)
require(rjags)
require(binom)
require(varhandle)
require(loo)
require(MCMCvis)

################################################################################
### Read in data and source functions
################################################################################

# read in model functions
source("Catalytic_model/catalyticModelfunctions.R")

# Reading in 5-year age grouped data
seroDatGrouped <- read_csv("Data/toxoSeroDatGrouped.csv")

# rename to Ellie's naming convention
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

###################
### define model
###################

jcode <- "model{ 
for (i in 1:length(n)){

n.pos[i] ~ dbin(seropos_est[i],n[i]) 

#reverse catalytic model
seropos_est[i] = (lambda / (lambda + delta)) * (1-exp(-(lambda+delta)*age[i])) 

## #calculate likelihood (used for WAIC)
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],n[i]) 
}
# Define priors
lambda ~ dunif(0,0.5)
delta ~ dunif(0,10)
}"


##### MALE #####
seroDatGrouped <- seroDatGrouped_all %>% filter(sex=="Male")

fiveYearAgeTotals <- seroDatGrouped$total

################################################################################
# Running the model
################################################################################

## Number of model iterations
mcmc.length=20000 

## Specify my data
jdat <- list(n.pos= seroDatGrouped$seropositive,
             age=seroDatGrouped$midpoint,
             n=seroDatGrouped$total)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=10, n.adapt=5000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","delta","loglik"), n.iter=mcmc.length)

plot(jpos) ## Check convergence of all chains

MCMCsummary(jpos, round = 2) ## Check ESS and Rhat
summary(jpos)

#convert mcmc.list to a matrix
mcmcMatrix <- as.matrix(jpos)

# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

#remove burn in
jpos <- window(jpos, start=10000)
plot(jpos)
MCMCsummary(jpos, round = 2)
summary(jpos)

# calculate the DIC 
dic.samples(jmod, n.iter = mcmc.length)

## extract log liklihood
logLik <- mcmcMatrix[,3:7]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)

# waic = 32.8, looic = 33.2

# Get point estimates
deltaPointEst_male <- mcmcMatrix[,"delta"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEst_male <- mcmcMatrix[,"lambda"] %>% quantile(probs=c(.5,.025,.975))
AR_male <- (1-exp(-lambdaPointEst_male))*100

################################################################################
## Get data in correct format for plotting
## Sample from mcmc chains for credible intervals
## Add sampling uncertainty
################################################################################

ageVector = 4:18
numberOfSamples = 1000
outDf <- matrix(NA,nrow=numberOfSamples, ncol = length(ageVector))
midpoints <- seroDatGrouped$midpoint

## 1. Sample from mcmc chain to get 95% credible intervals (model uncertainty)

randomlySampleMcmcChainModelUncertainty <- mcmcModelRandomSampler(1000,mcmcMatrix,ageVector)
# for each column in the matrix get quantiles by age
ageQuantilesModelUncertainty <- ageQuantiles(randomlySampleMcmcChainModelUncertainty)

## Create a df with model uncertainty
df_upperLower = data.frame(
  midpoint = ageVector,
  mean = (lambdaPointEst[1] / (lambdaPointEst[1]+deltaPointEst[1])) *(1 - exp(-ageVector*(lambdaPointEst[1]+deltaPointEst[1]))),
  upper = ageQuantilesModelUncertainty[,3],
  lower = ageQuantilesModelUncertainty[,2]
)

## 2. Sample uncertainty - accounts for the sample size of the underlying data
randomlySampleMcmcChain <- mcmcRandomSampler(1000,mcmcMatrix,midpoints,fiveYearAgeTotals)
ageQuantilesSamplingUncertainty <- ageQuantiles(randomlySampleMcmcChain)

## Create a df with sample uncertainty
df_sampling = data.frame(
  midpoint = seroDatGrouped$midpoint,
  mean = (lambdaPointEst[1] / (lambdaPointEst[1]+deltaPointEst[1])) *(1 - exp(-seroDatGrouped$midpoint*(lambdaPointEst[1]+deltaPointEst[1]))),
  upper = ageQuantilesSamplingUncertainty[,3],
  lower = ageQuantilesSamplingUncertainty[,2]
)

## Save the dataframes for plotting
saveRDS(df_sampling,"Catalytic_model/samplingUncertaintyRevCat_5_male.rds")
saveRDS(df_upperLower, "Catalytic_model/modelUncertaintyRevCat_5_male.rds")

##  Plot figures ####

## Read in dataframes for catalytic and reverse catalytic models (sampling and model uncertainty)
samplingUncertaintyCat <- readRDS("Catalytic_model/samplingUncertaintyRevCat_5_male.rds")
modelUncertaintyCat <- readRDS("Catalytic_model/modelUncertaintyRevCat_5_male.rds" )

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

lambdaEst_male  <- round(lambdaPointEst_male , 3)
lambdaEst_male <- format(lambdaEst_male , nsmall = 2) # to keep trailing zeroes

AR_male <- round(AR_male, 2)
AR_male<- format(AR_male, nsmall = 2) # to keep trailing zeroes

deltaPointEst_male <- round(deltaPointEst_male,3)
deltaPointEst_male <- format(deltaPointEst_male , nsmall = 2)

label_male <- paste0("FOI = ", lambdaEst_male[1], " (95% CrI: ", lambdaEst_male[2], ", ",lambdaEst_male[3],
                     ") \nAR = ", AR_male[1], "% (95% CrI: ", AR_male[2], "%, ",AR_male[3],"%)",
                     "\n",expression(omega)," = ", deltaPointEst_male[1], " (95% CrI: ", deltaPointEst_male[2],", ", deltaPointEst_male[3],")")

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

