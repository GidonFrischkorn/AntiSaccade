# Bayesian partial correlation between anti-saccade performance and WMC, controlling for prosaccade performance, in simulated data
rm(list=ls())
graphics.off()

# remotes::install_github("GidonFrischkorn/AntiSaccade")
library(AntiSaccade)
library(brms)
library(cmdstanr)
library(ppcor)

help(package="AntiSaccade")

########## Unsworth 2023 data ################

data <- Unsworth_2023_E2_binomial

nSims <- 1000

Correlations <- matrix(NA, nSims, 5)
colnames(Correlations) <- c("Anti.Pro", "Anti.WMC", "Pro.WMC", "Partial.Anti.WMC", "Partial.Pro.WMC")

for (sim in 1:nSims) {
  # simulate n_corr in both pro- and anti-saccade conditions based on the same latent ability (binding strength), partially determined by WMC
  bstrength <- 0.4*data$zWMC[data$Task=="Antisaccade"] + 0.6*rnorm(dim(data)[1]/2, 0, 1)
  data$bstrength[data$Task=="Prosaccade"] <- 0.7*bstrength + 0.3*rnorm(dim(data)[1]/2, 0, 1) + 2  # some of the ability is specific to the pro-/anti- condition. Pro is easier than anti, therefore add a larger constant
  data$bstrength[data$Task=="Antisaccade"] <- 0.7*bstrength + 0.3*rnorm(dim(data)[1]/2, 0, 1) + 0
  data$n_corr <- rbinom(dim(data)[1], size=data$n_trials, prob=inv_logit_scaled(data$bstrength, lb = 1/3))
  data$pC <- data$n_corr/data$n_trials
  df <- data.frame(data$pC[data$Task=="Antisaccade"], data$pC[data$Task=="Prosaccade"], data$zWMC[data$Task=="Prosaccade"])
  names(df) <- c("Anti", "Pro", "WMC")
  Correlations[sim,1] <- cor(df$Anti, df$Pro)
  Correlations[sim,2] <- cor(df$Anti, df$WMC)
  Correlations[sim,3] <- cor(df$Pro, df$WMC)
  PCorr <- pcor(df)
  Correlations[sim,4] <- PCorr$estimate[1,3]
  Correlations[sim,5] <- PCorr$estimate[2,3]
}
print(colMeans(Correlations))
apply(Correlations, 2, quantile, c(0.025, 0.975))

# Bayesian GLM with the last simulated data set

unsworth_formula <- bf(n_corr | trials(n_trials) ~  1 + Task * zWMC + (1 + Task | Sub),
                       family = binomial)

# set appropriate constrasts for the task variable
contrasts(data$Task) <- contr.treatment(2, base = 2)

# set priors
unsworth_prior <- prior(logistic(0,1), class = Intercept) +
  prior(logistic(0, 0.5), class = b)

# fit model including an intercept, main effects of task and WMC, and their interaction
simulation_fit <- brm(
  formula = unsworth_formula,
  data = data,
  family = binomial,
  prior = unsworth_prior,
  sample_prior = TRUE,
  iter = 8000,
  warmup = 1000,
  chains = 8,
  cores = 8,
  control = list(adapt_delta = 0.99),
  backend = "cmdstanr"
)

# print summary
summary(simulation_fit)
hyp_unsworth <- hypothesis(simulation_fit,
                           c(wmc_main = "zWMC = 0",
                             wmc_int = "Task1:zWMC = 0"))
hyp_unsworth

# Bayesian logistic analysis with the last simulated data set

options (mc.cores=parallel::detectCores()) # Run on multiple cores
nChains <- min(4, parallel::detectCores())
nIter <- 4000

fixefPrior <- set_prior("logistic(0,0.5)", class="b")
ranefPrior <- set_prior("gamma(1,0.01)", class="sd")

# models predicting the number of correct responses of each person from their latent ability on the logit scale
M.Anti <- brm(n_corr|trials(n_trials) ~ 1 + (1||Sub), data=subset(data, Task=="Antisaccade"), prior = c(ranefPrior), family=binomial(link="logit"),
            iter = nIter, warmup = 1000, chains = nChains, save_pars=save_pars(all=T), backend='cmdstanr')

M.Pro <- brm(n_corr|trials(n_trials) ~ 1 + (1||Sub), data=subset(data, Task=="Prosaccade"), prior = c(ranefPrior), family=binomial(link="logit"),
            iter = nIter, warmup = 1000, chains = nChains, save_pars=save_pars(all=T), backend='cmdstanr')

# Extract posterior samples of each person's latent ability on the logit scale (i.e., the random effect of the model)
post.Anti <- as_draws_df(M.Anti, variable="^r", regex=T)  # each row = one posterior draw; each column = one subject
post.Anti <- post.Anti[, 1:169] # cut off the last 3 variables, which are not posteriors of individual subjects

post.Pro <- as_draws_df(M.Pro, variable="^r", regex=T)
post.Pro <- post.Pro[, 1:169] # cut off the last 3 variables, which are not posteriors of individual subjects
WMC <- unlist(subset(data, Task=="Antisaccade")[, "zWMC"])

Corr.Unsw <- matrix(NA, dim(post.Anti)[1], 5)
colnames(Corr.Unsw) <- c("Pro.WMC", "Anti.WMC", "Pro.Anti", "Partial.Anti.WMC","Partial.Pro.WMC")

for (sample in 1:dim(post.Anti)[1]) {
  posterior.Pro <- as.numeric(post.Pro[sample,])  # select the vector of posterior estimates from this MCMC sample across the 169 subjects
  posterior.Anti <- as.numeric(post.Anti[sample,])
  regModel_anti <- lm(posterior.Anti ~ posterior.Pro)
  Residual.Anti <- regModel_anti$residuals
  regModel_pro <- lm(posterior.Pro ~ posterior.Anti)
  Residual.pro <- regModel_pro$residuals
  regModel_WMC_pro <- lm(WMC ~ posterior.Pro)
  Residual.WMC.pro <- regModel_WMC_pro$residuals
  regModel_WMC_anti <- lm(WMC ~ posterior.Anti)
  Residual.WMC.anti <- regModel_WMC_anti$residuals
  Corr.Unsw[sample, 1] <- cor(posterior.Pro, WMC)
  Corr.Unsw[sample, 2] <- cor(posterior.Anti, WMC)
  Corr.Unsw[sample, 3] <- cor(posterior.Pro, posterior.Anti)
  Corr.Unsw[sample, 4] <- cor(Residual.Anti, Residual.WMC.pro)
  Corr.Unsw[sample, 5] <- cor(Residual.pro, Residual.WMC.anti)
}

for (c in 1:5) {
  hist(Corr.Unsw[, c], main=colnames(Corr.Unsw)[c])
}

corr_means <- apply(Corr.Unsw, 2, mean)
qunatiles <- apply(Corr.Unsw, 2, quantile, probs=c(0.025, 0.975))
print(corr_means)
print(qunatiles)

