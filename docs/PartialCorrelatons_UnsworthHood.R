# Bayesian partial correlation between anti-saccade performance and WMC, controlling for prosaccade performance, in Unsworth et al (2023) and Hood et al (2022)

rm(list=ls())
graphics.off()

# remotes::install_github("GidonFrischkorn/AntiSaccade")
library(AntiSaccade)
library(brms)
library(cmdstanr)


########## Unsworth 2023 data ################

data <- Unsworth_2023_E2_binomial

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

x11()
layout(matrix(1:5,2,3, byrow=T))
for (c in 1:5) {
  hist(Corr.Unsw[, c], main=colnames(Corr.Unsw)[c])
}

corr_means <- apply(Corr.Unsw, 2, mean)
qunatiles <- apply(Corr.Unsw, 2, quantile, probs=c(0.025, 0.975))

# difference in correlation to WMC between anti- and prosaccade performance
quantile(Corr.Unsw[,2]-Corr.Unsw[,1], probs = c(0.025, 0.5, 0.975))

########## Hood et al data ################

data <- Hood_2022_behavioral
data <- subset(data, !is.na(data$z_ospan))

options (mc.cores=parallel::detectCores()) # Run on multiple cores
nChains <- min(4, parallel::detectCores())
nIter <- 4000

fixefPrior <- set_prior("logistic(0,0.5)", class="b")
ranefPrior <- set_prior("gamma(1,0.01)", class="sd")

M.Anti <- brm(n_corr|trials(n_trials) ~ 1 + (1||sub), data=subset(data, task=="anti"), prior = c(ranefPrior), family=binomial(link="logit"),
              iter = nIter, warmup = 1000, chains = nChains, save_pars=save_pars(all=T), backend='cmdstanr')

M.Pro <- brm(n_corr|trials(n_trials) ~ 1 + (1||sub), data=subset(data, task=="pro"), prior = c(ranefPrior), family=binomial(link="logit"),
             iter = nIter, warmup = 1000, chains = nChains, save_pars=save_pars(all=T), backend='cmdstanr')


post.Anti <- as_draws_df(M.Anti, variable="^r", regex=T)  # each row = one posterior draw; each column = one subject
post.Anti <- post.Anti[, 1:129] # cut off the last 3 variables, which are not posteriors of individual subjects

post.Pro <- as_draws_df(M.Pro, variable="^r", regex=T)
post.Pro <- post.Pro[, 1:129] # cut off the last 3 variables, which are not posteriors of individual subjects
WMC <- unlist(subset(data, task=="anti")[, "z_ospan"])

Corr.Hood <- matrix(NA, dim(post.Anti)[1], 5)
colnames(Corr.Hood) <- c("Pro.WMC", "Anti.WMC", "Pro.Anti", "Partial.Anti.WMC","Partial.Pro.WMC")

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
  Corr.Hood[sample, 1] <- cor(posterior.Pro, WMC)
  Corr.Hood[sample, 2] <- cor(posterior.Anti, WMC)
  Corr.Hood[sample, 3] <- cor(posterior.Pro, posterior.Anti)
  Corr.Hood[sample, 4] <- cor(Residual.Anti, Residual.WMC.pro)
  Corr.Hood[sample, 5] <- cor(Residual.pro, Residual.WMC.anti)
}

x11()
layout(matrix(1:4,2,2, byrow=T))
for (c in 1:5) {
  hist(Corr.Hood[, c], main=colnames(Corr.Hood)[c])
}

corr_means <- apply(Corr.Hood, 2, mean)
qunatiles <- apply(Corr.Hood, 2, quantile, probs=c(0.025, 0.975))
corr_means
qunatiles

# difference in correlation to WMC between anti- and prosaccade performance
quantile(Corr.Hood[,2]-Corr.Hood[,1], probs = c(0.025, 0.5, 0.975))
