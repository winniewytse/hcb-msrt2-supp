#### Reanalyze STAR for Planning A New Multisite Randomized Trial ####

# Setup --------------------------------------------------------------------

# Load required libraries
library(haven)
library(tidyverse)
library(lme4)
library(brms)
library(tidybayes)
library(bootmlm)
library(hcbr)

# Uncomment below to install bootmlm
# install.packages("remotes")
# remotes::install_github("marklhc/bootmlm")

# Data are available on 
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/10766

# Read data
dat <- read_spss(here::here(
  "analysis/exploration/STAR/data/STAR_Students.sav"
)) %>%
  as.data.frame() %>%
  # Combine regular class with and without aides
  mutate(g1class = recode(as.factor(g1classtype), 
                          `1` = 1, `2` = 2, `3` = 2))


# Helper functions

# Calculates expected power/assurance level with posterior draws
msrt2_pow_brm <- function(draws_delta, draws_rho, draws_omega, 
                          J, n, P = .5, K = 0, rsq1 = 0, rsq2 = 0, 
                          alpha = .05, pow = .8, goal = c("ep", "al")) {
  df <- J - 1
  cv <- qt(1 - alpha / 2, df)
  draws_ncp <- draws_delta * sqrt(
    P * (1 - P) * J * n /
      (draws_rho * draws_omega * (1 - rsq2) * P * (1 - P) * n +
         (1 - draws_rho) * (1 - rsq1))
  )
  draws_pow <- pt(cv, df = df, ncp = draws_ncp, lower.tail = FALSE) +
    pt(-cv, df = df, ncp = draws_ncp, lower.tail = TRUE)
  if (goal == "ep") {
    # expected power
    mean(draws_pow)
  } else if (goal == "al") {
    mean(draws_pow > pow)
  }
}

# Calculate assurance level of precision with posterior draws
msrt2_prec_brm <- function(draws_delta, draws_rho, draws_omega, 
                           J, n, P = .5, K = 0, rsq1 = 0, rsq2 = 0, 
                           alpha = .05, precision = .15) {
  df <- J - K - 1
  cv <- qt(1 - alpha / 2, df)
  precision_draws <- cv * sqrt(
    (draws_rho * draws_omega * (1 - rsq2) * 
       P * (1 - P) * n + (1 - draws_rho) * (1 - rsq1)) / 
      (P * (1 - P) * J * n)
  )
  # assurance level of preicison
  mean(precision_draws < precision)
}

# Classical Multilevel Analysis --------------------------------------------

mod_lme <- lmer(
  g1tmathss ~ g1class + (g1class | g1schid), 
  data = dat
)
sum_lme <- summary(mod_lme)

# Number and size of clusters
J <- sum_lme$ngrps
n <- sum_lme$devcomp$dims["N"] / J # average cluster size

# Random effects and residuals
vc <- as.data.frame(VarCorr(mod_lme))
tau_00 <- vc$sdcor[1]
tau_11 <- vc$sdcor[2]
sigma <- vc$sdcor[4]

# Estimated rho (intraclass correlation)
rho_est_c <- tau_00^2 / (tau_00^2 + sigma^2)
# Estimated standard error of rho
rho <- function(mod) {
  vc <- as.data.frame(VarCorr(mod))
  tau00 <- vc$sdcor[1]
  sigma <- vc$sdcor[4]
  return(tau00^2 / (tau00^2 + sigma^2))
}
boot_rho <- bootMer(mod_lme, rho, nsim = 1000)
rho_se_c <- sd(boot_rho$t)

# Estimated omega (treatment effect heterogeneity)
omega_est_c <- tau_11^2 / tau_00^2
# Estimated standard error of omega (with delta method)
omega <- function(mod) {
  vc <- as.data.frame(VarCorr(mod))
  tau00 <- vc$sdcor[1]
  tau11 <- vc$sdcor[2]
  return(tau11^2 / tau00^2)
}
boot_omega <- bootMer(mod_lme, omega, nsim = 1000)
omega_se_c <- sd(boot_omega$t)

# Estimated delta (treatment effect size)
gamma10 <- abs(fixef(mod_lme)[2])
d_est_c <- gamma10 / sqrt(tau_00^2 + sigma^2)
# Estimated SE of delta
d_se_c <- sqrt((rho_est_c * omega_est_c * n + 4 * (1 - rho_est_c)) / (J * n))


# Bayesian Multilevel Analysis ---------------------------------------------

# Run a random slope model in brms
mod_brm <- brm(
  g1tmathss ~ g1class + (g1class | g1schid),
  seed = 123,
  control = list(adapt_delta = .95),
  chains = 4,
  iter = 8000,
  data = dat
)
summary(mod_brm)

# Obtain posterior draws of tau's and sigma
sd_brm <- VarCorr(mod_brm, summary = FALSE)
draws_tau00 <- sd_brm$g1schid$sd[, "Intercept"] # tau00
draws_tau11 <- sd_brm$g1schid$sd[, "g1class"] # tau11
draws_tau10 <- sd_brm$g1schid$cor[, , "Intercept"][, "g1class"] # tau10
draws_sigma <- sd_brm$residual__$sd[, 1] # sigma

# Compute draws for rho (intraclass correlation)
draws_rho <- draws_tau00^2 / (draws_tau00^2 + draws_sigma^2)

# Compute draws for omega (treatment effect heterogeneity)
draws_omega <- draws_tau11^2 / draws_tau00^2

# Obtain posterior draws of gamma10 (treatment effect)
draws_gamma10 <- spread_draws(mod_brm, b_g1class)$b_g1class

# Compute draws for delta (treatment effect size)
draws_delta <- draws_gamma10 / sqrt(draws_tau00^2 + draws_sigma^2)

# Posterior mean and standard deviation of delta
(d_est_b <- abs(mean(draws_delta)))
(d_sd_b <- sd(draws_delta))

# Posterior mode and standard deviation of rho
dens_rho <- density(draws_rho)
(rho_est_b <- dens_rho$x[which.max(dens_rho$y)])
(rho_sd_b <- sd(draws_rho))

# Posterior mode and standard deviation of omega
(dens_omega <- density(draws_omega))
(omega_est_b <- dens_omega$x[which.max(dens_omega$y)])
(omega_sd_b = sd(draws_omega))


# Sample Size Planning ------------------------------------------------------

msrt2_pow_brm(draws_delta, draws_rho, draws_omega, J = 45, n = 20, goal = "al")
al_msrt2(J = 45, n = 20, d_est = d_est_b, d_sd = d_sd_b, 
         rho_est = rho_est_b, rho_sd = rho_sd_b, 
         omega_est = omega_est_b, omega_sd = omega_sd_b)

msrt2_prec_brm(draws_delta, draws_rho, draws_omega, J = 63, n = 20, precision = .12)
apr_msrt2(J = 63, n = 20, rho = rho_est_b, rho_sd = rho_sd_b, 
          omega = omega_est_b, omega_sd = omega_sd_b, precision = .12)

