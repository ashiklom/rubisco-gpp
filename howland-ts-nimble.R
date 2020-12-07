library(tidyverse)
library(nimble)
library(ggplot2)

source("funs.R")

dat <- read_csv("data/sites/site_US-Ho1/Ho1-flux.csv") %>%
  mutate(
    ts = as.Date(as.character(TIMESTAMP), format = "%Y%m%d"),
    Rsoil = if_else(Rsoil > 0, Rsoil, NA_real_),
    Rh = if_else(Rh > 0, Rh, NA_real_)
  ) %>%
  filter(!is.na(Rsoil)) %>%
  select(ts, NEE, Rsoil, Rh) %>%
  arrange(ts)

# Nimble
gpp_code <- nimbleCode({
  # Process model -- random walk
  gpp[1] ~ dnorm(0, 0.1)
  for (t in 2:nT) {
    gpp[t] ~ dnorm(gpp[t-1], tau_rw)
  }
  # Data model
  for (t in 1:nT) {
    r[t] <- r_chamber[t] * (1 + Rho / Beta - Rho)
    nep_pred[t] <- gpp_finality[t] - r[t]
    nep[t] ~ dnorm(nep_pred[t], tau_nep)
    gpp_finality[t] ~ dnorm(gpp[t], tau_finality)
  }
  # Priors
  tau_rw ~ dgamma(0.1, 0.1)
  ## tau_finality ~ dgamma(10, 10)
  ## Rho ~ dnorm(rho_mu, sd = rho_sd)
  ## Beta ~ dnorm(Beta_mu, sd = Beta_sd)
  Rho ~ dbeta(rho_a, rho_b)
  Beta ~ dnorm(Beta_a, Beta_b)
})
gpp_data <- list(
  r_chamber = dat$Rsoil,
  nep = dat$NEE,
  ## rho_mu = 0.56, rho_sd = 0.11,
  ## Beta_mu = 0.5, Beta_sd = 0.125,
  rho_a = 15, rho_b = 15,
  Beta_a = 15, Beta_b = 15,
  tau_nep = 1000, tau_finality = 1000
)
samples <- nimbleMCMC(
  code = gpp_code,
  data = gpp_data,
  constants = list(nT = nrow(dat)),
  monitors = c("gpp", "gpp_finality", "nep_pred",
               "Rho", "Beta",
               "tau_rw")
)

samples_gpp <- samples[,grep("gpp\\[", colnames(samples))]
samples_npp <- samples[,grep("nep_pred", colnames(samples))]
gpp_mean <- colMeans(samples_gpp)
nep_mean <- colMeans(samples_npp)
gpp_quant <- apply(samples_gpp, 2, quantile, c(0.025, 0.975))
dat2 <- dat %>%
  mutate(
    gpp_mean = gpp_mean,
    gpp_lo = gpp_quant[1,],
    gpp_hi = gpp_quant[2,]
  )

ggplot(dat2) +
  aes(x = ts, y = gpp_mean, ymin = gpp_lo, ymax = gpp_hi) +
  geom_ribbon(col = NA, fill = "deepskyblue") +
  geom_line()
  ## geom_line(aes(y = nep_mean, color = "NEP_pred")) +
  ## geom_line(aes(y = NEE, color = "NEP_obs"))

hist(samples[, "Beta"])
hist(samples[, "Rho"])

plot(dat$ts, gpp_mean, type = "n")
## polygon(c(dat$ts, rev(dat$ts)), c(gpp_quant[1,], rev(gpp_quant[2,])),
##         col = "deepskyblue")
lines(dat$ts, gpp_mean)

hist(samples[,"gpp_finality[2]"])
samples[, "gpp_finality[2]"]
hist(samples[,"tau_rw"])
