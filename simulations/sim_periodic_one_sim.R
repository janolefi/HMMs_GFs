####### Simulation experiment - additive GP trend #######

# Installing required packages
# Newest CRAN versions of RTMB (1.9) and LaMa (2.1.0) are required!
# install.packages(c("RTMB", "fmesher", "LaMa", "parallel"))


library(parallel) # for parallelising fits
library(LaMa)     # for HMM functions
library(fmesher)  # to construct mesh and SDPE matrices
library(RTMBdist) # for rgmrf()


### Colors
source("utils.R")


# Simulate one true field
nObs <- 1000
mesh <- fm_mesh_1d(1:nObs)
spde <- fm_fem(mesh)

sigma0 <- 0.03
rho0 <- 20
kappa0 <- sqrt(8)/rho0
tau0 <- 1/(sigma0 * sqrt(4*pi) * kappa0)
omega0 <- 0.99


# True parameters used for simulation
true_par <- list(
  mu = c(-1, 3),
  sigma = c(1, 2),
  beta0 = rep(qlogis(0.05), 2),
  tau = tau0,
  kappa = kappa0,
  omega = omega0
)


### Function that simulates one data set
sim_data <- function(p, # parameter list
                     nObs = 1000,
                     spde = spde
) {
  n <- nObs
  Gamma <- tpm(p$beta0)
  delta <- stationary(Gamma)

  # Sample Simple Markov chain
  s <- rep(NA, n)
  s[1] <- sample(1:2, size = 1, prob = delta) # first state
  for(t in 2:n) {
    s[t] <- sample(1:2, 1, prob = Gamma[s[t-1],])
  }

  # Sample latent observatios given chain
  x <- rnorm(n, p$mu[s], p$sigma[s])

  # Simulate field
  Q <- p$tau^2 * (p$kappa^4 * spde$c0 + 2 * cos(pi*p$omega) * p$kappa^2 * spde$g1 + spde$g2)
  true_field <- as.numeric(rgmrf(1, 0, Q))

  # Add true field
  y <- x + true_field

  data.frame(s = s, x = x, y = y, f = true_field)
}

# likelihood function
jnll <- function(par) {
  getAll(par, dat, warn = FALSE)
  REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  Gamma <- tpm(beta0)
  delta <- stationary(Gamma)
  x <- y - w # detrending observed sequence with field
  REPORT(w)
  lallprobs <- matrix(0, length(y), 2)
  for(j in 1:2) {
    lallprobs[,j] <- dnorm(x, mu[j], sigma[j], log = TRUE)
  }

  # HMM likelihood
  nll <- -forward(delta, Gamma, lallprobs, bw = bw, logspace = TRUE)

  # GMRF likelihood
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  cos_pi_omega <- 2*plogis(u)-1; omega <- acos(cos_pi_omega) / pi; REPORT(omega)
  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2*cos_pi_omega*kappa_sq * g1 + g2)
  nll <- nll - dgmrf(w, 0, Q, log = TRUE)

  nll
}

set.seed(1) # to get the same sample and same plot

# Simulate data
data <- sim_data(true_par, spde = spde)

dat <- list(
  y = data$y,
  c0 = spde$c0, g1 = spde$g1, g2 = spde$g2,
  bw = 15
)

# initialise with true values
par0 <- list(
  mu = true_par$mu,
  log_sigma = log(true_par$sigma),
  beta0 = true_par$beta0,
  w = rep(0, nrow(spde$c0)),
  log_tau_sq = log(true_par$tau^2),
  log_kappa_sq = log(true_par$kappa^2),
  u = qlogis((1 + cos(pi * true_par$omega)) / 2)
)

# fit model
obj <- MakeADFun(jnll, par0, random = "w")
opt <- nlminb(obj$par, obj$fn, obj$gr)
mod <- report(obj)
states <- viterbi(mod = mod)

# pdf("./figs/sim_p_sample.pdf", width = 7, height = 5)
plot(data$y, pch = 16, col = adjustcolor(c("#00000060", colorCB[3])[states], alpha.f = 0.4), bty = "n",
     ylab = "Y")
lines(data$f, lwd = 3, col = colorCB[6])
lines(mod$w, lwd = 3, col = colorCB[2], lty = 2)
legend("topright", legend = c("State 1", "State 2", "True trend", "Predicted Trend"),
       pch = c(16, 16, NA, NA), lwd = c(NA, NA, 3, 3), lty = c(NA, NA, 1, 2),
       bg = "white", box.col = "white",
       col = c(adjustcolor(c("#00000060", colorCB[3]), alpha.f = 0.4), colorCB[6], colorCB[2]))
# dev.off()
