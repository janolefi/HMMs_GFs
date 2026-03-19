####### Reproducing the analysis of Esquivel et al. (2025) #######
# https://iopscience.iop.org/article/10.3847/1538-4357/ad95f6/meta

# installing required packages
# install.packages(c("RTMB", "fmesher", "LaMa", "RTMBdist"))

# devtools::install_github("janolefi/LaMa")

library(LaMa)       # for HMM functions
library(RTMBdist)   # for ExGaussian distribution
library(fmesher)    # for mesh and FEM matrices
library(scales)     # for semi-transparent colors


### Colors for plotting
color <- c("#00000070", "red", "orange")


### Changing AD setting inside RTMB for faster calculations
TapeConfig(matmul = "plain")



# Data and preprocessing --------------------------------------------------

### Loading data
# Brightness measurements of TESS2018206045859 at 2 min cadence
data <- read.csv("./data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[,c("TIME","PDCSAP_FLUX")]
colnames(data) <- c("time", "y")

data <- data[-1, ] # exclude first observation because NA

# Interpolate time variable
time_na <- is.na(data$time)
data$time[time_na] <- approx(which(!time_na), data$time[!time_na], xout = which(time_na))$y

# Time series has NAs between these two indices
end_before_gap <- 9524
start_after_gap <- 10341

# Exclude this section and create trackID variable to handle tracks separately in likelihood
data <- data[-((end_before_gap + 1):(start_after_gap-1)), ]
data$trackID <- 2
data$trackID[1:end_before_gap] <- 1

# Centering data as in Esquivel et al. (2025)
data$y <- scale(data$y, scale = FALSE)



# Creating meshs and finite element matrices ------------------------------

mesh1 <- fm_mesh_1d(data$time[data$trackID == 1])
spde1 <- fm_fem(mesh1)

mesh2 <- fm_mesh_1d(data$time[data$trackID == 2])
spde2 <- fm_fem(mesh2)



# Simple model without state switching ------------------------------------

# likelihood function
nll0 <- function(par) {
  getAll(par, dat) # unpacks lists

  sigma <- exp(log_sigma); REPORT(sigma) # error variance
  f1 <- w1 + mu # adding mean to first field
  f2 <- w2 + mu # adding mean to second field
  f <- c(f1, f2); REPORT(f) # total smooth
  nll <- -sum(dnorm(y, f, sigma, log = TRUE), na.rm = TRUE) # observation likelihood

  ### GP likelihood
  # parameter transformations
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  cos_pi_omega <- 2 * plogis(u) - 1; omega <- acos(cos_pi_omega) / pi; REPORT(omega)

  # building precision matrices
  Q1 <- tau_sq * (kappa_sq*kappa_sq * spde1$c0 + 2 * cos_pi_omega * kappa_sq * spde1$g1 + spde1$g2)
  Q2 <- tau_sq * (kappa_sq*kappa_sq * spde2$c0 + 2 * cos_pi_omega * kappa_sq * spde2$g1 + spde2$g2)

  # evaluating densities
  nll <- nll - dgmrf(w1, 0, Q1, log = TRUE) # GP likelihood 1
  nll <- nll - dgmrf(w2, 0, Q2, log = TRUE) # GP likelihood 2

  nll
}

# Initial parameter list
par <- list(
  log_sigma = log(7),
  log_tau_sq = log(0.005^2),
  log_kappa_sq = log(30^2),
  u = -5,
  mu = -0.5,
  w1 = numeric(nrow(spde1$c0)),
  w2 = numeric(nrow(spde2$c0))
)

# Data list
dat <- list(
  y = data$y,
  spde1 = spde1,
  spde2 = spde2
)

# Constructing marginal likelihood via automatic Laplace
obj0 <- MakeADFun(nll0, par, random = c("w1", "w2"))

# Optimising marginal likelihood
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)

# Getting reported quantities
mod0 <- report(obj0)
mod0$tau
mod0$kappa
mod0$omega

# Sample from posterior distribution for parameter uncertainty
samples0 <- MCreport(obj0, report = TRUE)
# quantiles for smooth f
qf <- t(apply(simplify2array(samples0$f), 1,
              quantile, probs = c(0.025, 0.975)))

### visualising results
idx <- 6000:8000
plot(data$time[idx], data$y[idx], pch = 16, col = color[1],
     xlab = "Time (days)", ylab = "Flux", bty = "n")
lines(data$time[idx], mod0$f[idx], lwd = 2, col = "plum")
polygon(c(data$time[idx], rev(data$time[idx])),
        c(qf[idx, 1], rev(qf[idx,2])), col = alpha("plum", 0.6), border = NA)

# The predcited smooth is erratic whenever there is a flare -> 2 stage methods have problems



# HMM analysis ------------------------------------------------------------

# HMM likelihood function
jnll <- function(par) {
  getAll(par, dat)

  #### state process ####
  # restricted tpm
  Gamma <- diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  # 2 estimated initial distributions
  Delta <- cbind(1, exp(logit_Delta))
  Delta <- Delta / rowSums(Delta)

  #### state-dependent process ####
  # parameter transformations
  sigma <- exp(log_sigma); REPORT(sigma)
  r <- plogis(logit_r); REPORT(r)
  lambda <- exp(log_lambda); REPORT(lambda)

  # quasi-periodic trend
  f1 <- w1 + mu # smooth 1
  f2 <- w2 + mu # smooth 2
  f <- c(f1, f2); REPORT(f) # quasi-periodic smooth
  z <- y - f; REPORT(mu)

  # latent z's where observations are mising -> also integrated out via Laplace
  z[is.na(y)] <- z.star; REPORT(z)

  # evaluating state-dependent densities
  n <- length(z)
  idx <- list(2:(trackInd[2]-1), (trackInd[2]+1):n) # idx of the 2 time series
  lallprobs <- matrix(0, n, 3)
  for(i in 1:2) { # loop over time series
    ind <- idx[[i]]
    # quiet (Q):
    lallprobs[ind,1] <- dnorm(z[ind], 0, sigma, log = TRUE)
    # firing (F):
    lallprobs[ind,2] <- dexgauss(z[ind], z[ind-1], sigma, lambda, log = TRUE)
    # decaying (D):
    lallprobs[ind,3] <- dnorm(z[ind], r * z[ind-1], sigma, log = TRUE)
  }

  ### HMM likelihood
  nll <- -forward(Delta, Gamma, lallprobs,
                  trackID = trackID, # separate tracks
                  bw = bw,           # passing this enables sparsity in Hessian
                  logspace = TRUE)   # logspace calculations more stable

  ### GP likelihood
  # parameter transformations
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  cos_pi_omega <- 2 * plogis(u) - 1; omega <- acos(cos_pi_omega)/pi; REPORT(omega)

  # building precision matrices
  Q1 <- tau_sq * (kappa_sq*kappa_sq * spde1$c0 + 2 * cos_pi_omega * kappa_sq * spde1$g1 + spde1$g2)
  Q2 <- tau_sq * (kappa_sq*kappa_sq * spde2$c0 + 2 * cos_pi_omega * kappa_sq * spde2$g1 + spde2$g2)

  # evaluating densities
  nll <- nll - dgmrf(w1, 0, Q1, log = TRUE) # GP likelihood 1
  nll <- nll - dgmrf(w2, 0, Q2, log = TRUE) # GP likelihood 2

  nll
}


# Initial parameter list
par <- list(
  eta = rep(-2, 4),
  logit_Delta = matrix(0, 2, 2),
  log_sigma = log(7),
  logit_r = qlogis(0.8),
  log_lambda = log(0.025),
  log_tau_sq = log(0.005^2),
  log_kappa_sq = log(30^2),
  w1 = rep(0, nrow(spde1$c0)),
  w2 = rep(0, nrow(spde2$c0)),
  u = -5,
  mu = -0.5,
  z.star = rnorm(sum(is.na(data$y)), 0, 3)
)

# Data and hyperparameter list
dat <- list(
  y = data$y,
  spde1 = spde1, # contains finite element matrices for f1
  spde2 = spde2, # contains finite element matrices for f2
  bw = 20,       # bandwidth parameter enabling sparse Hessian
  trackID = data$trackID,
  trackInd = calc_trackInd(data$trackID)
)

t1 <- Sys.time()
# Constructing Laplace marginal likelihood
obj <- MakeADFun(jnll, par, random = c("w1", "w2", "z.star"))

# Optimising
system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr) # inital inner Hessian evaluation takes time
)
Sys.time() - t1

# Reporting quantities of interest
mod <- report(obj)
mod$kappa
mod$tau
mod$omega

# Sampling from posterior distribution for uncertainty
samples <- MCreport(obj, report = TRUE) # takes some time

# Computing 0.025 and 0.975 (pointwise) quantiles of f
f_quantiles <- t(apply(simplify2array(samples$f), 1,
                       quantile, probs = c(0.025, 0.975)))


### State decoding
# Hard decoding: Viterbi algorithm
states <- viterbi(mod = mod)
# Soft decoding: Pr(S_t = j | y_1, ..., y_T)
stateprobs <- stateprobs(mod = mod)

## Identify flare events
flare <- states == 2
flare_event <- states != 1
rle_f <- rle(flare_event)
rle_f$lengths[1:18*2] # lengths
round(mean(rle_f$lengths[1:18*2]) * 2, 1) # on average 17.4 minutes

rle(flare)
# 18 flaring events in total



# Plot result - main figure from the paper --------------------------------

# Choose which section to plot here
idx <- 7030:7350

# pdf("./figs/flare_result.pdf", width = 7, height = 4.5)

# Stacked barplot of state probabilities
layout(matrix(1:2, ncol = 1), heights = c(0.6, 1))

par(mar = c(1, 4, 1, 2))
barplot(t(stateprobs[idx,]), col = c("black", "red", "orange"), border = "white", space = 0,
        ylab = "State probability",
        main = "", yaxt = "n")
axis(2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1))

# Decoded time series
par(mar = c(5, 4, 0.5, 2), xpd = NA)
plot(data$time[idx], data$y[idx], col = color[states[idx]], pch = 16,
     xlab = "Time (days)", ylab = "Flux", bty = "n")
# 95 % CI
polygon(data$time[c(idx, rev(idx))],
        c(f_quantiles[idx,1], f_quantiles[rev(idx),2]),
        col = alpha("plum", 0.35), border = NA)
lines(data$time[idx], mod$f[idx], lwd = 2, col = "plum") # line width is wider thatn CI

# Legend
legend(x = 1335.415, y = 190,
       legend = c("Quiet", "Firing", "Decaying", "Trend"),
       pch = c(rep(16, 3), NA),
       lwd = c(rep(NA, 3), 2),
       col = c(color, "plum"), bty = "n")

# dev.off()



# Plot full time series (Appendix) ----------------------------------------

ch <- calc_trackInd(data$trackID)[2]

# pdf("./figs/flare_result_full.pdf", width = 8, height = 4)

# Decoded time series
par(mfrow = c(1,1), mar = c(5, 4, 0.5, 2))
plot(data$time, data$y, col = color[states], pch = 16, cex = 0.6,
     xlab = "Time (days)", ylab = "Flux", bty = "n", ylim = c(-30, 160))
lines(data$time[1:(ch-1)], mod$f[1:(ch-1)], lwd = 1.5, col = "plum")
lines(data$time[ch:nrow(data)], mod$f[ch:nrow(data)], lwd = 1.5, col = "plum")

rect(data$time[idx[1]-10], -30, data$time[idx[length(idx)]+10], 160,
     lty = 2, border = "gray")

# Legend
legend("topright",
       legend = c("Quiet", "Firing", "Decaying", "Trend"),
       pch = c(rep(16, 3), NA),
       lwd = c(rep(NA, 3), 1.5),
       cex = 0.8, bty = "n",
       col = c(color, "plum"), box.col = "gray")

# dev.off()

