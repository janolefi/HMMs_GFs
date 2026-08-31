####### Reproducing the analysis of Esquivel et al. (2025) #######
# https://iopscience.iop.org/article/10.3847/1538-4357/ad95f6/meta

# Installing required packages
# Versions used:
# - RTMB: 1.9
# - LaMa: 2.1.3
# - RTMBdist: 1.0.6
# - fmesher: 0.8.0
# install.packages(c("RTMB", "LaMa", "RTMBdist", "fmesher"))
# install.packages(c("viridis", "scales", "leaflet"))

library(LaMa)       # for HMM functions
library(RTMBdist)   # for ExGaussian distribution
library(fmesher)    # for mesh and finite element matrices
library(scales)     # for semi-transparent colors

source("utils.R")

### Colors for plotting
color <- c("#00000070", "tomato1", "orange")



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

# The predicted trend is erratic whenever there is a flare
# -> 2 stage methods are problematic

# Clear memory
rm(obj0); rm(samples0)



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
  z <- y - f; REPORT(mu) # integrate over w instead of z since it is a linear transformation

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
  bw = 15,       # bandwidth parameter enabling sparse Hessian
  trackID = data$trackID,
  trackInd = calc_trackInd(data$trackID)
)


# Constructing Laplace marginal likelihood
obj <- MakeADFun(jnll, par, random = c("w1", "w2", "z.star"))

# Profile log-likelihood for different bandwidth parameters
bw_chk <- bandwidth_check(obj)
# pdf("./figs/bw_check_flares.pdf", width = 7, height = 4)
plot(bw_chk$bw, bw_chk$llk, type = "l", bty = "n",
     xlab = "Bandwidth", ylab = "Log-likelihood",
     col = "darkblue", lwd = 2, main = "Bandwidth Check - Stellar flare case study")
# dev.off()
dat$bw <- 30 # suitable bandwidth

# clear memory before large fit
rm(obj)

t1 <- Sys.time()
# Reconstructing objective for accurate timing
obj <- MakeADFun(jnll, par, random = c("w1", "w2", "z.star"))
# Optimising
system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr) # inital inner Hessian evaluation takes time
)
Sys.time() - t1

# Check accurary of the Laplace approximation
l_chk <- laplace_check(obj)
l_chk
loo::psis(l_chk$logw)$diagnostics$pareto_k # is importance sampling reliable here? -> no

# Reporting quantities of interest
mod <- report(obj)
mod$kappa
mod$tau
mod$omega
round(mod$Gamma, 4)

# Sampling from posterior distribution for uncertainty
set.seed(123)
samples <- MCreport(obj, report = TRUE) # takes some time

# Computing 0.025 and 0.975 (pointwise) quantiles of f
f_quantiles <- t(apply(simplify2array(samples$f), 1,
                       quantile, probs = c(0.025, 0.975)))


### State decoding
# Hard decoding: Viterbi algorithm
states <- viterbi(mod = mod)
# Soft decoding: Pr(S_t = j | y_1, ..., y_T)
stateprobs <- stateprobs(mod = mod)


## Flare events across the full data set (firing + decaying merged)
rle_all <- rle(states != 1)                       # TRUE = firing/decaying
n_all   <- sum(rle_all$values)                    # number of flare events
dur_all <- rle_all$lengths[rle_all$values]        # durations, in observations

n_all                                             # 18 flares
round(mean(dur_all) * 2, 1)                       # mean duration in minutes (2-min cadence)

## --- Comparison window (first segment, before the gap) vs. Esquivel et al. ---
comp_idx <- 1:end_before_gap
rle_comp <- rle(states[comp_idx] != 1)
n_comp   <- sum(rle_comp$values)                  # number of flares in the window
dur_comp <- rle_comp$lengths[rle_comp$values]     # durations, in observations

n_comp                                            # 11 flares
round(mean(dur_comp), 1)                          # mean duration in observations




# Plot result - main figure in paper ------------------------------------

# Choose which section to plot here
idx <- 7030:7350

# pdf("./figs/flare_result_new.pdf", width = 7, height = 4.5)
cex_state <- c(0.7, 1.2, 0.9)     # quiet small, firing largest, decaying medium

# Stacked barplot of state probabilities
layout(matrix(1:2, ncol = 1), heights = c(0.6, 1))
par(mar = c(1, 4, 1, 2))
barplot(t(stateprobs[idx,c(1,3,2)]), col = c("black", "orange", "tomato1") , border = "white", space = 0,
        ylab = "State probability", main = "", yaxt = "n")
axis(2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1))

# Decoded time series
par(mar = c(5, 4, 0.5, 2), xpd = NA)
plot(data$time[idx], data$y[idx], cex = cex_state[states[idx]],
     col = color[states[idx]], pch = 16,
     xlab = "Time (days)", ylab = "Flux", bty = "n")
# 95 % CI
polygon(data$time[c(idx, rev(idx))],
        c(f_quantiles[idx,1], f_quantiles[rev(idx),2]),
        col = alpha("plum", 0.35), border = NA)
lines(data$time[idx], mod$f[idx], lwd = 2, col = "plum")
# Legend
legend(x = 1335.415, y = 183,
       legend = c("Quiet", "Firing", "Decaying", "Trend"),
       pch = 16, pt.cex = c(cex_state, NA),
       lwd = c(NA, NA, NA, 2),
       col = c(color, "plum"), bty = "n")
# dev.off()





# Plot full time series (Supplementary) ------------------------------------

ch <- calc_trackInd(data$trackID)[2]

# pdf("./figs/flare_result_full_new.pdf", width = 8, height = 4)

# Decoded time series
par(mfrow = c(1,1), mar = c(5, 4, 0.5, 2))
plot(data$time, data$y, col = color[states], pch = 16,
     xlab = "Time (days)", ylab = "Flux", bty = "n", ylim = c(-30, 160),
     cex = cex_state[states])
lines(data$time[1:(ch-1)], mod$f[1:(ch-1)], lwd = 1.5, col = "plum")
lines(data$time[ch:nrow(data)], mod$f[ch:nrow(data)], lwd = 1.5, col = "plum")

rect(data$time[idx[1]-10], -30, data$time[idx[length(idx)]+10], 160,
     lty = 2, border = "gray")

# Legend
legend("topright",
       legend = c("Quiet", "Firing", "Decaying", "Trend"),
       pch = 16, pt.cex = c(cex_state, NA),
       lwd = c(NA, NA, NA, 2),
       col = c(color, "plum"), bty = "n")

# dev.off()







###########################################################################
# Extra computations ------------------------------------------------------
###########################################################################


# Refitting the model to a subsection for more accurate Laplace check -----
data_small <- data[2:5000,] # first observation is NA, select subsection
mesh <- fm_mesh_1d(data_small$time)
spde <- fm_fem(mesh)

# HMM likelihood function
jnll_small <- function(par) {
  getAll(par, dat)

  #### state process ####
  # restricted tpm
  Gamma <- diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  # estimated initial distributions
  delta <- c(1, exp(logit_delta))
  delta <- delta / sum(delta)

  #### state-dependent process ####
  # parameter transformations
  sigma <- exp(log_sigma); REPORT(sigma)
  r <- plogis(logit_r); REPORT(r)
  lambda <- exp(log_lambda); REPORT(lambda)

  # quasi-periodic trend
  f <- w + mu # smooth 1
  REPORT(f) # quasi-periodic smooth
  z <- y - f; REPORT(mu)

  # latent z's where observations are mising -> also integrated out via Laplace
  z[is.na(y)] <- z.star; REPORT(z)

  # evaluating state-dependent densities
  n <- length(z); idx <- 2:n
  lallprobs <- matrix(0, n, 3)
  # regular measurement error
  lallprobs[idx,1] <- dnorm(z[idx], 0, sigma, log = TRUE)
  # firing
  lallprobs[idx,2] <- dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  # decaying
  lallprobs[idx,3] <- dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)

  ### HMM likelihood
  nll <- -forward(delta, Gamma, lallprobs,
                  bw = bw,           # passing this enables sparsity in Hessian
                  logspace = TRUE)   # logspace calculations more stable

  ### GP likelihood
  # parameter transformations
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  cos_pi_omega <- 2 * plogis(u) - 1; omega <- acos(cos_pi_omega)/pi; REPORT(omega)

  # building precision matrices
  Q <- tau_sq * (kappa_sq*kappa_sq * spde$c0 + 2 * cos_pi_omega * kappa_sq * spde$g1 + spde$g2)

  # evaluating densities
  nll <- nll - dgmrf(w, 0, Q, log = TRUE) # GP likelihood

  nll
}


# Initial parameter list
par <- list(
  eta = rep(-2, 4),
  logit_delta = c(0,0),
  log_sigma = log(7),
  logit_r = qlogis(0.8),
  log_lambda = log(0.025),
  log_tau_sq = log(0.005^2),
  log_kappa_sq = log(30^2),
  w = rep(0, nrow(spde$c0)),
  u = -5,
  mu = -0.5,
  z.star = rnorm(sum(is.na(data_small$y)), 0, 1)
)

# Data and hyperparameter list
dat <- list(
  y = data_small$y,
  spde = spde,   # contains finite element matrices for f
  bw = 30        # bandwidth parameter enabling sparse Hessian
)


# Fitting the model
obj_small <- MakeADFun(jnll_small, par, random = c("w", "z.star"))
opt_small <- nlminb(obj_small$par, obj_small$fn, obj_small$gr)

# Check accurary of the Laplace approximation
set.seed(123)
l_chk_small <- laplace_check(obj_small)
l_chk_small # -> reliable



# Track computation times -------------------------------------------------

# speed of the banded forward algorithm for different bandwidth parameters
ks <- 2:50
Speeds <- matrix(0, length(ks), 100)
for(j in 1:ncol(Speeds)) {
  speeds <- rep(0, length(ks))
  for(i in 1:length(ks)){
    s <- Sys.time()
    forward(mod$delta, mod$Gamma, mod$allprobs, trackID = mod$trackID,
            bw = ks[i], ad = TRUE) # set ad = TRUE to run R version (which has bw option)
    speeds[i] <- Sys.time()-s
  }
  Speeds[,j] <- speeds
}
# speed of the regular forward algorithm
speed0 <- rep(NA, 100)
for(i in 1:length(speed0)) {
  s <- Sys.time()
  forward(mod$delta, mod$Gamma, mod$allprobs, trackID = mod$trackID, ad = TRUE)
  speed0[i] <- Sys.time()-s
}


# pdf("./figs/computation_time.pdf", width = 6, height = 4)
plot(ks, rowMeans(Speeds), ylim = c(0, 0.1), pch = 20, xlim = c(0, 50),
     xlab = "Bandwidth parameter (k)", ylab = "Time (sec)", bty = "n")
segments(x0 = 0, x1 = 50, y0 = mean(speed0), col = "darkblue")
segments(x0 = 0, x1 = 50, y0 = 2*mean(speed0), col = "darkblue", lty = 2)
rect(0, 0, 5, 0.1, col = "#00000010", border = NA)
legend("topright", bty = "n",
       legend = c("Banded forward algorithm",
                  "Regular forward algorithm",
                  "2x regular forward algorithm"),
       lty = c(NA, 1, 2), pch = c(16, NA, NA), col = c("black", "darkblue", "darkblue"))
# dev.off()

# practically twice that of the regular forward algorithm
