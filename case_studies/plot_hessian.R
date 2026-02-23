library(LaMa)       # for HMM functions
library(RTMBdist)   # for ExGaussian distribution
library(fmesher)    # for mesh and FEM matrices
library(Matrix)     # for sparse matrices
library(scales)


### loading data
data <- read.csv("./data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[,c("TIME","PDCSAP_FLUX")]
colnames(data) <- c("time", "y")

data <- data[2:201,]

# linearly interpolating missing observations
data$y[is.na(data$y)] = approx(data$time, data$y, data$time[is.na(data$y)])$y

# centering data
data$y <- scale(data$y, scale = FALSE)


### HMM analysis
# HMM likelihood function
jnll <- function(par) {
  getAll(par, dat)

  ## state process ##
  # restricted tpm
  Gamma <- diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  # estimated initial distribution
  delta <- c(1, exp(logit_delta))
  delta <- delta / sum(delta)

  ## state-dependent process ##
  # parameter transformations
  sigma <- exp(log_sigma); REPORT(sigma)
  r <- plogis(logit_r); REPORT(r)
  lambda <- exp(log_lambda); REPORT(lambda)
  # state-dependent densities
  f <- f0 + mu; REPORT(f0); REPORT(f) # smooth function
  z <- y - f; REPORT(z); REPORT(mu)

  n <- length(z); idx = 2:n
  lallprobs <- matrix(0, n, 3)
  # regular measurement error
  lallprobs[idx,1] <- dnorm(z[idx], 0, sigma, log = TRUE)
  # firing
  lallprobs[idx,2] <- dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  # decaying
  lallprobs[idx,3] <- dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)

  # banded forward algorithm - HMM likelihood
  nll <- - forward(delta, Gamma, lallprobs, bw = k, logspace = TRUE)
  # lls <- forward_summands(delta, Gamma, lallprobs, logspace = TRUE)

  ### GP part ###
  # parameter transformations
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1
  # omega <- plogis(logit_omega); REPORT(omega)
  cos_pi_omega <- 2 * plogis(u) - 1
  omega <- acos(cos_pi_omega) / pi; REPORT(omega)

  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * cos_pi_omega * kappa_sq * g1 + g2); REPORT(Q)

  nll <- nll - dgmrf(f - mu, 0, Q, log = TRUE) # GP likelihood

  nll
}

### creating mesh and finite element matrices
mesh <- fm_mesh_1d(data$time)
spde <- fm_fem(mesh)

# initial parameter list
par <- list(
  eta = rep(-2, 4),
  logit_delta = rep(0, 2),
  log_sigma = log(7),
  logit_r = qlogis(0.8),
  log_lambda = log(0.025),
  log_tau_sq = log(0.005^2),
  log_kappa_sq = log(30^2),
  f0 = numeric(nrow(data)),
  u = -5,
  mu = -0.5
)

# data list
dat <- list(
  y = data$y,
  c0 = spde$c0, g1 = spde$g1, g2 = spde$g2,
  k = NULL
)

obj_dense <- MakeADFun(jnll, par, random = "f0")
H_dense <- obj_dense$env$spHess(random = TRUE)

dat$k <- 5
obj_banded <- MakeADFun(jnll, par, random = "f0")
H_banded <- obj_banded$env$spHess(random = TRUE)

idx <- 1:60

# left panel
pdf("./figs/denseH.pdf", width = 5, height = 5)
SparseM::image(H_dense[idx, idx])
dev.off()

# right panel
pdf("./figs/bandedH.pdf", width = 5, height = 5)
SparseM::image(H_banded[idx, idx])
dev.off()
