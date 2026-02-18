library(LaMa)
library(fmesher)
library(viridis)
# source("./functions/forward_banded.R")

set.seed(38475)

################
## Simulation ##
################
# True spatial function for mean step length
true_field <- function(x,y){
  (sin(2*pi*x/5) + cos(2*pi*y/5)) / 10
}
mu <- function(x, y, state) {
  beta0 <- log(c(0.05, 0.3))
  exp(beta0[state] + true_field(x,y))
}

# Simulate state process
tpm <- matrix(c(0.9,0.1,0.1,0.9), nrow = 2)
# n <- 10000
n <- 30000
state <- rep(1, n)
for(i in 2:n) {
  state[i] <- sample(1:2, size = 1, prob = tpm[state[i-1],])
}

# Simulate observation process
xy <- matrix(0, nrow = n, ncol = 2)
step <- step_mean <- rep(NA, n)
for(i in 2:n) {
  step_mean[i] <- mu(x = xy[i-1,1], y = xy[i-1,2], state = state[i])
  step_sd <- c(0.02, 0.1)[state[i]]
  step[i-1] <- rgamma2(n = 1, step_mean[i], step_sd)
  angle <- runif(n = 1, min = -pi, max = pi)
  xy[i,] <- xy[i-1,] + step[i-1] * c(cos(angle), sin(angle))
}
data <- data.frame(x = xy[,1], y = xy[,2], step = step, state = state)

plot(data$x, data$y, col = state)

plot(step_mean, col = state)


# Make a prediction grid
x_seq <- seq(min(data$x), max(data$x), length.out = 512)
y_seq <- seq(min(data$y), max(data$y), length.out = 512)
grid <- as.matrix(expand.grid(x_seq, y_seq))

z1 <- outer(x_seq, y_seq, true_field)

par(mfrow = c(1,1))
image(x_seq, y_seq, z1,
      xlab = "x", ylab = "y",
      col = viridis(100),
      main = "True spatial field", asp = 1, bty = "n")
points(data$x, data$y, pch = 20, cex = 0.3, col = "black")


## Spatial part part of model
loc <- cbind(data$x, data$y)  # Spatial coordinates
bnd1 <- fm_nonconvex_hull(loc, convex = -0.05)
bnd2 <- fm_nonconvex_hull(loc, convex = -0.3)

mesh <- fm_mesh_2d(
  loc=loc,
  boundary=list(bnd1, bnd2),
  # min.angle=21,
  max.edge=c(0.5, 3),
  cutoff=0.2,
  plot.delay=0.5
)

plot(mesh)

spde <- fm_fem(mesh) # Calculate the sparse matrices c0,g1, g2 need for precission matrix
dim(spde$c0)


jnll <- function(par) {
  getAll(par, dat)

  sigma <- exp(log_sigma); REPORT(sigma)

  Gamma <- tpm(eta)
  delta <- stationary(Gamma)

  ## compute means for all field positions
  logMu <- matrix(beta0, nrow(c0), 2, byrow = TRUE)
  logMu[,1] <- logMu[,1] + field
  logMu[,2] <- logMu[,2] + field
  Mu <- exp(logMu); REPORT(Mu)
  Mu_idx <- Mu[meshidxloc,]; REPORT(Mu_idx)

  lallprobs <- matrix(0, length(step), 2)
  idx <- which(!is.na(step))
  for(j in 1:2) {
    lallprobs[idx,j] <- dgamma2(step[idx], Mu_idx[idx,j], sigma[j], log = TRUE)
  }

  nll <- -forward(delta, Gamma, lallprobs, bw = 15, logspace = TRUE)
  # nll <- -forward(delta, Gamma, lallprobs, logspace = TRUE)

  ## Gaussian field part
  tau <- exp(log_tau); REPORT(tau)
  kappa <- exp(log_kappa); REPORT(kappa)

  Q <- tau^2*(kappa^4 * c0 + 2 * kappa^2 * g1 + g2)        ## GMRF prior
  nll <- nll - dgmrf(field, 0, Q, log=TRUE)        ## Negative log likelihood

  nll
}

range <- max(dist(loc)) / 5
kappa = sqrt(8)/range
sigma0 = 0.1
tau = 1/(sigma0 * sqrt(4*pi) * kappa)

par <- list(
  beta0 = log(c(0.05, 0.3)),
  log_sigma = log(c(0.02, 0.1)),
  eta = qlogis(rep(0.1, 2)),
  log_tau = log(tau),
  log_kappa = log(kappa),
  field = rep(0, nrow(spde$c0))
)

dat <- list(
  step = data$step,
  meshidxloc = mesh$idx$loc,
  c0 = spde$c0,
  g1 = spde$g1,
  g2 = spde$g2
)

# TapeConfig(matmul = 'plain') # faster

obj <- MakeADFun(jnll, par, random = "field")
# H <- obj$env$spHess(random=TRUE)
# par(mfrow=c(1,1))
# SparseM::image(H)

system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr)
)

mod <- obj$report()
mod$sigma
mod$tau
mod$kappa
mod$Gamma

field_est <- obj$env$last.par.best[names(obj$env$last.par.best) == "field"]


A_pred <- fm_basis(mesh, grid)       # projection matrix (grid × vertices)
field_grid <- as.vector(A_pred %*% field_est)

z <- matrix(field_grid, nrow = length(x_seq), ncol = length(y_seq))

par(mfrow = c(1,2))
image(x_seq, y_seq, z1,
      xlab = "x", ylab = "y",
      col = viridis(100),
      main = "True spatial field", bty = "n", asp = 1)
# lines(data$x, data$y, lwd = 1, col = "#00000050")
# points(data$x, data$y, pch = 20, cex = 0.3, col = "#00000010")

image(x_seq, y_seq, z,
      xlab = "x", ylab = "y",
      col = viridis(100),
      main = "Estimated spatial field", bty = "n", asp = 1)
# lines(data$x, data$y, lwd = 1, col = "#00000050")
# points(data$x, data$y, pch = 20, cex = 0.3, col = "#00000010")


par(mfrow =c(1,1))
image(x_seq, y_seq, z-z1,
      xlab = "x", ylab = "y",
      main = "Estimated spatial field", bty = "n", asp = 1)
lines(data$x, data$y, lwd = 1, col = "#00000050")


ind <- z != 0
cor(z1[ind], z[ind])
