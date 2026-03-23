####### Case study - lions in Kalahari #######

# Installing required packages
# Newest CRAN versions of RTMB (1.9), LaMa (2.1.0), and RTMBdist (1.0.2) are required!
# install.packages(c("RTMB", "fmesher", "LaMa", "RTMBdist"))

library(LaMa)       # for HMM functions
library(RTMBdist)   # for zero-inflated gamma distribution
library(fmesher)    # for mesh and finite element matrices
library(viridis)    # for field color palette
library(scales)     # for semi-transparent color
library(leaflet)    # for satellite images


### Colors for plotting and random init
source("utils.R")



# Data and satellite EDA --------------------------------------------------

### Loading data
data <- readRDS("./data/lions.rds")
head(data)
nrow(data)


### Plot data on satellite image
leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  fitBounds(lng1 = min(data$x, na.rm = TRUE),
            lat1 = min(data$y, na.rm = TRUE),
            lng2 = max(data$x, na.rm = TRUE),
            lat2 = max(data$y, na.rm = TRUE)) %>%
  addCircleMarkers(lng = data$x,
                   lat = data$y,
                   radius = 0.5,
                   color = "#00000040")




# Homogeneous 2-state model -----------------------------------------------

# likelihood function
nll0 <- function(par) {
  getAll(par, dat) # unpacks lists

  Gamma <- tpm(eta) # transition matrix from unconstrained (softmax)
  delta <- stationary(Gamma) # stationary dist

  # parameter transformations
  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  rho <- plogis(logit_rho); REPORT(rho); REPORT(mu_a)

  # computing state-dependent densities
  allprobs <- matrix(1, length(step), N)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:N) {
    allprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j]) *
      dwrpcauchy(angle[ind], mu_a[j], rho[j])
  }
  # HMM likelihood via forward algorithm (separated by tracks)
  -forward(delta, Gamma, allprobs, trackID = ID)
}


# Initial parameter list
par <- list(
  eta = rep(-2, 2),
  log_mu = log(c(0.01, 1)),
  log_sigma = log(c(0.01, 1)),
  logit_zprob = qlogis(c(0.01, 0.01)),
  logit_rho = qlogis(c(0.2, 0.3)),
  mu_a = c(pi, 0)
)

# Data and hyperparameter list
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2
)


# Constructing AD likelihood
map <- list(mu_a = factor(rep(NA, 2))) # fixing angle mean at init
obj0 <- MakeADFun(nll0, par, map = map)

# Optimising for many initial values
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)

# Getting reported quantities from best model
mod0 <- report(obj0)

# state-dependent parameters
mu <- mod0$mu
sigma <- mod0$sigma
(rho <- mod0$rho)
mod0$zprob

# state process parameters
mod0$Gamma
delta <- mod0$delta


### Plotting state-dependent distributions
par(mfrow = c(1,2))
# step length
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
# turning angle
hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30))
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
# marginal = sum
curve(delta[1] * dwrpcauchy(x, c(pi,0)[1], rho[1]) +
        delta[2] * dwrpcauchy(x, c(pi,0)[2], rho[2]), add = TRUE, lwd = 2, lty = 2)


### Pseudo-residuals
pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod0
)
plot(pres_step)
# no residuals for wrpcauchy possible


### State decoding via Viterbi algorithm
mod0$states <- viterbi(mod = mod0)


### Plotting locations with decoded states
par(mfrow = c(1,1))
idx <- 1:10000
plot(data$x[idx], data$y[idx], asp = 1,type = "l")
points(data$x[idx], data$y[idx], cex = 0.8,
       col = scales::alpha(colorCB[mod0$states[idx]], 0.4), pch = 20)



# Homogeneous 3-state model -----------------------------------------------

# Initial parameter list
par <- list(
  eta = rep(-2, 6),
  log_mu = log(c(0.01, 0.5, 1.5)),
  log_sigma = log(c(0.01, 0.5, 1.5)),
  logit_zprob = qlogis(rep(1e-3, 3)),
  logit_rho = qlogis(c(0.2, 0.3, 0.4)),
  mu_a = c(pi, pi, 0)
)

# Data and hyperparamter list
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 3 # now 3 states
)

obj0.3 <- MakeADFun(nll0, par)
opt0.3 <- nlminb(obj0.3$par, obj0.3$fn, obj0.3$gr)
mod0.3 <- report(obj0.3)

mu <- mod0.3$mu
sigma <- mod0.3$sigma
rho <- mod0.3$rho
delta <- mod0.3$delta


### Plotting state-dependent distributions
par(mfrow = c(1,2))
# step lengths
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2),
     main = "", xlab = "Step length")
for(j in 1:3) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
# turning angles
hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30),
     main = "", xlab = "Turning angle")
for(j in 1:3) {
  curve(delta[j] * dwrpcauchy(x, c(pi,pi,0)[j], rho[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho[1]) +
        delta[2] * dwrpcauchy(x, pi, rho[2]) +
        delta[3] * dwrpcauchy(x, 0, rho[3]), add = TRUE, lwd = 2, lty = 2)


### Pseudo-residuals
pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod0.3
)
plot(pres_step)


## State decoding
mod0.3$states <- viterbi(mod = mod0.3)

idx <- 1:10000
par(mfrow = c(1,1))
plot(data$x[idx], data$y[idx], asp = 1, type = "l", bty = "n",
     xlab = "Longitude", ylab = "Latitude")
points(data$x[idx], data$y[idx], cex = 0.8,
       col = scales::alpha(colorCB[mod0.3$states[idx]], 0.4), pch = 20)

# We are setteling for 2 states:
### Much overlap between states 2 and 3 (splitting movement into 2)
### but TPM: basically no persistence -> discrimination between 2 and 3 seems random
### Conclusion: probably we cannot infer different movement behaviours at this coarse temporal resolution.

rm(obj0.3)
gc()



# 2-state model with periodic variation in tpm ----------------------------
# (We know lions are nocturnal so we want to capture this)

# likelihood
nll1 <- function(par) {
  getAll(par, dat)

  Gamma <- tpm_g(Z, beta)
  Delta <- stationary_p(Gamma) # periodically stationary dist
  REPORT(Delta)
  delta <- Delta[hour[startInd], ] # initial dist for each track

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  rho <- plogis(logit_rho); REPORT(rho)

  allprobs <- matrix(1, length(step), N)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:N) {
    allprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j]) *
      dwrpcauchy(angle[ind], c(pi, 0)[j], rho[j])
  }
  -forward_g(delta, Gamma[,,hour], allprobs, trackID = ID)
}

# Constructing sine, cosine design matrix
Z <- cosinor(1:24, period = c(24, 12, 6)) # 3 sin/cos minimises AIC and BIC

# Initial parameters and data
par <- list(
  beta = cbind(mod0$par$eta, matrix(0, 2, ncol(Z))),
  log_mu = log(mod0$mu),
  log_sigma = log(mod0$sigma),
  logit_zprob = qlogis(mod0$zprob),
  logit_rho = qlogis(mod0$rho)
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2,
  Z = Z,
  hour = data$hour,
  startInd = calc_trackInd(data$ID)
)

obj1 <- MakeADFun(nll1, par)
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
mod1 <- report(obj1)


### Sampling from distribution of MLE
samples1 <- MCreport(obj1)
Delta2s <- sapply(samples1$beta, function(b){
  stationary_p(tpm_g(Z, b))[,2]
})
# Computing quantials of p-stationary distribution
Delta_quantiles <- t(apply(Delta2s, 1, quantile, prob = c(0.025, 0.975)))

(mu <- mod1$mu)
(sigma <- mod1$sigma)
(rho <- mod1$rho)
Gamma <- tpm_g(Z, mod1$beta)
Delta <- stationary_p(Gamma) # at MLE


### Plotting periodically stationary state distribution
par(mfrow = c(1,1))
plot(Delta[,2], type = "b", ylim = c(0,1), col = colorCB[2], lwd = 2,
     bty = "n", xlab = "Hour", ylab = "Pr(active)")
points(Delta_quantiles[,1], col = colorCB[2], pch = 20)
points(Delta_quantiles[,2], col = colorCB[2], pch = 20)

mod1$states <- viterbi_g(mod = mod1)
delta <- prop.table(table(mod1$states))


## Plotting state-dependent distributions
par(mfrow = c(1,2))
# step length
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2),
     main = "", xlab = "Step length")
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
# turning angle
hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30),
     main = "", xlab = "Turning angle")
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
# marginal
curve(delta[1] * dwrpcauchy(x, pi, rho[1]) +
        delta[2] * dwrpcauchy(x, 0, rho[2]), add = TRUE, lwd = 2, lty = 2)


### Pseudo-residuals
pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod1
)
plot(pres_step)



# Model with spatial field in the tpm -------------------------------------

# Constructing the mesh for GMRF approximation to spatial field
loc <- cbind(data$x_int, data$y_int) # coords with interpolated NAs
# Boundaries
bnd_inner <- fm_nonconvex_hull(loc, convex=0.04) # inner boundary
bnd_outer <- fm_nonconvex_hull(loc, convex=0.15) # outer boundary
# Actual mesh building
mesh <- fm_mesh_2d(
  loc=loc,
  boundary=list(bnd_inner, bnd_outer),
  min.angle=24,
  max.edge=c(0.03, 1),
  cutoff=0.01,
  plot.delay=0.5
)

# pdf("./figs/lion_mesh.pdf", width = 7, height = 5)
oldpar <- par(mfrow = c(1,1), mar = rep(0.5, 4))
plot(mesh)
points(data$x, data$y, pch = 4, col = alpha("darkblue", 0.2), cex = 0.5, lwd = 0.5)
par(oldpar)
# dev.off()


### Calculate finite element matrices c0,g1, g2 need for precision matrix Q
spde <- fm_fem(mesh)
dim(spde$c0)


### Prediction matrix for observed locations
X_p <- fm_basis(mesh, loc)


# likelihood
jnll <- function(par) {
  getAll(par, dat, warn = FALSE)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  rho_angle <- plogis(logit_rho_angle); REPORT(rho_angle)

  # State process tpm building
  Eta <- matrix(eta, length(step), 2, byrow = TRUE)
  Eta[,1] <- Eta[,1] + as.numeric(X_p %*% w) # field only for Pr(active -> resting)
  Gamma <- tpm_g(Eta = Eta)
  Delta <- stationary(Gamma[,,startInd]) # use this as initial

  # State-dependent densities on log-scale now for accuracy
  lallprobs <- matrix(0, length(step), 2)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:2) {
    lallprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j], log = TRUE) +
      dwrpcauchy(angle[ind], c(pi, 0)[j], rho_angle[j], log = TRUE)
  }

  ## HMM likelihood
  nll <- -forward_g(Delta, Gamma, lallprobs, trackID = ID,
                    bw = bw, # bandwidth parameter leading to banded Hessian
                    logspace = TRUE) # computations on log-scale

  ## GMRF likelihood
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # dist at which corr has dropped to 0.1

  # Precision matrix of the field
  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
  nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

  nll
}

# Initial parameter list
par <- list(
  eta = mod0$par$eta,
  log_mu = log(mod0$mu),
  log_sigma = log(mod0$sigma),
  logit_zprob = qlogis(mod0$zprob),
  logit_rho_angle = qlogis(mod0$rho),
  log_tau_sq = log(0.1^2),
  log_kappa_sq = log(10^2),
  w = rep(0, nrow(spde$c0))
)

# Data and hyperparameter list
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2,
  startInd = calc_trackInd(data$ID),
  X_p = X_p,
  c0 = spde$c0, # C
  g1 = spde$g1, # G
  g2 = spde$g2, # G C^{-1} G
  bw = 15
)

t1 <- Sys.time()
# Constructing marginal likelihood via automatic Laplace
obj_sp1 <- MakeADFun(jnll, par,
                     random = "w") # telling RTMB to integrate out w

# Optimising marginal likelihood
opt_sp1 <- nlminb(obj_sp1$par, obj_sp1$fn, obj_sp1$gr) # inital inner Hessian evaluation takes time
Sys.time() - t1

# Reporting at optimum
mod_sp1 <- report(obj_sp1)

# Sampling parameters from approx posterior distribution (including random effects)
samples_sp1 <- MCreport(obj_sp1)

mod_sp1$states <- viterbi_g(mod = mod_sp1)
delta <- prop.table(table(mod_sp1$states))

mu <- mod_sp1$mu
sigma <- mod_sp1$sigma
rho_angle <- mod_sp1$rho_angle


### Plotting state-dependent distributions
par(mfrow = c(1,2))
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2),
     main = "", xlab = "Step length")
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 1000)
}
hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30),
     main = "", xlab = "Step length")
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho_angle[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho_angle[1]) +
        delta[2] * dwrpcauchy(x, 0, rho_angle[2]), add = TRUE, lwd = 2, lty = 2)


### Visualising estimated spatial field
# fine (1024 x 1024) spatial grid
x_seq <- seq(min(data$x_int), max(data$x_int), length.out = 1024)
y_seq <- seq(min(data$y_int), max(data$y_int), length.out = 1024)
grid <- as.matrix(expand.grid(x_seq, y_seq))

# Projection matrix on grid
A <- fm_basis(mesh, grid)
field <- as.numeric(A %*% mod_sp1$par$w) + mod_sp1$par$eta[1]

z1 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))
g1 <- plogis(z1) # transition probability

par(mfrow = c(1,1))
image(x_seq, y_seq, g1,
      xlab = "x", ylab = "y",
      col = viridis(35),
      main = expression(Pr(active~"→"~resting)), bty = "n", asp = 1)


### Save memory
rm(obj_sp1)
gc()



# Model with spatial field and periodic variation in the tpm --------------

# likelihood
jnll2 <- function(par) {
  getAll(par, dat, warn = FALSE)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  rho_angle <- plogis(logit_rho_angle); REPORT(rho_angle)

  Eta <- cbind(1, Z) %*% t(beta) # periodic variation part, calculations for 1:24 only
  Eta <- Eta[hour, ] # repeat rows according to time of day
  Eta[,1] <- Eta[,1] + as.numeric(X_p %*% w) # field only for Pr(active -> resting)
  Gamma <- tpm_g(Eta = Eta)
  Delta <- stationary(Gamma[,,startInd])

  lallprobs <- matrix(0, length(step), 2)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:2) {
    lallprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j], log = TRUE) +
      dwrpcauchy(angle[ind], c(pi, 0)[j], rho_angle[j], log = TRUE)
  }

  ## HMM likelihood
  nll <- -forward_g(Delta, Gamma, lallprobs, trackID = ID,
                    bw = bw, logspace = TRUE)

  ## GMRF likelihood
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # dist at which corr has dropped to 0.1

  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
  nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

  nll
}

# Initial parameters and data
par <- list(
  beta = mod1$beta,
  log_mu = log(mod1$mu),
  log_sigma = log(mod1$sigma),
  logit_zprob = qlogis(mod1$zprob),
  logit_rho_angle = qlogis(mod1$rho),
  log_tau_sq = log(0.1^2),
  log_kappa_sq = log(10^2),
  w = rep(0, nrow(spde$c0))
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2,
  startInd = calc_trackInd(data$ID),
  Z = Z,
  hour = data$hour,
  X_p = X_p,
  c0 = spde$c0,
  g1 = spde$g1,
  g2 = spde$g2,
  bw = 15
)

t1 <- Sys.time()
# Constructing marginal likelihood via automatic Laplace
obj_sp2 <- MakeADFun(jnll2, par, random = "w")

# Optimising marginal likelihood
opt_sp2 <- nlminb(obj_sp2$par, obj_sp2$fn, obj_sp2$gr)
Sys.time() - t1

mod_sp2 <- report(obj_sp2)
samples_sp2 <- MCreport(obj_sp2)

mod_sp2$states <- viterbi_g(mod = mod_sp2)
delta <- prop.table(table(mod_sp2$states)) # Monte Carlo approx to state distribution

mu <- mod_sp2$mu
sigma <- mod_sp2$sigma
rho_angle <- mod_sp2$rho_angle


### Plotting state-dependent distributions
# pdf("./figs/lions_statedep.pdf", width = 7, height = 3.5)
par(mfrow = c(1,2), mar = c(5,4,4,2)+0.1, xpd = FALSE)
hist(data$step, prob = TRUE, bor = "white", breaks = 100,
     main = "", xlab = "Step length (km)",
     xlim = c(0,4), ylim = c(0,1))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 1000)
}
curve(delta[1] * dgamma2(x, mu[1], sigma[1])+
        delta[2] * dgamma2(x, mu[2], sigma[2]), add = TRUE, lwd = 2, lty = 2, n = 1000)

legend("topright", col = c(colorCB[1:2], "black"), lty = c(1,1,2), lwd = 2,
       legend = c("State 1", "State 2", "Marginal"), bty = "n")

hist(data$angle, prob = TRUE, bor = "white",
     breaks = seq(-pi, pi, length = 30),
     main = "", xlab = "Turning angle (radians)",
     xaxt = "n")
axis(1, at = seq(-pi, pi, by = pi/2),
     labels = c(expression(-pi), expression(-pi/2), expression(0), expression(pi/2), expression(pi)))
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho_angle[j]), add = TRUE, col = colorCB[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho_angle[1]) +
        delta[2] * dwrpcauchy(x, 0, rho_angle[2]), add = TRUE, lwd = 2, lty = 2, n = 500)
# dev.off()


### Visualising estimated field in this model

field <- as.numeric(A %*% mod_sp2$par$w) + mod1$par$beta[1,1]
z2 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))
g2 <- plogis(z2)


# cairo_pdf("./figs/lions_spatial_field.pdf", width = 7.365, height = 4.95)
oldpar <- par(mfrow = c(1,1), mar = c(5,4,4,4), xpd = NA)
image(x_seq, y_seq, g2,
      xlab = "Longitude", ylab = "Latitude",
      # col = viridis(300),
      col = hcl.colors(35),
      main = expression(Pr(active~"→"~resting)), bty = "n", asp = 1,
      useRaster = TRUE)

# Move legend outside plot (to the right)
usr <- par("usr")
legend_x <- c(usr[2] + 0.025 * diff(usr[1:2]),
              usr[2] + 0.045 * diff(usr[1:2]))
legend_y <- c(usr[3]+0.03,
              usr[4]-0.03)
legend_colors <- as.raster(matrix(rev(hcl.colors(100)), ncol = 1))

rasterImage(legend_colors,
            legend_x[1], legend_y[1],
            legend_x[2], legend_y[2])

legend_values <- round(seq(min(g2), max(g2), length.out = 6), 1)
legend_positions <- seq(legend_y[1], legend_y[2], length.out = 6)

text(legend_x[2] + 0.02 * diff(usr[1:2]),
     legend_positions,
     labels = legend_values,
     adj = 0,
     cex = 0.9)
par(oldpar)
# dev.off()


# Comparing field in both models
par(mfrow = c(1,2))
image(x_seq, y_seq, z1,
      xlab = "x", ylab = "y",
      col = viridis(15),
      main = "spatial model", bty = "n", asp = 1)
image(x_seq, y_seq, z2,
      xlab = "x", ylab = "y",
      col = viridis(15),
      main = "spatial model with periodic var.", bty = "n", asp = 1)


### Plotting periodically stationary state distribution
betas <- samples_sp2$beta
B <- length(betas)
tseq <- seq(0, 24, length = 200)
Deltas <- array(dim = c(length(tseq), 2, B))
Delta <- matrix(NA, length(tseq), 2)

# Compute periodically stationary (and quantiles) on fine grid
for(t in seq_along(tseq)) { # takes some time
  ts <- tseq[t] + 0:23
  ts <- ts %% 24
  Z <- cosinor(ts, period = c(24,12,6))
  Delta[t,] <-  stationary_p(tpm_g(Z, mod_sp2$par$beta), t = 1)
  for(i in 1:B) {
    Gamma <- tpm_g(Z, betas[[i]])
    Deltas[t, ,i] <- stationary_p(Gamma, t = 1)
  }
}
Delta_q <- apply(Deltas, 1:2, quantile, probs = c(0.025, 0.975))


# pdf("./figs/lions_pstationary.pdf", width = 6, height = 4)

par(mfrow = c(1,1))
plot(NA, bty = "n", ylim = c(0,0.8), xlim = c(0,24), ylab = "Pr(active)", xlab = "Time of day", xaxt = "n")
for(t in 0:47){
  polygon(x = c(t/2, (t+1)/2, (t+1)/2, t/2), y = c(-0.03, -0.03, 0, 0), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
}
polygon(x = c(0, 6.5, 6.5, 0), y = c(0.01, 0.01, 0.8, 0.8), col = "gray90", border = FALSE)
polygon(x = c(18.5, 24, 24, 18.5), y = c(0.01, 0.01, 0.8, 0.8), col = "gray90", border = FALSE)

polygon(c(tseq, rev(tseq)), c(Delta_q[1, ,2], rev(Delta_q[2, ,2])),
        col = alpha("black", 0.25), border = FALSE)
lines(tseq, Delta[,2], lwd = 3, col = alpha("black", 0.7))

axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# dev.off()


### Model comparison in terms of information criteria
AIC(mod0, mod1, mod_sp1, mod_sp2)
BIC(mod0, mod1, mod_sp1, mod_sp2)
# both prefer spatial model with periodic variation


### Residual comparison in the two models
pres1 <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod_sp1
)
pres2 <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod_sp2
)
par(mfrow = c(2,3))
plot(pres1)
plot(pres2)
