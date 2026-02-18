library(LaMa)
library(RTMBdist)
library(fmesher)
library(viridis)

color <- c("orange", "deepskyblue", "seagreen2")

data <- readRDS("./data/lions.rds")

head(data)
nrow(data)

hist(data$step, prob = TRUE)
hist(data$angle, prob = TRUE)

nll0 <- function(par) {
  getAll(par, dat)

  Gamma <- tpm(eta)
  delta <- stationary(Gamma)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  # kappa <- exp(log_kappa); REPORT(kappa)

  allprobs <- matrix(1, length(step), N)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:N) {
    allprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j]) # *
      # dvm(angle[ind], c(pi, 0, 0)[j], kappa[j])
  }
  -forward(delta, Gamma, allprobs, trackID = ID)
}

## 2 state model

par <- list(
  eta = rep(-2, 2),
  log_mu = log(c(0.2, 2)),
  log_sigma = log(c(0.2, 1.5)),
  logit_zprob = qlogis(c(0.01, 0.01))#,
  # log_kappa = log(c(0.5, 1))
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2
)

obj0 <- MakeADFun(nll0, par)
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)

mod0 <- report(obj0)
mu <- mod0$mu
sigma <- mod0$sigma
# kappa <- mod0$kappa
mod0$zprob
mod0$Gamma
delta <- mod0$delta

hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}

# hist(data$angle, prob = TRUE, bor = "white")
# for(j in 1:2) {
#   curve(delta[j] * dvm(x, c(pi,0)[j], kappa[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
# }
# curve(delta[1]*dvm(x, pi, kappa[1])+delta[2]*dvm(x, 0, kappa[2]), add = TRUE, lwd = 2, lty = 2)

pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod0
)
# pres_angle <- pseudo_res(
#   data$angle,
#   "vm",
#   list(mu = c(pi, 0), kappa = kappa),
#   mod = mod0
# )
plot(pres_step)
# plot(pres_angle)

mod0$states <- viterbi(mod = mod0)

idx <- 1:1000
plot(data$step[idx], col = color[mod0$states[idx]], type = "h")
plot(data$x[idx], data$y[idx], col = color[mod0$states[idx]], asp = 1, pch = 16)



# Model with time of day --------------------------------------------------

Z <- cosinor(1:24, period = c(24, 12))

nll1 <- function(par) {
  getAll(par, dat)

  Gamma <- tpm_g(Z, beta)
  Gamma <- Gamma[,, tod]
  Delta <- matrix(0, length(startInd), 2)
  for(i in 1:length(startInd)){
    Delta[i, ] <- stationary_p(Gamma[,,startInd[i]+0:23], t = 1)
  }

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  # kappa <- exp(log_kappa); REPORT(kappa)

  allprobs <- matrix(1, length(step), N)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:N) {
    allprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j]) # *
      # dvm(angle[ind], c(pi, 0, 0)[j], kappa[j])
  }
  -forward_g(Delta, Gamma, allprobs, trackID = ID)#, trackID = ID)
}

## 2 state model

par <- list(
  beta = cbind(rep(-2, 2), matrix(0, 2, ncol(Z))),
  log_mu = log(c(0.01, 1)),
  log_sigma = log(c(0.01, 2)),
  logit_zprob = qlogis(c(0.001, 0.001))# ,
  # log_kappa = log(c(0.5, 1))
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2,
  Z = Z,
  tod = data$tod,
  startInd = calc_trackInd(data$ID)
)

obj1 <- MakeADFun(nll1, par)
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

mod1 <- report(obj1)
mu <- mod1$mu
sigma <- mod1$sigma
# kappa <- mod1$kappa
Gamma <- tpm_g(Z, mod1$beta)
Delta <- stationary_p(Gamma)

plot(Delta[,2], type = "b", ylim = c(0,1), col = color[2], lwd = 2)

mod1$states <- viterbi_g(mod = mod1)
delta <- prop.table(table(mod1$states))

hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}

# hist(data$angle, prob = TRUE, bor = "white")
# for(j in 1:2) {
#   curve(delta[j] * dvm(x, c(pi,0)[j], kappa[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
# }
# curve(delta[1]*dvm(x, pi, kappa[1])+delta[2]*dvm(x, 0, kappa[2]), add = TRUE, lwd = 2, lty = 2)

pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod1
)
# pres_angle <- pseudo_res(
#   data$angle,
#   "vm",
#   list(mu = c(pi, 0), kappa = kappa),
#   mod = mod1
# )
plot(pres_step)
# plot(pres_angle)


# Model with spatial field ------------------------------------------------


## Spatial part part of model
loc <- cbind(data$x_int, data$y_int)  # Spatial coordinates

bnd_inner <- fm_nonconvex_hull(loc, convex=0.04)
bnd_outer <- fm_nonconvex_hull(loc, convex=0.15)

mesh <- fm_mesh_2d(
  loc=loc,
  boundary=list(bnd_inner, bnd_outer),
  min.angle=24,
  max.edge=c(0.03, 1),
  cutoff=0.01,
  plot.delay=0.5
)

plot(mesh)
points(data$x, data$y)

spde <- fm_fem(mesh) # Calculate the sparse matrices c0,g1, g2 need for precission matrix
dim(spde$c0)

X_p <- fm_basis(mesh, loc) # prediction matrix for observed locations

jnll <- function(par) {
  getAll(par, dat, warn = FALSE)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  # kappa_angle <- exp(log_kappa_angle); REPORT(kappa_angle)

  # Eta <- as.matrix(X_p %*% t(w)) + matrix(eta, length(step), 2, byrow = TRUE)
  Eta <- matrix(eta, length(step), 2, byrow = TRUE)
  Eta[,1] <- Eta[,1] + as.numeric(X_p %*% w) # field only for Pr(active -> resting)

  Gamma <- tpm_g(Eta = Eta)
  Delta <- stationary(Gamma[,,startInd])

  lallprobs <- matrix(0, length(step), 2)
  ind <- which(!is.na(step))# & !is.na(angle))
  for(j in 1:2) {
    lallprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j], log = TRUE)# +
     # dvm(step[ind], mu_angle[j], kappa_angle[j], log = TRUE)
  }

  ## HMM likelihood
  nll <- -forward_g(Delta, Gamma, lallprobs, trackID = ID, bw = bw, logspace = TRUE)

  ## GMRF likelihood
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # dist at which corr has dropped to 0.1

  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
  nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

  nll
}


# initial field parameters
sigma0 <- 0.1 # marginal sd, initialise small
range <- diff(range(data$x, na.rm = TRUE)) * sqrt(2) / 5 # largest distance devided by 5 -> smooth field
kappa0 <- sqrt(8) / range
tau0 <- 1 / (sigma0 * sqrt(4*pi) * kappa0)

par <- list(
  eta = qlogis(c(0.1, 0.1)),
  log_mu = log(c(0.01, 1)),
  log_sigma = log(c(0.01, 1)),
  logit_zprob = qlogis(rep(1e-4, 2)),
  # log_kappa_angle = log(c(0.4, 0.5)),
  # mu_angle = c(pi, 0),
  log_tau_sq = log(0.1^2),
  log_kappa_sq = log(10^2),
  w = rep(0, nrow(spde$c0))
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2,
  ID = data$ID,
  startInd = calc_trackInd(data$ID),
  X_p = X_p,
  c0 = spde$c0,
  g1 = spde$g1,
  g2 = spde$g2,
  bw = 10
)

t1 <- Sys.time()
obj_sp1 <- MakeADFun(jnll, par, random = "w")
opt_sp1 <- nlminb(obj_sp1$par, obj_sp1$fn, obj_sp1$gr)
Sys.time() - t1

mod_sp1 <- report(obj_sp1)

mod_sp1$states <- viterbi_g(mod = mod_sp1)
delta <- prop.table(table(mod_sp1$states))

mu <- mod_sp1$mu
sigma <- mod_sp1$sigma
# kappa <- mod$kappa_angle

hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 1000)
}

# hist(data$angle, prob = TRUE, bor = "white")
# for(j in 1:2) {
#   curve(delta[j] * dvm(x, c(pi,0)[j], kappa[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
# }
# curve(delta[1]*dvm(x, pi, kappa[1])+delta[2]*dvm(x, 0, kappa[2]), add = TRUE, lwd = 2, lty = 2)


## estimated fields

w <- obj_sp1$env$last.par.best[names(obj_sp1$env$last.par.best) == "w"]
# w <- matrix(w, 2, nrow(spde$c0))

x_seq <- seq(min(data$x_int), max(data$x_int), length.out = 1024)
y_seq <- seq(min(data$y_int), max(data$y_int), length.out = 1024)
grid <- as.matrix(expand.grid(x_seq, y_seq))

A <- fm_basis(mesh, grid)       # projection matrix (grid × vertices)
field <- as.numeric(A %*% w)

z21 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))
# z12 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))

par(mfrow = c(1,1))
image(x_seq, y_seq, z21,
      xlab = "x", ylab = "y",
      col = viridis(100),
      main = "logit(Pr(active -> resting))", bty = "n", asp = 1)

# idx <- which(mod_sp1$states == 1)
# points(data$x[-idx], data$y[-idx], col = "#00000040")




# Spatial field and periodic variation ------------------------------------

jnll2 <- function(par) {
  getAll(par, dat, warn = FALSE)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  # kappa_angle <- exp(log_kappa_angle); REPORT(kappa_angle)

  Eta <- cbind(1, Z) %*% t(beta); REPORT(beta) # periodic variation part, calculations for 1:24 only
  Eta <- Eta[tod, ] # rep according to time of day
  Eta[,1] <- Eta[,1] + as.numeric(X_p %*% w) # field only for Pr(active -> resting)
  Gamma <- tpm_g(Eta = Eta)
  Delta <- stationary(Gamma[,,startInd])

  lallprobs <- matrix(0, length(step), 2)
  ind <- which(!is.na(step))# & !is.na(angle))
  for(j in 1:2) {
    lallprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j], log = TRUE) # +
     # dvm(step[ind], c(pi, 0)[j], kappa_angle[j], log = TRUE)
  }

  ## HMM likelihood
  nll <- -forward_g(Delta, Gamma, lallprobs, trackID = ID, bw = bw, logspace = TRUE)

  ## GMRF likelihood
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # dist at which corr has dropped to 0.1

  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
  nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

  nll
}

# initial field parameters
# sigma0 <- 0.1 # marginal sd, initialise small
# range <- diff(range(data$x, na.rm = TRUE)) * sqrt(2) / 5 # largest distance devided by 5 -> smooth field
# kappa0 <- sqrt(8) / range
# tau0 <- 1 / (sigma0 * sqrt(4*pi) * kappa0)

par <- list(
  beta = cbind(rep(-2, 2), matrix(0, 2, ncol(Z))),
  log_mu = log(c(0.01, 1)),
  log_sigma = log(c(0.01, 1)),
  logit_zprob = qlogis(rep(1e-4, 2)),
  # log_kappa_angle = log(c(0.4, 0.5)),
  log_tau_sq = log(0.1^2),
  log_kappa_sq = log(10^2),
  w = rep(0, nrow(spde$c0))
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2,
  ID = data$ID,
  startInd = calc_trackInd(data$ID),
  Z = Z,
  tod = data$tod,
  X_p = X_p,
  c0 = spde$c0,
  g1 = spde$g1,
  g2 = spde$g2,
  bw = 10
)

t1 <- Sys.time()
obj_sp2 <- MakeADFun(jnll2, par, random = "w")
opt_sp2 <- nlminb(obj_sp2$par, obj_sp2$fn, obj_sp2$gr)
Sys.time() - t1

mod_sp2 <- report(obj_sp2)

mod_sp2$states <- viterbi_g(mod = mod_sp2)
delta <- prop.table(table(mod_sp2$states))

mu <- mod_sp2$mu
sigma <- mod_sp2$sigma
# kappa <- mod$kappa_angle

hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 1000)
}

# hist(data$angle, prob = TRUE, bor = "white")
# for(j in 1:2) {
#   curve(delta[j] * dvm(x, c(pi,0)[j], kappa[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
# }
# curve(delta[1]*dvm(x, pi, kappa[1])+delta[2]*dvm(x, 0, kappa[2]), add = TRUE, lwd = 2, lty = 2)


## estimated fields

w <- obj_sp2$env$last.par.best[names(obj_sp2$env$last.par.best) == "w"]
# w <- matrix(w, 2, nrow(spde$c0))

x_seq <- seq(min(data$x_int), max(data$x_int), length.out = 1024)
y_seq <- seq(min(data$y_int), max(data$y_int), length.out = 1024)
grid <- as.matrix(expand.grid(x_seq, y_seq))

A <- fm_basis(mesh, grid)       # projection matrix (grid × vertices)
field <- as.numeric(A %*% w)

z21 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))
# z12 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))

par(mfrow = c(1,1))
image(x_seq, y_seq, z21,
      xlab = "x", ylab = "y",
      col = viridis(100),
      main = "logit(Pr(active -> resting))", bty = "n", asp = 1)

# idx <- which(mod_sp2$states == 1)
# points(data$x[-idx], data$y[-idx], col = "#00000040")


Gamma <- tpm_g(Z, mod_sp$beta)
Delta <- stationary_p(Gamma)

plot(Delta[,2], type = "b", ylim = c(0,1), col = color[2], lwd = 2)



## model comparison
AIC(mod0, mod1, mod_sp1, mod_sp2)
BIC(mod0, mod1, mod_sp1, mod_sp2)


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
