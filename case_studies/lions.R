### installing dev version of packages
# devtools::install_github("janolefi/LaMa")
# devtools::install_github("kaskr/RTMB", subdir = "RTMB")
# devtools::install_github("janolefi/RTMBdist")

library(LaMa)
library(RTMBdist)
library(fmesher)
library(viridis)
library(leaflet)
library(scales)

color <- c("orange", "deepskyblue", "seagreen2")

# lood lion GPS data
data <- readRDS("./data/lions.rds")
head(data)
nrow(data)


## plot data on satellite image
# leaflet() %>%
#   addProviderTiles(providers$Esri.WorldImagery) %>%
#   fitBounds(lng1 = min(data$x, na.rm = TRUE),
#             lat1 = min(data$y, na.rm = TRUE),
#             lng2 = max(data$x, na.rm = TRUE),
#             lat2 = max(data$y, na.rm = TRUE)) %>%
#   addCircleMarkers(lng = data$x,
#                    lat = data$y,
#                    radius = 0.5,
#                    color = "#00000040")


# Homogeneous model -------------------------------------------------------

nll0 <- function(par) {
  getAll(par, dat)

  Gamma <- tpm(eta)
  delta <- stationary(Gamma)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  rho <- plogis(logit_rho); REPORT(rho); REPORT(mu_a)

  allprobs <- matrix(1, length(step), N)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:N) {
    allprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j]) *
      dwrpcauchy(angle[ind], mu_a[j], rho[j])
  }
  -forward(delta, Gamma, allprobs, trackID = ID)
}

## 2-state model

par <- list(
  eta = rep(-2, 2),
  log_mu = log(c(0.01, 1)),
  log_sigma = log(c(0.01, 1)),
  logit_zprob = qlogis(c(0.01, 0.01)),
  logit_rho = qlogis(c(0.2, 0.3)),
  mu_a = c(pi, 0)
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 2
)

obj0 <- MakeADFun(nll0, par, map = list(mu_a = factor(rep(NA, 2))))
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)

mod0 <- report(obj0)
mu <- mod0$mu
sigma <- mod0$sigma
(rho <- mod0$rho)

mod0$zprob
mod0$Gamma
delta <- mod0$delta

par(mfrow = c(1,2))
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}

hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30))
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, c(pi,0)[1], rho[1]) +
        delta[2] * dwrpcauchy(x, c(pi,0)[2], rho[2]), add = TRUE, lwd = 2, lty = 2)

pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod0
)
plot(pres_step)

# no residuals for wrpcauchy possible


mod0$states <- viterbi(mod = mod0)

par(mfrow = c(1,1))
idx <- 1:10000
# plot(data$step[idx], col = color[mod0$states[idx]], type = "h")
plot(data$x[idx], data$y[idx], asp = 1,type = "l")
points(data$x[idx], data$y[idx], cex = 0.8,
       col = scales::alpha(color[mod0$states[idx]], 0.4), pch = 20)


## 3-state model

par <- list(
  eta = rep(-2, 6),
  log_mu = log(c(0.01, 0.5, 1.5)),
  log_sigma = log(c(0.01, 0.5, 1.5)),
  logit_zprob = qlogis(rep(1e-3, 3)),
  logit_rho = qlogis(c(0.2, 0.3, 0.4)),
  mu_a = c(pi, pi, 0)
)
dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  N = 3
)

obj0.3 <- MakeADFun(nll0, par)
opt0.3 <- nlminb(obj0.3$par, obj0.3$fn, obj0.3$gr)

mod0.3 <- report(obj0.3)

mu <- mod0.3$mu
sigma <- mod0.3$sigma
rho <- mod0.3$rho
delta <- mod0.3$delta


par(mfrow = c(1,2))
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:3) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}

hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30))
for(j in 1:3) {
  curve(delta[j] * dwrpcauchy(x, c(pi,pi,0)[j], rho[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho[1]) +
        delta[2] * dwrpcauchy(x, pi, rho[2]) +
        delta[3] * dwrpcauchy(x, 0, rho[3]), add = TRUE, lwd = 2, lty = 2)

pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod0.3
)
plot(pres_step)

mod0.3$states <- viterbi(mod = mod0.3)

idx <- 1:10000
# plot(data$step[idx], col = color[mod0$states[idx]], type = "h")
plot(data$x[idx], data$y[idx], asp = 1,type = "l")
points(data$x[idx], data$y[idx], cex = 0.8,
       col = scales::alpha(color[mod0.3$states[idx]], 0.4), pch = 20)


# We are setteling for 2 states:
### Much overlap between states 2 and 3 (splitting movement into 2)
### but TPM: basically no persistence -> discrimination between 2 and 3 seems random
### Conclusion: probably we cannot infer different movement behaviours at this coarse temporal resolution.

rm(obj0.3)
gc()


# Model with periodic variation in the tpm --------------------------------

# We know lions are nocturnal so we want to capture this.

nll1 <- function(par) {
  getAll(par, dat)

  Gamma <- tpm_g(Z, beta)
  Gamma <- Gamma[,, hour]
  Delta <- matrix(0, length(startInd), 2)
  for(i in 1:length(startInd)){
    Delta[i, ] <- stationary_p(Gamma[,,startInd[i]+0:23], t = 1)
  }

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
  -forward_g(Delta, Gamma, allprobs, trackID = ID)
}

## 2 state model

Z <- cosinor(1:24, period = c(24, 12, 6)) # 3 sin/cos minimises AIC and BIC

par <- list(
  beta = cbind(rep(-2, 2), matrix(0, 2, ncol(Z))),
  log_mu = log(c(0.01, 1)),
  log_sigma = log(c(0.01, 2)),
  logit_zprob = qlogis(c(0.001, 0.001)),
  logit_rho = qlogis(c(0.2, 0.3))
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

(mu <- mod1$mu)
(sigma <- mod1$sigma)
(rho <- mod1$rho)
Gamma <- tpm_g(Z, mod1$beta)
Delta <- stationary_p(Gamma)

par(mfrow = c(1,1))
plot(Delta[,2], type = "b", ylim = c(0,1), col = color[2], lwd = 2)

mod1$states <- viterbi_g(mod = mod1)
delta <- prop.table(table(mod1$states))


par(mfrow = c(1,2))
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}

hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30))
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho[1]) +
        delta[2] * dwrpcauchy(x, 0, rho[2]), add = TRUE, lwd = 2, lty = 2)

pres_step <- pseudo_res(
  data$step,
  "gamma2",
  list(mean = mu, sd = sigma),
  mod = mod1
)
plot(pres_step)



# Model with spatial field in the tpm -------------------------------------

# constructing the mesh for GMRF approximation to spatial field
loc <- cbind(data$x_int, data$y_int) # coords with interpolated NAs

bnd_inner <- fm_nonconvex_hull(loc, convex=0.04) # inner boundary
bnd_outer <- fm_nonconvex_hull(loc, convex=0.15) # outer boundary

mesh <- fm_mesh_2d(
  loc=loc,
  boundary=list(bnd_inner, bnd_outer),
  min.angle=24,
  max.edge=c(0.03, 1),
  cutoff=0.01,
  plot.delay=0.5
)

par(mfrow = c(1,1))
plot(mesh)
# points(data$x, data$y)

# Calculate finite element matrices c0,g1, g2 need for precision matrix Q
spde <- fm_fem(mesh)
dim(spde$c0)

# prediction matrix for observed locations
X_p <- fm_basis(mesh, loc)

jnll <- function(par) {
  getAll(par, dat, warn = FALSE)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  rho_angle <- plogis(logit_rho_angle); REPORT(rho_angle)

  Eta <- matrix(eta, length(step), 2, byrow = TRUE); REPORT(eta)
  Eta[,1] <- Eta[,1] + as.numeric(X_p %*% w); REPORT(w) # field only for Pr(active -> resting)

  Gamma <- tpm_g(Eta = Eta)
  Delta <- stationary(Gamma[,,startInd])

  lallprobs <- matrix(0, length(step), 2)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:2) {
    lallprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j], log = TRUE) +
      dwrpcauchy(angle[ind], c(pi, 0)[j], rho_angle[j], log = TRUE)
  }

  ## HMM likelihood
  nll <- -forward_g(Delta, Gamma, lallprobs,
                     trackID = ID, bw = bw, logspace = TRUE, ad = T)

  ## GMRF likelihood
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # dist at which corr has dropped to 0.1

  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
  nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

  nll
}

par <- list(
  eta = qlogis(c(0.1, 0.1)),
  log_mu = log(mod0$mu),
  log_sigma = log(mod0$sigma),
  logit_zprob = qlogis(mod0$zprob),
  logit_rho_angle = qlogis(mod0$rho),
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
  c0 = spde$c0, # C
  g1 = spde$g1, # G
  g2 = spde$g2, # G C^{-1} G
  bw = 15
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
rho_angle <- mod_sp1$rho_angle

par(mfrow = c(1,2))
hist(data$step, prob = TRUE, bor = "white", breaks = 150, xlim = c(0,4), ylim = c(0,2))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 1000)
}

hist(data$angle, prob = TRUE, bor = "white", breaks = seq(-pi, pi, length = 30))
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho_angle[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho_angle[1]) +
        delta[2] * dwrpcauchy(x, 0, rho_angle[2]), add = TRUE, lwd = 2, lty = 2)


## visualise estimated field

x_seq <- seq(min(data$x_int), max(data$x_int), length.out = 1024)
y_seq <- seq(min(data$y_int), max(data$y_int), length.out = 1024)
grid <- as.matrix(expand.grid(x_seq, y_seq))

# Projection matrix on grid
A <- fm_basis(mesh, grid)
field <- as.numeric(A %*% mod_sp1$w) + mod_sp1$eta[1]

z1 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))
g1 <- plogis(z1) # transition probability

par(mfrow = c(1,1))
image(x_seq, y_seq, g1,
      xlab = "x", ylab = "y",
      col = viridis(35),
      main = expression(Pr(active~"→"~resting)), bty = "n", asp = 1)

# idx <- which(mod_sp1$states == 1)
# points(data$x[-idx], data$y[-idx], col = "#00000040")

rm(obj_sp1)
gc()


# Model with spatial field and periodic variation in the tpm --------------

jnll2 <- function(par) {
  getAll(par, dat, warn = FALSE)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  zprob <- plogis(logit_zprob); REPORT(zprob)
  rho_angle <- plogis(logit_rho_angle); REPORT(rho_angle)

  Eta <- cbind(1, Z) %*% t(beta); REPORT(beta) # periodic variation part, calculations for 1:24 only
  Eta <- Eta[hour, ] # rep according to time of day
  Eta[,1] <- Eta[,1] + as.numeric(X_p %*% w); REPORT(w) # field only for Pr(active -> resting)
  Gamma <- tpm_g(Eta = Eta)
  Delta <- stationary(Gamma[,,startInd])

  lallprobs <- matrix(0, length(step), 2)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:2) {
    lallprobs[ind,j] <- dzigamma2(step[ind], mu[j], sigma[j], zprob[j], log = TRUE) +
      dwrpcauchy(angle[ind], c(pi, 0)[j], rho_angle[j], log = TRUE)
  }

  ## HMM likelihood
  nll <- -forward_g(Delta, Gamma, lallprobs,
                    trackID = ID, bw = bw, logspace = TRUE)

  ## GMRF likelihood
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # dist at which corr has dropped to 0.1

  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
  nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

  nll
}

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
  ID = data$ID,
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
obj_sp2 <- MakeADFun(jnll2, par, random = "w")
opt_sp2 <- nlminb(obj_sp2$par, obj_sp2$fn, obj_sp2$gr)
Sys.time() - t1

mod_sp2 <- report(obj_sp2)

sdr <- sdreport(obj_sp2, getJointPrecision = TRUE)
mod_sp2$sdr <- sdr

mod_sp2$states <- viterbi_g(mod = mod_sp2)
delta <- prop.table(table(mod_sp2$states))

mu <- mod_sp2$mu
sigma <- mod_sp2$sigma
rho_angle <- mod_sp2$rho_angle


# pdf("./figs/lions_statedep.pdf", width = 7, height = 3.5)

par(mfrow = c(1,2), mar = c(5,4,4,2)+0.1, xpd = FALSE)

hist(data$step, prob = TRUE, bor = "white", breaks = 100,
     main = "", xlab = "Step length (km)",
     xlim = c(0,4), ylim = c(0,1))
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2, n = 1000)
}
curve(delta[1] * dgamma2(x, mu[1], sigma[1])+
        delta[2] * dgamma2(x, mu[2], sigma[2]), add = TRUE, lwd = 2, lty = 2, n = 1000)

legend("topright", col = c("orange", "deepskyblue", "black"), lty = c(1,1,2), lwd = 2,
       legend = c("State 1", "State 2", "Marginal"), bty = "n")

hist(data$angle, prob = TRUE, bor = "white",
     breaks = seq(-pi, pi, length = 30),
     main = "", xlab = "Turning angle (radians)",
     xaxt = "n")
axis(1, at = seq(-pi, pi, by = pi/2),
     labels = c(expression(-pi), expression(-pi/2), expression(0), expression(pi/2), expression(pi)))
for(j in 1:2) {
  curve(delta[j] * dwrpcauchy(x, c(pi,0)[j], rho_angle[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho_angle[1]) +
        delta[2] * dwrpcauchy(x, 0, rho_angle[2]), add = TRUE, lwd = 2, lty = 2, n = 500)

# dev.off()


## visualise estimated field

field <- as.numeric(A %*% mod_sp2$w) + mod1$beta[1,1]

z2 <- matrix(field, nrow = length(x_seq), ncol = length(y_seq))
g2 <- plogis(z2)


# cairo_pdf("./figs/lions_spatial_field.pdf", width = 7.365, height = 4.95)
par(mfrow = c(1,1), mar = c(5,4,4,4), xpd = NA)
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
# dev.off()


# both fields
par(mfrow = c(1,2))
image(x_seq, y_seq, z1,
      xlab = "x", ylab = "y",
      col = viridis(15),
      main = "spatial model", bty = "n", asp = 1)
image(x_seq, y_seq, z2,
      xlab = "x", ylab = "y",
      col = viridis(15),
      main = "spatial model with periodic var.", bty = "n", asp = 1)


### plot periodically stationary state distribution

Sigma <- mod_sp2$sdr$cov.fixed
B <- 1000
thetas <- mvtnorm::rmvnorm(B, mod_sp2$sdr$par.fixed, Sigma)
betas <- array(dim = c(2, 7, B))
for(i in 1:B) betas[,,i] <- matrix(thetas[i,1:14], 2, 7)

tseq <- seq(0, 24, length = 200)
Deltas <- array(dim = c(length(tseq), 2, B))
Delta <- matrix(NA, length(tseq), 2)

for(t in seq_along(tseq)) {
  print(t)
  ts <- tseq[t] + 0:23
  ts <- ts %% 24
  Z <- cosinor(ts, period = c(24,12,6))
  Delta[t,] <-  stationary_p(tpm_g(Z, mod_sp2$beta), t = 1)
  for(i in 1:B) {
    Gamma <- tpm_g(Z, betas[,,i])
    Deltas[t, ,i] <- stationary_p(Gamma, t = 1)
  }
}
Delta_q <- apply(Deltas, 1:2, quantile, probs = c(0.025, 0.975))


# pdf("./figs/lions_pstationary.pdf", width = 6, height = 4)

par(mfrow = c(1,1))
plot(tseq, Delta[,2], type = "l", lwd = 3,
     ylim = c(0,0.8), col = color[2], bty = "n",
     ylab = "Pr(active)", xlab = "Time of day", xaxt = "n")
polygon(c(tseq, rev(tseq)), c(Delta_q[1, ,2], rev(Delta_q[2, ,2])),
        col = scales::alpha(color[2], 0.25), border = NA)
axis(1, at = seq(0, 24, by = 6), labels = seq(0, 24, by = 6))

# dev.off()


pdf("./figs/lions_pstationary.pdf", width = 6, height = 4)

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
# lines(tseq[-ind], Delta[-ind,2], lwd = 2, col = alpha("black", 0.5))

axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

dev.off()



## model comparison
AIC(mod0, mod1, mod_sp1, mod_sp2)
BIC(mod0, mod1, mod_sp1, mod_sp2)


## residual comparison
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


