####### Simulation experiment - spatial field in TPM #######

# Installing required packages
# Versions used:
# - RTMB: 1.9
# - LaMa: 2.1.3
# - fmesher: 0.8.0
# install.packages(c("RTMB", "fmesher", "LaMa", "viridis", "sf"))


library(parallel) # parallelising fits
library(LaMa)     # HMM functions


### Colors
source("utils.R")


# True spatial function for mean step length
true_field <- function(x,y){
  2 * (sin(2*pi*x / 40) + cos(2*pi*y / 40))
}

# True parameters used for simulation
true_par <- list(
  delta = c(0.5, 0.5),
  mu = c(0.2, 5),
  sigma = c(0.5, 3),
  beta0 = c(-2, -2)
)

# Plot true state-dependent distributions
curve(dgamma2(x, true_par$mu[1], true_par$sigma[1]), xlim = c(0, 10), ylim = c(0,1),
      n = 500, bty = "n", ylab = "Density", xlab = "Step length", lwd = 3, col = colorCB[1])
curve(dgamma2(x, true_par$mu[2], true_par$sigma[2]), add = TRUE, n = 500,
      lwd = 3, col = colorCB[2])



### Function that simulates one data set
sim_data <- function(n, # length of time series
                     p, # parameter list
                     kappa_pull = 0.3 # concentration for biased random walk
                     ) {
  s <- rep(NA, n)
  s[1] <- sample(1:2, size = 1, prob = p$delta) # first state
  loc <- matrix(0, n, 2) # initialise location matrix
  step <- rep(NA, n) # initialise step length vector
  phi <- 0 # initialise absolute heading
  field_val <- rep(NA, n)
  for(t in 1:(n-1)) {
    # simulate step
    step[t] <- rgamma2(1, p$mu[s[t]], p$sigma[s[t]])
    # angle: add a pull towards c(0,0)
    pull_angle <- atan2(-loc[t,2], -loc[t,1]) - phi
    # draw absolute angle (independent of state)
    angle <- rvm(1, pull_angle, kappa_pull)
    # update heading
    phi <- phi + angle
    # next location
    loc[t+1,1] <- loc[t,1] + step[t] * cos(phi)
    loc[t+1,2] <- loc[t,2] + step[t] * sin(phi)
    # tpm evalutated at next location
    eta <- p$beta0
    field_val[t+1] <- true_field(loc[t+1,1], loc[t+1,2])
    eta[1] <- eta[1] + field_val[t+1]
    Gamma <- tpm(eta)
    # print(Gamma[2,1])
    # sample next state based on that tpm
    s[t+1] <- sample(1:2, 1, prob = Gamma[s[t], ])
  }
  data.frame(step = step, x = loc[,1], y = loc[,2], state = s, field_val = field_val)
}



### Function that simulates data set and fits model for a given specification
one_fit <- function(dummy,
                    nObs = 5000,
                    par = true_par,
                    bw = 10) {
  set.seed(dummy) # for reproducibility of each iteration

  library(RTMB)
  library(LaMa) # also loads RTMB

  ### Changing AD setting inside RTMB for faster calculations
  TapeConfig(matmul = "plain")

  # Simulate data
  data <- sim_data(nObs, par)

  # Create a 512x512 prediction grid
  x_grid <- seq(min(data$x), max(data$x), length.out = 512)
  y_grid <- seq(min(data$y), max(data$y), length.out = 512)
  grid_data <- expand.grid(x = x_grid, y = y_grid)
  grid <- as.matrix(grid_data)
  z21 <- outer(x_grid, y_grid, true_field)

  # likelihood function
  jnll <- function(par) {
    getAll(par, dat, warn = FALSE)

    mu <- exp(log_mu); REPORT(mu)
    sigma <- exp(log_sigma); REPORT(sigma)

    beta <- cbind(beta0, w)
    Gamma <- tpm_g(X, beta)
    delta <- c(0.5, 0.5)

    lallprobs <- matrix(0, length(step), 2)
    ind <- which(!is.na(step))
    for(j in 1:2) {
      lallprobs[ind,j] <- dgamma2(step[ind], mu[j], sigma[j], log = TRUE)
    }

    # HMM likelihood
    nll <- -forward_g(delta, Gamma, lallprobs, bw = bw, logspace = TRUE)

    # GMRF likelihood
    tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
    kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
    Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
    nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

    nll
  }
  environment(jnll) <- environment()

  # Constructing GMRF approximation to spatial field
  loc <- cbind(data$x, data$y)  # Spatial coordinates
  bnd1 <- fmesher::fm_nonconvex_hull(loc, convex = -0.03)
  bnd2 <- fmesher::fm_nonconvex_hull(loc, convex = -0.2)
  mesh <- fmesher::fm_mesh_2d(
    loc=loc,
    boundary=list(bnd1, bnd2),
    max.edge=c(5, 50),
    cutoff=1,
    plot.delay=0.5
  )
  spde <- fmesher::fm_fem(mesh)
  X <- fmesher::fm_basis(mesh, loc)

  dat <- list(
    step = data$step,
    X = X,
    c0 = spde$c0, g1 = spde$g1, g2 = spde$g2,
    bw = bw
  )

  # initial parameter list
  range <- max(dist(loc)) / 5
  kappa0 <- sqrt(8)/range
  sigma0 <- 3 # initial marginal variance of the field
  tau0 <- 1 / (sigma0 * sqrt(4*pi) * kappa0)

  # initialise with true values
  par0 <- list(
    log_mu = log(par$mu),
    log_sigma = log(par$sigma),
    beta0 = par$beta0,
    w = matrix(0, 2, nrow(spde$c0)),
    log_tau_sq = log(tau0^2),
    log_kappa_sq = log(kappa0^2)
  )

  # fit model
  t1 <- Sys.time()
  obj <- MakeADFun(jnll, par0, random = "w")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  time <- Sys.time() - t1

  mod <- report(obj)
  mod$fitting_time <- time
  mod$beta0 <- mod$beta[,1]

  #### Compute quantities of interest ####
  # linear predictors at observed locations
  Eta <- as.matrix(cbind(1, X) %*% t(mod$beta))
  eta21_est <- Eta[-1,1]
  eta12_est <- Eta[-1,2]
  eta21 <- par$beta0[1] + data$field_val[-1] # true linear predictor in simulation
  eta21_diff <- mean(eta21_est - eta21) # mean difference
  eta12_diff <- mean(eta12_est - par$beta0[2]) # mean difference

  # linear predictors on grid
  X_grid <- fmesher::fm_basis(mesh, grid)
  Eta_grid <- as.matrix(cbind(1, X_grid) %*% t(mod$beta))
  eta21_est_grid <- matrix(Eta_grid[,1], nrow = length(x_grid), ncol = length(y_grid))
  eta12_est_grid <- matrix(Eta_grid[,2], nrow = length(x_grid), ncol = length(y_grid))
  eta21_grid <- par$beta0[1] + z21

  # construct indicator which data points are inside bnd1
  locs_sf <- sf::st_as_sf(grid_data, coords = c("x", "y"), crs = sf::st_crs(bnd1))
  inside <- sf::st_within(locs_sf, bnd1, sparse = FALSE)[,1]

  # only those inside non-convex hull
  eta21_inside <- as.numeric(eta21_grid)[inside]
  eta21_est_inside <- as.numeric(eta21_est_grid)[inside]
  eta12_est_inside <- as.numeric(eta12_est_grid)[inside]

  # transform from logit scale to transition probability scale
  g21_grid <- plogis(eta21_grid)
  g21_est_grid <- plogis(eta21_est_grid)
  g12_est_grid <- plogis(eta12_est_grid)
  g21_inside <- plogis(eta21_inside)
  g21_est_inside <- plogis(eta21_est_inside)
  g12_est_inside <- plogis(eta12_est_inside)

  # compute correlation and MSE inside non-convex hull
  cor21 <- cor(eta21_inside, eta21_est_inside)
  cor21_g <- cor(g21_inside, g21_est_inside)
  me21 <- mean(eta21_est_inside - eta21_inside)
  me21_g <- mean(g21_est_inside - g21_inside)
  mse21 <- mean((eta21_est_inside - eta21_inside)^2)
  mse21_g <- mean((g21_est_inside - g21_inside)^2)
  me12 <- mean(eta12_est_inside - par$beta[2])
  me12_g <- mean(g12_est_inside - plogis(par$beta[2]))
  mse12 <- mean((eta12_est_inside - par$beta[2])^2)
  mse12_g <- mean((g12_est_inside - plogis(par$beta[2]))^2)

  # construct return object
  ret <- mod[c("mu", "sigma", "beta0",
               "tau", "kappa",
               "fitting_time")]

  ret$eta21_diff <- eta21_diff
  ret$eta12_diff <- eta12_diff
  ret$cor21 <- cor21
  ret$cor21_g <- cor21_g
  ret$me21 <- me21
  ret$me21_g <- me21_g
  ret$mse21 <- mse21
  ret$mse21_g <- mse21_g
  ret$me12 <- me12
  ret$me12_g <- me12_g
  ret$mse12 <- mse12
  ret$mse12_g <- mse12_g

  # convergence information
  ret$opt_convergence <- opt$convergence
  ret$opt_msg <- opt$message

  # Coverage via the approximate joint posterior
  nSamples <- 1000
  sam <- MCreport(obj, nSamples=nSamples)

  # pointwise coverage of the linear predictors at the observed locations
  Xobs <- cbind(1, X)                                  # n x (1 + nNodes)
  b1 <- sapply(1:nSamples, function(i) c(sam$beta0[[i]][1], sam$w[[i]][1,]))
  b2 <- sapply(1:nSamples, function(i) c(sam$beta0[[i]][2], sam$w[[i]][2,]))
  eta21_s <- as.matrix(Xobs[-1, ] %*% b1)              # (n-1) x nSamples
  eta12_s <- as.matrix(Xobs[-1, ] %*% b2)
  ci21 <- t(apply(eta21_s, 1, quantile, probs = c(0.025, 0.975)))
  ci12 <- t(apply(eta12_s, 1, quantile, probs = c(0.025, 0.975)))
  ret$cover21 <- mean(eta21 >= ci21[, 1] & eta21 <= ci21[, 2])
  ret$cover12 <- mean(par$beta0[2] >= ci12[, 1] & par$beta0[2] <= ci12[, 2])

  # coverage of the state-dependent parameters (0/1 per state, this run)
  mu_s <- sapply(sam$log_mu, exp)
  sigma_s <- sapply(sam$log_sigma, exp)
  ci_mu <- apply(mu_s, 1, quantile, probs = c(0.025, 0.975))
  ci_sigma <- apply(sigma_s, 1, quantile, probs = c(0.025, 0.975))
  ret$cover_mu <- as.numeric(par$mu >= ci_mu[1, ] & par$mu <= ci_mu[2, ])
  ret$cover_sigma <- as.numeric(par$sigma >= ci_sigma[1, ] & par$sigma <= ci_sigma[2, ])

  # cleanup
  rm(sam)
  rm(obj)
  gc()

  return(ret)
}

### Safe wrapper to avoid failing inside mclapply()
one_fit_safe <- function(dummy, ...){
  tryCatch(
    one_fit(dummy, ...),
    error = function(e) {
      message("Error in one_fit: ", conditionMessage(e))
      return("failed")
    }
  )
}

### Call one_fit in a loop over bws for the SAME data set
one_rep <- function(dummy, bws, ...) {
  res <- list()
  # loop over bws for same dummy -> same data set
  for(bw in bws) {
    nm <- paste0("bw", bw)
    res[[nm]] <- one_fit_safe(dummy, bw = bw, ...)
  }
  return(res)
}



##### Run simulation #####
nCores <- 1 # number of cores to use
nSim <- 200 # number of data sets to simulate

nObs <- 5000 # time series length (has to be set to 5000, 10000 manually)
bws <- c(2, 5, 10, 15) # bandwidths explored


# parallelise over data sets
# remove commenting to run
# res <- mclapply(1:nSim, one_rep,
#                 bws = bws,
#                 nObs = nObs,
#                 par = true_par,
#                 mc.cores = nCores) # number of cores to parallelise on
#
# # Save results
# nm <- paste0("./simulations/results/results_cover_nObs", nObs, ".rds")
# saveRDS(res, nm)
#
# gc() # global cleanup




#### Results --------------------------------------------------------------
library(tidyr); library(dplyr); library(patchwork); library(ggplot2)
source("utils.R")
col_vio <- c("2" = colorCB[1], "5" = colorCB[2], "10" = colorCB[3], "15" = colorCB[4])

nSim     <- 200
nObs_vec <- c(5000, 10000)
bws      <- c(2, 5, 10, 15)

# safe accessor: NA for failed fits / missing fields
get1 <- function(r, nm, i = 1) {
  if (identical(r, "failed") || is.null(r[[nm]])) return(NA_real_)
  r[[nm]][i]
}

extract_one <- function(nObs) {
  res <- readRDS(paste0("./simulations/results/results_cover_nObs", nObs, ".rds"))
  do.call(rbind, lapply(seq_along(res), function(s)
    do.call(rbind, lapply(bws, function(bw) {
      r <- res[[s]][[paste0("bw", bw)]]
      data.frame(
        nObs = nObs, run = s, bw = bw,
        eta21_diff = get1(r, "eta21_diff"), eta12_diff = get1(r, "eta12_diff"),
        beta0_21 = get1(r, "beta0", 1),     beta0_12 = get1(r, "beta0", 2),
        cor21   = get1(r, "cor21"),   cor21_g = get1(r, "cor21_g"),
        me21    = get1(r, "me21"),    me21_g  = get1(r, "me21_g"),
        mse21   = get1(r, "mse21"),   mse21_g = get1(r, "mse21_g"),
        me12    = get1(r, "me12"),    me12_g  = get1(r, "me12_g"),
        mse12   = get1(r, "mse12"),   mse12_g = get1(r, "mse12_g"),
        mu1     = get1(r, "mu", 1),   mu2     = get1(r, "mu", 2),
        sigma1  = get1(r, "sigma", 1), sigma2 = get1(r, "sigma", 2),
        cover21 = get1(r, "cover21"), cover12 = get1(r, "cover12"),
        cover_mu1 = get1(r, "cover_mu", 1),    cover_mu2 = get1(r, "cover_mu", 2),
        cover_sigma1 = get1(r, "cover_sigma", 1), cover_sigma2 = get1(r, "cover_sigma", 2),
        converged = !identical(r, "failed") && isTRUE(r$opt_convergence == 0)
      )
    }))))
}

res_all <- do.call(rbind, lapply(nObs_vec, extract_one))
res_all <- res_all[res_all$converged, ]          # drop failed / non-converged fits

# root mean squared errors
res_all$rmse21   <- sqrt(res_all$mse21)
res_all$rmse12   <- sqrt(res_all$mse12)
res_all$rmse21_g <- sqrt(res_all$mse21_g)
res_all$rmse12_g <- sqrt(res_all$mse12_g)


#### Visualise results ----------------------------------------------------

# install.packages(c("tidyr", "dplyr", "patchwork", "ggplot2"))

##### Plot for eta_12 and eta_21 #####

eta12 <- ggplot(res_all, aes(x = factor(bw), y = eta12_diff, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ nObs, nrow = 1) +
  theme_light() +
  labs(x = "", y = expression(hat(eta)[12] - eta[12]), fill = "Bandwidth") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = col_vio)

eta21 <- ggplot(res_all, aes(x = factor(bw), y = eta21_diff, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ nObs, nrow = 1) +
  theme_light() +
  labs(x = "Bandwidth", y = expression(hat(eta)[21] - eta[21]),
       fill = "Bandwidth")+
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = col_vio)

# pdf("./figs/sim_eta_tall.pdf", width = 8, height = 5)
(eta12 / eta21) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
# dev.off()




##### Plot for correlation and MSE #####

cor21 <- ggplot(res_all, aes(x = factor(bw), y = cor21_g, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ nObs, nrow = 1) +
  theme_light() +
  labs(x = "", y = expression(Cor(gamma[21], hat(gamma)[21])),
       fill = "Bandwidth")+
  scale_fill_manual(values = col_vio)

rmse21 <- ggplot(res_all, aes(x = factor(bw), y = rmse21_g, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ nObs, nrow = 1) +
  theme_light() +
  labs(x = "Bandwidth", y = expression(RMSE(hat(gamma)[21])),
       fill = "Bandwidth")+
  scale_fill_manual(values = col_vio)

rmse12 <- ggplot(res_all, aes(x = factor(bw), y = rmse12_g, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ nObs, nrow = 1) +
  theme_light() +
  labs(x = "", y = expression(RMSE(hat(gamma)[12])),
       fill = "Bandwidth")+
  scale_fill_manual(values = col_vio)


# pdf("./figs/sim_cor_rmse.pdf", width = 8, height = 6)
(cor21 / rmse12 / rmse21) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
# dev.off()



##### Plot for mu and sigma #####

mu_long <- res_all %>%
  select(bw, nObs, mu1, mu2) %>%
  pivot_longer(cols = c(mu1, mu2),
               names_to = "state",
               values_to = "value")
mu_long$state <- recode(mu_long$state, mu1 = "mu[1]", mu2 = "mu[2]")
true_mu <- data.frame(
  state = c("mu[1]", "mu[1]", "mu[2]", "mu[2]"),
  nObs  = c(5000, 10000, 5000, 10000),
  true  = c(0.2, 0.2, 5, 5)
)

mu_plot <- ggplot(mu_long, aes(x = factor(bw), y = value, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_grid(state ~ nObs,
             scales = "free_y",
             labeller = labeller(state = label_parsed)) +
  geom_hline(data = true_mu,
             aes(yintercept = true),
             linetype = "dashed",
             inherit.aes = FALSE) +
  theme_light() +
  labs(x = "",
       y = expression(mu),
       fill = "Bandwidth") +
  scale_fill_manual(values = col_vio)


sigma_long <- res_all %>%
  select(bw, nObs, sigma1, sigma2) %>%
  pivot_longer(cols = c(sigma1, sigma2),
               names_to = "state",
               values_to = "value")
sigma_long$state <- recode(sigma_long$state, sigma1 = "sigma[1]", sigma2 = "sigma[2]")


true_sigma <- data.frame(
  state = c("sigma[1]", "sigma[1]", "sigma[2]", "sigma[2]"),
  nObs  = c(5000, 10000, 5000, 10000),
  true  = c(0.5,0.5, 3,3)
)

sigma_plot <- ggplot(sigma_long, aes(x = factor(bw), y = value, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_grid(state ~ nObs,
             scales = "free_y",
             labeller = labeller(state = label_parsed)) +
  geom_hline(data = true_sigma,
             aes(yintercept = true),
             linetype = "dashed",
             inherit.aes = FALSE) +
  theme_light() +
  labs(x = "Bandwidth",
       y = expression(sigma),
       fill = "Bandwidth") +
  scale_fill_manual(values = col_vio)

# pdf("./figs/sim_mu_sigma.pdf", width = 8, height = 6)
(mu_plot / sigma_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
# dev.off()


coverage_summary <- res_all %>%
  group_by(nObs, bw) %>%
  summarise(across(c(cover21, cover12, cover_mu1, cover_mu2,
                     cover_sigma1, cover_sigma2), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")
print(as.data.frame(round(coverage_summary, 2)))
