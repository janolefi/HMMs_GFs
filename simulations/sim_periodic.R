####### Simulation experiment - additive GP trend #######

# Installing required packages
# Versions used:
# - RTMB: 1.9
# - LaMa: 2.1.3
# - RTMBdist: 1.0.6
# - fmesher: 0.8.0
# install.packages(c("RTMB", "LaMa", "RTMBdist", "fmesher"))


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
Q <- tau0^2 * (kappa0^4 * spde$c0 + 2 * cos(pi*omega0) * kappa0^2 * spde$g1 + spde$g2)

true_field <- as.numeric(rgmrf(1, 0, Q))
plot(true_field, type = "l", lwd = 2, bty = "n")


# True parameters used for simulation
true_par <- list(
  mu = c(-1, 3),
  sigma = c(1, 2),
  beta0 = rep(qlogis(0.05), 2),
  tau = tau0,
  kappa = kappa0,
  omega = omega0
)

# Plot true state-dependent distributions
curve(dnorm(x, true_par$mu[1], true_par$sigma[1]), xlim = c(-4, 8),
      n = 500, bty = "n", ylab = "Density", xlab = "X", lwd = 3, col = colorCB[1])
curve(dnorm(x, true_par$mu[2], true_par$sigma[2]), add = TRUE, n = 500,
      lwd = 3, col = colorCB[2])



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

### Function that simulates data set and fits model for a given specification
one_fit <- function(dummy,
                    par = true_par,
                    bw = 10,
                    spde = spde) {
  set.seed(dummy) # for reproducibility of each iteration

  library(LaMa) # also loads RTMB

  # Simulate data
  data <- sim_data(par, spde = spde)

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
  environment(jnll) <- environment()

  dat <- list(
    y = data$y,
    c0 = spde$c0, g1 = spde$g1, g2 = spde$g2,
    bw = bw
  )

  # initialise with true values
  par0 <- list(
    mu = par$mu,
    log_sigma = log(par$sigma),
    beta0 = par$beta0,
    w = rep(0, nrow(spde$c0)),
    log_tau_sq = log(par$tau^2),
    log_kappa_sq = log(par$kappa^2),
    u = qlogis((1 + cos(pi * par$omega)) / 2)
  )

  # fit model
  t1 <- Sys.time()
  obj <- MakeADFun(jnll, par0, random = "w")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  time <- Sys.time() - t1

  mod <- report(obj)
  mod$fitting_time <- time

  # construct return object
  ret <- mod[c("mu", "sigma", "Gamma",
               "tau", "kappa", "omega",
               "fitting_time")]

  # convergence information
  ret$opt_convergence <- opt$convergence
  ret$opt_msg <- opt$message

  # point estimate of the trend
  f_hat <- mod$w
  true_field <- data$f

  # --- field (trend) recovery ---
  ret$cor_field <- cor(f_hat, true_field)                 # offset-invariant
  ret$me_field  <- mean(f_hat - true_field)               # level confounded with mu (see note)
  ret$mse_field <- mean((f_hat - true_field)^2)

  # --- coverage via the approximate joint posterior ---
  nSamples <- 1000
  sam <- MCreport(obj, nSamples = nSamples)

  # field: pointwise coverage of the trend
  f_draws  <- sapply(sam$w, identity)                    # nObs x nSamples
  ci_field <- t(apply(f_draws, 1, quantile, probs = c(0.025, 0.975)))
  ret$cover_field <- mean(true_field >= ci_field[, 1] & true_field <= ci_field[, 2])

  # state-dependent parameters (0/1 per state)
  mu_s    <- sapply(sam$mu, identity)                    # 2 x nSamples
  sigma_s <- sapply(sam$log_sigma, exp)
  ci_mu    <- apply(mu_s,    1, quantile, probs = c(0.025, 0.975))
  ci_sigma <- apply(sigma_s, 1, quantile, probs = c(0.025, 0.975))
  ret$cover_mu    <- as.numeric(par$mu    >= ci_mu[1, ]    & par$mu    <= ci_mu[2, ])
  ret$cover_sigma <- as.numeric(par$sigma >= ci_sigma[1, ] & par$sigma <= ci_sigma[2, ])

  # field hyperparameters (scalars, 0/1)
  tau_s   <- sqrt(exp(unlist(sam$log_tau_sq)))
  kappa_s <- sqrt(exp(unlist(sam$log_kappa_sq)))
  omega_s <- acos(2 * plogis(unlist(sam$u)) - 1) / pi
  ci_tau   <- quantile(tau_s,   probs = c(0.025, 0.975))
  ci_kappa <- quantile(kappa_s, probs = c(0.025, 0.975))
  ci_omega <- quantile(omega_s, probs = c(0.025, 0.975))
  ret$cover_tau   <- as.numeric(par$tau   >= ci_tau[1]   & par$tau   <= ci_tau[2])
  ret$cover_kappa <- as.numeric(par$kappa >= ci_kappa[1] & par$kappa <= ci_kappa[2])
  ret$cover_omega <- as.numeric(par$omega >= ci_omega[1] & par$omega <= ci_omega[2])

  # transition probabilities (off-diagonal, 0/1 per row)
  gamma_s  <- sapply(sam$beta0, plogis)                  # 2 x nSamples
  ci_gamma <- apply(gamma_s, 1, quantile, probs = c(0.025, 0.975))
  ret$cover_gamma <- as.numeric(plogis(par$beta0) >= ci_gamma[1, ] &
                                  plogis(par$beta0) <= ci_gamma[2, ])

  rm(sam); rm(obj); gc()
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
nCores <- 4 # number of cores to use
nSim <- 500 # number of data sets to simulate

bws <- c(2, 5, 10, 15) # bandwidths explored

# parallelise over data sets
# remove commenting to run
# res <- mclapply(1:nSim, one_rep,
#                 bws = bws,
#                 par = true_par,
#                 spde = spde,
#                 mc.cores = nCores) # number of cores to parallelise on
#
# # Save results
# saveRDS(res, "./simulations/results/periodic_sim2_500.rds")
res <- readRDS("./simulations/results/periodic_sim2_500.rds")

##### Extract & plot results #####
library(tidyr); library(dplyr); library(ggplot2); library(patchwork)

# res is already in memory from the run above (or: res <- readRDS("./simulations/results/periodic_sim.rds"))

# helper: pull scalar / indexed element, NA for failed fits
pull <- function(r, nm, i = NULL) {
  if (identical(r, "failed") || is.null(r[[nm]])) return(NA_real_)
  if (is.null(i)) r[[nm]] else r[[nm]][i]
}

res_all <- do.call(rbind, lapply(seq_along(res), function(s)
  do.call(rbind, lapply(bws, function(bw) {
    r <- res[[s]][[paste0("bw", bw)]]
    data.frame(
      run = s, bw = bw,
      mu1 = pull(r,"mu",1),       mu2 = pull(r,"mu",2),
      sigma1 = pull(r,"sigma",1), sigma2 = pull(r,"sigma",2),
      tau = pull(r,"tau"), kappa = pull(r,"kappa"), omega = pull(r,"omega"),
      gamma12 = if (identical(r,"failed")) NA_real_ else r$Gamma[1,2],
      gamma21 = if (identical(r,"failed")) NA_real_ else r$Gamma[2,1],
      cor_field  = pull(r,"cor_field"),
      rmse_field = sqrt(pull(r,"mse_field")),
      me_field   = pull(r,"me_field"),
      cover_field = pull(r,"cover_field"),
      cover_mu1 = pull(r,"cover_mu",1), cover_mu2 = pull(r,"cover_mu",2),
      cover_sigma1 = pull(r,"cover_sigma",1), cover_sigma2 = pull(r,"cover_sigma",2),
      cover_tau = pull(r,"cover_tau"), cover_kappa = pull(r,"cover_kappa"),
      cover_omega = pull(r,"cover_omega"),
      cover_gamma1 = pull(r,"cover_gamma",1), cover_gamma2 = pull(r,"cover_gamma",2),
      converged = !identical(r,"failed") && r$opt_convergence == 0
    )
  }))))
res_all <- res_all[res_all$converged, ]

col_vio <- c("2"=colorCB[1], "5"=colorCB[2], "10"=colorCB[3], "15"=colorCB[4])

### Coverage summary (mean over runs, by bandwidth)
coverage_summary <- res_all %>%
  group_by(bw) %>%
  summarise(across(starts_with("cover_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
print(as.data.frame(round(coverage_summary, 3)))

### mu and sigma
mu_long <- res_all %>% select(bw, mu1, mu2) %>%
  pivot_longer(c(mu1, mu2), names_to = "state", values_to = "value") %>%
  mutate(state = recode(state, mu1 = "mu[1]", mu2 = "mu[2]"))
true_mu <- data.frame(state = c("mu[1]","mu[2]"), true = true_par$mu)

mu_plot <- ggplot(mu_long, aes(factor(bw), value, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ state, scales = "free_y", labeller = label_parsed) +
  geom_hline(data = true_mu, aes(yintercept = true), linetype = "dashed", inherit.aes = FALSE) +
  theme_light() + labs(x = "", y = expression(mu), fill = "Bandwidth") +
  scale_fill_manual(values = col_vio)

sigma_long <- res_all %>% select(bw, sigma1, sigma2) %>%
  pivot_longer(c(sigma1, sigma2), names_to = "state", values_to = "value") %>%
  mutate(state = recode(state, sigma1 = "sigma[1]", sigma2 = "sigma[2]"))
true_sigma <- data.frame(state = c("sigma[1]","sigma[2]"), true = true_par$sigma)

sigma_plot <- ggplot(sigma_long, aes(factor(bw), value, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ state, scales = "free_y", labeller = label_parsed) +
  geom_hline(data = true_sigma, aes(yintercept = true), linetype = "dashed", inherit.aes = FALSE) +
  theme_light() + labs(x = "Bandwidth", y = expression(sigma), fill = "Bandwidth") +
  scale_fill_manual(values = col_vio)

# pdf("./figs/sim_p_mu_sigma.pdf", width = 8, height = 6)
(mu_plot / sigma_plot) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# dev.off()


### Field hyperparameters + transition probabilities
lab_map <- c(gamma12 = "gamma[12]", gamma21 = "gamma[21]",
             kappa = "kappa", tau = "tau", omega = "omega")
lab_levels <- c("gamma[12]", "gamma[21]", "kappa", "tau", "omega")

scal_long <- res_all %>%
  select(bw, gamma12, gamma21, kappa, tau, omega) %>%
  pivot_longer(-bw, names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(lab_map[parameter], levels = lab_levels))
true_scal <- data.frame(
  parameter = factor(lab_levels, levels = lab_levels),
  true = c(plogis(true_par$beta0[1]), plogis(true_par$beta0[2]),
           true_par$kappa, true_par$tau, true_par$omega))

scal_plot <- ggplot(scal_long, aes(factor(bw), value, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ parameter, scales = "free_y", labeller = label_parsed, nrow = 1) +
  geom_hline(data = true_scal, aes(yintercept = true), linetype = "dashed", inherit.aes = FALSE) +
  theme_light() + labs(x = "Bandwidth", y = "", fill = "Bandwidth") +
  scale_fill_manual(values = col_vio)

# pdf("./figs/sim_p_gamma_spde.pdf", width = 9, height = 3.5)
scal_plot + theme(legend.position = "bottom")
# dev.off()


### Field (trend) recovery
field_long <- res_all %>%
  select(bw, cor_field, rmse_field) %>%
  pivot_longer(-bw, names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         cor_field  = "Cor(u, hat(u))",
                         rmse_field = "RMSE(hat(u))"))

field_plot <- ggplot(field_long, aes(factor(bw), value, fill = factor(bw))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y", labeller = label_parsed) +
  theme_light() + labs(x = "Bandwidth", y = "", fill = "Bandwidth") +
  scale_fill_manual(values = col_vio)

# pdf("./figs/sim_p_field.pdf", width = 8, height = 3.5)
field_plot + theme(legend.position = "bottom")
# dev.off()


# average corrlation and RMSE
mean(res_all$cor_field[res_all$bw == 15],  na.rm = TRUE)
mean(res_all$rmse_field[res_all$bw == 15], na.rm = TRUE)



