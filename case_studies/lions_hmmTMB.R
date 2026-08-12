####### Case study - lions in Kalahari; fitting with hmmTMB #######

# Installing and loading required packages
# install.packages(c("hmmTMB", "viridis", "scales"))
library(hmmTMB)
library(viridis)    # for field color palette
library(scales)     # for semi-transparent color


# Data and satellite EDA --------------------------------------------------

# Loading data
data <- readRDS("./data/lions.rds")


# Homogeneous 2-state model -----------------------------------------------

# Specify Markov chain
hid0 <- MarkovChain$new(
  data = data,
  n_states = 2,
  initial_state = "stationary"
)

# Specify observation model
obs0 <- Observation$new(
  data = data,
  dists = list(
    step = "zigamma2",
    angle = "wrpcauchy"
  ),
  par = list(
    step = list(mean = c(0.05, 1), sd = c(0.05, 1), z = c(0.001, 0.001)),
    angle = list(mu = c(3.14, 0), rho = c(0.2, 0.3))
  )
)

# Initialise simple HMM
hmm0 <- HMM$new(hid = hid0, obs = obs0)

# Fit model
hmm0$fit()
hmm0



# Spatial model -----------------------------------------------------------

# Formula for Transition active -> resting
fml <- "~ s(x_int, y_int, bs = 'gp', k = 70)"
# Matrix with transition-specific formulas
form <- matrix(c(
  ".", "~ 1",
  fml, "."
), 2, byrow = TRUE)

# Specify Markov chain
hid <- MarkovChain$new(
  formula = form,
  data = data,
  n_states = 2,
  initial_state = "stationary"
)
hid$update_lambda(5) # larger initial smoothing parameter

# Specify observation model
obs <- Observation$new(
  data = data,
  dists = list(
    step = "zigamma2",
    angle = "wrpcauchy"
  ),
  par = list(
    step = list(mean = c(0.007, 0.9), sd = c(0.007, 1), z = c(0.002, 0.005)),
    angle = list(mu = c(3.141, 0), rho = c(0.2, 0.3))
  )
)

# Initialise HMM
hmm <- HMM$new(hid = hid, obs = obs)

# Fit model
s <- Sys.time()
hmm$fit()
(time <- Sys.time() - s)

# create 256 by 256 prediction grid
x_seq <- seq(min(data$x_int), max(data$x_int), length.out = 256)
y_seq <- seq(min(data$y_int), max(data$y_int), length.out = 256)
grid <- expand.grid(x_int = x_seq, y_int = y_seq)

# Predict spatial field
Ghat <- hmm$predict(
  what = "tpm",
  newdata = grid
)

# Reshape to matrix
g21 <- Ghat[2,1,]
g21 <- matrix(g21, nrow = length(x_seq), ncol = length(y_seq))

# Plot spatial field
# cairo_pdf("./figs/lions_spatial_field_hmmTMB.pdf", width = 7.365, height = 4.95)
image(x_seq, y_seq, g21,
      xlab = "x", ylab = "y",
      col = viridis(35),
      main = expression(Pr(active~"→"~resting)), bty = "n", asp = 1)
# dev.off()



