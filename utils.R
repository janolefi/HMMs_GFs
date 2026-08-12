#### Colors used for plotting

# color-blind friendly palette
colorCB <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7", "#F0E442", "#0072B2", "#D55E00")

# colors of the sun cycle
sun_cycle_colors <- c(
  "#0b0d3e",  # 00:00 - Midnight (Night)
  "#0d1046",  # 00:30
  "#0f124e",  # 01:00
  "#121557",  # 01:30
  "#14185f",  # 02:00
  "#161a67",  # 02:30
  "#191d6f",  # 03:00
  "#1b1f78",  # 03:30
  "#1d227f",  # 04:00
  "#1f2587",  # 04:30
  "#22288f",  # 05:00
  "#564d7b",  # 05:30
  "#877060",  # 06:00
  "#e67e22",  # 06:30 - Sunrise (Vibrant Orange)
  "#f28c32",  # 07:00
  "#f89b42",  # 07:30
  "#fba554",  # 08:00 - Morning
  "#fcaf67",  # 08:30
  "#fcc978",  # 09:00
  "#fcd68d",  # 09:30
  "#fce3a2",  # 10:00
  "#fde1b7",  # 10:30
  "#feeacb",  # 11:00
  "#d8eef5",  # 11:30 - Gradual shift from morning to day
  "#b7e9f9",  # 12:00 - Noon (Light Blue)
  "#a2e2f9",  # 12:30
  "#8cdcf9",  # 13:00
  "#75d5fa",  # 13:30
  "#5ecefa",  # 14:00
  "#47c7fa",  # 14:30
  "#30c0fb",  # 15:00
  "#2ca9e2",  # 15:30
  "#2593cb",  # 16:00 - Afternoon (Deeper Blue)
  "#1f7db4",  # 16:30
  "#19679d",  # 17:00
  "#7f6d63",  # 17:30
  "#e57328",  # 18:00 - Sunset (Vibrant Orange)
  "#cc6925",  # 18:30
  "#b36022",  # 19:00
  "#9a571f",  # 19:30
  "#7c4550",  # 20:00 - Nightfall
  "#603a59",  # 20:30
  "#4c2f55",  # 21:00
  "#38264b",  # 21:30
  "#291f3f",  # 22:00 - Nighttime
  "#1c1732",  # 22:30
  "#140f28",  # 23:00 - Late Night
  "#0e0a1e"   # 23:30
)


# Computes log-likelihood at current parameter value for different bandwidths and plots the result
# this aids choosing an appropriate bandwidth
# bw should be chosen sufficiently large that the log-likelihood value has stabilised
bandwidth_check <- function(obj, bws = 10:40) {
  mod <- obj$report() # run the reporting to extract required quantities

  # informative error if obj is not a likelihood built based on forward()
  if(is.null(mod$delta) | is.null(mod$Gamma) | is.null(mod$allprobs)) {
    stop("obj does not seem to be an HMM objective function; have you called forward() in the likelihood?")
  }

  # loop over bandwidths and store log-likelihood value
  llks <- rep(NA, length(bws))
  for(i in 1:length(bws)) {
    llks[i] <- forward(mod$delta, mod$Gamma, mod$allprobs, mod$trackID,
                       ad = TRUE, bw = bws[i])
  }

  # plot the profile
  par(mfrow = c(1,1))
  plot(bws, llks, type = "l", bty = "n",
       xlab = "Bandwidth", ylab = "Log-likelihood",
       col = "darkblue", lwd = 2)
  invisible(data.frame(bw = bws, llk = llks))
}
