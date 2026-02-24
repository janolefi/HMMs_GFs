r <- function(d, omega = 0.98, kappa = 5) {
  1 / (2 * sin(pi * omega) * kappa^2) *
    exp(-kappa * cos(pi * omega / 2) * d) *
    sin(pi * omega / 2 + kappa * sin(pi * omega / 2) * d)
}

pdf("./figs/osc_cov.pdf", width = 7, height = 4)
par(mfrow = c(1,1))
curve(r(x, 0.98, 7), xlim = c(0, 10), lwd = 2, n = 1000,
      bty = "n", xlab = "d = |u - v|", ylab = "r(u, v)")
dev.off()
