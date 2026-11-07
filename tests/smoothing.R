library(radialFundSelect)

assert_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}

# A symmetric sample with one exactly zero residual is well-defined for
# positive smoothing; the legacy algorithm must retain its old failure rule.
v1 <- c(1, 2, 3, 4); v2 <- c(-2, 1, 4, 3); v3 <- c(2, 4, -1, 1)
special <- rbind(rep(0, 4), v1, -v1, v2, -v2, v3, -v3)
smooth <- select_funds(special)
legacy <- select_funds(special, smoothing = 0)
stopifnot(smooth$fit_converged, smooth$fit$min_radius_ratio == 0,
          all(is.finite(smooth$fit$z)), all(smooth$fit$z == 0),
          !legacy$fit_converged, !any(legacy$selected))

# Verify the general variance against a direct QR/projection calculation.
# This would reject reuse of the special inverse-norm identity at lambda>0.
panel <- simulate_panel(distribution = "t3", rho = .9, delta = .48, replicate = 9)
result <- select_funds(panel$Y, panel$factors)
old <- select_funds(panel$Y, panel$factors, smoothing = 0)
fit <- result$fit
stopifnot(result$fit_converged, !old$fit_converged,
          fit$smoothing == .5, fit$implementation == "smoothed_inverse_norm",
          sum(result$selected) == 20, all(panel$alpha[result$selected] > 0))
Q <- qr.Q(qr(panel$factors))
h <- 1 - as.vector(Q %*% colSums(Q))
Z <- panel$Y - Q %*% crossprod(Q, panel$Y)
e <- sweep(Z - tcrossprod(h, fit$theta), 2, fit$d, "/")
r2 <- rowSums(e^2)
a <- 1 / (r2 + fit$lambda^2)
pi <- as.vector(Q %*% crossprod(Q, h * a))
A <- mean(h^2 * a)
B <- mean(r2 * (h * a - pi)^2)
ols <- lm.fit(cbind(1, panel$factors), panel$Y)
pilot_r <- sqrt(rowSums(sweep(ols$residuals, 2, fit$d, "/")^2))
stopifnot(abs(A - fit$A) < 1e-12, abs(B - fit$B) < 1e-12,
          abs(fit$lambda - .5 * median(pilot_r)) < 1e-12,
          abs(B - (A + mean(r2 * pi^2))) > 1e-5,
          fit$location_residual < 1e-8)
three <- fit_three(panel$Y, panel$factors, capture_errors = TRUE)
stopifnot(identical(three$IN$z, fit$z), identical(three$IN$selected, result$selected))

# Defaults and explicit strengths propagate through the public simulation
# entry point. Independently formed panels recover the operational averages.
cases <- scenario_inventory(); cases <- cases[c(1, 27, 53), ]
default <- simulate_study(reps = 2, workers = 1, scenarios = cases, progress = FALSE)
explicit <- simulate_study(reps = 2, workers = 1, scenarios = cases,
                           smoothing = .5, progress = FALSE)
control <- simulate_study(reps = 2, workers = 1, scenarios = cases,
                          smoothing = 0, progress = FALSE)
stopifnot(identical(default$summary, explicit$summary),
          identical(default$paired_differences, explicit$paired_differences),
          default$config$IN_smoothing == .5,
          control$config$IN_implementation == "legacy_inverse_norm")
for (i in seq_len(nrow(cases))) {
  vals <- vapply(1:2, function(r) {
    p <- simulate_panel(T = cases$T[i], distribution = cases$distribution[i],
                        rho = cases$rho[i], delta = cases$delta[i], replicate = r)
    a <- select_funds(p$Y, p$factors)
    c(FDR = sum(a$selected & p$alpha <= 0) / max(sum(a$selected), 1),
      TDP = sum(a$selected & p$alpha > 0) / sum(p$alpha > 0),
      failed = !a$fit_converged)
  }, numeric(3))
  row <- default$summary[default$summary$case_id == cases$case_id[i] &
                          default$summary$method == "IN", ]
  stopifnot(abs(row$FDR - mean(vals[1, ])) < 1e-12,
            abs(row$TDP - mean(vals[2, ])) < 1e-12,
            row$failed == sum(vals[3, ]))
}
for (bad in list(-.1, NA_real_, Inf, TRUE, c(.1,.5))) {
  assert_error(select_funds(special, smoothing = bad))
  assert_error(fit_three(special, smoothing = bad))
  assert_error(simulate_study(reps = 1, workers = 1, smoothing = bad))
}
cat("Smoothed-default, variance, zero-radius, legacy-control, and simulation checks passed.\n")
