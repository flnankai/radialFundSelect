library(radialFundSelect)

assert_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}

# Six pairwise distances have middle values 3 and 4. The inverse empirical
# CDF must choose 3, not their interpolated median 3.5.
Y_small <- outer(c(0, 1, 3, 7), c(1, 2, 3))
scale_small <- pqr_scale(Y_small)
stopifnot(identical(scale_small$order_statistic, 3L),
          scale_small$pair_count == 6,
          isTRUE(all.equal(scale_small$quantiles, c(3, 6, 9))),
          abs(mean(scale_small$d^2) - 1) < 1e-14)

# BH is a step-up rule: the first ordered p-value may fail its own
# threshold yet still be selected when the second passes. Ties are kept.
stopifnot(identical(bh(c(.035, .04, .5, .9)), c(TRUE, TRUE, FALSE, FALSE)),
          identical(bh(c(.02, .02, .5, .9)), c(TRUE, TRUE, FALSE, FALSE)),
          identical(bh(c(.025, .8, .9, 1)), c(TRUE, FALSE, FALSE, FALSE)),
          identical(bh(rep(1, 4)), rep(FALSE, 4)),
          identical(bh(numeric()), logical()))
assert_error(bh(c(.2, NA_real_)))
assert_error(bh(c(.2, 1.01)))
assert_error(bh(.1, q = 0))

set.seed(827)
T <- 60L
N <- 40L
F <- matrix(rnorm(T * 3L), T, 3L)
Y <- matrix(rnorm(T * N), T, N)
Y <- Y + F %*% matrix(rnorm(3L * N), 3L, N) +
  outer(rep(1, T), seq(-.2, .3, length.out = N))

# Partial regression must give precisely the ordinary OLS intercept and
# the original RSS/T statistic, not RSS/(T-p-1).
ls <- select_funds(Y, F, method = "LS")
ols <- lm.fit(cbind(1, F), Y)
alpha_ols <- ols$coefficients[1, ]
Q <- qr.Q(qr(F))
h <- 1 - as.vector(Q %*% colSums(Q))
sigma_ols <- sqrt(colMeans(ols$residuals^2))
z_ols <- sqrt(sum(h^2)) * alpha_ols / sigma_ols
stopifnot(ls$fit_converged,
          max(abs(ls$fit$alpha_first_order - alpha_ols)) < 1e-12,
          max(abs(ls$fit$d - sigma_ols)) < 1e-12,
          max(abs(ls$fit$z - z_ols)) < 1e-11,
          max(abs(ls$pvalue - pnorm(z_ols, lower.tail = FALSE))) < 1e-12)

# PQR is exactly invariant to adding an arbitrary intercept vector and
# an arbitrary observed-factor component (up to numerical projection error).
Y_shift <- Y + outer(rep(1, T), seq_len(N) / 7) +
  F %*% matrix(seq_len(3L * N) / 60, 3L, N)
scale_original <- pqr_scale(Y, F)
scale_shift <- pqr_scale(Y_shift, F)
stopifnot(max(abs(scale_original$quantiles - scale_shift$quantiles)) < 1e-12,
          max(abs(scale_original$d - scale_shift$d)) < 1e-12)

# Calibration and convergence criteria are invariant to the return unit.
in_fit <- select_funds(Y, F, method = "IN")
tiny_fit <- select_funds(Y * 1e-8, F, method = "IN")
stopifnot(in_fit$fit_converged, tiny_fit$fit_converged,
          max(abs(in_fit$fit$z - tiny_fit$fit$z)) < 1e-9,
          identical(in_fit$selected, tiny_fit$selected),
          in_fit$fit$residual_df == T - 4L,
          abs(in_fit$fit$df_multiplier - sqrt((T - 4L) / T)) < 1e-15,
          max(abs(in_fit$fit$z - in_fit$fit$z_uncalibrated *
                    sqrt((T - 4L) / T))) < 1e-13,
          in_fit$fit$location_residual < 1e-8)

fits <- fit_three(Y, F, capture_errors = TRUE)
stopifnot(identical(names(fits), c("LS", "SS", "IN")),
          all(!vapply(fits, inherits, logical(1), "error")),
          max(abs(fits$SS$alpha_first_order - fits$SS$theta / fits$SS$omega)) == 0)

# Failed fits remain explicit and make no selections; they are not
# treated as successful fits or silently omitted from simulation means.
for (method in c("LS", "SS", "IN")) {
  failed <- select_funds(matrix(0, 6, 5), method = method)
  stopifnot(!failed$fit_converged, is.null(failed$fit),
            all(failed$pvalue == 1), !any(failed$selected),
            identical(failed$diagnostics$failure_policy, "no rejections"))
}
for (method in c("SS", "IN")) {
  unfinished <- select_funds(Y, F, method = method, max_iter = 1L)
  stopifnot(!unfinished$fit_converged, !any(unfinished$selected),
            all(unfinished$pvalue == 1))
}
failed_three <- fit_three(matrix(0, 6, 5), capture_errors = TRUE)
stopifnot(all(vapply(failed_three, inherits, logical(1), "error")))
assert_error(select_funds(Y, cbind(1, F)))
assert_error(select_funds(Y, cbind(F, F[, 1])))
assert_error(select_funds(Y, F, max_iter = 1.5))
assert_error(select_funds(Y, F, tol = 0))
assert_error(pqr_scale(Y, F, quantile = 1))

cat("Algorithm and failure-policy checks passed.\n")
