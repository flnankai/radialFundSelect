# Regularized inverse-norm implementation. No simulation data are saved.

.rfs_smooth_pilot <- function(Z, projection, d) {
  Z <- .finite_matrix(Z, "Z")
  h <- projection$h
  if (nrow(Z) != length(h) || length(d) != ncol(Z) ||
      any(!is.finite(d)) || any(d <= 0))
    stop("Invalid projection or PQR scale", call. = FALSE)
  h2 <- sum(h * h)
  if (!is.finite(h2) || h2 <= 0)
    .fit_error("IN: intercept is unidentified")
  alpha <- as.vector(crossprod(h, Z)) / h2
  e <- sweep(Z - tcrossprod(h, alpha), 2L, d, "/")
  radii <- sqrt(rowSums(e * e))
  median_radius <- stats::median(radii)
  if (any(!is.finite(radii)) || !is.finite(median_radius) || median_radius <= 0)
    .fit_error("IN: OLS pilot median radius is not positive and finite")
  list(alpha = alpha, radii = radii, median_radius = median_radius)
}

.rfs_fit_smooth_in <- function(Z, projection, d, lambda, q = 0.1,
                                     tol = 1e-8, max_iter = 1000L,
                                     initial_alpha = NULL, keep_trace = FALSE) {
  .fit_controls(q, tol, max_iter)
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
    stop("lambda must be finite and nonnegative", call. = FALSE)
  if (!is.logical(keep_trace) || length(keep_trace) != 1L || is.na(keep_trace))
    stop("keep_trace must be TRUE or FALSE", call. = FALSE)
  Z <- .finite_matrix(Z, "Z")
  T <- nrow(Z)
  N <- ncol(Z)
  h <- projection$h
  if (length(d) != N || any(!is.finite(d)) || any(d <= 0))
    stop("d must be a finite positive vector of length N", call. = FALSE)
  if (length(h) != T || any(!is.finite(h)))
    stop("Projection must have a finite h value per date", call. = FALSE)
  h2 <- sum(h * h)
  if (!is.finite(h2) || h2 <= 0)
    .fit_error("IN: intercept is unidentified")

  # The zero-smoothing control is exactly the original operational method,
  # including its original radius-collapse and failure rules.
  if (lambda == 0) {
    fit <- fit_in_projected(Z, projection, d, q = q,
                           initial_alpha = initial_alpha, tol = tol,
                           max_iter = max_iter)
    e <- sweep(Z - tcrossprod(h, fit$theta), 2L, d, "/")
    r2 <- rowSums(e * e)
    date_weights <- h * h / r2
    fit$alpha <- fit$theta
    fit$A_h <- fit$A
    fit$B_c <- fit$B
    fit$score <- fit$location_residual
    fit$lambda <- lambda
    fit$converged <- TRUE
    fit$failure_reason <- NA_character_
    fit$max_date_weight <- max(date_weights) / sum(date_weights)
    fit$trace <- NULL
    return(fit)
  }

  alpha <- if (is.null(initial_alpha))
    as.vector(crossprod(h, Z)) / h2 else as.numeric(initial_alpha)
  if (length(alpha) != N || any(!is.finite(alpha)))
    stop("initial_alpha must be a finite vector of length N", call. = FALSE)
  lambda2 <- lambda * lambda
  if (!is.finite(lambda2) || lambda2 <= 0)
    stop("lambda squared must be positive and finite", call. = FALSE)
  trace <- if (keep_trace) matrix(NA_real_, max_iter, 6L,
      dimnames = list(NULL, c("iteration", "step", "min_radius_ratio",
                             "normalized_raw_score", "objective", "max_date_weight"))) else NULL
  trace_length <- 0L
  initial_objective <- NA_real_
  previous_objective <- NA_real_
  max_objective_increase <- 0
  fail <- function(reason, iteration, radii = NULL, score = NULL) {
    detail <- list(reason = reason, lambda = lambda, iteration = iteration,
                   omega = projection$omega)
    if (keep_trace) detail$trace <- trace[seq_len(trace_length), , drop = FALSE]
    if (!is.null(radii)) {
      med <- stats::median(radii)
      detail$min_radius_ratio <- if (is.finite(med) && med > 0) min(radii) / med else 0
    }
    if (!is.null(score)) detail$normalized_raw_score <- score
    .fit_error(paste0("IN: ", reason), detail)
  }
  # Positive lambda defines e/(r^2+lambda^2) continuously at e=0. A single
  # zero radius is legitimate and is never converted into a fit failure.
  evaluate <- function(alpha) {
    e <- sweep(Z - tcrossprod(h, alpha), 2L, d, "/")
    r2 <- rowSums(e * e)
    if (any(!is.finite(r2))) return(NULL)
    r <- sqrt(r2)
    median_r <- stats::median(r)
    inv <- 1 / (r2 + lambda2)
    score_vector <- colMeans(e * (h * inv))
    # The ordinary median-radius criterion is retained. In a degenerate
    # repeated-data case with median radius zero, use lambda as its positive,
    # unit-equivariant extension rather than multiplying the score by zero.
    score_scale <- if (median_r > 0) median_r else lambda
    list(e = e, r = r, r2 = r2, median_r = median_r, inv = inv,
         score_scale = score_scale,
         score = sqrt(N) * score_scale * max(abs(score_vector)),
         objective = 0.5 * sum(log(r2 + lambda2)))
  }
  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    state <- evaluate(alpha)
    if (is.null(state)) fail("nonfinite residual radius", iteration)
    if (iteration == 1L) initial_objective <- state$objective
    if (is.finite(previous_objective))
      max_objective_increase <- max(max_objective_increase,
                                    state$objective - previous_objective)
    previous_objective <- state$objective
    date_weights <- h * h * state$inv
    denominator <- sum(date_weights)
    if (!is.finite(denominator) || denominator <= 0)
      fail("nonpositive or nonfinite curvature", iteration, state$r, state$score)
    alpha_new <- colSums(Z * (h * state$inv / denominator))
    if (any(!is.finite(alpha_new)))
      fail("nonfinite alpha update", iteration, state$r, state$score)
    loc_change <- max(abs((alpha_new - alpha) / d)) /
      (state$score_scale / sqrt(N))
    if (keep_trace) {
      trace[iteration, ] <- c(iteration, loc_change,
          if (state$median_r > 0) min(state$r) / state$median_r else 0,
          state$score, state$objective, max(date_weights) / denominator)
      trace_length <- iteration
    }
    alpha <- alpha_new
    if (loc_change < tol) {
      check <- evaluate(alpha)
      if (is.null(check)) fail("nonfinite stopping residual", iteration)
      if (is.finite(check$score) && check$score < tol) {
        converged <- TRUE
        break
      }
    }
  }
  if (!converged)
    fail(paste0("did not converge in ", max_iter, " iterations"),
         max_iter, state$r, state$score)
  state <- evaluate(alpha)
  max_objective_increase <- max(max_objective_increase,
                                state$objective - previous_objective)
  A_h <- mean(h * h * state$inv)
  pi <- if (!ncol(projection$Q)) rep(0, T) else
    as.vector(projection$Q %*% crossprod(projection$Q, h * state$inv))
  # General radial-weight studentization. The inverse-norm-specific identity
  # B_c = A_h + mean(r^2*pi^2) is not valid for positive lambda.
  coefficients <- state$r * (h * state$inv - pi)
  B_c <- mean(coefficients * coefficients)
  residual_df <- projection$residual_df
  if (length(residual_df) != 1L || !is.finite(residual_df) || residual_df <= 0)
    stop("IN requires positive residual degrees of freedom", call. = FALSE)
  if (!is.finite(A_h) || A_h <= 0 || !is.finite(B_c) || B_c <= 0)
    fail("nonpositive or nonfinite fitted variance", iteration, state$r, state$score)
  z_uncalibrated <- sqrt(T * N) * A_h * alpha / (d * sqrt(B_c))
  df_multiplier <- sqrt(residual_df / T)
  z <- df_multiplier * z_uncalibrated
  if (any(!is.finite(z)))
    fail("nonfinite statistic", iteration, state$r, state$score)
  date_weights <- h * h * state$inv
  fit <- .new_fit("IN", alpha, d, z, q, projection$omega, A_h, B_c,
                  iteration, state$score, NaN,
                  if (state$median_r > 0) min(state$r) / state$median_r else 0,
                  df_multiplier, residual_df, z_uncalibrated)
  fit$method <- "IN"
  fit$alpha <- alpha
  fit$A_h <- A_h
  fit$B_c <- B_c
  fit$score <- state$score
  fit$lambda <- lambda
  fit$converged <- TRUE
  fit$failure_reason <- NA_character_
  fit$max_date_weight <- max(date_weights) / sum(date_weights)
  fit$initial_objective <- initial_objective
  fit$final_objective <- state$objective
  fit$max_objective_increase <- max_objective_increase
  fit$trace <- if (keep_trace) trace[seq_len(trace_length), , drop = FALSE] else NULL
  fit
}


.rfs_smoothing <- function(smoothing) {
  if (!is.numeric(smoothing) || length(smoothing) != 1L ||
      !is.finite(smoothing) || smoothing < 0)
    stop("smoothing must be a finite nonnegative number", call. = FALSE)
  smoothing
}

fit_regularized_in_projected <- function(Z, projection, d, q = 0.1,
                                         tol = 1e-8, max_iter = 1000L,
                                         smoothing = 0.5, pilot_radius = NULL) {
  smoothing <- .rfs_smoothing(smoothing)
  if (smoothing == 0) {
    fit <- fit_in_projected(Z, projection, d, q = q, tol = tol, max_iter = max_iter)
    fit$lambda <- 0
    fit$pilot_radius <- NA_real_
    fit$implementation <- "legacy_inverse_norm"
  } else {
    if (inherits(pilot_radius, "error"))
      .fit_error(paste0("IN: pilot failed: ", conditionMessage(pilot_radius)))
    if (is.null(pilot_radius)) pilot_radius <- .rfs_smooth_pilot(Z, projection, d)$median_radius
    if (length(pilot_radius) != 1L || !is.finite(pilot_radius) || pilot_radius <= 0)
      .fit_error("IN: pilot radius must be positive and finite")
    fit <- .rfs_fit_smooth_in(Z, projection, d, smoothing * pilot_radius,
                             q = q, tol = tol, max_iter = max_iter)
    fit$pilot_radius <- pilot_radius
    fit$implementation <- "smoothed_inverse_norm"
  }
  fit$smoothing <- smoothing
  fit
}
