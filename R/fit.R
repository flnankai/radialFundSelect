.new_fit <- function(method, theta, d, z, q, omega, A = NaN, B = NaN,
                     iterations = 0L, location_residual = NaN,
                     scale_residual = NaN, min_radius_ratio = NaN,
                     df_multiplier = 1, residual_df = NULL,
                     z_uncalibrated = NULL) {
  pvalue <- stats::pnorm(z, lower.tail = FALSE)
  structure(list(method = method, theta = theta, d = d, z = z,
    pvalue = pvalue, selected = bh(pvalue, q), omega = omega, kappa = 1 - omega,
    A = A, B = B, iterations = iterations, location_residual = location_residual,
    scale_residual = scale_residual, min_radius_ratio = min_radius_ratio,
    df_multiplier = df_multiplier, residual_df = residual_df,
    z_uncalibrated = z_uncalibrated,
    alpha_first_order = if (method == "SS") theta / omega else theta,
    scale_estimator = switch(method, LS = "OLS_RSS_over_T",
                             SS = "joint_spatial_sign", IN = "fixed_PQR"),
    fit_converged = TRUE), class = "radial_fit")
}

fit_ls_projected <- function(Z, projection, q = 0.1) {
  Z <- .finite_matrix(Z, "Z")
  T <- nrow(Z)
  h <- projection$h
  if (T != length(h)) stop("Invalid projected observation matrix", call. = FALSE)
  h2 <- sum(h * h)
  if (h2 <= 0) .fit_error("LS: intercept is unidentified")
  alpha <- as.vector(crossprod(h, Z)) / h2
  residual <- Z - tcrossprod(h, alpha)
  variance <- colMeans(residual * residual)
  if (any(variance <= 0) || any(!is.finite(variance)))
    .fit_error("LS: nonpositive or nonfinite marginal OLS residual variance")
  d <- sqrt(variance)
  z <- sqrt(h2) * alpha / d
  .new_fit("LS", alpha, d, z, q, projection$omega,
           location_residual = max(abs(crossprod(h, residual))) / (T * max(d)))
}

fit_ls <- function(Y, factors = NULL, q = 0.1) {
  Y <- .finite_matrix(Y, "Y")
  projection <- observed_projection(nrow(Y), factors)
  fit_ls_projected(project_observations(Y, projection), projection, q)
}

radial_calibration <- function(r, method = "SS", kappa) {
  if (method != "SS") stop("This calibration belongs to the original SS fit only")
  A <- mean(1 / r)
  B <- mean((1 - kappa * A * r)^2) +
    kappa * (1 - kappa) * A * A * mean(r * r)
  if (!is.finite(A + B) || A <= 0 || B <= 0)
    .fit_error("Nonpositive or nonfinite radial calibration")
  list(A = A, B = B)
}

fit_projected <- function(Z, omega, method = "SS", q = 0.1,
                          initial = NULL, tol = 1e-8, max_iter = 1000L) {
  .fit_controls(q, tol, max_iter)
  if (method != "SS") stop("IN requires fit_in_projected with the fixed PQR scale")
  Z <- .finite_matrix(Z, "Z")
  if (length(omega) != 1L || !is.finite(omega) || omega <= 0 || omega > 1)
    stop("omega must lie in (0,1]", call. = FALSE)
  T <- nrow(Z)
  N <- ncol(Z)
  kappa <- 1 - omega
  if (is.null(initial)) {
    theta <- colMeans(Z)
    centered <- sweep(Z, 2L, theta, "-")
    d <- sqrt(colSums(centered * centered) / (T - 1))
  } else {
    theta <- initial[[1L]]
    d <- initial[[2L]]
    if (length(theta) != N || length(d) != N || any(!is.finite(theta)))
      stop("initial must contain finite location and scale vectors", call. = FALSE)
  }
  if (any(d <= 0) || any(!is.finite(d)))
    .fit_error("Initial coordinate scale is degenerate")
  d <- d / exp(mean(log(d)))
  radius_tol <- 100 * .Machine$double.eps
  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    e <- sweep(sweep(Z, 2L, theta, "-"), 2L, d, "/")
    r <- sqrt(rowSums(e * e))
    median_r <- stats::median(r)
    if (any(!is.finite(r)) || min(r) <= radius_tol * median_r)
      .fit_error("SS: collapsed or nonfinite residual radius")
    U <- e / r
    inverse_r <- 1 / r
    theta_new <- colSums(Z * (inverse_r / sum(inverse_r)))
    d_new <- d * sqrt(N * colMeans(U * U))
    if (any(d_new <= 0) || any(!is.finite(d_new)))
      .fit_error("SS: degenerate scale update")
    d_new <- d_new / exp(mean(log(d_new)))
    loc_change <- max(abs((theta_new - theta) / d)) / (median_r / sqrt(N))
    scale_change <- max(abs(log(d_new / d)))
    theta <- theta_new
    d <- d_new
    if (max(loc_change, scale_change) < tol) {
      check_e <- sweep(sweep(Z, 2L, theta, "-"), 2L, d, "/")
      check_r <- sqrt(rowSums(check_e * check_e))
      if (any(!is.finite(check_r)) ||
          min(check_r) <= radius_tol * stats::median(check_r))
        .fit_error("SS: residual radius collapsed at stopping candidate")
      check_U <- check_e / check_r
      location_equation <- sqrt(N) * max(abs(colMeans(check_U)))
      scale_equation <- max(abs(N * colMeans(check_U * check_U) - 1))
      if (max(location_equation, scale_equation) < tol) {
        converged <- TRUE
        break
      }
    }
  }
  if (!converged)
    .fit_error(paste0("SS: did not converge in ", max_iter, " iterations"))
  e <- sweep(sweep(Z, 2L, theta, "-"), 2L, d, "/")
  r <- sqrt(rowSums(e * e))
  if (min(r) <= radius_tol * stats::median(r))
    .fit_error("SS: final residual radius collapsed")
  U <- e / r
  calibration <- radial_calibration(r, "SS", kappa)
  A <- calibration$A
  B <- calibration$B
  z <- sqrt(T * N) * A * theta / (d * sqrt(B))
  .new_fit("SS", theta, d, z, q, omega, A, B, iteration,
           max(abs(colMeans(U))), max(abs(N * colMeans(U * U) - 1)),
           min(r) / stats::median(r))
}

fit_in_projected <- function(Z, projection, d, q = 0.1, initial_alpha = NULL,
                             tol = 1e-8, max_iter = 1000L) {
  .fit_controls(q, tol, max_iter)
  Z <- .finite_matrix(Z, "Z")
  T <- nrow(Z)
  N <- ncol(Z)
  if (length(d) != N || any(!is.finite(d)) || any(d <= 0))
    stop("d must be a finite positive scale vector of length N", call. = FALSE)
  h <- projection$h
  if (length(h) != T || any(!is.finite(h)))
    stop("Projection must have one finite h value per date", call. = FALSE)
  h2 <- sum(h * h)
  if (h2 <= 0) .fit_error("IN: intercept is unidentified")
  alpha <- if (is.null(initial_alpha)) as.vector(crossprod(h, Z)) / h2 else initial_alpha
  if (length(alpha) != N || any(!is.finite(alpha)))
    stop("initial_alpha must be a finite vector of length N", call. = FALSE)
  radius_tol <- 100 * .Machine$double.eps
  # This trace remains in memory and is attached only to a failed fit.
  trace <- matrix(NA_real_, max_iter, 4L,
                  dimnames = list(NULL, c("iteration", "step", "min_radius_ratio",
                                          "normalized_raw_score")))
  trace_length <- 0L
  fail <- function(message, iteration, radii = NULL, raw_score = NULL) {
    detail <- list(iteration = iteration,
                   trace = trace[seq_len(trace_length), , drop = FALSE],
                   h_min = min(h), h_max = max(h), omega = projection$omega)
    if (!is.null(radii)) detail$min_radius_ratio <- min(radii) / stats::median(radii)
    if (!is.null(raw_score)) detail$normalized_raw_score <- raw_score
    .fit_error(message, detail)
  }
  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    e <- sweep(Z - tcrossprod(h, alpha), 2L, d, "/")
    r <- sqrt(rowSums(e * e))
    median_r <- stats::median(r)
    if (any(!is.finite(r)) || min(r) <= radius_tol * median_r)
      fail("IN: collapsed or nonfinite residual radius", iteration, r)
    inverse_r2 <- 1 / (r * r)
    current_score <- colMeans(e * (h * inverse_r2))
    raw_equation <- sqrt(N) * median_r * max(abs(current_score))
    denominator <- sum(h * h * inverse_r2)
    if (!is.finite(denominator) || denominator <= 0)
      fail("IN: nonpositive or nonfinite h-weighted curvature", iteration, r, raw_equation)
    alpha_new <- colSums(Z * (h * inverse_r2 / denominator))
    loc_change <- max(abs((alpha_new - alpha) / d)) / (median_r / sqrt(N))
    trace[iteration, ] <- c(iteration, loc_change, min(r) / median_r, raw_equation)
    trace_length <- iteration
    alpha <- alpha_new
    if (loc_change < tol) {
      check_e <- sweep(Z - tcrossprod(h, alpha), 2L, d, "/")
      check_r <- sqrt(rowSums(check_e * check_e))
      if (any(!is.finite(check_r)) ||
          min(check_r) <= radius_tol * stats::median(check_r))
        fail("IN: residual radius collapsed at stopping candidate", iteration, check_r)
      check_inverse <- 1 / (check_r * check_r)
      check_score <- colMeans(check_e * (h * check_inverse))
      equation <- sqrt(N) * stats::median(check_r) * max(abs(check_score))
      # Do not divide this raw-score criterion by the exploding A_h.
      if (equation < tol) {
        converged <- TRUE
        break
      }
    }
  }
  if (!converged)
    fail(paste0("IN: did not converge in ", max_iter, " iterations"),
         max_iter, r, raw_equation)
  e <- sweep(Z - tcrossprod(h, alpha), 2L, d, "/")
  r <- sqrt(rowSums(e * e))
  inverse_r2 <- 1 / (r * r)
  A_h <- mean(h * h * inverse_r2)
  score <- colMeans(e * (h * inverse_r2))
  score_residual <- sqrt(N) * stats::median(r) * max(abs(score))
  pi <- if (!ncol(projection$Q)) rep(0, T) else
    as.vector(projection$Q %*% crossprod(projection$Q, h * inverse_r2))
  B_c <- A_h + mean((r * pi)^2)
  z_uncalibrated <- sqrt(T * N) * A_h * alpha / (d * sqrt(B_c))
  residual_df <- projection$residual_df
  if (residual_df <= 0) stop("IN requires positive residual degrees of freedom")
  df_multiplier <- sqrt(residual_df / T)
  z <- df_multiplier * z_uncalibrated
  .new_fit("IN", alpha, d, z, q, projection$omega, A_h, B_c, iteration,
           score_residual, NaN, min(r) / stats::median(r),
           df_multiplier, residual_df, z_uncalibrated)
}

fit_three_projected <- function(Z, projection, q = 0.1, tol = 1e-8,
                                max_iter = 1000L, capture_errors = FALSE,
                                pqr_scale = NULL, smoothing = 0.5, in_pilot_radius = NULL) {
  smoothing <- .rfs_smoothing(smoothing)
  .fit_controls(q, tol, max_iter)
  if (is.null(pqr_scale))
    pqr_scale <- tryCatch(pqr_scale_projected(Z, projection), error = identity)
  if (inherits(pqr_scale, "error") && !capture_errors) stop(pqr_scale)
  out <- stats::setNames(vector("list", 3L), c("LS", "SS", "IN"))
  for (method in names(out)) {
    result <- tryCatch({
      switch(method,
        LS = fit_ls_projected(Z, projection, q),
        SS = fit_projected(Z, projection$omega, "SS", q = q, tol = tol,
                           max_iter = max_iter),
        IN = {
          if (inherits(pqr_scale, "error"))
            .fit_error(paste0("IN: PQR estimation failed: ", conditionMessage(pqr_scale)))
          fit_regularized_in_projected(Z, projection, pqr_scale$d, q, tol = tol,
                           max_iter = max_iter, smoothing = smoothing,
                           pilot_radius = in_pilot_radius)
        })
    }, error = identity)
    if (inherits(result, "error") && !capture_errors) stop(result)
    out[[method]] <- result
  }
  out
}

fit_three <- function(Y, factors = NULL, q = 0.1, tol = 1e-8,
                      max_iter = 1000L, capture_errors = FALSE, smoothing = 0.5) {
  Y <- .finite_matrix(Y, "Y")
  projection <- observed_projection(nrow(Y), factors)
  fit_three_projected(project_observations(Y, projection), projection,
                      q, tol, max_iter, capture_errors, smoothing = smoothing)
}

select_funds <- function(Y, factors = NULL, method = c("IN", "SS", "LS"),
                         q = 0.1, tol = 1e-8, max_iter = 1000L, smoothing = 0.5) {
  smoothing <- .rfs_smoothing(smoothing)
  method <- match.arg(method)
  .fit_controls(q, tol, max_iter)
  Y <- .finite_matrix(Y, "Y")
  projection <- observed_projection(nrow(Y), factors)
  Z <- project_observations(Y, projection)
  fit <- tryCatch(switch(method,
    LS = fit_ls_projected(Z, projection, q),
    SS = fit_projected(Z, projection$omega, "SS", q = q, tol = tol, max_iter = max_iter),
    IN = fit_regularized_in_projected(Z, projection, pqr_scale_projected(Z, projection)$d,
                          q = q, tol = tol, max_iter = max_iter, smoothing = smoothing)),
    radial_fit_error = identity)
  if (inherits(fit, "radial_fit_error"))
    return(structure(list(method = method, pvalue = rep(1, ncol(Y)),
      selected = rep(FALSE, ncol(Y)), fit_converged = FALSE, fit = NULL,
      diagnostics = c(list(error = conditionMessage(fit), failure_policy = "no rejections",
                           smoothing = if (method == "IN") smoothing else NA_real_),
                      fit$diagnostics)), class = "radial_selection"))
  structure(list(method = method, pvalue = fit$pvalue, selected = fit$selected,
    fit_converged = TRUE, fit = fit,
    diagnostics = fit[c("iterations", "location_residual", "min_radius_ratio",
                        "scale_estimator", "df_multiplier", "residual_df",
                        "implementation", "smoothing", "lambda", "pilot_radius",
                        "max_date_weight")]),
    class = "radial_selection")
}

discovery_metrics <- function(selected, alpha) {
  if (length(selected) != length(alpha)) stop("selected and alpha must have equal length")
  positive <- alpha > 0
  R <- sum(selected)
  V <- sum(selected & !positive)
  S <- sum(selected & positive)
  n_positive <- sum(positive)
  c(R = R, V = V, S = S, n_positive = n_positive,
    FDP = V / max(R, 1), TDP = if (n_positive) S / n_positive else NaN)
}
