.pair_indices <- local({
  cache <- list()
  function(T) {
    key <- as.character(T)
    if (is.null(cache[[key]])) {
      indices <- which(upper.tri(matrix(FALSE, T, T)), arr.ind = TRUE)
      if (length(cache) >= 8L) cache <<- cache[-1L]
      cache[[key]] <<- list(left = indices[, 1L], right = indices[, 2L])
    }
    cache[[key]]
  }
})

pqr_scale_from_residuals <- function(residual, quantile = 0.5) {
  residual <- .finite_matrix(residual, "residual")
  if (length(quantile) != 1L || !is.finite(quantile) ||
      quantile <= 0 || quantile >= 1)
    stop("PQR quantile must lie in (0,1)", call. = FALSE)
  indices <- .pair_indices(nrow(residual))
  M <- length(indices$left)
  order <- as.integer(ceiling(M * quantile))
  quantiles <- vapply(seq_len(ncol(residual)), function(i) {
    x <- residual[, i]
    distances <- abs(x[indices$right] - x[indices$left])
    sort.int(distances, partial = order)[order]
  }, numeric(1))
  if (any(quantiles <= 0) || any(!is.finite(quantiles)))
    .fit_error("PQR: a coordinate pairwise quantile is zero or nonfinite")
  normalizer <- sqrt(mean(quantiles * quantiles))
  if (!is.finite(normalizer) || normalizer <= 0)
    .fit_error("PQR: trace normalizer is degenerate")
  list(d = quantiles / normalizer, quantiles = quantiles,
       quantile = quantile, pair_count = M, order_statistic = order)
}

pqr_scale_projected <- function(Z, projection, quantile = 0.5) {
  Z <- .finite_matrix(Z, "Z")
  h <- projection$h
  if (nrow(Z) != length(h))
    stop("Invalid projected observation matrix", call. = FALSE)
  alpha <- as.vector(crossprod(h, Z)) / sum(h * h)
  unrestricted <- Z - tcrossprod(h, alpha)
  pqr_scale_from_residuals(unrestricted, quantile)
}

pqr_scale <- function(Y, factors = NULL, quantile = 0.5) {
  Y <- .finite_matrix(Y, "Y")
  projection <- observed_projection(nrow(Y), factors)
  pqr_scale_projected(project_observations(Y, projection), projection, quantile)
}
