# Observed factors exclude the intercept. These helpers retain the original
# factor-only projection used by all three published benchmark procedures.
.finite_matrix <- function(x, name, minimum = 3L) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L ||
      min(dim(x)) < minimum || any(!is.finite(x)))
    stop(name, " must be a finite numeric T by N matrix, with T,N >= ",
         minimum, call. = FALSE)
  storage.mode(x) <- "double"
  x
}

.fit_error <- function(message, diagnostics = list()) {
  stop(structure(list(message = message, call = NULL, diagnostics = diagnostics),
                 class = c("radial_fit_error", "error", "condition")))
}

.fit_controls <- function(q, tol, max_iter) {
  if (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1 ||
      length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter))
    stop("Require 0 < q < 1, tol > 0, and an integer max_iter >= 1",
         call. = FALSE)
}

observed_projection <- function(T, factors = NULL) {
  if (length(T) != 1L || !is.finite(T) || T < 3 || T != as.integer(T))
    stop("Require integer T >= 3", call. = FALSE)
  if (is.null(factors) || !length(factors))
    return(list(Q = matrix(numeric(), T, 0L), h = rep(1, T), omega = 1,
                rank_X = 1L, residual_df = T - 1L))
  F <- as.matrix(factors)
  if (!is.numeric(F) || nrow(F) != T || any(!is.finite(F)))
    stop("factors must be finite with one row per date", call. = FALSE)
  p <- ncol(F)
  if (p >= T - 1L)
    stop("factors must have full column rank and p < T-1", call. = FALSE)
  singular_values <- svd(F, nu = 0L, nv = 0L)$d
  rank_F <- sum(singular_values > max(dim(F)) * .Machine$double.eps *
                 max(singular_values))
  if (rank_F != p)
    stop("factors must have full column rank and p < T-1", call. = FALSE)
  # LAPACK QR is used for the projection; pivoting does not change its span.
  Q <- qr.Q(qr(F, LAPACK = TRUE), complete = FALSE)
  h <- 1 - as.vector(Q %*% colSums(Q))
  omega <- sum(h * h) / T
  if (!is.finite(omega) || omega <= 1e-10 || omega > 1 + 1e-12)
    stop("Intercept is unidentified: do not include an intercept in factors",
         call. = FALSE)
  list(Q = Q, h = h, omega = min(omega, 1), rank_X = p + 1L,
       residual_df = T - p - 1L)
}

project_observations <- function(Y, projection) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || nrow(Y) != length(projection$h) ||
      any(!is.finite(Y)))
    stop("Y must be finite with one row per date", call. = FALSE)
  if (!ncol(projection$Q)) return(Y)
  Y - projection$Q %*% crossprod(projection$Q, Y)
}

observed_factor_residuals <- function(Y, factors = NULL) {
  Y <- .finite_matrix(Y, "Y")
  projection <- observed_projection(nrow(Y), factors)
  list(Z = project_observations(Y, projection), omega = projection$omega,
       kappa = 1 - projection$omega)
}

bh <- function(pvalue, q = 0.1) {
  if (!is.numeric(pvalue) || !is.null(dim(pvalue)) ||
      any(!is.finite(pvalue)))
    stop("pvalue must be a finite numeric vector", call. = FALSE)
  if (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1 ||
      any(pvalue < 0 | pvalue > 1))
    stop("Require 0 < q < 1 and p-values in [0,1]", call. = FALSE)
  N <- length(pvalue)
  if (!N) return(logical())
  ordered <- sort.int(pvalue)
  passing <- which(ordered <= q * seq_len(N) / N)
  if (!length(passing)) return(rep(FALSE, N))
  pvalue <= ordered[passing[length(passing)]]
}
