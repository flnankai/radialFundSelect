\name{select_in_bh}
\alias{select_in_bh}
\alias{select_ss_bh}
\alias{select_ls_bh}
\title{Select funds with positive alpha using a specified method}
\description{Apply the inverse-norm (IN-BH), spatial-sign (SS-BH), or least-squares
(LS-BH) procedure to monthly fund excess returns. Each method uses upper-tail
p-values and the Benjamini-Hochberg rule.}
\usage{
select_in_bh(Y, factors = NULL, q = 0.1, smoothing = 0.5,
             tol = 1e-8, max_iter = 1000L)
select_ss_bh(Y, factors = NULL, q = 0.1, tol = 1e-8, max_iter = 1000L)
select_ls_bh(Y, factors = NULL, q = 0.1)
}
\arguments{
\item{Y}{Finite numeric matrix of excess returns. Rows are dates, columns are funds.}
\item{factors}{Numeric factor matrix aligned to the dates of \code{Y}, without
an intercept. \code{NULL} fits an intercept-only model.}
\item{q}{Nominal false discovery rate between zero and one. The empirical
application uses 0.05.}
\item{smoothing}{Nonnegative smoothing strength. A value of 0.5 gives the
implemented IN-BH rule in the paper. Zero gives the unregularized control.}
\item{tol}{Positive convergence tolerance.}
\item{max_iter}{Maximum number of location or scale updates.}}
\value{A \code{radial_selection} list. \code{selected} flags funds in column
order, \code{pvalue} gives the one-sided p-values, and \code{fit_converged}
indicates success. \code{fit} contains the fitted quantities, including
\code{alpha_first_order} and \code{z}. See \code{\link{select_funds}} for
failure handling, calibration and diagnostics.}
\examples{
panel <- simulate_panel(N = 30L)
in_fit <- select_in_bh(panel$Y, panel$factors)
ss_fit <- select_ss_bh(panel$Y, panel$factors)
ls_fit <- select_ls_bh(panel$Y, panel$factors)
which(in_fit$selected)
}
