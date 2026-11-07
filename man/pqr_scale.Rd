\name{pqr_scale}
\alias{pqr_scale}
\title{Pairwise-quantile scales from unrestricted regression residuals}
\description{Removes both the intercept and observed factors, computes a
marginal quantile of absolute pairwise residual differences, and normalizes
the resulting scales so their squared mean is one.}
\usage{pqr_scale(Y, factors = NULL, quantile = 0.5)}
\arguments{
\item{Y}{Finite numeric matrix of dates by funds; both dimensions at least three.}
\item{factors}{Observed factors without an intercept column, or \code{NULL}.}
\item{quantile}{Quantile probability in the open interval \eqn{(0,1)}.}
}
\details{With \eqn{M=T(T-1)/2} pairs, the estimator uses the order statistic
of one-based rank \eqn{\lceil M\,\mathrm{quantile}\rceil}. Thus a probability
of one half uses the lower middle observation when \eqn{M} is even, rather
than an interpolated median. This matches the manuscript's empirical
inverse-CDF convention. The regression removes any column-specific alpha
and the observed factor exposures before differences are formed.}
\value{A list with positive scales \code{d}, unnormalized \code{quantiles},
the probability \code{quantile}, \code{pair_count}, and \code{order_statistic}.
A zero or nonfinite scale raises a numerical fit error.}
\examples{
set.seed(3)
scale <- pqr_scale(matrix(rnorm(30 * 10), 30, 10))
mean(scale$d^2)
}
