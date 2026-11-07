\name{select_funds}
\alias{select_funds}
\alias{fit_three}
\title{Select positive-alpha funds using LS-BH, SS-BH, or IN-BH}
\description{Fits the original least-squares and spatial-sign benchmarks or
IN-BH with a fixed pairwise-quantile scale, projected
intercept weights, and residual degrees-of-freedom calibration. Smoothing is
an implementation parameter; the method name remains IN-BH.}
\usage{
select_funds(Y, factors = NULL, method = c("IN", "SS", "LS"),
             q = 0.1, tol = 1e-8, max_iter = 1000L, smoothing = 0.5)
fit_three(Y, factors = NULL, q = 0.1, tol = 1e-8,
          max_iter = 1000L, capture_errors = FALSE, smoothing = 0.5)
}
\arguments{
\item{Y}{Finite numeric matrix: rows are dates and columns are funds. Both
dimensions must be at least three.}
\item{factors}{Observed factor matrix, with one row per date and no intercept
column. The columns must have full rank and leave an identifiable intercept.
\code{NULL} uses an intercept-only model.}
\item{method}{\code{"LS"}, \code{"SS"}, or \code{"IN"}.}
\item{q}{Nominal Benjamini--Hochberg level in the open interval \eqn{(0,1)}.}
\item{tol}{Positive convergence tolerance, applied to both updates and
estimating equations.}
\item{max_iter}{Positive integer iteration limit.}
\item{smoothing}{Nonnegative relative smoothing strength for IN only.
The default is 0.5. Zero explicitly requests the legacy unregularized IN
implementation from package version 0.1.0. LS and SS are unaffected.}
\item{capture_errors}{Whether \code{fit_three} retains each failed fit as an
error object instead of stopping.}
}
\details{LS uses unrestricted OLS residual variance with divisor \eqn{T}.
SS jointly updates location and geometrically normalized diagonal scale,
using the same old residuals for both updates. IN estimates a trace-normalized
PQR scale from unrestricted residuals and holds it fixed. For positive
\code{smoothing}, let \eqn{r_{t,0}} be the norm of the PQR-standardized
unrestricted OLS residual, and set
\eqn{\lambda=\mathtt{smoothing}\times\mathrm{median}_t(r_{t,0})}.
This pilot does not use the true alpha and is invariant to an intercept or
observed-factor shift. Lambda is fixed throughout the location iterations.
The estimating equation is \eqn{\sum_t h_t e_t/(r_t^2+\lambda^2)=0},
where \eqn{e_t=D^{-1/2}(Z_t-h_t\alpha)}. It retains the
factor-projection term in its studentizer and multiplies the statistic by
\eqn{\sqrt{(T-\mathrm{rank}([1,F]))/T}}. All methods use upper normal tails.
The positive-lambda variance uses the general radial expression
\eqn{A=\mathrm{mean}_t\{h_t^2/(r_t^2+\lambda^2)\}} and
\eqn{B=\mathrm{mean}_t\{r_t^2[h_t/(r_t^2+\lambda^2)-\pi_t]^2\}},
with \eqn{\pi=P_F\{h/(r^2+\lambda^2)\}}. The special inverse-norm
variance identity is used only when smoothing is zero.
No radii are clipped and no alternate roots are substituted. A zero residual
radius is valid for positive lambda; it is a failure under the legacy algorithm.
The IN raw score is checked without division by its sensitivity, to avoid
false convergence at a collapsing residual radius.

The positive fixed smoothing strength defines a regularized implementation,
not the exact inverse-norm weight that appears in an optimality theorem.
The pilot-dependent rule is empirically evaluated; the numerical results do
not establish a new FDR theorem for this rule. In particular, the finite-sample
FDR need not be below the nominal level in every simulation setting.}
\value{\code{select_funds} returns a list with \code{method}, \code{pvalue},
\code{selected} (logical vector), \code{fit_converged}, \code{fit}, and
\code{diagnostics}. A numerical fit failure yields p-values all equal to one,
no selections, a false convergence flag, and the failure diagnostics. Invalid
arguments raise errors.
\code{fit_three} returns a list named \code{LS}, \code{SS}, and \code{IN}.
A successful fit contains \code{theta}, \code{alpha_first_order}, \code{d},
\code{z}, \code{pvalue}, \code{selected}, calibration quantities, and
convergence diagnostics. IN additionally records \code{implementation},
\code{smoothing}, \code{lambda}, and \code{pilot_radius}; positive smoothing
also records \code{max_date_weight}. Under \code{capture_errors = TRUE}, each failed
fit is an error condition and must count as zero rejections in simulation
summaries. Neither function writes data to disk.}
\examples{
set.seed(7)
Y <- matrix(rnorm(60 * 20), 60, 20)
ans <- select_funds(Y, method = "LS")
sum(ans$selected)
fits <- fit_three(Y, capture_errors = TRUE)
vapply(fits, inherits, logical(1), what = "error")
}
