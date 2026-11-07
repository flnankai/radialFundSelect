\name{simulate_panel}
\alias{simulate_panel}
\title{Generate One Observed-Factor Return Panel}
\usage{simulate_panel(T = 60L, N = 200L, distribution = "normal",
  rho = 0.5, delta = 0.06, seed = 20261105L, replicate = 1L)}
\arguments{
\item{T}{Number of dates, either 60 or 120.}
\item{N}{Number of funds, at least 10; the manuscript uses 200.}
\item{distribution}{One of \code{"normal"}, \code{"t3"}, or \code{"mixture"}.}
\item{rho}{Cross-sectional AR(1) correlation, between zero and one, excluding one.}
\item{delta}{Nonnegative increment added to each positive alpha.}
\item{seed}{Nonnegative integer master seed.}
\item{replicate}{One-based repetition index.}}
\description{Generates a panel from the manuscript design using independent R random
streams. No files are written, and the caller's random-number state is restored.}
\details{There are three observed Gaussian factors with means drawn independently
from uniform(0,1) and covariance entries \eqn{0.5^{|j-k|}}. The loading ranges are
(0.2,2), (-1,1.5), and (-1.5,1.5). A random ten percent of funds have alpha uniform
(0.15,0.30) plus \code{delta}, another ten percent have alpha uniform(-0.30,-0.20),
and the rest have zero alpha. Nuisance parameters are redrawn every repetition.

Gaussian errors have cross-sectional covariance \eqn{\rho^{|i-j|}}. Student errors
divide the correlated Gaussian vector by the square root of one chi-squared(3)
draw per date. Mixture errors multiply that vector by 1 or 3 with probabilities
0.8 and 0.2, then divide by the square root of 2.6. Both scales are shared across
funds at a date; all three laws have the same covariance matrix.

L'Ecuyer-CMRG streams separate sample sizes and nuisance/error sources; repetitions
use independent substreams. For fixed sample size and repetition, nuisance draws
are shared across distributions, and errors are shared across rho and delta.
The RNG reproduces the distributions and pairing design, not the numerical draws
of Python/NumPy. The seed is not numerically interchangeable with the Python seed.}
\value{A list containing \code{Y}, \code{factors}, \code{alpha}, \code{beta},
\code{errors}, positive and negative index sets, factor means, and draw metadata.}
\examples{panel <- simulate_panel()
dim(panel$Y)
table(sign(panel$alpha))}
\seealso{\code{\link{simulate_study}}, \code{\link{scenario_inventory}}}
