\name{scenario_inventory}
\alias{scenario_inventory}
\title{The 78 Distinct Monte Carlo Scenarios}
\usage{scenario_inventory()}
\description{Constructs the complete design of the manuscript without running a simulation.}
\value{A data frame with 78 rows and columns \code{case_id}, \code{distribution},
\code{N}, \code{T}, \code{rho}, \code{delta}, \code{signal_grid}, and
\code{rho_sensitivity}. The shared scenario belongs to both displays but occurs once.}
\details{The three errors are Gaussian, covariance-standardized Student t with three
degrees of freedom, and the symmetric normal scale mixture. Signal grids use
\eqn{T=60,120}, \eqn{\rho=0.5}, and nine increments at each sample size. The
correlation grid uses \eqn{T=60}, \eqn{\delta=0.48}, and \eqn{\rho=0.1,\ldots,0.9}.
Every scenario has 200 funds, including 20 positive and 20 negative alphas.}
\examples{cases <- scenario_inventory()
stopifnot(nrow(cases) == 78L, sum(cases$signal_grid) == 54L,
          sum(cases$rho_sensitivity) == 27L)}
