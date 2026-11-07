\name{simulate_smooth_study}
\alias{simulate_smooth_study}
\title{Paired Comparison of Smoothed and Unregularized IN-BH}
\usage{
simulate_smooth_study(reps = 2000L, workers = 1L, seed = 20261105L,
  c_grid = 0.5, q = 0.1, tol = 1e-8, max_iter = 1000L,
  block_size = 5L, scenarios = scenario_inventory(), progress = TRUE)
}
\arguments{
\item{reps}{Number of independent repetitions for each setting.}
\item{workers}{Number of workers, from 1 to 25.}
\item{seed}{Integer master seed for the Monte Carlo streams.}
\item{c_grid}{Distinct positive relative smoothing strengths. Each is compared with zero.}
\item{q}{Nominal level for the BH rule.}
\item{tol}{Tolerance for numerical convergence.}
\item{max_iter}{Maximum number of updates.}
\item{block_size}{Repetitions per scheduling block.}
\item{scenarios}{Data frame containing one or more case identifiers from \code{scenario_inventory()}.}
\item{progress}{Logical. Print progress at batch completion.}
}
\description{Fits unregularized and smoothed IN-BH on the same panels and estimates
the paired change in failure rate, FDR, and TDP. For group-level intervals the
code first averages the settings within each repetition, preserving dependence
from the common draws. All repetitions enter the denominator.}
\value{A list containing per-setting \code{summary} and \code{paired} data frames,
\code{group_summary} and \code{group_paired} estimates, \code{failure_reasons},
the scenario design and configuration. Rates and differences are in proportion
units. The paired frame reports conditional results for rescued fits and fits
on which both procedures succeeded.}
\examples{
cases <- scenario_inventory()[1L, , drop = FALSE]
fit <- simulate_smooth_study(reps = 2L, scenarios = cases, progress = FALSE)
fit$group_paired[fit$group_paired$group == "Overall", ]
}
\seealso{\code{\link{simulate_study}}, \code{\link{reproduce_simulations}}}
