\name{reproduce_simulations}
\alias{reproduce_simulations}
\alias{plot_simulation_figures}
\title{Reproduce the Paper's Simulation Results and Figures}
\usage{
reproduce_simulations(reps = 2000L, workers = 1L, output_dir = NULL,
  progress = TRUE)

plot_simulation_figures(study, output_dir)
}
\arguments{
\item{reps}{Repetitions in each of the 78 distinct settings. The paper uses 2000.}
\item{workers}{Number of workers between 1 and 25. Worker count affects runtime.}
\item{output_dir}{A directory for three PDF figures. With \code{NULL},
\code{reproduce_simulations} returns the numerical results in memory.}
\item{progress}{Print Monte Carlo progress.}
\item{study}{Result of \code{simulate_study()} using all 78 scenarios.}
}
\description{Runs the paper's three-method comparison and paired smoothing
comparison using the original common random streams, seed 20261105, and the
paper's model and tuning parameters. The curves in the figures are drawn
from the setting-level estimates returned
by the simulation, with pointwise 95 percent Monte Carlo intervals.}
\details{In the signal experiment the panels compare three error distributions
at nine signal strengths for each sample size. The correlation experiment
uses nine residual correlations. Rates are displayed as percentages.
The \code{simulation_facts} and \code{smoothing_facts} entries report proportions.
Multiply by 100 for percentage points, including contributions to the mean
change in TDP. The full run requires substantial computation.}
\value{\code{reproduce_simulations} returns \code{study}, \code{smoothing},
\code{simulation_facts}, \code{smoothing_facts}, and the paths to any requested
\code{figures}. \code{plot_simulation_figures} invisibly returns the three PDF paths.}
\examples{
\dontrun{
paper <- reproduce_simulations(workers = 4L, output_dir = "figures")
paper$smoothing_facts$overall_paired
}
}
\seealso{\code{\link{simulate_study}}, \code{\link{simulate_smooth_study}}}
