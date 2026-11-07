\name{simulate_study}
\alias{simulate_study}
\title{Run the Paired Study With In-Memory Aggregation}
\usage{simulate_study(reps = 2000L, workers = 25L, seed = 20261105L,
  progress = TRUE, q = 0.10, tol = 1e-8, max_iter = 1000L,
  block_size = 5L, scenarios = scenario_inventory(), smoothing = 0.5)}
\arguments{
\item{reps}{Repetitions per distinct scenario; 2000 reproduces the study size.}
\item{workers}{Maximum worker count from 1 to 25. The default uses 25 PSOCK workers
when there are at least 25 work blocks. One requests serial execution.}
\item{seed}{Nonnegative integer master seed for R's independent random streams.}
\item{progress}{Whether to print progress after each batch of blocks.}
\item{q}{Common nominal BH level.}
\item{tol}{Location and scale convergence tolerance passed to the estimators.}
\item{max_iter}{Maximum iterations passed to the estimators.}
\item{smoothing}{IN relative smoothing strength, default 0.5. Zero selects
legacy unregularized IN. See \code{\link{select_funds}}.}
\item{block_size}{Number of repetitions aggregated per task. Changes scheduling
and the rounding order of sums, but not the generated panels.}
\item{scenarios}{A data frame identifying a nonempty subset of
\code{scenario_inventory()} by \code{case_id}. The default runs all 78 cases;
all design columns are resolved from the canonical inventory.}}
\description{Evaluates LS-BH, SS-BH, and IN-BH on paired panels, aggregating sums,
squared sums, paired differences, and failure counts in memory. IN is smoothed
with strength 0.5 by default; zero smoothing is an explicit legacy control. It returns only
final consolidated summaries and writes no files.}
\details{All three methods receive the same panel in each repetition. Error
sources are shared across rho and delta; nuisance draws are also shared across
the three distributions. The rho=0.5, T=60, delta=0.48 overlap is simulated once
and marked for both displays. The unrestricted PQR scale is computed once per
rho and base draw and held fixed across alpha shifts. The OLS/PQR pilot radius
is also computed from that base draw. Its invariance to the projected alpha
shift makes this the same pilot as from the actual panel. Location fitting
always uses the actual shifted panel and its own OLS intercept initial value.
The IN variance uses the general radial expression for positive smoothing;
it does not use the special unsmoothed inverse-norm identity.

The full default run has 78 times 2000 = 156000 panels and 468000 method fits.
It is never run automatically by loading the package. PSOCK workers are stopped
on normal completion or error. At most 25 workers are started, and numerical
library thread limits are set to one before workers start. The caller's RNG and
environment variables are restored. Installed-package and sourced-code execution
are supported.

Each failed fit makes no rejections and contributes zero FDP and TDP. Failures
remain in all attempted-repetition and paired-difference denominators; they are
not replaced by fresh draws. FDR and TDP standard errors are sample standard
deviations divided by the square root of the number of repetitions, calculated
from sums and squared sums. Paired standard errors use within-panel differences,
not independent-method approximations. Standard errors are missing for one
repetition. Confidence intervals are pointwise normal Monte Carlo intervals.

Random streams depend on sample size, source, and one-based repetition index,
not worker count, load balancing, scenario subset, or block size. Aggregation
with different block sizes can differ by floating-point roundoff. R's
L'Ecuyer-CMRG/Inversion/Rejection generator reproduces the Python design in
distribution; it does not reproduce NumPy SeedSequence/PCG64 values or historical
Monte Carlo results bit for bit.

Only per-scenario additive statistics are retained after each block. No
per-repetition records, block files, RDS files, CSV dumps, checkpoints, or seed
files are created. Save or format the returned final summaries explicitly if
desired.}
\value{A \code{radial_study} list with:
\item{summary}{One row per scenario and method, including FDR, TDP, mean discovery
counts, Monte Carlo standard errors, successful and failed fit counts, and
failure-sensitivity bounds. Operational summaries include failed fits.}
\item{paired_differences}{Within-panel SS-LS, IN-LS, and IN-SS differences for FDR,
TDP, and discovery counts, with their Monte Carlo standard errors and intervals.}
\item{scenarios}{The distinct scenarios actually requested.}
\item{config}{Design, RNG, fitting, parallelism, runtime, and provenance metadata.}}
\examples{cases <- scenario_inventory()
\dontrun{
# Full requested study. This may require substantial computing time.
result <- simulate_study(reps = 2000, workers = 25)
# Small numerical check, without saving any intermediate output.
check <- simulate_study(reps = 2, workers = 2, scenarios = cases[c(1, 27), ],
                        progress = FALSE)
}}
\seealso{\code{\link{simulate_panel}}, \code{\link{scenario_inventory}}}
