\name{reproduce_paper}
\alias{reproduce_paper}
\title{Recompute the simulation and empirical results in the paper}
\usage{
reproduce_paper(data, rf_file = NULL, output_dir = "paper-results",
                reps = 2000L, workers = 1L, progress = TRUE)
}
\description{Runs the Section 5 Monte Carlo experiments and the Section 6
rolling mutual fund analysis. Writes Figures 1 through 4 as vector PDFs.}
\arguments{
\item{data}{The authors' original mutual fund MAT file or the directory with
the paper's aligned CSV inputs, as specified in \code{\link{reproduce_empirical}}.}
\item{rf_file}{Optional path to the French risk-free CSV or ZIP. Required for
MAT input. With a data directory the default is its \file{rf.csv}.}
\item{output_dir}{Destination for the four named PDFs, created if necessary.
Existing figures with the same names are replaced.}
\item{reps}{Monte Carlo repetitions per setting. The paper uses 2000.}
\item{workers}{Number of parallel R processes, from 1 through 25.}
\item{progress}{Display simulation progress.}
}
\value{A list with \code{simulations} from \code{\link{reproduce_simulations}},
\code{empirical} from \code{\link{reproduce_empirical}} and the four PDF paths
in \code{figures}. Numerical summaries are returned in memory.}
\examples{
\dontrun{
run <- reproduce_paper("CaseStudy/mutual_fund.mat", "F-F_Research_Data_Factors_CSV.zip",
                       workers = 4L)
run$figures
run$empirical$summary
}
}
