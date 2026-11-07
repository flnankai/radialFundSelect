\name{reproduce_empirical}
\alias{reproduce_empirical}
\alias{plot_empirical_wealth}
\title{Reproduce the U.S. equity fund application and portfolio wealth figure}
\description{
Runs the paper's rolling four-factor analysis on local copies of the
Gao--Peng--Yan mutual fund panel and the monthly risk-free series. Calculates
the full-period and common-endpoint results and draws Figure 4 to a PDF.
}
\usage{
reproduce_empirical(data, rf_file = NULL, figure = NULL)

plot_empirical_wealth(result, file)
}
\arguments{
\item{data}{Path to the authors' original \file{CaseStudy/mutual_fund.mat}
(requires the optional package \pkg{R.matlab}), or a directory containing
the aligned, headerless \file{fund.csv} (726 funds by 368 dates) and
\file{factor.csv} (368 dates by four factors). Returns and factors are
decimals. See Details for the months, factor order and missing-value code.}
\item{rf_file}{Path to the monthly U.S. Fama/French 3 Factors CSV or its
original ZIP archive, or to a headerless file with 368 aligned risk-free
returns in percent (one per line). For the directory input the default is
\file{rf.csv} within that directory. Required for the MAT input.}
\item{figure}{Optional output path for the three-panel vector PDF.}
\item{result}{Object returned by \code{reproduce_empirical}.}
\item{file}{Output PDF path in an existing directory.}
}
\details{
The original MAT archive contains monthly matrices \code{R} (668 by 726)
and \code{F} (668 by 4). Rows 301--668 of both matrices span February 1987
through September 2017. The four archived factors are market excess return,
SMB, HML and momentum. The French file supplies only the risk-free series.

In each 69-month estimation window a zero in a raw fund series codes a
missing observation and excludes that fund from the candidate set. Subtract
the monthly risk-free rate after screening. \code{\link{fit_three}} fits
LS-BH, SS-BH and IN-BH at \eqn{q=0.05}, smoothing 0.5, tolerance
\eqn{10^{-8}} and at most 1,000 updates. Fitting errors stop the calculation.
The next complete 3, 6 or 9 months form the holding period. The next
training window advances by the same number of months.

Equal dollar investments compound raw fund returns with weights drifting
within the holding period. An empty selection earns the risk-free rate.
On its first missing return, a position earns the risk-free rate from that
month through the end of the period. Under stress, the affected position
loses all value at its first missing return. The baseline and stress
calculations use the same selections. No trading costs are charged.
}
\value{
\code{reproduce_empirical} returns four data frames in memory.
\code{summary} reports dates, numbers of months and formations, eligible
fund ranges, empty selection counts, affected fund positions, terminal
wealth, compound annual growth rate (\code{cagr}), peak-to-trough loss
(\code{max_drawdown}, a positive fraction), \code{stress_cagr}, and
the common-endpoint growth rate and drawdown through April 2017. Rates
and drawdowns are decimals. \code{periods} records historical and holding
dates, candidate counts, selection counts and affected positions for each
formation and method. \code{wealth} and \code{stress_wealth} report
monthly wealth for each design, including the initial unit investment at
end-October 1992. \code{plot_empirical_wealth} returns its file path invisibly.
}
\references{
Gao, J., B. Peng, and Y. Yan (2024). Higher-order expansions and inference for
panel data models. \emph{Journal of the American Statistical Association}
119, 2760--2771.
}
\seealso{
The source links and instructions are installed at
\code{system.file("replication", "empirical-data.md", package = "radialFundSelect")}.
}
\examples{
\dontrun{
result <- reproduce_empirical("CaseStudy/mutual_fund.mat",
                             "F-F_Research_Data_Factors_CSV.zip",
                             figure = "empirical_T69.pdf")
result$summary
}
}
