\name{bh}
\alias{bh}
\title{Benjamini--Hochberg step-up rule}
\description{Selects all p-values at or below the largest ordered p-value
passing its Benjamini--Hochberg threshold.}
\usage{bh(pvalue, q = 0.1)}
\arguments{
\item{pvalue}{Finite numeric vector of p-values in \eqn{[0,1]}.}
\item{q}{Nominal level in \eqn{(0,1)}.}
}
\value{A logical vector of selections, with the same length as \code{pvalue}.
If no ordered p-value passes, every entry is false.}
\examples{bh(c(0.001, 0.02, 0.2, 0.9), q = 0.1)}
