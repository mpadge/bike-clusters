% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/get-num-clusts.R
\name{get.num.clusts}
\alias{get.num.clusts}
\title{get.num.clusts}
\usage{
get.num.clusts(city = "nyc", method = "complete")
}
\arguments{
\item{city}{nyc, washington, chicago, boston, london (case insensitive)}

\item{method}{= (ward, k-means, complete)}
}
\value{
data frame
}
\description{
Just loads up the output of calc.pnc(), and returns the final numbers of
clusters corresponding to minimal joint probabilities of peak spacings and
depths.
}

