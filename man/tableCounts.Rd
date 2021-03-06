% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tableCounts.R
\name{tableCounts}
\alias{tableCounts}
\title{tableCounts}
\usage{
tableCounts(gammalist, nobs.a, nobs.b, n.cores)
}
\arguments{
\item{gammalist}{A list of objects produced by gammaKpar, gammaCK2par, or
gammaCKpar.}

\item{nobs.a}{number of observations in dataset 1}

\item{nobs.b}{number of observations in dataset 2}

\item{dedupe}{logical - indicator for internal linkage / deduplication of a single dataset}

\item{n.cores}{Number of cores to parallelize over. Default is NULL.}
}
\value{
\code{tableCounts} returns counts of all unique mathching patterns, which can be
fed directly into \code{emlinkMAR} to get posterior matching probabilities for each unique pattern.
}
\description{
Count pairs with the same pattern in the cross product between two datasets.
}
\examples{
\dontrun{
## Calculate gammas
g1 <- gammaCKpar(dfA$firstname, dfB$firstname)
g2 <- gammaCKpar(dfA$middlename, dfB$middlename)
g3 <- gammaCKpar(dfA$lastname, dfB$lastname)
g4 <- gammaKpar(dfA$birthyear, dfB$birthyear)

## Run tableCounts
tc <- tableCounts(list(g1, g2, g3, g4), nobs.a = nrow(dfA), nobs.b = nrow(dfB))
}
}
\author{
Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
}
