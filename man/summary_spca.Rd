% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{summary_spca}
\alias{summary_spca}
\title{\strong{ Prints summaries from an spca object}}
\usage{
summary_spca(obj, rtn = FALSE, prn = TRUE)
}
\arguments{
\item{obj}{An spca object.}

\item{rtn}{Logical: should the summary matrix of summaries be returneded?}

\item{prn}{Logical: should anything be printed? Takes priority on prnload.}
}
\value{
If rtn = TRUE, a numerical matrix with the summaries.
}
\description{
Prints summaries and comparisons with the full PCA solutions for a set of LS
SPCA loadings.
}
\details{
The summaries are printed as formatted text, if rtn = TRUE, the value
returned is a numerical matrix.

the values printed are
VEXP is the percentage variance explained\cr
CVEXP is the percentage cumulative variance explained\cr
RCVEXP  is the percentage cumulative variance explained relative to that of the corresponding principal components\cr
Card  is the cardinality, that is the number of non zero loadings\cr
}
\seealso{
Examples in \code{\link{lsspca}}.
}
