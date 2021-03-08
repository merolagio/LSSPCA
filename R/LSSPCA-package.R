#'
#' LSSPCA: A Function to compute Least Squares Sparse Principal Components.
#'
#' This package is an appendix to the paper \emph{my paper} and comes with
#' only function.
#' Sparse principal components are combinations of only few of the observed
#' variables.
#' The LS SPCA components give the best possible data approximation to the
#' data under sparsity constraints.
#' The \code{lsspca} function takes a data matrix and a goodness of fit
#' specification, (percent variance explained) or a cardinality.
#' Optimal orthogonal components are computed choosing \code{method = "u"},
#' suboptimal correlated components with \code{method = "c"} and less
#' computationally demanding ones with \code{method = "p"}.
#'
#' Subsets of variables can be chosen with different search algorithms and
#' variables can be foced in or out from these subsets.
#'
#'
#' @name LSSPCA-package
#' @docType package
## #' @importFrom leaps summary.regsubsets
#' @importFrom stats cor coef
#' @references Giovanni M. Merola. 2014. \emph{Least Squares Sparse Principal
#' Component Analysis: a Backward Elimination approach to attain large
#' loadings.} Austr.&NZ Jou. Stats. 57, pp 391-429\cr\cr
#' Giovanni M. Merola and Gemai Chen. 2019. \emph{Sparse Principal Component Analysis: an
#' efficient Least Squares approach.} Jou. Multiv. Analysis 173, pp 366--382 \url{http://arxiv.org/abs/1406.1381}
#' @keywords package
#' @seealso \code{\link{lsspca}} for usage examples.
## #' @useDynLib LSSPCA
NULL
