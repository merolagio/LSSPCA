

#' Utility for computing the variance explained by a set of components
#'
#' @usage makevexp(scores, X, unc = FALSE)
#'
#' @param scores, a matrix of components scores
#' @param X the data matrix.
#' @param unc logical, are the components orthogonal.
#'
#' @return a vector with the variance explained by each component
#'
#' @examples{
#'  if(FALSE){
#'    A = matrix(c(rep(1, 16), rep(c(-1, 0.5, 1, 0.5), 4)), 16, 2)
#'    sc = hitters %*% A
#'    makevexp(sc, hitters, FALSE)
#'
#'  ## this is wrong
#'    makevexp(sc, hitters, TRUE)
#'
#'  ## because the components are not orthogonal
#'    cor(sc)
#'  }
#' }
#' @author Giovanni Merola
#'
#' @export
makevexp = function(scores, X, unc = FALSE){
  n = ncol(scores)
  vexp = rep(0, n)
  if(unc){
    Q <- sweep(scores, 2, sqrt(colSums(scores^2)), "/")
    vexp = colSums((crossprod(X, Q))^2)
  }
  else{
    vexp[1] = sum((crossprod(X, scores[, 1, drop = FALSE])^2))/sum(scores[, 1]^2)
    for (j in 2:n){
      Q = qr.Q(qr(scores[, 1:j]))
      vexp[j] = colSums((crossprod(X, Q[, j, drop = FALSE]))^2)
    }
  }
  return(vexp)
}
