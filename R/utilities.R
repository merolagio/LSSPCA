# # # # print_spca
#' \strong{Prints the sparse loadings from an spca object}
#'
#' Prints sparse loadings omitting the zero ones and giving the cumulative
#' variance explained.
#'
#'
#' @param obj An spca object.
#' @param only.nonzero  Logical: if = TRUE only the nonzero loadings are printed.
#' otherwise all loadings are printed.
#' @param contributions Logical: should the loadings be standardised to unit \eqn{L_1}
#' norm (and printed as percentage contributions)?
#' @param digits Integer: number of decimal figures.
#' @param threshold Value below which loadings are considered zero and not
#' printed.
#' @param rtn Logical: should the formatted (text) table be returned?
#' @param  prn Logical: should the table be printed?
#' @return If rtn = TRUE, it returns a text table formatted as specified by the
#' arguments.
#' @seealso Examples in \link{lsspca}.
#' @export
print_spca = function(obj, contributions = TRUE, only.nonzero = TRUE,
                      digits = 1, threshold = 0.001,
                      prn = TRUE, rtn = FALSE){
  if(is.list(obj)){
    if (contributions){
      if (!is.null(obj$contributions)){
        A = obj$contributions
      }
      else
       if (!is.null(obj$loadings)){
         A = obj$loadings
         A = apply(A, 2, function(x) x/sum(abs(x)))
       }
       else{
         stop("need to pass an object created withlsspca or a matrix of loadings")
       }
    }
    else
      if (!is.null(obj$loadings)){
        A = obj$loadings
      }
      else
          stop("need to pass an object created withlsspca or a matrix of loadings")
  }
  else
    if (is.matrix(obj)){
      if (contributions)
        A = apply(A, 2, function(x) x/sum(abs(x)))
      else
        A = apply(A, 2, function(x) x/sqrt(sum(x^2)))
    }
    else
      stop("need to pass an object created with lsspca or a matrix of loadings")
  if(only.nonzero){
    A = A[!apply(A, 1, function(x, a) all(abs(x) < a), a = threshold), ]
  }
  A = round(A*100, digits)
  A[abs(A) < threshold] = " "
  A <- format(A, drop0trailing = TRUE, justify = "centre")
  if (prn)
    print(A, quote = FALSE)
  if (rtn)
    return(A)
}

#' \strong{ Prints summaries from an spca object}
#'
#' Prints summaries and comparisons with the full PCA solutions for a set of LS
#' SPCA loadings.
#'
#' The summaries are printed as formatted text, if rtn = TRUE, the value
#' returned is a numerical matrix.
#'
#'
#' @param obj An spca object.
#' @param rtn Logical: should the summary matrix of summaries be returneded?
#' @param prn Logical: should anything be printed? Takes priority on prnload.
#' @return If rtn = TRUE, a numerical matrix with the summaries.
#' @details the values printed are
#' VEXP is the percentage variance explained\cr
#' CVEXP is the percentage cumulative variance explained\cr
#' RCVEXP  is the percentage cumulative variance explained relative to that of the corresponding principal components\cr
#' Card  is the cardinality, that is the number of non zero loadings\cr
#' @seealso Examples in \code{\link{lsspca}}.
#' @export
summary_spca = function(obj, rtn = FALSE, prn = TRUE){
  if (any(c(is.null(obj$vexp), is.null(obj$cvexp), is.null(obj$rcvexp),
    is.null(obj$cardinality))))
      stop("Need to pass an object created with lsspca")
  else{
    out = with(obj, rbind(round(rbind(vexp, cvexp, rcvexp)*100, 1),
                            cardinality))
    dimnames(out) = list(c("Vexp", "Cvexp", "Rcvexp", "card"),
                         paste("Comp", 1:obj$ncomps))
  }
  if (prn)
    print(out)
  if (rtn)
    return(out)
}
