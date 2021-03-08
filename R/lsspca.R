
##------------------------------------------------------------##
## lsspca FUNCTION FOR TUTORIAL PAPER
## Author Giovanni Merola please acknowledge my work
## email giovanni.merola@xjtlu.edu.cn
##------------------------------------------------------------##



#' @title Computes LS SPCA components using different variable selection algorithms
#'
#'
#' @description For each component, the variables are selected so as to explain
#' a percentage \emph{alpha} of the vexp by the corresponding principal component.
#' \emph{ind_blocks} is a list containing the indices for each component,
#'
#' @usage lsspca(X, alpha = 0.95, maxcard = 0, ncomps = 4,
#' spcaMethod = c("u", "c", "p"), scalex = FALSE,
#' subsetSelection = c("exhaustive", "seqrep", "backward", "forward", "lasso"),
#' really.big = FALSE, force.in = NULL, force.out = NULL, selectfromthese = NULL,
#' lsspca_forLasso = TRUE, lasso_penalty = 0.5)
#'
#' @param X The data matrix or data.frame.
#' @param alpha Real in [0,1]. percentage of variance of the PCs explained by the sparse component.
#' @param maxcard an integer vector or an integer. Missing values filled with last value.
#' @param ncomps number of components to compute
#' @param spcaMethod how LS SPCA components are computed:
#' 'u' for uncorrelated, 'c' for correlated and 'p' for projection
#' @param scalex Logical, if TRUE variables are scaled to unit variance.default  FALSE
#'    Variables are automatically centered to zero if they aren't already.
#' @param subsetSelection how the variables for each component are selected
#' 'exhaustive' all subsets, 'seqrep' stepwise, 'backward', 'forward', 'lasso'
#' @param really.big Must be TRUE to perform exhaustive search on more than 50 variables.
#' @param force.in NULL or list of indices that must be in component. not for lasso. [NULL]
#' @param force.out NULL or list of indices cannot be in component. [NULL]
#' @param selectfromthese NULL or list of indices from which model chosen. [NULL]
#' @param lsspca_forLasso use lsspca with indices selected with lasso or just the lasso regression
#' @param lasso_penalty real between 0 and 1; 0 -> ridge regression, 1 -> lasso
#'
#' @details  for USPCA, \code{maxcard} cannot be smaller than the order of the components
#'    computed, so \code{maxcard = c(1, 1, 1)} will be automatically changed to
#'    \code{maxcard = c(1, 2, 3)}. Exaustive search can be slow for matrices with
#'    30 or more variables. See the documentation for \code{regsubset} in the package
#'    \code{leaps} for the option \code{really.big}.
#' @return a list
## #' \loadmathjax
#' \describe{
## #' \item{loadings}{Matrix with the loadings scaled to unit \mjeqn{L_2}{ASCII representation} norm.}
## #' \item{contributions}{Matrix of loadings scaled to unit \eqn{L_1}{ASCII representation} norm.}
#' \item{loadings}{Matrix with the loadings scaled to unit \eqn{L_2} norm.}
#' \item{contributions}{Matrix of loadings scaled to unit \eqn{L_1} norm.}
#' \item{ncomps}{integer number of components computed. Default is 4.}
#' \item{cardinality}{Vector with the cardinalities of each loadings.}
#' \item{ind}{List with the indices of the non-zero loadings for each component.}
#' \item{loadlist}{A list with only the nonzero ladings for each component.}
#' \item{vexp}{Vector with the \% variance explained by each component.}
#' \item{vexpPC}{Vector with the \% variance explained by each principal component.}
#' \item{cvexp}{Vector with the \% cumulative variance explained by each component.}
#' \item{rcvexp}{Vector with the \% proportion of cumulative variance explained by each component to that explained by the PCs.}
#' \item{scores}{the SPCs scores.}
#' \item{PCloadings}{Matrix with the PCs loadings scaled to unit \eqn{L_2} norm. }
#' \item{PCs scores}{the PCs scores.}
#' \item{method}{method used to compute the loadings}
#' \item{corComp}{Matrix of correlations among the sparse components. Only if ncomps > 1.}
#' \item{Call}{The called with its arguments.}
#' }
#' @references Giovanni M. Merola. 2014. \emph{Least Squares Sparse Principal
#' Component Analysis: a Backward Elimination approach to attain large
#' loadings.} Austr.&NZ Jou. Stats. 57, pp 391-429\cr\cr
#' Giovanni M. Merola and Gemai Chen. 2019. \emph{Sparse Principal Component Analysis: an
#' efficient Least Squares approach.} Jou. Multiv. Analysis 173, pp 366--382
#' \url{http://arxiv.org/abs/1406.1381}

#' @examples
#' \dontrun{
#' library(LSSPCA)
#' data(hitters)
#'
#' dim(hitters)
#' ## USPCA 95
#' hit_uspca95 = lsspca(X = hitters, alpha = 0.95, ncomps = 4,
#'                      spcaMethod = "u", subsectSelection = "e")
#' #> Warning message:
#' #>  In log(vr) : NaNs produced
#' ## the warnings come from the variable selection, don't worry
#'
#' ## summaries
#' t(data.frame(card = hit_uspca95$cardinality,
#'              cvexp = round(hit_uspca95$cvexp, 2),
#'              rcvexp = round(hit_uspca95$rcvexp, 2)))
#'
#' ## print loadings individually
#' lapply(hit_uspca95$loadlist, function(x) round(x, 2))
#' ## print contributions individually
#' lapply(hit_uspca95$loadlist, function(x) round(x/sum(abs(x)), 2))
#'
#' ## plot PC and USPC loadings
#' par(mfrow = c(1, 2))
#' barplot(-hit_uspca95$PCloadings[, 1], main = "PCA")
#' barplot(-hit_uspca95$loadings[, 1], main = "USPCA")
#' par(mfrow = c(1,1))
#'
#' ## Holzinger data
#' data(holzinger)
#' dim(holzinger)
#'
#' ## CSPCA
#' hol_cspca95 = lsspca(X = holzinger, alpha = 0.95, ncomps = 4,
#'                      spcaMethod = "c", subsectSelection = "e")
#'
#' ## summaries
#' t(data.frame(card = hol_cspca95$cardinality,
#'              cvexp = round(hol_cspca95$cvexp, 2),
#'              rcvexp = round(hol_cspca95$rcvexp, 2)))
#'
#' ## print loadings
#' lapply(hol_cspca95$loadlist, function(x) round(x, 2))
#' ## print contributions
#' lapply(hol_cspca95$loadlist, function(x) round(x/sum(abs(x)), 2))
#'
#' ## correlation between SPCs
#' round(hol_cspca95$corComp, 2)
#'
#' ## plot contributions
#' barplot(-hol_cspca95$contributions[, 1])
#'
#' ## SPCs scores against PC scores
#' plot(hol_cspca95$scores[, 1], hol_cspca95$PCscores[, 1], pch = 16)
#' regline = lm(hol_cspca95$PCscores[, 1] ~ hol_cspca95$scores[, 1]- 1)$coef
#' abline(a = 0, b = regline, col = 2)
#'
#'
#' ## SPCA on each ability separately
#' h_groups = lapply(seq(1, 10, 3), function(x) x:(x + 2))
#'
#' ## projection SPCA
#' hol_block_spca95 = lsspca(X = holzinger, alpha = 0.95, ncomps = 4,
#'                      spcaMethod = "p", subsectSelection = "e",
#'                      selectfromthese = h_groups)
#'
#' ## summaries
#' t(data.frame(card = hol_block_spca95$cardinality,
#'              cvexp = round(hol_block_spca95$cvexp, 2),
#'              rcvexp = round(hol_block_spca95$rcvexp, 2)))
#'
#' ## print loadings
#' lapply(hol_block_spca95$loadlist, function(x) round(x, 2))
#'
#' ## print contributions
#' lapply(hol_block_spca95$loadlist, function(x) round(x/sum(abs(x)), 2))
#'
#' ## correlation between SPCs
#' round(hol_block_spca95$corComp, 2)
#'
#' ## plot the contributions for each SPC
#' par(mfrow = c(2, 2))
#' for(k in 1:4){
#'   barplot(-hol_block_spca95$contributions[, k])
#' }
#' par(mfrow = c(1, 1))
#' }
#'
#' @author Giovanni Merola
#'
#' @export
lsspca <- function(X, alpha = 0.95, maxcard = 0, ncomps = 4,
                   spcaMethod = c("u", "c", "p"), scalex = FALSE,
                   subsetSelection = c("exhaustive", "seqrep", "backward", "forward", "lasso"),
                   really.big = FALSE,
                   force.in = NULL, force.out = NULL, selectfromthese = NULL,
                   lsspca_forLasso = TRUE, lasso_penalty = 0.5) {

  ##----------------------------------------##
  ## Creating the environment and error checking
  ##----------------------------------------##

  p <- ncol(X)
  n <- nrow(X)
  if (is.null(colnames(X)))
    colnames(X) <- paste0("V", 1:p)
  ## convert selectfromthese to force out
  if (is.list(selectfromthese)) {
    fx <- function(x, m) (1:m)[-x]
    force.out <- lapply(selectfromthese, fx, m = p)
  }
  if (is.list(force.in) | is.list(force.out)) {
    ncomps <- max(length(force.in), length(force.out))
    if (length(force.in) < ncomps)
      force.in <- c(force.in, vector("list", ncomps - length(force.in)))
    if (length(force.out) < ncomps)
      force.out <- c(force.out, vector("list", ncomps - length(force.out)))
  }
  else if (ncomps == 0) {
    ncomps <- p
  }
  if (!any(spcaMethod == c("u", "c", "p")))
    stop("only one of these options allowed for spcaMethod: 'u', 'c', 'p'")
  if (length(maxcard) > 1)
    ncomps <- length(maxcard)
  if (is.vector(maxcard)) {
    maxcard[maxcard == 0] <- p
    if (length(maxcard < ncomps)) {
      lm <- length(maxcard)
      maxcard <- c(maxcard, rep(maxcard[lm], ncomps - lm))
    }
  } else {
    stop("must pass a vector or an integer as maxcard")
  }
  subsetSelection <- switch(stringr::str_sub(subsetSelection[1], 1, 1), e = "exhaustive", b = "backward", f = "forward",
                             s = "seqrep", l = "lasso")
  if (is.null(subsetSelection))
    stop("please pass a valid variable selection option")
  if((ncol(X)> 30) & (stringr::str_sub(subsetSelection[1], 1, 1) == "e"))
    if(really.big == FALSE){
      warning("exhaustive search with so many variables requires really.big = TRUE")
      really.big == TRUE
  }
  spcaMethod <- spcaMethod[1]

  if (!requireNamespace("geigen", quietly = TRUE)) {
    stop("Package \"geigen\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package \"stringr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (subsetSelection == "lasso")
    if (!requireNamespace("elasticnet", quietly = TRUE)) {
      stop("Package \"elasticnet\" needed for lasso selection to work. Please install it.",
           call. = FALSE)
    }
  else
    if (!requireNamespace("leaps", quietly = TRUE)) {
      stop("Package \"leaps\" needed for variable selection to work. Please install it.",
           call. = FALSE)
    }
  ##-------------------------------------------------------------------##
  ## create objects for output
  ##-------------------------------------------------------------------##
  card <- rep(0, ncomps)
  ind <- as.list(card)
  vexp <- rep(0, ncomps)
  cvexp <- rep(0, ncomps)
  A <- matrix(0, p, ifelse(ncomps == 0, p, ncomps))
  contributions <- A
  loadlist <- as.list(1:p)
  scores <- matrix(0, n, ncomps)
  if (spcaMethod[1] == "u")
    R <- matrix(0, p, ncomps)

  if (is.null(colnames(X)))
    namx <- paste("Var", 1:p) else namx <- colnames(X)

  ##-------------------------------------------------------------------##
  ## center and scale the variables
  ##-------------------------------------------------------------------##
  if (scalex == TRUE)
    X <- scale(X)
  else
    if (any(abs(colMeans(X)) > 10^-6)) {
      message("Centering the columns of X to zero mean")
      X <- scale(X, scale = FALSE)
    }
  ##-------------------------------------------------------------------##
  ## compute of the components
  ##-------------------------------------------------------------------##
  K <- X
  nc <- 0
  j <- 1
  stopComp <- FALSE
  while (stopComp == FALSE) {
    D <- crossprod(K)
    if (j == 1)
      Cmat <- D
    r_ee <- eigen(D)
    pc <- K %*% r_ee$vec[, 1]
    lambda <- r_ee$val[1]
    if (j == 1) {
      vexpPC <- r_ee$val
      totvar <- sum(r_ee$val)
      PCloadings = r_ee$vec[, 1:ncomps]
      PCscores = X %*% PCloadings
    }

    ## variable selection ==========================
    if (stringr::str_sub(subsetSelection[1], 1, 1) == "l") {
      if ((length(force.out) > 0) | (length(force.in) > 0))
        stop("lasso not implemented with force.out or force.in")
      else g_fit <- glmnet::glmnet(X, pc, intercept = FALSE, alpha = lasso_penalty)
      mod_ind <- ifelse(any(g_fit$dev.ratio > alpha), min(which(g_fit$dev.ratio > alpha)), length(g_fit$dev.ratio))
      coe <- as.numeric(stats::coef(g_fit, s = g_fit$lambda[mod_ind]))[-1]
      ind[[j]] <- which(coe != 0)
      card[j] <- length(ind[[j]])
      a <- coe[ind[[j]]]
    }
    else {
      ssr <- leaps::regsubsets(x = X, y = pc, spcaMethod = subsetSelection[1],
                               nbest = 1, force.in = force.in[[j]],
                               force.out = force.out[[j]], nvmax = maxcard[j],
                               intercept = FALSE, really.big = really.big)

      aa <- summary(ssr)
      aa$which <- aa$which[, order(match(colnames(aa$which), colnames(X)))]
      mrsq <- any(aa$rsq >= alpha)
      if (mrsq == TRUE)
        if (spcaMethod == "u")
          iord <- which((aa$rsq >= alpha) & ((1:length(aa$rsq)) >= j))  ## rsquared
      else iord <- which(aa$rsq >= alpha)  ## rsquared
      else {
        iord <- length(aa$rsq)
        if (maxcard[j] > min(iord))
          message(paste("warning: lsspca component ", j, "could not reach Rsquared", alpha))
      }

      indall <- 1:p
      ind[[j]] <- indall[aa$which[iord[1], ]]
      card[j] <- length(ind[[j]])
    }
    ##  compute LSSPCs =====================
    if ((subsetSelection != "l") | ((subsetSelection == "l") & (lsspca_forLasso == TRUE))) {
      if (length(ind[[j]]) > 1) {
        Xd <- X[, ind[[j]]]
        Sd <- crossprod(Xd)
        if (spcaMethod == "u")
        {
          if (j == 1) {
            M <- tcrossprod(crossprod(Xd, X))
            ga <- geigen::geigen(M, Sd, symmetric = TRUE)
            a <- ga$vec[, card[j]]
          }
          else {
            xxd <- crossprod(X, Xd)
            H <- crossprod(A[, 1:(j - 1), drop = FALSE], xxd)
            Sm <- solve(Sd)
            E <- H %*% Sm
            G <- diag(card[j]) - crossprod(H, solve(tcrossprod(E, H))) %*% E
            M <- crossprod(tcrossprod(xxd, G))
            ga <- geigen::geigen(M, Sd)
            a <- ga$vec[, card[j]]
          }
        }
        else {
          if (spcaMethod == "c") {
            M <- crossprod(crossprod(K, Xd))
            ga <- geigen::geigen(M, Sd, symmetric = TRUE)
            a <- ga$vec[, card[j]]
          }
          if (spcaMethod == "p") {
            alm <- lm(pc ~ Xd - 1)
            a <- alm$coefficients
          }
        }
        ## save results for j-th SPC ==========================

        a <- a/sqrt(sum(a^2))
        if (all(a <= 0))
          a <- -a
        scores[, j] <- Xd %*% a
        bb <- stats::cor(pc, scores[, j])
        if (bb < 0) {
          a <- -a
          scores[, j] <- -scores[, j]
        }
        names(a) <- namx[ind[[j]]]
        A[ind[[j]], j] <- a
        contributions[ind[[j]], j] <- a/sum(abs(a))
        loadlist[[j]] <- a
      }
      else {
        a <- 1
        scores[, j] <- X[, ind[[j]]]
        bb <- cor(pc, scores[, j])
        if (bb < 0) {
          a <- -a
          scores[, j] <- -scores[, j]
        }
        names(a) <- namx[ind[[j]]]
        A[ind[[j]], j] <- a
        contributions[ind[[j]], j] <- a
        loadlist[[j]] <- a
      }
    }
    ## deflate X
    if (j <= ncomps)
      K <- K - tcrossprod(scores[, j]) %*% K/sum(scores[, j]^2)  #spcaTutoPack:::deflXCforR(c(scores[, j]), K)
    cvexp[j] <- totvar - sum(K^2)
    vexp[j] <- ifelse(j > 1, cvexp[j] - cvexp[j - 1], cvexp[j])
    nc <- nc + 1
    if (j < ncomps)
      j <- j + 1 else stopComp <- TRUE
  }
  ##-------------------------------------------------------------------##
  ## create output
  ##-------------------------------------------------------------------##
  ncomps <- nc
  rownames(contributions) <- namx
  rownames(A) <- namx
  if (ncomps > 1) {
    dimnames(PCloadings) = dimnames(A)
    dimnames(PCscores) = dimnames(scores)
  }
  else{
    names(PCloadings) = rownames(A)
  }

  out <- list(loadings = A[, 1:ncomps], contributions = contributions,
              ncomps = ncomps, cardinality = card[1:ncomps], indices = ind,
              loadlist = loadlist[1:ncomps],
              vexp = vexp[1:ncomps]/totvar, vexpPC = vexpPC[1:ncomps]/totvar,
              cvexp = cvexp[1:ncomps]/totvar, rcvexp = cvexp[1:ncomps]/cumsum(vexpPC[1:ncomps]),
              scores = scores[, 1:ncomps],  PCloadings = PCloadings,
              PCscores = PCscores, spcaMethod = spcaMethod)
  if (ncomps > 1) {
    if (spcaMethod[1] != "u") {
      out$corComp <- cor(scores)
    }
  }
  out$Call <- match.call()
  return(out)
}

