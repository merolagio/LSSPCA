
##------------------------------------------------------------##
## lsspca FUNCTION FOR TUTORIAL PAPER
## Author Giovanni Merola please acknowledge my work
##------------------------------------------------------------##



#' @title Computes LS SPCA components using different variable selection algorithms
#'
#' @description  For each component, the variables are selected so as to explain
#' a percentage \emph{alpha} of the vexp by the corresponding principal component.
#' \emph{ind_blocks} is a list containing the indeces for each component,
#'
#' @param X The data matrix.
#' @param alpha Real in [0,1]. percentage of variance of the PCs explained by the sparse component.
#' @param ncomps number of components to compute
#' @param spcaMethod character vector how LS SPCA components are computed:
#' "u" for uncorrelated, "c" for correlated and "p" for projection.
#' If only one value, the same method is used for all components.
#' @param variableSelection how the variables for each components are selected
#' 'seqrep' stepwise, 'exhaustive' all subsets 'backward', 'forward', 'lasso'
#' @param scalex = FALSE, whether to scale the variables to unit variance. Variables are
#' scaled to zero mean (if needed) even if scaleX = FALSE
#' @param maxcard a vector or an integer. Missing values filled with last value.
#' @param force.in NULL or list of indeces that must be in component. not for lasso. [NULL]
#' @param force.out NULL or list of indeces cannot be in component. [NULL]
#' @param selectfromthese NULL or list of indeces from which model chosen. [NULL]
#' @param lsspca_forLasso use lsspca with indeces selected with lasso or just the lasso regression
#' @param lasso_penalty real between 0 and 1, , 0-> ridge regression, 1 -> lasso
#' @return a list  with vif for each component
#'
#' @author Giovanni Merola
#'
#' @export
lsspca <- function(X, alpha = 0.95,  ncomps = 4,
                   spcaMethod = "u",
                   variableSelection = c("exhaustive","seqrep", "backward", "forward", "lasso"),
                   scalex = FALSE, maxcard = 0, force.in = NULL, force.out = NULL, selectfromthese = NULL,
                   lsspca_forLasso = TRUE, lasso_penalty = 0.5){

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
  me = unique(spcaMethod)
  lme = length(me)
  if (lme == 1)
    if (!any(spcaMethod == c("u", "c", "p")))
      stop("only one of these options allowed for method: 'u', 'c', 'p'")
  else
    for(m in 1:lme){
      if (!any(me[m] == c("u", "c", "p")))
        stop("only one of these options allowed for method: 'u', 'c', 'p'")
    }
  lm = length(spcaMethod)
  if (lm < ncomps){
    spcaMethod = c(spcaMethod, rep(spcaMethod[lm], ncomps - lm))
  }

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
  variableSelection <- switch(stringr::str_sub(variableSelection[1], 1, 1), e = "exhaustive", b = "backward", f = "forward",
                             s = "seqrep", l = "lasso")
  if (is.null(variableSelection))
    stop("need to pass a valid search direction")


  stopifnot(requireNamespace("geigen", quietly = TRUE))
  stopifnot(requireNamespace("stringr", quietly = TRUE))

  if (variableSelection == "lasso")
    stopifnot(requireNamespace("elasticnet", quietly = TRUE))
  else
    stopifnot(requireNamespace("leaps", quietly = TRUE))


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


  if (is.null(colnames(X)))
    namx <- paste("Var", 1:p) else namx <- colnames(X)

  ##-------------------------------------------------------------------##
  ## center and scale the variables
  ##-------------------------------------------------------------------##
  if (scalex == TRUE)
    X <- scale(X) else if (any(abs(colMeans(X)) > 10^-6)) {
      message("You need to center the columns of X to zero mean, I'll do it for you")
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
    D <- crossprod(K)  #SPCA3::ataC(K)
    if (j == 1)
      Cmat <- D
    r_ee <- eigen(D)  #SPCA3::EigenC(D)
    pc <- K %*% r_ee$vec[, 1]
    lambda <- r_ee$val[1]
    if (j == 1) {
      vexpPC <- r_ee$val
      totvar <- sum(r_ee$val)
      PCloadings = r_ee$vec[, 1:ncomps, drop = FALSE]
      PCscores = X %*% PCloadings
    }

    ## variable selection ==========================
    if (stringr::str_sub(variableSelection[1], 1, 1) == "l") {
      if ((length(force.out) > 0) | (length(force.in) > 0))
        stop("lasso not implemented with force.out or force.in")
      else g_fit <- glmnet::glmnet(X, pc, intercept = FALSE, alpha = lasso_penalty)
      mod_ind <- ifelse(any(g_fit$dev.ratio > alpha), min(which(g_fit$dev.ratio > alpha)), length(g_fit$dev.ratio))
      coe <- as.numeric(coef(g_fit, s = g_fit$lambda[mod_ind]))[-1]
      ind[[j]] <- which(coe != 0)
      card[j] <- length(ind[[j]])
      a <- coe[ind[[j]]]
    }
    else {
      ssr <- leaps::regsubsets(x = X, y = pc, method = variableSelection[1], nbest = 1, force.in = force.in[[j]],
                               force.out = force.out[[j]], nvmax = maxcard[j], intercept = FALSE)

      aa <- summary(ssr)
      aa$which <- aa$which[, order(match(colnames(aa$which), colnames(X)))]
      mrsq <- any(aa$rsq >= alpha)
      if (mrsq == TRUE)
        if (spcaMethod[j] == "u")
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
    if ((variableSelection != "l") | ((variableSelection == "l") & (lsspca_forLasso == TRUE))) {
      if (length(ind[[j]]) > 1) {
        Xd <- X[, ind[[j]]]
        Sd <- crossprod(Xd)
        if (spcaMethod[j] == "u")
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
          if (spcaMethod[j] == "c") {
            M <- crossprod(crossprod(K, Xd))
            ga <- geigen::geigen(M, Sd, symmetric = TRUE)  ##SPCA3::GenEigenC(M, Sd)
            a <- ga$vec[, card[j]]
          }
          if (spcaMethod[j] == "p") {
            alm <- lm(pc ~ Xd - 1)
            a <- alm$coefficients
          }
        }
        ## save results for j-th SPC ==========================

        a <- a/sqrt(sum(a^2))
        if (all(a <= 0))
          a <- -a
        scores[, j] <- Xd %*% a
        bb <- cor(pc, scores[, j])
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
    ## deflate K
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
  dimnames(PCloadings) = dimnames(A)
  dimnames(PCscores) = dimnames(scores)


  out <- list(loadings = A[, 1:ncomps], contributions = contributions, vexp = vexp[1:ncomps]/totvar, vexpPC = vexpPC[1:ncomps]/totvar,
              cvexp = cvexp[1:ncomps]/totvar, rcvexp = cvexp[1:ncomps]/cumsum(vexpPC[1:ncomps]),
              scores = scores[, 1:ncomps], ncomps = ncomps, indices = ind, cardinality = card[1:ncomps], loadlist = loadlist[1:ncomps],
              PCloadings = PCloadings, PCscores = PCscores, spcaMethod = spcaMethod)
  if (ncomps > 1) {
    if (any(spcaMethod != "u")) {
      out$corComp <- cor(scores)
    }
  }
  out$Call <- match.call()
  return(out)
}
