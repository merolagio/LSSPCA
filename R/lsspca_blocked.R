##------------------------------------------------------------##
## Author Giovanni Merola please acknowledge my work
## email giovanni.merola@xjtlu.edu.cn
##------------------------------------------------------------##

#' @title Computes LS SPCA components on selected groups of variables
#'
#' @description Each component, the variables are selected so as to explain
#' a percentage \emph{alpha} of the variance explained by the corresponding principal component.
#' \emph{blocks_list} is a list containing the indeies for each component,
#' Subsets can be overlapping and need not be exhaustive.
#'
#' @usage lsspca_blocked(X, alpha = 0.95, blocks_list = list(),
#' ncomps_per_block = 1, blocks_names = NA, maxcard = 0,
#' spcaMethod = c("u", "c", "p"), scalex = FALSE,
#' variableSelection = c("exhaustive", "seqrep", "backward", "forward", "lasso"),
#' lsspca_forLasso = FALSE, lasso_penalty = 0.5, rtn_all_spca = FALSE)
#'
#' @param X the data matrix.
#' @param alpha Real (or vector) in [0,1]. percentage of variance of the PCs explained by the sparse component.
#' @param blocks_list a list of indices or a vector or factor of codes for each block
#' @param ncomps_per_block number of components per block, integer or vector of. Default = 1.
#' @param blocks_names names of each block, if possible names are taken for the names of blocks_list
#' @param maxcard a vector or an integer of maximal cardinality for each block, if USPCA selected
#'   maxcard[j] cannot be < j
#' @param spcaMethod char how lsspca is computed "u" = USPCA (default), "c" = CSPCA, "p" = PSPCA
#' @param scalex Logical, if TRUE variables are scaled to unit variance.default  FALSE
#'    Variables are automatically centered to zero if they aren't already.
#' @param variableSelection how the variables in each components are selected
#'  "exhaustive" = all subsets, "seqrep" = stepwise, "backward", "forward", "lasso".
#'  See documentation in packages \code{leaps}, and \code{glmnet} for lasso.
#' @param lasso_penalty real between 0 and 1, , 0-> ridge regression, 1 -> lasso
#' @param lsspca_forLasso use lsspca with indeies selected with lasso, otherwise
#' just the lasso regression
#' @param rtn_all_spca logical, should the lsspca results for each blovck be returned?
#'
#' @details Differently from the lsspca function, the components are computed so as
#' to explain a given percentage of the variance explained by the PCs of the blocks
#' of variables under consideration.
#'
#' @return a list
#' \describe{
#' \item{loadings}{Matrix with the loadings scaled to unit \eqn{L_2} norm.}
#' \item{contributions}{Matrix of loadings scaled to unit \eqn{L_1} norm.}
#' \item{ncomps}{integer number of components computed. Default is 4.}
#' \item{cardinality}{Vector with the cardinalities of each loadings.}
#' \item{ind}{List with the indices of the non-zero loadings for each component.}
#' \item{loadingslist}{A list with only the nonzero ladings for each component.}
#' \item{vexp}{Vector with the \% variance explained by each component.}
#' \item{vexpPC}{Vector with the \% variance explained by each principal component.}
#' \item{cvexp}{Vector with the \% cumulative variance explained by each component.}
#' \item{rcvexp}{Vector with the \% proportion of cumulative variance explained by each component to that explained by the PCs.}
#' \item{scores}{the SPCs scores.}
#' \item{PCloadings}{Matrix with the PCs loadings scaled to unit \eqn{L_2} norm. }
#' \item{PCscores}{the PCs scores.}
#' \item{method}{method used to compute the loadings}
#' \item{corComp}{Matrix of correlations among the sparse components. Only if ncomps > 1.}
#' \item{Call}{The called with its arguments.}
#' }
#'
#' @references Giovanni M. Merola. 2014. \emph{Least Squares Sparse Principal
#' Component Analysis: a Backward Elimination approach to attain large
#' loadings.} Austr.&NZ Jou. Stats. 57, pp 391-429\cr\cr
#' Giovanni M. Merola and Gemai Chen. 2019. \emph{Sparse Principal Component Analysis: an
#' efficient Least Squares approach.} Jou. Multiv. Analysis 173, pp 366--382
#' \url{http://arxiv.org/abs/1406.1381}
#'
#' @seealso \code{\link{lsspca}}
#' @export
lsspca_blocked = function(X, alpha = 0.95, blocks_list = list(),
                          ncomps_per_block = 1, blocks_names = NA,
                          maxcard = 0, spcaMethod = c("u", "c", "p"), scalex = FALSE,
                          variableSelection = c("exhaustive", "seqrep", "backward", "forward", "lasso"),
                          lsspca_forLasso = FALSE, lasso_penalty = 0.5, rtn_all_spca = FALSE){
  p = ncol(X)
  n = nrow(X)
  ##-------------------------------------------##
  ## Creating the environment and error checking
  ##-------------------------------------------##
  nblocks = length(blocks_list)
  if(nblocks <= 1){
    stop("must pass at least of 2 sets of indeies\ otherwise, use lsspca")
  }
  if(!all(sapply(blocks_list, function(x) all(is.integer(x))))){
      if(is.vector(blocks_list) | is.factor(blocks_list))
        blocks_list <- tapply(1:length(blocks_list), blocks_list, function(x) x)
      else
        stop("please pass a list of indeies or a vector or factor of codes as block_list")
  }

  if ((ncomps_per_block[1] == 0) | any(ncomps_per_block > sapply(blocks_list, length))){
    stop("ncomps_per_block cannot be zero or larger than block size")
  }
  else{
    ncomps_per_block = makevec(vec = ncomps_per_block, n = nblocks)
  }
  if (is.na(ncomps_per_block[1]))
    ncomps_per_block = rep(1, nblocks)
  else{
    ncomps_per_block = makevec(vec = ncomps_per_block, n = nblocks)
  }
  if(is.na(blocks_names[1])){
    if (!is.null(names(blocks_list)))
      blocks_names = rep(names(blocks_list), times = ncomps_per_block)
    else
      blocks_names = paste0("Block", rep(1:nblocks, times = ncomps_per_block))
  }
  else
    if (length(blocks_names) == length(blocks_list))
      blocks_names = rep(blocks_names, times = ncomps_per_block)
  else{
    warning("something wring with blocks_names")
    blocks_names = paste0("Block", rep(1:nblocks, times = ncomps_per_block))

  }

    maxcard = makevec(maxcard, nblocks)
  if (any(maxcard == 0)){
    maxcard[maxcard == 0] = sapply(blocks_list[maxcard == 0], length)
   # message("0 maxcard set to maximum")
  }

  alpha = makevec(alpha, nblocks)

  variableSelection = switch(stringr::str_sub(variableSelection[1], 1, 1), "e" = "exhaustive",
                     "b"= "backward", "f"= "forward", "s"= "seqrep", "l" = "lasso")
  if (is.null(variableSelection))
    stop("need to pass a valid search direction")
  spcaMethod = spcaMethod[1]
  ##-------------------------------------------------------------------##
  ## create objects for output
  ##-------------------------------------------------------------------##
  ncomps <- sum(ncomps_per_block)
  card = rep(0, ncomps)
  ind = as.list(ncomps)
  vexp = rep(0, ncomps)
  cvexp = rep(0, ncomps)
  loadings = matrix(0, p, ncomps)
  contributions = loadings
  loadingslist = as.list(1:ncomps)
  scores = matrix(0, n, ncomps)

  if (is.null(colnames(X)))
    namx = paste("Var", 1:p)
  else
    namx = colnames(X)

  if (scalex == TRUE)
    X = scale(X)
  else
    if (any(abs(colMeans(X)) > 10^-6)){
      message("Centering columns of X to zero mean")
      X = scale(X, scale = FALSE)
    }

  spca_list = list()
  ind_list = list()
  loadingslist = list()
  ncompsdone = 0
  ##-------------------------------------------------------------------##
  ## compute of the components
  ##-------------------------------------------------------------------##
  for (j in 1:nblocks){
    spca_list[[j]] <- lsspca(X[,  blocks_list[[j]]], alpha = alpha[j],
                             maxcard = maxcard[j], ncomps = ncomps_per_block[j],
                             spcaMethod = spcaMethod, scalex = FALSE,
                             variableSelection = variableSelection,
                             force.in = NULL, force.out = NULL, selectfromthese = NULL,
                             lsspca_forLasso = lsspca_forLasso, lasso_penalty = lasso_penalty)

    ind = lapply(spca_list[[j]]$ind, function(x, ii) ii[x], ii = blocks_list[[j]])
    ind_list = c(ind_list, ind)
    loadingslist = c(loadingslist, spca_list[[j]]$loadingslist)
    scores[, (ncompsdone + 1):(ncompsdone + ncomps_per_block[j])] = spca_list[[j]]$scores
    loadings[blocks_list[[j]], (ncompsdone + 1):(ncompsdone + ncomps_per_block[j])] =
      spca_list[[j]]$loadings
    ncompsdone =  ncompsdone + ncomps_per_block[j]
  }
  ##-------------------------------------------------------------------##
  ## create output
  ##-------------------------------------------------------------------##

  contributions = sweep(loadings, 2, colSums(abs(loadings)), "/")
  rownames(contributions) = namx
  colnames(contributions) = blocks_names
  rownames(loadings) = namx
  colnames(loadings) = blocks_names
  names(ind_list) = blocks_names
  names(loadingslist) = blocks_names
  Rx = crossprod(X)
  r_ee = eigen(Rx)
  PCscores = X %*% r_ee$vec[, 1:ncomps]
  vexpPC = r_ee$val
  totvar = sum(r_ee$val)
  vexpPC = r_ee$val[1: ncomps]/totvar
  vexp = makevexp(scores, X)/totvar
  cvexp = cumsum(vexp)

  out = list(loadings = loadings, contributions = contributions,
             vexp = vexp, vexpPC = vexpPC[1:ncompsdone],
             cvexp = cvexp,
             rcvexp = cvexp/cumsum(vexpPC[1:ncompsdone]),
             scores = scores,
             ncomps = ncomps,
             ind = ind_list, cardinality = sapply(ind_list, length),
             loadingslist = loadingslist, spcaMethod = spcaMethod)

    out$corComp = cor(scores)
    if (rtn_all_spca)
      out$all_spca = spca_list
    out$Call = match.call()

  return(out)
}


makevec = function(vec, n){
  le = length(vec)
  if(le == 0)
    stop("need to pass a vector of lenght > 0")
  if(le < n){
    m = n - le
    v = c(vec, rep(vec[le], m))
  }
  if (le == n)
    v = vec
  return(v)
}
