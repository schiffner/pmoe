#' @title Penalized Mixtures of Experts
#'
#' @description
#' Train a penalized mixture of experts model by IPOPT.
#'
#' @param formula A \code{formula} with two right-hand sides of the form \code{y ~ expert model | gating model}.
#'  The normal formula syntax applies (see \code{\link[Formula]{Formula}}).
#'  In order to use all available columns one or both right-hand sides can be dots (\code{y ~ . | .}, \code{y ~ . | gating model}).
#'  Both right-hand sides normally contain an intercept (which is not penalized) that can be removed by adding \code{0} or \code{- 1}.
#'  Offsets are allowed, too.
#  In the binary case numeric vector as long as the number of training observations.
#  In the multi-class case matrix or vector? FIXME
#'  The left-hand side should be a vector indicating the class membership, preferably a \code{factor}.
#' @param data A \code{data.frame} from which variables specified in \code{formula} are to be taken.
#' @param X (Required if no \code{formula} is given as principal argument.) A \code{matrix} or \code{data.frame} or \code{Matrix} containing
#'   the explanatory variables (must not contain an intercept column).
#' @param y (Required if no \code{formula} is given as principal argument.) A \code{factor} specifying
#'   the class membership for each observation.
#' @param colsGating Names or indices of columns in \code{X} to be used for the gating model. Default is all columns.
#' @param colsExperts Names or indices of columns in \code{X} to be used for the expert models. Default is all columns.
#' @param interceptGating Logical. Does the gating model include an intercept? If \code{TRUE}, an intercept column is added to \code{X}.
#'   Defaults to \code{TRUE}.
#' @param interceptExperts Logical. Does the expert model include an intercept? If \code{TRUE}, an intercept column is added to \code{X}.
#'   Defaults to \code{TRUE}.
#' @param offsetGating Offset term for the gating model.
#   A numeric vector. FIXME: matrix?
#' @param offsetExperts Offset term for the expert model.
#   A numeric vector. FIXME: matrix?
#' @param J The number of experts / mixture components. Defaults to 2.
#' @param lambda Penalty parameter. Can be a scalar or a vector of length \code{1+J} with different components for the gating and the
#'   \code{J} expert models. All components must be >= 0.
#' @param alpha Mixing parameter for the elastic net penalty. Can be a scalar or a vector of length \code{1+J} with different components
#'   for the gating and the \code{J} expert models. All components must be in [0,1].
#'   Defaults to 1.
#' @param type.multinomial \code{"grouped"}, \code{"ungrouped"}.
#' @param standardize Logical. Should the columns of \code{X} be standardized prior to fitting the model?
#'   Defaults to \code{FALSE}. If \code{TRUE} the coefficients are returned on the original scale.
#   when necessary? return on original scale? macht das sinn fuer spase X?
#' @param genetic Logical.
#' @param ipopt.max.iter The maximum number of IPOPT iterations.
#' @param ipopt.tol Tolerance for IPOPT convergence.
#' @param subset A subset...
#' @param na.action ...
#' @param penalty ...
#' @param model ...
#' @param \dots Further arguments.
#'
#' @return An object of class \code{pmoe}.
#  A \code{list} with entries:
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#   \item{}{.}
#'
#' @references
#' Waechter, A. and Biegler, L. T. (2006),
#' On the Implementation of an Interior-Point Filter Line-Search Algorithm for Large-Scale Nonlinear Programming,
#' \emph{Mathematical Programming}, \strong{106}, 25-57.
#'
#' @rdname pmoe
#'
#' @export

pmoe = function(X, ...) {
  UseMethod("pmoe")
}



#' @rdname pmoe
#'
#' @export

pmoe.formula = function(formula, data, subset, na.action, ...) {
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf = mf[c(1,m)]
  f = Formula(formula)
  mf[[1]] = as.name("model.frame")
  mf$formula = f
  mf = eval(mf, parent.frame())
  y = model.response(mf)
  ## offsets are not in model.matrix
  offsetGating = model.offset(model.part(f, data = mf, rhs = 2, terms = TRUE))   ## either a numerical vector/matrix or NULL
  offsetExperts = model.offset(model.part(f, data = mf, rhs = 1, terms = TRUE))  ## not included in model.matrix
  ## intercepts
  TermsGating = terms(f, data = data, rhs = 2)
  TermsExperts = terms(f, data = data, rhs = 1)
  interceptGating = as.logical(attr(TermsGating, "intercept"))
  interceptExperts = as.logical(attr(TermsExperts, "intercept"))
  ## build model matrix (w/o intercept column)
  Terms = terms(f, data = data, rhs = NULL)
  attr(Terms, "intercept") = 0
#  l = vector(attr(Terms, "terms.labels"), )
  X = model.matrix(Terms, mf, rhs = NULL) ## contrast
  colsGating = attr(TermsGating, "term.labels")
  colsExperts = attr(TermsExperts, "term.labels")
  res = pmoe.default(X, y, colsGating, colsExperts, interceptGating, interceptExperts, offsetGating, offsetExperts, ...)
  res$terms = Terms  # FIXME: all what is needed to reconstruct model.frame from data
  res$contrasts = attr(X, "contrasts")
  res$xlevels = .getXlevels(Terms, mf)
  res
}


#' @rdname pmoe
#'
#' @export

pmoe.default = function(X, y, colsGating = 1:ncol(X), colsExperts = 1:ncol(X), interceptGating = TRUE, interceptExperts = TRUE,
  offsetGating = NULL, offsetExperts = NULL, J = 2, lambda, alpha = 1, penalty = c("ungrouped", "grouped"), type.multinomial = c("ungrouped", "grouped"), model = c("binomial", "multinomial"), standardize = FALSE, genetic = FALSE, ipopt.max.iter = 500, ipopt.tol = 1e-6, ...) {

  ## checks
  model = match.arg(model)
  type.multinomial = match.arg(type.multinomial)
  penalty = match.arg(penalty)
  # if (K == 2 && type.multinomial == "grouped") {
    # warning("")
    # type.multinomial = "ungrouped"  
  # }
  # if (J == 2)
  ## FIXME: man braucht type.multinomial eigentlich fuer gating und experts separat, vektor?
  
  ## check X
  if (is.null(dim(X))) 
    stop("'X' is not a matrix")
  if (!is.matrix(X) & !inherits(X, "Matrix"))
    X = as.matrix(X)
  ## FIXME: Matrix
  # if (any(!is.finite(X)))
    # stop("infinite, NA or NaN values in 'X'")
  
  ## check that X is numeric
  # FIXME if (!is.numeric(X@x))
  # if (mode(X) != "numeric") funktioniert nicht fuer sparse
  #  stop("'X' is not 'numeric'")
  ## FIXME: remove not needed columns (neither in colsGating nor colsExperts from X)
  ## this might change colsGating, colsExperts
  ## scan for constant columns  
  ## check if intercept column is there after all
  ## FIXME: Matrix kann zu groß werden
  const = (apply(X, 2, sd) < 1e-08) ## FIXME: < tolerance parameter
  if (any(const)) {
    warning(sprintf(ngettext(sum(const), "Column %s appears to be constant", 
      "columns %s appear to be constant"), paste(colnames(X)[const], collapse = ", ")), domain = NA)
  }
  # FIXME: entfernen? in scrime nachschauen?
  ## standardize numeric columns in X (except intercept)
  cen = sc = NULL
  if (standardize) {
    X = scale(X)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
  }
  ## FIXME: backscaling, scaling for sparse X

  ## check colsGating, colsExperts
  V = ncol(X)
  VGating = length(colsGating)
  VExperts = length(colsExperts)
  if (is.character(colsGating)) {
    ## correct names
    ## coerce to numeric?
    if (interceptGating)
      colsGating = c("(Intercept)", colsGating)
  } else if (is.numeric(colsGating)) {
    ## 
    if (any(colsGating > V))
      stop("'colsGating' contains undefined columns")
    if (interceptGating)
      colsGating = c(1, colsGating + 1)
    else if (interceptExperts)
      colsGating = colsGating + 1
  } else 
    stop("'colsGating' is neither a character nor a numeric vector")  
  if (is.character(colsExperts)) {
    ## correct names
    ## coerce to numeric?
    if (interceptExperts)
      colsExperts = c("(Intercept)", colsExperts)
  } else if (is.numeric(colsExperts)) {
    ## 
    if (any(colsExperts > V))
      stop("'colsExperts' contains undefined columns")
    if (interceptExperts)
      colsExperts = c(1, colsExperts + 1)
    else if (interceptGating)
      colsExperts = colsExperts + 1
  } else 
    stop("'colsExperts' is neither a character nor a numeric vector")
  
  ## X and intercept column
  if (interceptGating || interceptExperts)
    X = cbind(`(Intercept)` = 1, X) ## geht das schlauer?
  ## Important: for glmnet initialization no intercept column is needed
  
  ## check for intercept column, add if it is not there?
  ## FIXME: requirements on X and intercepts
  # should X contain intercept column, what indices are in colsGating, colsExperts? 
  
  ## check offset
  # FIXME: binomial, multinomial, matrix
  N = nrow(X)
  # if (N != length(weights)) 
    # stop("nrow(X) and length(weights) are different")
  # if (any(weights < 0)) 
    # stop("weights have to be larger or equal to zero")
  # if (all(weights == 0)) 
    # stop("all weights are zero")
  # names(weights) = rownames(x)

  ## check y
  if (N != length(y))
    stop("'nrow(X)' and 'length(y)' are different")
  if (!is.factor(y)) 
    warning("'y' was coerced to a factor")
  y = as.factor(y)
  lev = lev1 = levels(y)
  counts = as.vector(table(y))
  if (any(counts == 0)) {
    empty = lev[counts == 0]
    warning(sprintf(ngettext(length(empty), "group %s is empty", 
      "groups %s are empty"), paste(empty, collapse = ", ")), 
      domain = NA)
    lev1 = lev[counts > 0]
    y = factor(y, levels = lev1)
    counts = as.vector(table(y))
  }
  if (length(lev1) < 2L) 
    stop("training data from only one group given")
  names(counts) = lev1
  prior = counts/N
  ## check for groups with small numbers of training observations
  if (any(counts == 1))
    stop("One or more class has only 1 training observation")
  # if (any(counts < 8))
    # warning("One or more class has small number of training observations (< 8)")  ## FIXME: welche toleranz ist sinnvoll? glmnet?
  K = length(lev1)

  ## J
  if (J < 2)
    stop("'J' must be at least 2")

  ## lambda
  if (any(lambda < 0))
    stop("'lambda' must be >= 0")
  if (length(lambda) == 1) {
    lambda = rep(lambda, 2 + J)
  } else if (length(lambda) != 2 + J) {
    stop("'length(lambda)' is not 2 + J")
  }

  ## alpha
  if (any(alpha < 0 | alpha > 1))
    stop("'alpha' must be in [0,1]")
  if (length(alpha) == 1) {
    alpha = rep(alpha, 1 + J)
  } else if (length(alpha) != 1 + J) {
    stop("'length(alpha)' is not 1 + J")
  }
  
  ## set up the model and choose correct ipopt helper functions
  # VGating = length(colsGating) - interceptGating    # w/o intercept
  # VExperts = length(colsExperts) - interceptExperts    # w/o intercept
  ## FIXME: clean up, make this more practical
  indGating = indExperts = integer(0)      ## indices of intercept coefficients
  if (J == 2 && model == "binomial") {
    if (K == 2) {
      lenW = (VGating + interceptGating) + 2 * (VExperts + interceptExperts)
      lenU = VGating + 2 * VExperts
      indWGating = 1:(VGating + interceptGating)
      indWExperts = matrix((VGating + interceptGating) + 1:(2 * (VExperts + interceptExperts)), nrow = 2, byrow = TRUE)  ## 2 x (VExperts + interceptExperts) matrix
      indUGating = 1:VGating
      indUExperts = matrix(VGating + 1:(2 * VExperts), nrow = 2, byrow = TRUE)  ## 2 x VExperts matrix
      if (interceptGating)
        indGating = 1
      if (interceptExperts)
        indExperts = seq((VGating + interceptGating) + 1, lenW, VExperts+1) ## indWExperts[,1]
        ## FIXME: 1. Spalte der Matrix indWExperts nehmen?
# print(indGating)
# print(indExperts)
      eval_f = eval_f_BinaryBinary
      eval_grad_f = eval_grad_f_BinaryBinary
    }
    else if (K > 2) {
      lenW = (VGating + interceptGating) + 2 * K * (VExperts + interceptExperts)
      lenU = VGating + 2 * K * VExperts
      indWGating = 1:(VGating + interceptGating)
      indWExperts = matrix((VGating + interceptGating) + 1:(2 * K * (VExperts + interceptExperts)), nrow = 2, byrow = TRUE)  ## 2 x K*(VExperts + interceptExperts) matrix
      indUGating = 1:VGating
      indUExperts = matrix(VGating + 1:(2 * K * VExperts), nrow = 2, byrow = TRUE)  ## 2 x K*VExperts matrix
      if (interceptGating)
        indGating = 1
      if (interceptExperts)
        indExperts = seq((VGating + interceptGating) + 1, lenW, VExperts+1)
      eval_f = eval_f_BinaryMulti
      eval_grad_f = eval_grad_f_BinaryMulti
    }
    else {
      stop("incorrect specification of 'K'")
    }
  }
  else if (J > 2 || model == "multinomial") {
    if (K == 2 && model == "binomial") {
      lenW = J * (VGating + interceptGating) + J * (VExperts + interceptExperts)
      lenU = J * VGating + J * VExperts
      indWGating = 1:(J * (VGating + interceptGating))
      indWExperts = matrix(J * (VGating + interceptGating) + 1:(J * (VExperts + interceptExperts)), nrow = J, byrow = TRUE)  ## J x (VExperts + interceptExperts) matrix
      indUGating = 1:(J * VGating)
      indUExperts = matrix(J * VGating + 1:(J * VExperts), nrow = J, byrow = TRUE)    ## J x VExperts matrix    
      if (interceptGating)
        indGating = seq(1, J * (VGating + interceptGating), VGating + 1)
      if (interceptExperts)
        indExperts = seq(J * (VGating + interceptGating) + 1, lenW, VExperts+1)
      eval_f = eval_f_MultiBinary
      eval_grad_f = eval_grad_f_MultiBinary
    }
    else if (K > 2 || model == "multinomial") {
      lenW = J * (VGating + interceptGating) + J * K * (VExperts + interceptExperts)
      lenU = J * VGating + J * K * VExperts
      indWGating = 1:(J * (VGating + interceptGating))
      indWExperts = matrix(J * (VGating + interceptGating) + 1:(J * K * (VExperts + interceptExperts)), nrow = J, byrow = TRUE)  ## J x K*(VExperts + interceptExperts) matrix
      indUGating = 1:(J * VGating)
      indUExperts = matrix(J * VGating + 1:(J * K * VExperts), nrow = J, byrow = TRUE)    ## J x K*VExperts matrix
      if (interceptGating)
        indGating = seq(1, J * (VGating + 1), VGating + 1)
      if (interceptExperts)
        indExperts = seq(J * (VGating + interceptGating) + 1, lenW, VExperts+1)
      eval_f = eval_f_MultiMulti
      eval_grad_f = eval_grad_f_MultiMulti
    }
    else {
      stop("incorrect specification of 'K'")
    }
  } else {
    stop("incorrect specification of 'J'")
  }
  indIntercept = c(indGating, indExperts)

  # f = function(x, auxdata) 100*eval_f(x, auxdata)
  # grad_f = function(x, auxdata) 100*eval_grad_f(x, auxdata)

  ## initialization by an EM step
  # # initPmoe = function(X, y, J, K, alpha, lambda, lenW, lenU, colsGating, interceptGating, offsetGating, indWGating, colsExperts, interceptExperts, offsetExperts, indWExperts, indIntercept, standardize) {
    # ## y should be a factor
    # ## FIXME: important: X without intercept
    # x0 = numeric(lenW + lenU)  ## sparse?
    # cluster = kmeans(x = X[,colsGating][,-1], centers = J)$cluster  ## fuer welche Dimensionen geht das noch? Clustern mit Dim-Reduktion? FIXME: ordinal data FIXME: OHNE INTERCEPT
    # ## FIXME: kmedoids/other method for discrete data???
    # cluster = diag(J)[cluster,]    # J column matrix, FIXME: softness/overlap parameter statt 0.9? an hme orientieren, 
    # cluster[cluster == 1] = 0.9
    # cluster[cluster == 0] = (1 - 0.9)/(J - 1)
    # ## calculate models for lambda sequence, extract coefficients and fitted values for the user-defined lambda
    # ## even if there is no intercept in the model, there is an intercept row in the coef matrix that needs to be excluded
    # if (J == 2) {
      # gating = glmnet::glmnet(as.matrix(X[,colsGating][,-1]), cluster, family = "binomial", alpha = alpha[1], intercept = interceptGating, offset = offsetGating, standardize = standardize)
# ## FIXME: max und min lambda merken
      # if (interceptGating)
        # x0[indWGating] = glmnet::coef.glmnet(gating, s = lambda[1])[,1]    ## sparse lassen? FIXME
      # else
        # x0[indWGating] = glmnet::coef.glmnet(gating, s = lambda[1])[-1,1]  ## sparse lassen? FIXME      
    # } else {
      # gating = glmnet::glmnet(X[,colsGating][,-1], cluster, family = "multinomial", alpha = alpha[1], intercept = interceptGating, offset = offsetGating, standardize = standardize)
      # #, type.multinomial = type.multinomial[1])
      # ## FIXME: glmnet docu says offset has to be a N x K matrix for multinomial
      # if (interceptGating)
        # x0[indWGating] = as.vector(sapply(1:J, function(j) glmnet::coef.glmnet(gating, s = lambda[1])[[j]][, 1])) ## FIXME
      # else
        # x0[indWGating] = as.vector(sapply(1:J, function(j) glmnet::coef.glmnet(gating, s = lambda[1])[[j]][-1, 1])) ## FIXME      
    # }
    # # weights = getS3method("predict", class(gating))(gating, newx = X[,colsGating][,-1], s = lambda[1], type = "response", offset = offsetGating)  ## probabilities
    # weights = predict(gating, newx = X[,colsGating][,-1], s = lambda[1], type = "response", offset = offsetGating)  ## probabilities
# print(summary(weights))    
    # if (J == 2) {
      # weights = cbind(1 - weights, weights)
    # } else {
      # weights = weights[,,1]
    # }
    # if (K == 2) {
      # for (j in 1:J) {
        # expert = glmnet::glmnet(as.matrix(X[,colsExperts][,-1]), y, weights = weights[,j], family = "binomial", alpha = alpha[1+j], intercept = interceptExperts,
          # offset = offsetExperts, standardize = standardize)
        # if (interceptExperts)
          # x0[indWExperts[j,]] = glmnet::coef.glmnet(expert, s = lambda[1+j])[,1] ## FIXME
        # else
          # x0[indWExperts[j,]] = glmnet::coef.glmnet(expert, s = lambda[1+j])[-1,1] ## FIXME        
      # }
    # } else {
      # for (j in 1:J) {
        # expert = glmnet::glmnet(X[,colsExperts][,-1], y, weights = weights[,j], family = "multinomial", alpha = alpha[1+j], intercept = interceptExperts,
          # offset = offsetExperts, standardize = standardize) #, type.multinomial = type.multinomial[1+j])
        # if (interceptExperts)
          # x0[indWExperts[j,]] = as.vector(sapply(1:K, function(k) glmnet::coef.glmnet(expert, s = lambda[1+j])[[k]][, 1]))
        # else
          # x0[indWExperts[j,]] = as.vector(sapply(1:K, function(k) glmnet::coef.glmnet(expert, s = lambda[1+j])[[k]][-1, 1]))          
      # }
    # }
    # ## initialize u
    # if (length(indIntercept)) {
      # x0[lenW + 1:lenU] = abs(x0[1:lenW][-indIntercept]) + 0.2
    # } else {
      # x0[lenW + 1:lenU] = abs(x0[1:lenW]) + 0.2
    # }
    # # # return(x0)
  # # # }
  
  ## sparsity structure
  # eval_jac_g_structure = lapply(rep(1:lenU, 2), function(x) {
    # if (length(indIntercept)) {
      # c((1:lenW)[-indIntercept][x], lenW + x)
    # } else {
      # c((1:lenW)[x], lenW + x)
    # }  
  # })
  # if (length(indIntercept)) {
    # eval_jac_g_structure = lapply(rep(1:lenU, 2), function(x) c((1:lenW)[-indIntercept][x], lenW + x))
  # } else {
    # eval_jac_g_structure = lapply(rep(1:lenU, 2), function(x) c((1:lenW)[x], lenW + x))
  # }
  ## habe ich nicht schon Vektoren, die so aussehen?

  eval_jac_g_structure = rep(1:lenU, 2) + lenW
  if (length(indIntercept)) {
    m1 = rep((1:lenW)[-indIntercept], 2)
  } else {
    m1 = rep((1:lenW), 2)
  }
  eval_jac_g_structure = cbind(m1, eval_jac_g_structure, deparse.level = 0) 
  eval_jac_g_structure = lapply(1:(2*lenU), function(i) eval_jac_g_structure[i,])
  # eval_jac_g_structure2 = lapply(1:lenU, function(x) {
    # lenW + x
  # })
  # c eval_jac_g_structure[[2*lenU + 1]] = lenW + indUGating

  # The constraint functions are bounded from below by zero
  constraint_lb = rep(0.0, 2 * lenU)
  constraint_ub = rep(Inf, 2 * lenU)
  # c constraint_lb = c(rep(0.0, 2 * lenU), tol)
  # c constraint_ub = rep(Inf, 2 * lenU + 1)

  y = unclass(y)
  ## FIXME: check if entries are 1:K

  ## auxiliary data
  # auxdata = new.env()
  auxdata = list()
  auxdata$X = X
  auxdata$y = y
  auxdata$lenW = lenW
  auxdata$lenU = lenU
  auxdata$colsGating = colsGating
  auxdata$colsExperts = colsExperts
  auxdata$indWGating = indWGating
  auxdata$indWExperts = indWExperts
  auxdata$indUGating = indUGating
  auxdata$indUExperts = indUExperts
  auxdata$indIntercept = indIntercept
  auxdata$lambda = lambda
  auxdata$alpha = alpha
  auxdata$N = N
  auxdata$J = J
  auxdata$K = K
  auxdata$type.multinomial = type.multinomial
  auxdata$penalty = penalty
  # auxdata$tol = tol
  ## rm(list = names(auxdata))

  # # x0 = initPmoe(X, y, K, J, alpha, lambda, lenW, lenU, colsGating, interceptGating, offsetGating, indWGating, colsExperts, interceptExperts, offsetExperts, indWExperts, indIntercept, standardize)
  if (genetic) {
#    gen = GenSA(par = NULL, fn = eval_f, lower = rep(c(-5,0),c(lenW,lenU)), upper = rep(5, lenW + lenU), control = list(max.call = max.call, verbose = TRUE), auxdata = auxdata)
    gen = GenSA(par = NULL, fn = function(x, ...) eval_f(c(x, abs(x[-indIntercept])), ...), lower = rep(-5,lenW), upper = rep(5, lenW), control = list(max.time = 60, verbose = TRUE), auxdata = auxdata)
# print(gen)
    x0 = c(gen$par, abs(gen$par[-indIntercept]) + 0.1)
  } else {
    # x0Experts = rnorm(length(indWExperts[1,]), sd = 1)
    # x0 = c(rnorm(length(indWGating), sd = 1), x0Experts, -x0Experts, rep(1,lenU))
    x0 = c(rnorm(lenW, sd = 1), rep(1,lenU))
    # x0 = rep(0:1, c(lenW, lenU))
    # x0[4] = 1
    # x0[indWExperts[1,]][2] = 1
    # x0[indWExperts[2,]][2] = -1
  }

  ## ipopt options
  ipoptr_opts = list(
    "jac_d_constant" = 'yes', # Indicates whether all inequality constraints are linear
    "hessian_approximation" = 'limited-memory',
    #"derivative_test" = 'first-order',
    "obj_scaling_factor" = 100,
    # "point_perturbation_radius" = 0.1,
    # "derivative_test_print_all" = 'yes',
    #"hessian_approximation_space" = 'all-variables',
    #"hessian_constant" = 'yes',  ## FIXME: BFGS?
    "bound_relax_factor" = 0,
    "mu_strategy" = 'adaptive',
    "max_iter" = ipopt.max.iter,
    "tol" = ipopt.tol#,
    #acceptable_tol = 1e-10
  )


# print(address(`auxdata$X`))
# library(numDeriv)
# print("numDeriv")
# print(grad(eval_f, x0, auxdata = auxdata))
# print("formula")
# print(eval_grad_f(x0, auxdata))

  ## call ipopt
  res = ipoptr(
    x0 = x0, 
    eval_f = eval_f, 
    eval_grad_f = eval_grad_f, 
    eval_g = eval_g, 
    eval_jac_g = eval_jac_g,
    eval_jac_g_structure = eval_jac_g_structure,
    constraint_lb = constraint_lb, 
    constraint_ub = constraint_ub,
    opts = ipoptr_opts,
    auxdata = auxdata
  )
  
  ## FIXME: names x0
  ## FIXME: active 
  
  # res$loglik = eval_f(res$solution, auxdata)
  # res$penGating = pen*()
  # res$penExperts = pen*()

  res$lev = lev
  res$lev1 = lev1
  res$prior = prior
  res$counts = counts
  res$model = model
  # VGating
  # VExperts
  
  ## extract coefficients  
  res$coefs = res$solution[1:lenW]
  ## für path besser alles untereinanderschreiben und nur vernünftig benenen
  ## schwellwert wählen
  ## coefs sparse machen
  names(res$coefs) = c(
    paste("gating", if (J == 2 && model == "binomial") 2 else rep(1:J, each = (VGating + interceptGating)), colnames(X)[colsGating], sep = "."), 
    paste("expert", if (K == 2 && model == "binomial") rep(1:J, each = (VExperts + interceptExperts)) else rep(1:J, each = K * (VExperts+interceptExperts)), if (K == 2 && model == "binomial") lev1[2] else rep(lev1, each = (VExperts+interceptExperts)), colnames(X)[colsExperts], sep = ".")
  )

  if (J == 2 && model == "binomial") {
    ## length(colsGating) x 1 matrix
    res$coefs.gating = as.matrix(res$coefs[indWGating])
    rownames(res$coefs.gating) = colnames(X)[colsGating]
    colnames(res$coefs.gating) = "2"
  } else if (J > 2 || model == "multinomial") {
    ## length(colsGating) x J matrix
    res$coefs.gating = matrix(res$coefs[indWGating], ncol = J)
    rownames(res$coefs.gating) = colnames(X)[colsGating]
    colnames(res$coefs.gating) = paste("gate", 1:J, sep = ".")
  }
  if (K == 2 && model == "binomial") {
    ## length(colsExperts) x 1 x J array
    res$coefs.experts = array(res$coefs[t(indWExperts)], dim = c(VExperts+interceptExperts, 1, J),
      dimnames = list(colnames(X)[colsExperts], lev1[2], 1:J))
  } else if (K > 2 || model == "multinomial") {
    ## length(colsExperts) x K x J array
    res$coefs.experts = array(res$coefs[t(indWExperts)], dim = c(VExperts+interceptExperts, K, J),
      dimnames = list(colnames(X)[colsExperts], lev1, 1:J))
  }
  ## coefs sparse machen
  
  ## log-likelihood and penalties
  attr = eval_f(res$solution, auxdata)
  res$loglik = attr(attr, "loglik")
  res$pen = attr(attr, "pen")
  res$penGating = attr(attr, "penGating")
  res$penExperts = attr(attr, "penExperts")
  res$cen = cen
  res$sc = sc
  
  res = c(res, auxdata[-c(1:2)])
  ## loglik
  ## penalties zurückgeben
  
  ## offset  
  ## check if some experts are the same and remove them?
  ## extract results, i.e. parameters, log-likelihood, penalties
  
  ## sparse matrices gating, experts
  ## FIXME: path?
  # g = Matrix(res$solution[indWGating])
  # e = Matrix(res$solution[indWExperts], ncol = J)

  class(res) = "pmoe" ## "binary", "multi"
  res
}



#' @title Predict Penalized Mixtures of Experts
#'
#' @description Predict penalized mixtures of experts.
#'
# @details ...
#'
#' @param object An object of class \code{"pmoe"}.
#' @param newdata New data to be predicted.
#'   A \code{data.frame} of cases to be classified or, if \code{object} has a
#'          \code{formula}, a \code{data.frame} with columns of the same names as the
#'          variables used.  A vector will be interpreted as a row
#'          vector.
# If \code{newdata} is missing, an attempt will be made to
#         retrieve the data used to fit the \code{wlda} object.
#' @param \dots Further arguments. Currently unused.
#'
#' @return A \code{list} with components:
#'   \item{class}{The predicted class labels (a \code{factor}).}
#'   \item{posterior}{Matrix of class posterior probabilities.}
#'   \item{gating}{Matrix of gating probabilities.}
#'   \item{experts}{Matrix of class posterior probabilities from individual experts.}
#' 
#' @rdname predict.pmoe
#'
#' @export
#'
# @examples
#
#

## FIXME: also return linear predictors for gating and experts?
## FIXME: retrieve data if newdata is missing?
predict.pmoe = function(object, newdata, ...) {
  if (!inherits(object, "pmoe"))
    stop("'object' not of class 'pmoe'")
  ## newdata -> X  FIXME
  if (!is.null(Terms <- object$terms)) {
    Terms = delete.response(Terms)
    newdata = model.frame(Terms, newdata, na.action = na.pass, xlev = object$xlevels) ##?? save xlevels
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, newdata)
    X = model.matrix(Terms, newdata, contrasts = object$contrasts, rhs = NULL) ## save contrast???
  } else {
    if (is.null(dim(newdata)))  # a vector is treated as one row
      dim(newdata) = c(1, length(newdata))
      X = as.matrix(newdata)
    }
  if (!is.null(object$cen))
    X = scale(X, object$cen, object$sc)
  if (length(object$indIntercept))
    X = cbind(`(Intercept)` = 1, X)
  # if (ncol(X) != length(union(object$colsGating, object$colsExperts)))
        # stop("wrong number of variables")
    # if (length(colnames(x)) > 0L && any(colnames(x) != dimnames(object$means)[[2L]]))
        # warning("variable names in 'newdata' do not match those in 'object'")
  object$X = X
  object$N = nrow(X)
  # w = coef(object)
  w = object$solution[1:object$lenW]  ## FIXME: better extraction method
  ## calculate class posterior probabilities
  ## FIXME: gating and experts functions are defined differently, environmant, X?
  if (object$J == 2 && object$model == "binomial") {
    if (object$K == 2) {
      gating = gatingBinary(w, object)                # N x 1 matrix (group 2)
      experts = expertsBinary(w, object)              # N x 2 matrix (class 2)
      post = rowSums(cbind(1 - gating, gating) * experts)    # N x 1 matrix (class 2)
      post = cbind(1 - post, post)                    # N x 2 matrix
    } else if (object$K > 2) {
      gating = gatingBinary(w, object)                # N x 1 matrix (group 2)
      experts = expertsMulti(w, object)               # N x K x 2 array
      post = (1 - gating[,1]) * experts[,,1] + gating[,1] * experts[,,2]  # N x K matrix , FIXME: make gating a vector?
    } else {
      stop("something is wrong with 'object$K'")
    }
  } else if (object$J > 2 || object$model == "multinomial") {
    if (object$K == 2 && object$model == "binomial") {
      gating = gatingMulti(w, object)            # N x J matrix
      experts = expertsBinary(w, object)         # N x J matrix (class 2)
      post = rowSums(gating * experts)           # N x 1 matrix (class 2)
      post = cbind(1 - post, post)               # N x 2
    } else if (object$K > 2 || object$model == "multinomial") {
      gating = gatingMulti(w, object)            # N x J matrix
      experts = expertsMulti(w, object)          # N x K x J array
      post = matrix(rowSums(sapply(1:object$J, function(z) gating[,z] * experts[,,z])), ncol = object$K)  # N x K matrix 
    } else {
      stop("something is wrong with 'object$K'")
    }
  } else {
    stop("something is wrong with 'object$J'")
  }
  rownames(post) = rownames(X)
  colnames(post) = object$lev1
  cl = factor(object$lev1[max.col(post)], levels = object$lev)
  return(list(class = cl, posterior = post, gating = gating, experts = experts))
}
