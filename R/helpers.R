# @title Logistic Function
# @param w [\code{numeric}]\cr
#   Parameter vector (only relevant entries).
# @param ind [\code{character} | \code{integer}]\cr
#   Relevant columns of X.
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
# @return Returns probabilities of class 2 (assuming class labels 1 and 2), an nrow(X) x 1 matrix.
#
#' @noRd
logistic = function(w, ind, auxdata) {
  eta = as.matrix(auxdata$X[, ind, drop = FALSE] %*% w)    ## X: N x C, w: C x 1, eta: N x 1, as.matrix necessary if auxdata$X is sparse
#  probs = 1 / (1 + exp(-eta))
  probs = binomial()$linkinv(eta)
  return(probs)
}
# 1/(1 + exp(-eta))
# if (netinput < -45) netoutput = 0;
# else if (netinput > 45) netoutput = 1;
# else netoutput = 1 / (1+exp(-netinput));



# @title Softmax Function
# @param w [\code{numeric}]\cr
#   Parameter vector (only relevant entries).
# @param ind [\code{character} | \code{integer}]\cr
#   Relevant columns of X.
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
# @return Returns probabilities of all groups/classes, an nrow(X) x no. of groups (J or K) matrix.
#
#' @noRd
softmax = function(w, ind, auxdata) {
  eta = as.matrix(auxdata$X[, ind, drop = FALSE] %*% matrix(w, nrow = length(ind)))  ## X: N x C, w: C x J, eta: N x J, as.matrix necessary if auxdata$X is sparse
# print(exp(eta)/rowSums(exp(eta)))
  maxeta = apply(eta, 1, max)
  eta = eta - maxeta
  probs = exp(eta)
  probs = probs/rowSums(probs)
# print(probs)
  return(probs)
}



## Note: In the binary case the group lasso reduces to ordinary lasso, therefore type.multinomial is ignored
# @param u [\code{numeric}]\cr
#   Vector of slack variables (only relevant entries).
# @param ind [\code{integer(1)}]\cr
#   Which penalty parameters (\code{lambda}, \code{alpha}) should be used? (1 <= ind <= 1+J).
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
#' @noRd
penBinary = function(u, ind, auxdata) {
  penL2 = 0.5 * (1 - auxdata$alpha[ind]) * sum(u^2)
  penL1 = auxdata$alpha[ind] * sum(u)
  if (auxdata$penalty == "ungrouped") {
    pen = auxdata$lambda[ind] * (penL2 + penL1)
  } else if (auxdata$penalty == "grouped"){
    pen = auxdata$lambda[ind] * 0.5 * (penL2 + penL1)^2
  }
  return(list(pen = pen, penL1 = penL1, penL2 = penL2))
}



## Note: In the binary case the group lasso reduces to ordinary lasso, therefore type.multinomial is ignored
# @param u [\code{numeric}]\cr
#   Vector of slack variables (only relevant entries).
# @param ind [\code{integer(1)}]\cr
#   Which penalty parameters (\code{lambda}, \code{alpha}) should be used? (1 <= ind <= 1+J).
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
#' @noRd
gradPenBinary = function(u, ind, auxdata) {
  gradPen = auxdata$lambda[ind] * ((1 - auxdata$alpha[ind]) * u + auxdata$alpha[ind])
  if (auxdata$penalty == "grouped"){
    penL2 = 0.5 * (1 - auxdata$alpha[ind]) * sum(u^2)
    penL1 = auxdata$alpha[ind] * sum(u)
    gradPen = gradPen * (penL2 + penL1)
  }
  return(gradPen)
}



# @param u [\code{numeric}]\cr
#   Vector of slack variables (only relevant entries).
# @param ind [\code{integer(1)}]\cr
#   Which penalty parameters should be used? (1 <= ind <= 1+J).
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
# @param size
#   The group size (either number of classes K or mixture components J).
# @param V
#   The number of groups.
#
#' @noRd
penMulti = function(u, ind, auxdata, size, V) { ## FIXME: size and V: one of them can be sufficient: V = length(u)/size
  penL2 = 0.5 * (1 - auxdata$alpha[ind]) * sum(u^2)
  if (auxdata$type.multinomial == "ungrouped") {
    penL1 = auxdata$alpha[ind] * sum(u)
  } else {
    # V = length(u)/size #???
    penL1 = auxdata$alpha[ind] * sum(sapply(seq_len(V), function(v) sqrt(sum(u[seq(v, size * V, V)]^2))))
## FIXME: check this
  }
  pen = auxdata$lambda[ind] * (penL2 + penL1)
  return(list(pen = pen, penL1 = penL1, penL2 = penL2))
}
## FIXME andere penalty option



# @param w [\code{numeric}]\cr
#   Parameter vector.
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
# @return Returns an nrow(X) x 1 matrix with probabilities of group 2.
#' @noRd
gatingBinary = function(w, auxdata) {
## FIXME: pass only relevant entries of w to gatingBinary?
  gating = logistic(w[auxdata$indWGating], ind = auxdata$colsGating, auxdata)        # N x 1 matrix
  # FIXME: 
  # gating = gating / (gating + 1 - gating)
#print(rowSums(gating))
  return(gating)
}



# @param w [\code{numeric}]\cr
#   Parameter vector.
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
# @return Returns an nrow(X) x # groups matrix
#' @noRd
gatingMulti = function(w, auxdata) {
## FIXME: pass only relevant entries of w to gatingBinary?
  gating = softmax(w[auxdata$indWGating], ind = auxdata$colsGating, auxdata)        # N x J matrix
  return(gating)
}



# @param w [\code{numeric}]\cr
#   Parameter vector (only relevant entries???).
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
# @return An N x J matrix with one column per mixture component with probabilities for class 2 (assuming class labels 1 and 2)
#' @noRd
expertsBinary = function(w, auxdata) {
  experts = matrix(0, nrow = auxdata$N, ncol = auxdata$J)                  # J matrices of dimension N x 1
  experts[,1:(auxdata$J)] = sapply(1:(auxdata$J),
    function(z) logistic(w[auxdata$indWExperts[z,]], ind = auxdata$colsExperts, auxdata))
  return(experts)
}



# @param w [\code{numeric}]\cr
#   Parameter vector (only relevant entries???).
# @param auxdata [\code{environment} | \code{list}]\cr
#   Relevant data and parameters.
#
# @return An array with dimensions N, K, J.
#' @noRd
expertsMulti = function(w, auxdata) {
  experts = array(0, dim = c(auxdata$N, auxdata$K, auxdata$J))              # J matrices of dimension N x K
  experts[,,1:auxdata$J] = sapply(1:auxdata$J, function(z) softmax(w[auxdata$indWExperts[z,]], ind = auxdata$colsExperts, auxdata))
  return(experts)
}



## target functions
## CONVENTION: in the binary case we model class 2 (assuming class labels 1, 2) and group 2

## FIXME: Numerik log
## FIXME: check signs of loglik and penalties


#' @noRd
eval_f_BinaryBinary = function(x, auxdata) {

  ## split x into two parts
  w = x[1:(auxdata$lenW)]                          # parameters
  u = x[auxdata$lenW + 1:(auxdata$lenU)]                  # slack variables (intercept excluded)

  ## gating
  gating = gatingBinary(w, auxdata)                    # N x 1 matrix, probabilities for group 2
#print(head(gating))
  gating = cbind(1 - gating, gating)
  ## penalty gating
  penGating = penBinary(u[auxdata$indUGating], ind = 1, auxdata)      # scalar, positive
#print(penGating)

  ## experts
  experts = expertsBinary(w, auxdata)                    # N x 2 matrix, probabilities for y = 2
#print(head(experts))
  ind = (auxdata$y == 1)                          # indices of first class
  experts[ind,] = 1 - experts[ind,]                    # N x J matrix with probabilities of actual y
  ## penalty experts, FIXME: Extra-Funktion?
  penExperts = lapply(1:(auxdata$J),
    function(z) penBinary(u[auxdata$indUExperts[z,]], ind = 1+z, auxdata))  # J scalars, positive
# print(penExperts)  

  ## calculate negative penalized log-likelihood
  loglik = log(rowSums(gating * experts))                  # vector of length N
# print(sum(loglik))
  loglik = - sum(loglik)/(auxdata$N)                    # scalar, negative log-likelihood therefore > 0
# print(loglik)

  ## penalty
  pen = penGating$pen + sum(sapply(penExperts, function(p) p$pen))
  if (auxdata$penalty == "grouped") {
    if (!length(auxdata$indIntercept)) {
      penDiff = - auxdata$lambda[4] * 0.5 * sum((w[auxdata$indWExperts[1,]] - w[auxdata$indWExperts[2,]])^2)
    } else {
      penDiff = - auxdata$lambda[4] * 0.5 * sum((w[-auxdata$indIntercept][auxdata$indUExperts[1,]] - w[-auxdata$indIntercept][auxdata$indUExperts[2,]])^2)
    }
# print(penDiff)
    penExperts$penDiff = penDiff
# print(penExperts)
    pen = pen + penDiff
  }

  res = loglik + pen
  attr(res, "loglik") = loglik
  attr(res, "pen") = pen
  attr(res, "penGating") = penGating
  attr(res, "penExperts") = penExperts
  return(res)  
}
## FIXME: auch einzelne Terme ausgeben?


#' @noRd
eval_f_BinaryMulti = function(x, auxdata) {
  
  ## split x into two parts
  w = x[1:(auxdata$lenW)]                          # parameters
  u = x[auxdata$lenW + 1:(auxdata$lenU)]                  # slack variables (intercept excluded)

  ## gating
  gating = gatingBinary(w, auxdata)                    # N x 1 matrix, probabilities for group 2
  gating = cbind(1 - gating, gating)
  ## penalty gating
  penGating = penBinary(u[auxdata$indUGating], ind = 1, auxdata)      # scalar, positive

  ## experts
  experts = expertsMulti(w, auxdata)                    # N x K x J array
  gr = cbind(rep(1:auxdata$N, auxdata$J), rep(auxdata$y, auxdata$J), rep(1:auxdata$J, each = auxdata$N))    # indices of actual class
  experts = matrix(experts[gr], ncol = auxdata$J)              # N x J matrix with probabilities of actual y
  ## penalty experts, FIXME eigene Funktion
  penExperts = sapply(1:auxdata$J, function(z) {
      penMulti(u[auxdata$indUExperts[z,]], 1+z, auxdata, auxdata$K, auxdata$VExperts) ## FIXME NOW
    }
  )                                    # vector of length J, positive
  # penExperts = sapply(1:auxdata$J, function(z) {
      # penMulti(u[indUExperts[z,]], lambda[1+z], alpha[1+z], type.multinomial, K, VExperts) ## FIXME NOW
    # }
  # )                                    # vector of length J
  
  ## calculate negative penalized log-likelihood
  loglik = log(rowSums(gating * experts))                  # vector of length N
  loglik = - sum(loglik)/auxdata$N                      # scalar, negative log-likelihood therefore > 0

  return(loglik + penGating + sum(penExperts))
}



#' @noRd
eval_f_MultiBinary = function(x, auxdata) {
  
  ## split x into two parts
  w = x[1:auxdata$lenW]                          # parameters
  u = x[auxdata$lenW + 1:auxdata$lenU]                    # slack variables (intercept excluded)

  ## gating
  gating = gatingMulti(w, auxdata)                      # N x J matrix
  ## penalty gating
  penGating = penMulti(u[auxdata$indUGating], 1, auxdata, auxdata$J, auxdata$VGating)  # scalar, FIXME NOW

  ## experts
  experts = expertsBinary(w, auxdata)                    # N x 2 matrix, probabilities for y = 2
  ind = (auxdata$y == 1)                          # indices of first class
  experts[ind,] = 1 - experts[ind,]                    # N x J matrix with probabilities of actual y
  ## penalty experts, FIXME Extra-Funktion?
  penExperts = sapply(1:auxdata$J, function(z) penBinary(u[auxdata$indUExperts[z,]], 1+z, auxdata))  # vector of length J
  
  ## calculate negative penalized log-likelihood
  loglik = log(rowSums(gating * experts))                  # vector of length N
  loglik = - sum(loglik)/auxdata$N                      # scalar, negative log-likelihood therefore > 0

  return(loglik + penGating + sum(penExperts))  
}



#' @noRd
eval_f_MultiMulti = function(x, auxdata) {
  
  ## split x into two parts
  w = x[1:auxdata$lenW]                            # parameters
  u = x[auxdata$lenW + 1:auxdata$lenU]                      # slack variables (intercept excluded)
#print(matrix(w, nrow = length(auxdata$indWGating)))
  ## gating
  gating = gatingMulti(w, auxdata)                        # N x J matrix
#print(head(gating))
  ## penalty gating
  penGating = penMulti(u[auxdata$indUGating], 1, auxdata, auxdata$J, auxdata$VGating)  # scalar
  # penGating = penMulti(u[indUGating], lambda[1], alpha[1], type.multinomial, J, VGating)  # scalar
#print(penGating)

  ## experts
  experts = expertsMulti(w, auxdata)                      # N x K x J array
#print(experts)
  gr = cbind(rep(1:auxdata$N, auxdata$J), rep(auxdata$y, auxdata$J), rep(1:auxdata$J, each = auxdata$N))    # indices of actual class
#print(auxdata$y)
  experts = matrix(experts[gr], ncol = auxdata$J)                # N x J matrix with probabilities of y
#print(experts)
  ## penalty experts
  penExperts = lapply(1:auxdata$J, function(z) penMulti(u[auxdata$indUExperts[z,]], 1+z, auxdata, auxdata$K, auxdata$VExperts))  # vector of length J
  # penExperts = sapply(1:auxdata$J, function(z) penMulti(u[indUExperts[z,]], lambda[1+z], alpha[1+z], type.multinomial, K, VExperts))  # vector of   
#print(penExperts)

  ## calculate negative penalized log-likelihood
  loglik = log(rowSums(gating * experts))                    # vector of length N
  loglik = - sum(loglik)/auxdata$N                        # scalar, negative log-likelihood therefore > 0

  ## penalty
  pen = penGating$pen + sum(sapply(penExperts, function(p) p$pen))
  res = loglik + pen
  attr(res, "loglik") = loglik
  attr(res, "pen") = pen
  attr(res, "penGating") = penGating
  attr(res, "penExperts") = penExperts
  
  return(res)  
}


## gradient

#' @noRd
eval_grad_f_BinaryBinary = function(x, auxdata) {
#eval_grad_f = function(x, auxdata) {
#print(address(`auxdata$X`))
#print(str(auxdata$X))  

  ## split x in two parts
  w = x[1:(auxdata$lenW)]                          # parameters
  u = x[auxdata$lenW + 1:(auxdata$lenU)]                  # later on: absolute value of all pars (intercept excluded)
# print(x)

  ## calculate probabilities
  gating = gatingBinary(w, auxdata)                    # N x 1 matrix, probabilities of group 2
  experts = expertsBinary(w, auxdata)                    # N x 2 matrix, probability of y = 2 for groups 1 and 2
# print(summary(gating))
# print(summary(experts))
  ## same code as in eval_f
  ind = (auxdata$y == 1)                          # indices of class 1
  expertsy = experts
  expertsy[ind,] = 1 - expertsy[ind,]                    # N x 2 matrix with probabilities of actual y for groups 1 and 2  
  ge = cbind(1 - gating, gating) * expertsy
  sge = rowSums(ge)
#print("log sge")
#print(sum(log(sge)))
  gesge = ge/sge

  ## gradient gating
  B = gating - gesge[,2]
  # gradGating = 1/(auxdata$N) * as.vector(t(B) %*% auxdata$X[,auxdata$colsGating])
  # gradGating = 1/(auxdata$N) * as.vector(t(auxdata$X[,auxdata$colsGating]) %*% B)
  gradGating = 1/(auxdata$N) * colSums(B[,1] * auxdata$X[,auxdata$colsGating, drop = FALSE])
  # system.time(replicate(1000,1/(auxdata$N) * as.vector(t(B) %*% auxdata$X[,auxdata$colsGating])))
  # system.time(replicate(1000,1/(auxdata$N) * as.vector(t(auxdata$X[,auxdata$colsGating]) %*% B)))
  # system.time(replicate(1000,1/auxdata$N * colSums(B[,1] * auxdata$X[,auxdata$colsGating])))
  # FIXME: try for large V, small N
#print(gradGating)

  ## gradient experts
  # B = (2 * auxdata$y - 3) * gesge * (1 - expertsy)
# print(head(2 * auxdata$y - 3))
  # B = ifelse(auxdata$y == 2, 1, -1) * gesge * (1 - expertsy)
  B = ifelse(auxdata$y == 2, 1, -1)/sge * cbind(1 - gating, gating) * experts * (1 - experts)
  
  # system.time(replicate(1000,as.vector(t(auxdata$X[,auxdata$colsExperts]) %*% B)))
  # system.time(replicate(1000,as.vector(sapply(1:(auxdata$J), function(z) colSums(B[,z] * auxdata$X[,auxdata$colsExperts])))))
  # FIXME: try for large V, small N
  gradExperts = as.vector(sapply(1:(auxdata$J), function(z) -1/(auxdata$N) * colSums(B[,z] * auxdata$X[,auxdata$colsExperts, drop = FALSE])))
  # gradExperts = -1/(auxdata$N) * as.vector(t(auxdata$X[,auxdata$colsExperts]) %*% B)
#print(gradExperts)
  
  ## gradient penalty
  gradPenGating = gradPenBinary(u[auxdata$indUGating], ind = 1, auxdata)
  gradPenExperts = as.vector(sapply(1:(auxdata$J), function(z) gradPenBinary(u[auxdata$indUExperts[z,]], ind = 1+z, auxdata)))

  if (auxdata$penalty == "grouped") {
    ## FIXME: bug in intercept
#    penDiff = - auxdata$lambda[1] * 0.5 * sum((w[-auxdata$indIntercept][auxdata$indUExperts[1,]] - w[-auxdata$indIntercept][auxdata$indUExperts[2,]])^2)
    gradPenDiff = auxdata$lambda[4] * (w[-auxdata$indIntercept][auxdata$indUExperts[1,]] - w[-auxdata$indIntercept][auxdata$indUExperts[2,]])
    gradPenDiff = c(-gradPenDiff, gradPenDiff)
# print(gradPenDiff)
    gradExperts[-(auxdata$indIntercept[-1]-length(auxdata$indWGating))] = gradExperts[-(auxdata$indIntercept[-1]-length(auxdata$indWGating))] + gradPenDiff
  }
  
  # gradPenGating = auxdata$lambda[1] * ((1 - auxdata$alpha[1]) * u[auxdata$indUGating] + auxdata$alpha[1])
  # gradPenExperts = as.vector(sapply(1:(auxdata$J), function(z) auxdata$lambda[1+z] * ((1 - auxdata$alpha[1+z]) * u[auxdata$indUExperts[z,]] + auxdata$alpha[1+z])))

  return(c(gradGating, gradExperts, gradPenGating, gradPenExperts))  
}



# @noRd
# eval_h_BinaryBinary = function(x, auxdata) {
  # ## split x in two parts
  # w = x[1:(auxdata$lenW)]                          # parameters
  # u = x[auxdata$lenW + 1:(auxdata$lenU)]                  # later on: absolute value of all pars (intercept excluded)

  # ## calculate probabilities
  # gating = gatingBinary(w, auxdata)                    # N x 1 matrix, probabilities of group 2
  # experts = expertsBinary(w, auxdata)                    # N x 2 matrix, probability of y = 2 for groups 1 and 2
  
  # ind = (auxdata$y == 1)                          # indices of class 1
  # expertsy = experts
  # expertsy[ind,] = 1 - expertsy[ind,]                    # N x 2 matrix with probabilities of actual y for groups 1 and 2  
  # ge = cbind(1 - gating, gating) * expertsy
  # sge = rowSums(ge)
  # gesge = ge/sge
  
  # ## gating
  # B = gating * (1 - gating) - sge[,1] * sge[,2]
  # hGatingGating = 1/(auxdata$N) * as.vector(t(auxdata$X[,auxdata$colsGating]) %*% diag(B) %*% auxdata$X[,auxdata$colsGating])
  # # colSums(B[,1] * auxdata$X[,auxdata$colsGating, drop = FALSE] * FIXME)

  # ## experts
  
  # ## mixed

  # ## penalty gating
  # hPenGatingGating = 1 - auxdata$alpha[1]
  
  
  
  
# }



#' @noRd
eval_grad_f_BinaryMulti = function(x, auxdata) {
  
  ## split x in two parts
  w = x[1:auxdata$lenW]                        # parameters
  u = x[auxdata$lenW + 1:auxdata$lenU]                # later on: absolute value of all pars (intercept excluded)

  ## calculate probabilities
  gating = gatingBinary(w, auxdata)                  # N x 1 matrix
  experts = expertsMulti(w, auxdata)                  # N x J x K array

  ## same code as in eval_f
  gr = cbind(rep(1:auxdata$N, auxdata$J), rep(auxdata$y, auxdata$J), rep(1:auxdata$J, each = auxdata$N))
  expertsy = matrix(experts[gr], ncol = auxdata$J)          # N x J matrix with probabilities of y
  ge = cbind(1 - gating, gating) * expertsy              # N x J matrix
  sge = rowSums(ge)
  gesge = ge/sge                            # N x J matrix, with J = 2

  ## gradient gating
  B = gating - gesge[,2]
  gradGating = 1/(auxdata$N) * colSums(B[,1] * auxdata$X[,auxdata$colsGating])

  ## gradient experts
  experts[gr] = experts[gr] - 1
  B = matrix(sapply(1:auxdata$J, function(z) gesge[,z] * experts[,,z]), ncol = auxdata$J*auxdata$K)
  gradExperts = as.vector(sapply(1:(auxdata$J*auxdata$K), function(z) 1/auxdata$N*colSums(B[,z] * auxdata$X[,auxdata$colsExperts])))

  ## gradient penalty  
  gradPenGating = auxdata$lambda[1] * ((1 - auxdata$alpha[1]) * u[auxdata$indUGating] + auxdata$alpha[1])
  gradPenExperts = as.vector(sapply(1:auxdata$J, function(z) auxdata$lambda[1+z] * ((1 - auxdata$alpha[1+z]) * u[auxdata$indUExperts[z,]] + auxdata$alpha[1+z])))
  
  return(c(gradGating, gradExperts, gradPenGating, gradPenExperts))  
}



#' @noRd
eval_grad_f_MultiBinary = function(x, auxdata) {
  
  ## split x in two parts
  w = x[1:auxdata$lenW]                        # parameters
  u = x[auxdata$lenW + 1:auxdata$lenU]                # later on: absolute value of all pars (intercept excluded)

  ## calculate probabilities
  gating = gatingMulti(w, auxdata)                  # N x J matrix
  experts = expertsBinary(w, auxdata)                  # N x J matrix

  ## same code as in eval_f
  ind = (auxdata$y == 1)                        # indices of class 0
  expertsy = experts
  expertsy[ind,] = 1 - expertsy[ind,]                  # N x J matrix with probabilities of y  
  ge = gating * expertsy
  sge = rowSums(ge)
  gesge = ge/sge

  ## gradient gating
  B = gating - gesge
  gradGating = as.vector(sapply(1:auxdata$J, function(z) 1/auxdata$N * colSums(B[,z] * auxdata$X[,auxdata$colsGating])))

  ## gradient experts
  B = -(2 * (auxdata$y-1) - 1) * gesge * (1 - expertsy)
  gradExperts = as.vector(sapply(1:auxdata$J, function(z) 1/auxdata$N * colSums(B[,z] * auxdata$X[,auxdata$colsExperts])))
  
  ## gradient penalty  
  gradPenGating = auxdata$lambda[1] * ((1 - auxdata$alpha[1]) * u[auxdata$indUGating] + auxdata$alpha[1])
  gradPenExperts = as.vector(sapply(1:auxdata$J, function(z) auxdata$lambda[1+z] * ((1 - auxdata$alpha[1+z]) * u[auxdata$indUExperts[z,]] + auxdata$alpha[1+z])))
  
  return(c(gradGating, gradExperts, gradPenGating, gradPenExperts))  
}



#' @noRd
eval_grad_f_MultiMulti = function(x, auxdata) {
  
  ## split x in two parts
  w = x[1:auxdata$lenW]                        # parameters
  u = x[auxdata$lenW + 1:auxdata$lenU]                # later on: absolute value of all pars (intercept excluded)
  
  ## calculate probabilities
  gating = gatingMulti(w, auxdata)                  # N x J matrix
  experts = expertsMulti(w, auxdata)

  ## same code as in eval_f
  gr = cbind(rep(1:auxdata$N, auxdata$J), rep(auxdata$y, auxdata$J), rep(1:auxdata$J, each = auxdata$N))
  expertsy = matrix(experts[gr], ncol = auxdata$J)          # N x J matrix with probabilities of y
  ge = gating * expertsy                        # N x J matrix
  sge = rowSums(ge)
  gesge = ge/sge                            # N x J matrix

  ## gradient gating
  B = gating - gesge
  gradGating = as.vector(sapply(1:auxdata$J, function(z) 1/auxdata$N * colSums(B[,z] * auxdata$X[,auxdata$colsGating])))

  ## gradient experts
  experts[gr] = experts[gr] - 1
  B = matrix(sapply(1:auxdata$J, function(z) gesge[,z] * experts[,,z]), ncol = auxdata$J*auxdata$K)
  gradExperts = as.vector(sapply(1:(auxdata$J*auxdata$K), function(z) 1/auxdata$N * colSums(B[,z] * auxdata$X[,auxdata$colsExperts])))

  ## gradient of penalty (q = 1)
  ## FIXME: type.multinomial = "grouped"
  gradPenGating = auxdata$lambda[1] * ((1 - auxdata$alpha[1]) * u[auxdata$indUGating] + auxdata$alpha[1])
  gradPenExperts = as.vector(sapply(1:auxdata$J, function(z) auxdata$lambda[1+z] * ((1 - auxdata$alpha[1+z]) * u[auxdata$indUExperts[z,]] + auxdata$alpha[1+z])))

  return(c(gradGating, gradExperts, gradPenGating, gradPenExperts))
}


# eval_grad_f <- function(x) {
  # # separate x in two parts
  # w <- x[1:mw]                            # parameters
  # u <- x[mw + 1:mu]                          # later on: absolute value of all pars except intercept
  # ## gradient gating
  # gating <- softmax(X[,colsGating], w[1:(J*(VGating+interceptGating))])    # N x J matrix
  # experts <- array(0, dim = c(N, K, J))                # J matrices of dimension N x K
  # experts[,,1:J] <- sapply(1:J, function(z) {
      # softmax(X[,colsExperts],
        # w[J*(VGating+interceptGating) + ((z-1)*K*(VExperts+interceptExperts) + 1):(z*K*(VExperts+interceptExperts))])
    # }
  # )
  # gr <- cbind(rep(1:N, J), rep(y, J), rep(1:J, each = N))
  # expertsy <- matrix(experts[gr], ncol = J)              # N x J matrix with probabilities of y
  # ge <- gating * expertsy            # N x J matrix
  # sge <- rowSums(ge)
  # gesge <- ge/sge                # N x J matrix
  # B <- gating - gesge
  # gradGating <- as.vector(sapply(1:J, function(z) 1/N*colSums(B[,z] * X[,colsGating])))
  # ## gradient experts
  # experts[gr] <- experts[gr] - 1
  # B <- matrix(sapply(1:J, function(z) gesge[,z] * experts[,,z]), ncol = J*K)
  # gradExperts <- as.vector(sapply(1:(J*K), function(z) 1/N*colSums(B[,z] * X[,colsExperts])))
  # ## gradient of penalty (q = 1)
  # gradPen <- lambda * ((1 - alpha) * u + alpha)
  # return(c(gradGating, gradExperts, gradPen))
# }



## second derivative






## sparsity structure Hessian





# non-linear constraints
#' @noRd
# eval_g = function(x) {
  # ## split x into two parts
  # w = x[1:lenW]            # parameters
  # u = x[lenW + 1:lenU]          # later on: absolute value of all pars (intercept excluded)
  # if (lenW > lenU) {
    # w = w[-indIntercept]      # remove intercept pars from w
  # }
  # return(c(u + w, u - w))
# }
# eval_g = function(x, auxdata) {
  # ## split x into two parts
  # w = x[1:auxdata$lenW]            # parameters
  # u = x[auxdata$lenW + 1:auxdata$lenU]    # later on: absolute value of all pars (intercept excluded)
  # if (auxdata$lenW > auxdata$lenU) {
    # w = w[-c(auxdata$indIntercept)]      # remove intercept pars from w
  # }
  # return(c(u + w, u - w))
# }
eval_g = function(x, auxdata) {
  ## split x into two parts
  w = x[1:auxdata$lenW]            # parameters
  u = x[auxdata$lenW + 1:auxdata$lenU]      # later on: absolute value of all pars (intercept excluded)
  if (auxdata$lenW > auxdata$lenU) {
    w = w[-c(auxdata$indIntercept)]      # remove intercept pars from w
  }
  # c return(c(u + w, u - w, sum(u[auxdata$indUGating])))  ## FIXME: change constraints
  return(c(u + w, u - w))
}

## 




# Jacobian of g
# dim: (lenU + lenU) x (lenW + lenU)
#
# d (u+w) | d (u+w)
# d w     | d u
# ----------------- 
# d (u-w) | d (u-w)
# d w     | d u
# -----------------
# d(sum(u)) | d(sum(u))
# d w       | d u
#
#
# X_mu,mw  | I_mu
# ---------------
# -X_mu,mw | I_mu
# ---------------
# 0        | I_mu
#
# X_mu,mw is I_mu interspersed with 0 columns

## Jacobian of g:
#  J = [  I  I
#        -I  I 
#         0  I ],
# where I is and identity matrix of size m

#' @noRd
# eval_jac_g = function(x) {
    # # return a vector of 1 and minus 1, since those are the values of the non-zero elements
    # return(c(rep(1, 2 * lenU), rep(c(-1, 1), lenU)))
# }

eval_jac_g = function(x, auxdata) {
    # return a vector of 1 and minus 1, since those are the values of the non-zero elements
    return(c(rep(1, 2 * auxdata$lenU), rep(c(-1, 1), auxdata$lenU)))
}
# c eval_jac_g = function(x, auxdata) {
# c    # return a vector of 1 and minus 1, since those are the values of the non-zero elements
# c    return(c(rep(1, 2 * auxdata$lenU), rep(c(-1, 1), auxdata$lenU), rep(1, length(auxdata$indUGating))))
# c}



# ## sparsity structure
# eval_jac_g_structure = lapply(rep(1:lenU, 2), function(x) {
  # if (length(indIntercept)) {
    # c((1:lenW)[-indIntercept][x], lenW + x)
  # } else {
    # c((1:lenW)[x], lenW + x)
  # }  
   # # return(c((1:lenW)[-indIntercept][x], lenW + x))
# })



# # The constraint functions are bounded from below by zero.
# constraint_lb = rep(0.0, 2 * lenU)
# constraint_ub = rep(Inf, 2 * lenU)
