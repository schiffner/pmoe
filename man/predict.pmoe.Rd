% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmoe.R
\name{predict.pmoe}
\alias{predict.pmoe}
\title{Predict Penalized Mixtures of Experts}
\usage{
\method{predict}{pmoe}(object, newdata, ...)
}
\arguments{
\item{object}{An object of class \code{"pmoe"}.}

\item{newdata}{New data to be predicted.
A \code{data.frame} of cases to be classified or, if \code{object} has a
       \code{formula}, a \code{data.frame} with columns of the same names as the
       variables used.  A vector will be interpreted as a row
       vector.}

\item{\dots}{Further arguments. Currently unused.}
}
\value{
A \code{list} with components:
  \item{class}{The predicted class labels (a \code{factor}).}
  \item{posterior}{Matrix of class posterior probabilities.}
  \item{gating}{Matrix of gating probabilities.}
  \item{experts}{Matrix of class posterior probabilities from individual experts.}
}
\description{
Predict penalized mixtures of experts.
}

