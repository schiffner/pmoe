#' @noRd
#'
#' @export

coef.pmoe = function(object, ...) {
	object$coefs
}



#' @title Extract Interactions
#'
#' @description Extract interactions.
#'
#' @param object An object of class \code{"pmoe"}.
#' @param ... Currently unused.
#'
#' @export

extractInteractions = function(object, ...) {
	if (object$K == 2) {
		if (object$J == 2) {
			## remove intercept
			d = object$coefs.experts[,,2] - object$coefs.experts[,,1]
			inter = outer(d, object$coefs.gating)			
			# inter = sapply(1:ncol(coefs$coefs.experts), function(z) outer(d[,z], cbind(1 - coefs$coefs.gating, coefs$coefs.gating)[,z]))
			# rn = rownames(coefsExperts$Comp.1[[1]])
			# lrn = length(rn)
			# rownames(inter) = paste(rep(rn, lrn), rep(rn, each = lrn), sep = ":")
		} else if (object$J > 2) {
			# pairwise differences
			stop("not yet")
		}
	} else if (object$K > 2) {
		stop("not yet")
	} else
		stop("something is wrong with the number of classes")
	return(list(coefs = object$coefs, inter = inter))
}
