# Folding test of unimodality
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

#' Computes the folding ratio of the input data
#'
#' @param X nxd matrix (n observations, d dimensions)
#'
#' @return the folding ratio
#'
#' @examples
#' X = matrix(runif(n = 1000, min = 0., max = 1.), ncol = 1)
#' phi = folding.statistics(X)
#'
#' @export
folding.ratio <- function(X) {
	n = dim(X)[1]   # number of observations
	Sigma = cov(X)  # covariance matrix
	D = sum(diag(Sigma)) # trace of the covariance matrix
	X_squared_norm = matrix(rowSums(X**2),ncol = 1) # |X|²
	C = cov(X,X_squared_norm) # cov(X,|X|²)
	s2star = 0.5*solve(Sigma,C)
	S = t(matrix(rep(s2star,n),ncol=n)) # repetition of s*
	Ynorm = sqrt(matrix(rowSums((X-S)**2),ncol = 1)) # |X-s*|
	return(as.numeric(var(Ynorm)/D))
}


#' Computes the pivot \eqn{s_2} (approximate pivot)
#'
#' @param X nxd matrix (n observations, d dimensions)
#'
#' @return the approximate pivot
#'
#' @examples
#' library(MASS)
#' mu = c(0,0)
#' Sigma = matrix(c(1,0.5,1,0.5), ncol = 2)
#' X = mvrnorm(n = 5000, mu = mu, Sigma = Sigma)
#' Phi = pivot.approx(X)
#'
#' @export
pivot.approx <- function(X) {
	n = dim(X)[1]   # number of observations
	Sigma = cov(X)  # covariance matrix
	D = sum(diag(Sigma)) # trace of the covariance matrix
	X_squared_norm = matrix(rowSums(X**2),ncol = 1) # |X|²
	C = cov(X,X_squared_norm) # cov(X,|X|²)
	return(as.numeric(0.5*solve(Sigma,C)))
}


#' Computes the folding statistics of the input data
#'
#' @param X nxd matrix (n observations, d dimensions)
#'
#' @return the folding statsistics
#'
#' @examples
#' library(MASS)
#' mu = c(0,0)
#' Sigma = matrix(c(1,0.5,1,0.5), ncol = 2)
#' X = mvrnorm(n = 5000, mu = mu, Sigma = Sigma)
#' Phi = folding.statistics(X)
#'
#' @export
folding.statistics <- function(X) {
	d = dim(X)[2] # dimension
	return(folding.ratio(X)*(d+1)**2)
}


#' Computes the confidence bound for the significance level p
#'
#' @param n sample size
#' @param d dimension
#' @param p significance level (between 0 and 1, the lower, the more significant)
#'
#' @return the confidence bound q (the bounds are 1-q and 1+q)
#'
#' @examples
#' n = 2000  	# number of observations
#' d = 2 	  	# 2 dimensional data
#' p = 0.05  	# we want the bound at the level 0.05 (classical p-value)
#' q = folding.test.bound(n,d,p)
#'
#' @export
folding.test.bound <- function(n,d,p) {
	a = 0.4785
	b = 0.1946
	c = 2.0287
	return( (a*(1-p) - b*log(p)) * (log(d)+c) / sqrt(n) )
}


#' Computes the p-value of the folding test
#'
#' @param Phi the folding statistics
#' @param n sample size
#' @param d dimension
#'
#' @return the p-value (the lower, the more significant)
#'
#' @examples
#' library(MASS)
#' n = 5000
#' d = 2
#' mu = c(0,0)
#' Sigma = matrix(c(1,0.5,1,0.5), ncol = d)
#' X = mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' Phi = folding.statistics(X)
#' p = folding.test.pvalue(Phi,n,d)
#'
#' @export
folding.test.pvalue <- function(Phi,n,d) {
	a = 0.4785
	b = 0.1946
	c = 2.0287
	e = abs(Phi-1)*sqrt(n)/(log(d)+c)
	f = function(x) {a*(1-x) - b*log(x) - e}
	return (uniroot(f, lower = 0.0, upper = 1.0)$root)
}


#' Perform the folding test of unimodality
#'
#' @param X $nxd$ matrix (n observations, d dimensions)
#'
#' @return 1 if unimodal, 0 if multimodal
#'
#' @examples
#' library(MASS)
#' n = 10000
#' d = 3
#' mu = c(0,0,0)
#' Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1), ncol = d)
#' X = mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' m = folding.test(X)
#'
#' @export
folding.test <- function(X) {
	n = nrow(X)
	d = ncol(X)
	Phi = folding.statistics(X)
	p = folding.test.pvalue(Phi,n,d)
	if (Phi>1) {
		y <- list(unimodal = TRUE, pvalue = p, Phi = Phi)
		return(y)
	}
	else {
		y <- list(unimodal = FALSE, pvalue = p, Phi = Phi)
		return(y)
	}
}
