#' Perform Weighted Random Sampling
#' 
#' @param A A vector/list of indices to sample
#' @param W The weight of each corresponding index in A
#' @param k The number of indices to sample
#' @return \code{k} indices sampled according to weights
#' @examples
#' sample(1:100, rep(1/100, 100), k = 10)
#' @export
sample <- function(A, W, k = 1) {
	if (length(A) != length(W)) {
		stop("input A and W must be of same length")
	}

	#' @useDynLib Rsampl sample_impl
	idx <- .Call(sample_impl, W, k)
	return(idx)
}
