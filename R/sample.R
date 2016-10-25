#' @export
sample <- function(A, W, k = 1) {
	if (length(A) != length(W)) {
		stop("input A and W must be of same length")
	}

	#' @useDynLib Rsampl sample_impl
	idx <- .Call(sample_impl, W, k)
	return(idx)
}
