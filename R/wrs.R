#' Perform Weighted Random Sampling
#' 
#' @param A A vector/list of indices to sample. This value is not needed if
#'        struct has been initialized and passed.
#' @param W The weight of each corresponding index in A. This value is not
#'        needed if struct has been initialized and passed.
#' @param k The number of indices to sample.
#' @param method Must be either 'binary' or 'rstree'. This determines the
#'               sampling method used.
#' @param struct a preprocessed chunk of data.
#' @param draws The number of samples to draw. By default this method returns
#'              a vector of indices. If draws is specified, the function
#'              returns a matrix with dimensions (k, draws).
#' @return \code{k} indices sampled according to weights
#' @examples
#' wrs_sample(1:100, rep(1/100, 100), k = 10, method='binary')
#'
#' wrs_sample(1:100, rep(1/100, 100), k = 10, method='rstree')
#'
#' # To use preprocessing do the following: 
#' struct <- wrs_preprocess(1:100, rep(1/100, 100), k = 10)
#' samples <- wrs_sample(k=10, struct=struct)
#'
#' # Do many draws in a single sampling operation, same as running a for loop
#' # and sampling 5 times, but faster. Returns a matrix.
#' struct <- wrs_preprocess(1:100, rep(1/100, 100), k = 10)
#' samples <- wrs_sample(k=10, struct=struct, draws=5)
#'
#' @export
wrs_sample <- function(A = NULL, W = NULL, k = 1, method = 'binary', struct, draws=NULL) {
	check_preconditions_for_wrs(A, W, k, method, struct)

	if (missing(struct)) {
		struct = NULL
	}

	if (!is.null(draws)) {
		draws <- as.integer(draws)
	}

	#' @useDynLib Rsampl r_wrs_sample
	idx <- .Call(r_wrs_sample, W, k, method, struct, draws)

	# Select the indexes and remove zero indexing
	if (!is.null(draws)) {
		# This seems to be the fastest method.
		dims = dim(idx)
	}

	if (!is.null(struct)) {
		# Use the indexes we have saved from the preprocessing.
		idx <- attributes(struct)$A[idx + 1]
	} else {
		idx <- A[idx + 1]
	}

	if (!is.null(draws)) {
		dim(idx) <- dims
	}

	return(idx)
}

#' Perform Weighted Random Sampling Preprocessing
#' 
#' @param A A vector/list of indices to sample
#' @param W The weight of each corresponding index in A
#' @param k The number of indices to sample
#' @param method Must be either 'binary' or 'rstree'. This determines the
#'               sampling method used.
#' @return \code{S} Prepared structure for sampling.
#' @export
wrs_preprocess <- function(A, W, k = 1, method = 'binary') {
	check_preconditions_for_wrs(A, W, k, method)

	#' @useDynLib Rsampl r_wrs_preprocess
	str <- .Call(r_wrs_preprocess, W, k, method)

	attributes(str) <- list(A = A)

	return(str)
}

check_preconditions_for_wrs <- function(A, W, k = 1, method = 'binary', struct = NULL) {
	if (!missing(A)) {
		if (length(A) != length(W)) {
			stop("input A and W must be of same length")
		}
	}

	if (missing(struct)) {
		if (missing(A) | missing(W)) {
			stop("missing either struct or both A and W")
		}
	}
}

#' @export
rstree_draw <- function(filename, struct) {
	#' @useDynLib Rsampl print_rstree
	.Call(print_rstree, "out.dot", struct)
}
