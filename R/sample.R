check_preconditions_for_wrs <- function(A, W, k = 1, method = 'binary', struct) {
	if (!missing(A)) {
		if (length(A) != length(W)) {
			stop("input A and W must be of same length")
		}
	}

	if (missing(struct)) {
		if (missing(A) | missing(W)) {
			error("missing either struct or both A and W")
		}
	}

	#methods <- c('binary', 'rstree')

	#if (!method %in% methods) {
	#	stop(paste("Parameter method must be in: (", paste(methods, collapse=", "), ")", sep=""))
	#}
}

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
		idx <- apply(idx, 1:2, function(x,i ) {return(i[x+1])}, A)
	} else {
		idx <- A[idx]
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
	check_preconditions_for_wrs(A, W, k, method, struct)

	#' @useDynLib Rsampl r_wrs_preprocess
	str <- .Call(r_wrs_preprocess, W, k, method)

	return(str)
}
