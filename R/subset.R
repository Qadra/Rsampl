#' Perform subset sampling
#'
#' @param W A vector of weights to sample from. Not needed if preprocessing has
#'          been performed.
#' @param struct The preprocessed structured
#' @param draws Meant to specify the amount of draws you want, currently does
#'              nothing.
#' @return A vector of indices, where each index \code{i} is sampled with
#'         probability \code{W[i]}.
#' @examples
#'
#' # First create some data in (0,1]
#' W <- runif(100)
#'
#' # Sample indices from it
#' idxs <- sample_subset(W)
#'
#' # .. or use preprocessing
#' structure <- preprocess_sample_subset(W)
#' idxs <- sample_subset(struct = structure)
#'
#' @export
sample_subset <- function(W = NULL, struct = NULL, draws = NULL) {
	check_preconditions(W, struct, draws)

	if (!is.null(draws)) {
		draws <- as.integer(draws)
	}

	#' @useDynLib Rsampl r_subset_sample
	idx <- .Call(r_subset_sample, W, struct, draws)

	return(idx);
}

#' Preprocess weights for subset sampling.
#'
#' @param W Input weights where \code{W[i]} should be between 0 and 1.
#'
#' @return A datastructure that allows faster weighted random sampling through
#'         use of the preprocessing data structure.
#' @examples
#'
#' # First create some data in (0,1]
#' W <- runif(100)
#'
#' # Sample with preprocessing
#' structure <- preprocess_sample_subset(W)
#' idxs <- sample_subset(struct = structure)
#'
#' @export
preprocess_sample_subset <- function(W) {
	#' @useDynLib Rsampl r_subset_preprocess
	str <- .Call(r_subset_preprocess, W)

	return(str);
}

check_preconditions <- function(W, struct, draws) {
	if (is.null(W)) {
		if (is.null(struct)) {
			stop("Either W or struct must be given as argument")
		}
	}
}
