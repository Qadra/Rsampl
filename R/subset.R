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
