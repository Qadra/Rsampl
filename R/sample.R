sample <- function(A, k = 1) {
	n <- length(A)
	return(.Call(sample_, A, n, k))
}
