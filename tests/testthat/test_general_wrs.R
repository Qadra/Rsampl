library(Rsampl)

context("general tests for wrs")

						   # DEFINE TEST FUNCTIONS #

# {{{ Test Functions without preprocessing

test_invalid_k <- function(method) {
	test_that("k > n is invalid", {
				  S <- 1:100
				  W <- runif(100)
				  W <- W/sum(W)
				  k <- 101

				  expect_error(Rsampl::wrs_sample(S, W, k=k, method=method))
	})
}

test_for_same_length_input <- function(method) {
	test_that("indices and weights must have the same length", {
				  S <- 1:99
				  W <- runif(100)
				  W <- W/sum(W)
				  k <- 10

				  expect_error(Rsampl::wrs_sample(S, W, k=k, method=method))
	})
}

test_real_samples <- function(method) {
	for (n in c(32, 45, 50, 64, 95, 100)) {
		test_that(paste("sampling should work with n =", n, "and k = 10"), {
					  S <- 1:n
					  W <- seq(1/n, n)
					  k <- 10

					  # Set seed to fixed value to make randomness repeatable
					  set.seed(142)
					  res <- Rsampl::wrs_sample(S, W, k=k, method=method)

					  expect_equal(length(res), k)

					  # Simple test
					  for (i in 1:k) {
						  expect_true(res[i] >= 0 && res[i] <= n)
					  }
		})
	}
}

test_uniform_input <- function(method) {
	test_that("should handle uniform input", {
				  S <- 1:100
				  W <- rep(1/100, 100)
				  k <- 10

				  res <- Rsampl::wrs_sample(S, W, k=k, method=method)
				  expect_equal(length(res), k)
		})
}

test_single_input <- function(method) {
	test_that("should handle single input", {
				  S <- c(1)
				  W <- c(1)
				  k <- 1

				  res <- Rsampl::wrs_sample(S, W, k=k, method=method)
				  expect_equal(length(res), k)
		})
}

#}}}

# {{{ Test functions with preprocessing

test_that_preprocessing_and_sampling_works <- function(method)  {
	n <- 100
	I <- 1:n
	W <- runif(n)
	W <- W/sum(W)
	k <- n

	ptr <- wrs_preprocess(I, W, k=k, method=method)

	sample <- wrs_sample(k=k, struct=ptr, method=method)
	sample <- wrs_sample(k=k, struct=ptr, method=method)

	rm(ptr)
}

# }}}

# PERFORM ACTUAL TESTING #

test_that("wrs_sample rejects invalid methods",  {
			  S <- 1:100
			  W <- rep(1/100, 100)
			  k <- 10

			  expect_error(Rsampl::wrs_sample(k=k, method="unknown"))
})

methods <- c('binary', 'rstree')

for (method in methods) {
	test_invalid_k(method)
	test_for_same_length_input(method)
	test_real_samples(method)
	test_uniform_input(method)

	test_that_preprocessing_and_sampling_works(method)
}
