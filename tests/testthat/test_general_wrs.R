library(Rsampl)

context("general tests for wrs")

						   # DEFINE TEST FUNCTIONS #

# Load functions
source("generic_wrs_test_functions.R")

# PERFORM ACTUAL TESTING #

test_that("wrs_sample rejects invalid methods",  {
			  S <- 1:100
			  W <- rep(1/100, 100)
			  k <- 10

			  expect_error(Rsampl::wrs_sample(S, W, k=k, method="unknown"))
})

# Old method takes 12 seconds for 500k samples

test_that("distributions should be similar", {
			  set.seed(124124)
			  n <- 20
			  n_samples <- 30000
			  A <- 1:n
			  W <- runif(n)
			  W <- W/sum(W)
			  k <- 10


			  #ptm <- proc.time()

			  bin_struct <- wrs_preprocess(A, W, k, method="binary")
			  rstree_struct <- wrs_preprocess(A, W, k, method="rstree")

			  bin    <- vector(mode='integer', length=n)
			  rstree <- vector(mode='integer', length=n)

			  sample <- wrs_sample(A, W, k=k, struct=bin_struct, draws=n_samples)
			  sample <- wrs_sample(A, W, k=k, struct=rstree_struct, draws=n_samples)

			  #print(proc.time() - ptm)


			  sum_of_squared_errors <- sqrt(sum((bin/n_samples - rstree/n_samples)^2))

			  #print(sum_of_squared_errors)

			  #expect_equal(0, sum_of_squared_errors)
})
