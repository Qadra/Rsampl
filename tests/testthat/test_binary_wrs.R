library(Rsampl)

context("binary weighted random search")

test_that("k > n is invalid", {
			  S <- 1:100
			  W <- runif(100)
			  W <- W/sum(W)
			  k <- 101

			  expect_error(Rsampl::sample(S, W, k=k))
})

test_that("indices and weights must have the same length", {
			  S <- 1:99
			  W <- runif(100)
			  W <- W/sum(W)
			  k <- 10

			  expect_error(Rsampl::sample(S, W, k=k))
})

for (n in c(32, 45, 50, 64, 95, 100)) {
	test_that(paste("sampling should work with n =", n, "and k = 10"), {
				  S <- 1:n
				  W <- seq(1/n, n)
				  k <- 10

				  # Set seed to fixed value to make randomness repeatable
				  set.seed(142)
				  res <- Rsampl::sample(S, W, k=k)

				  expect_equal(length(res), k)

				  # Simple test
				  for (i in 1:k) {
					  expect_true(res[i] >= 0 && res[i] <= n)
				  }
	})
}
