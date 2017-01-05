library(Rsampl)

context("Testing Subset Sampling")

test_that_sampling_works <- function() {
	test_that("Subsetsampling works", {
		n <- 100
		W <- runif(100)

		idxs <- sample_subset(W)

		expect_true(!is.null(idxs))
	})
}

test_that_sampling_with_preprocessing_works <- function() {
	test_that("Subset sampling with preprocessing works", {
		n <- 100
		W <- runif(100)

		str = preprocess_sample_subset(W)

		for (i in 1:100) {
			idxs <- sample_subset(struct=str)
			expect_true(!is.null(idxs))
		}

		rm(str)
	})
}

test_that_we_must_supply_an_argument <- function() {
	test_that("Subset sampling warns when element is missing", {
		expect_error(sample_subset())
	})
}

test_that_multiple_samples_work <- function() {
	test_that("Multiple draws with subset sampling works", {
		n <- 100
		W <- runif(100)

		str = preprocess_sample_subset(W)

		idxs <- sample_subset(struct=str, draws=100)
		expect_true(!is.null(idxs))

		expect_equal(length(idxs), 100)

		rm(str)
	})
}

test_that_sampling_works()
test_that_sampling_with_preprocessing_works()
test_that_we_must_supply_an_argument()
test_that_multiple_samples_work()
