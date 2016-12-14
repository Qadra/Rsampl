library(Rsampl)

context("Testing SSS")

test_that_sampling_works <- function() {
	n <- 100
	W <- runif(100)

	idxs <- sample_subset(W)

	expect_true(!is.null(idxs))
}

test_that_sampling_with_preprocessing_works <- function() {
	n <- 100
	W <- runif(100)

	str = preprocess_sample_subset(W)

	for (i in 1:100) {
		idxs <- sample_subset(struct=str)
		expect_true(!is.null(idxs))
	}

	rm(str)
}

test_that_we_must_supply_an_argument <- function() {
	expect_error(sample_subset())
}

test_that_sampling_works()
test_that_sampling_with_preprocessing_works()
test_that_we_must_supply_an_argument()
