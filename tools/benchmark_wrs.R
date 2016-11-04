#sample_int_R(n, size, prob)
#sample_int_ccrank(n, size, prob)
#sample_int_crank(n, size, prob)
#sample_int_expj(n, size, prob)
#sample_int_expjs(n, size, prob)
#sample_int_rank(n, size, prob)
#sample_int_rej(n, size, prob)

library(ggplot2)
library(Rsampl)
library(wrswoR)
library(microbenchmark)

benchmark <- function(n = 1000, k = 250, times=1000) {
	B <- rexp(n)
	B <- B/sum(B)

	bin_pre <- wrs_preprocess(1:n, B, k=k, method="binary")
	rstree_pre <- wrs_preprocess(1:n, B, k=k, method="rstree")

	bm <- microbenchmark(times=times,
						 BIN=Rsampl::wrs_sample(1:n, B, k, method="binary"),
						 BINP=Rsampl::wrs_sample(k=k, struct=bin_pre),
						 RSTREE=Rsampl::wrs_sample(1:n, B, k, method="rstree"),
						 RSTREEP=Rsampl::wrs_sample(k=k, struct=rstree_pre),
						 CCRANK=sample_int_ccrank(n, k, B),
						 CRANK=sample_int_crank(n, k, B),
						 EXPJ=sample_int_expj(n, k, B),
						 #EXPJS=sample_int_expjs(n, k, B),
						 #RANK=sample_int_rank(n, k, B),
						 REJ=sample_int_rej(n, k, B)
						 #R=sample_int_R(n, k, B)
						 )

	return(bm)
}
