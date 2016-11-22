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
						 BINP=Rsampl::wrs_sample(1:n, B, k=k, struct=bin_pre),
						 RSTREE=Rsampl::wrs_sample(1:n, B, k, method="rstree"),
						 RSTREEP=Rsampl::wrs_sample(1:n, B, k=k, struct=rstree_pre),
						 SRMETHOD=Rsampl::wrs_sample(1:n, B, k=k, method="srmethod"),
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

do_sample_fast_bin <- function(A, W, k, s) {
	bin_pre <- wrs_preprocess(A, W, k=250, method='binary')

	samples = wrs_sample(A=A, W=W, k=k, struct=bin_pre, draws=s)
}

do_sample_slow_bin <- function(A, W, k, s) {
	bin_pre <- wrs_preprocess(A, W, k=250, method='binary')

	samples = matrix(nrow=k, ncol=0)
	for (i in 1:s) {
		sample = wrs_sample(A=A, W=W, k=k, struct=bin_pre)
		samples <- cbind(samples, sample)
	}
}

benchmark_multisample <- function(n = 1000, k = 250, s = 1000, times = 100) {
	W <- rexp(n)
	W <- W/sum(W)
	A <- 1:n

	bm <- microbenchmark(times = times,
						 FAST=do_sample_fast_bin(A, W, k, s),
						 SLOW=do_sample_slow_bin(A, W, k, s)
						 )

	return(bm)
}

benchmark_normal_distribution <- function(n = 10000, times=100, sd=250, limit = 1e9) {
	library(data.table)

	D <- dnorm(1:n, sd=sd, mean=as.integer(n/2))
	#D <- runif(1:n); D <- D/sum(D)

	I <- 1:n
	k <- 10

	out <- data.table()

	rstree <- wrs_preprocess(I, D, k, method='rstree')
	srmethod <- wrs_preprocess(I, D, k, method='srmethod')
	bin <- wrs_preprocess(I, D, k, method='binary')

	repeat {
		bm <- microbenchmark(times=times,
							 RSTree=  wrs_sample(I,D,k, struct=rstree,   method='rstree'),
							 SRMethod=wrs_sample(I,D,k, struct=srmethod, method='srmethod'),
							 Binary=  wrs_sample(I,D,k, struct=bin,      method='binary'),
							 unit='ns',
							 control=list(warmup=8)
							 )

		# Calculate the average of the amount of samples
		dt <- data.table(bm)[, lapply(.SD, mean), by=(expr)]
		# Add a column with k
		dt <- dt[, k := k]

		out <- rbind(out, dt)

		if (dt[which(dt$expr == 'SRMethod')]$time > limit) {
			break
		}

		if (k > 1000) {
			break
		}

		k <- k + 10

	}

	names(out) <- c('Method', 'ns', 'k')

	return(out)
}
