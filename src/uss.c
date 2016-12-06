#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "defs.h"

#include "sample.h"
#include "uss.h"
#include "xutil.h"

#include <inttypes.h>

typedef union {
	uint64_t i;
	double   f;
} T;

static double pow_lookup[64] = {
	1.00000000000000000000e+00, 5.00000000000000000000e-01,
	2.50000000000000000000e-01, 1.25000000000000000000e-01,
	6.25000000000000000000e-02, 3.12500000000000000000e-02,
	1.56250000000000000000e-02, 7.81250000000000000000e-03,
	3.90625000000000000000e-03, 1.95312500000000000000e-03,
	9.76562500000000000000e-04, 4.88281250000000000000e-04,
	2.44140625000000000000e-04, 1.22070312500000000000e-04,
	6.10351562500000000000e-05, 3.05175781250000000000e-05,
	1.52587890625000000000e-05, 7.62939453125000000000e-06,
	3.81469726562500000000e-06, 1.90734863281250000000e-06,
	9.53674316406250000000e-07, 4.76837158203125000000e-07,
	2.38418579101562500000e-07, 1.19209289550781250000e-07,
	5.96046447753906250000e-08, 2.98023223876953125000e-08,
	1.49011611938476562500e-08, 7.45058059692382812500e-09,
	3.72529029846191406250e-09, 1.86264514923095703125e-09,
	9.31322574615478515625e-10, 4.65661287307739257812e-10,
	2.32830643653869628906e-10, 1.16415321826934814453e-10,
	5.82076609134674072266e-11, 2.91038304567337036133e-11,
	1.45519152283668518066e-11, 7.27595761418342590332e-12,
	3.63797880709171295166e-12, 1.81898940354585647583e-12,
	9.09494701772928237915e-13, 4.54747350886464118958e-13,
	2.27373675443232059479e-13, 1.13686837721616029739e-13,
	5.68434188608080148697e-14, 2.84217094304040074348e-14,
	1.42108547152020037174e-14, 7.10542735760100185871e-15,
	3.55271367880050092936e-15, 1.77635683940025046468e-15,
	8.88178419700125232339e-16, 4.44089209850062616169e-16,
	2.22044604925031308085e-16, 1.11022302462515654042e-16,
	5.55111512312578270212e-17, 2.77555756156289135106e-17,
	1.38777878078144567553e-17, 6.93889390390722837765e-18,
	3.46944695195361418882e-18, 1.73472347597680709441e-18,
	8.67361737988403547206e-19, 4.33680868994201773603e-19,
	2.16840434497100886801e-19, 1.08420217248550443401e-19,
};

static int64_t mask = ~0xFFFFFFFFFFFFF,
		mask2 = 0x7FF;
		//invmask2 = ~0x7FF;

static inline int bucket_index(double pi, int64_t L) {
	int64_t exp;
	T rep = {
		.f = pi
	};
	
	//int preexp;
	//frexp(pi, &preexp);

	exp = (rep.i & mask) >> 52;
	//int sign = (exp & invmask2) >> 11;

	exp &= mask2;
	exp -= 1022;

	if (exp == 1) {
		exp = 0;
	}

	exp = -exp;

	if (exp > L) {
		exp = L;
	}

	//printf("%d == %"PRId64", %f\n", preexp, exp, pi);

	return exp;
}

/* TODO:
 * - Finding the bucket index requires another round down, as 2^{-k-1}
 *   should be rounded down to 2^k.
 */
uss_t *uss_preprocess(double *weights, int n) {
	int i, idx, offset,
		L          = xceil_log2(n), // ceil(log2(n));
	   *index      = XALLOC(n, int),
	   *Bcnt       = XALLOC(L + 1, int),
	   *prefix_sum = XALLOC(L + 1, int);

	double mu      = 0.,
		  *w       = weights,
	      *p_hat   = XALLOC(n, double);

	uss_t *u       = XALLOC(1, uss_t);

	// Counts elements in bucket to calculate prefix sum, then used for
	// storing indices into each bucket.
	memset(Bcnt,  0, sizeof(int)*(L + 1));

	for (i = 0; i < n; i++) {
		mu += w[i];

		// TODO: We need to round the k^{-k-1} elements down another notch.
		idx = bucket_index(w[i], L);
		assert(idx >= 0 && idx <= L);
		Bcnt[idx]++;
	}

	prefix_sum[0] = 0;
	for (i = 1; i <= L; i++) {
		// TODO SIMD instructions?
		prefix_sum[i] = prefix_sum[i-1] + Bcnt[i-1];
	}

	memset(Bcnt,  0, sizeof(int)*(L + 1));

	for (i = 0; i < n; i++) {
		idx    = bucket_index(w[i], L);
		assert(idx >= 0 && idx <= L);
		offset = prefix_sum[idx] + Bcnt[idx];

		assert(offset < n);

		//p_hat[offset] = pow(2, -idx);
		p_hat[offset] = pow_lookup[idx];
		index[offset] = i;
		Bcnt[idx]++;
	}

	u->index = index;
	u->p_hat = p_hat;
	u->w     = w;
	u->n     = n;
	u->mu    = 2 * mu + 1;
	u->samples = XALLOC(u->mu, int);

	/* Pass on the next method to use in the sampling */
	// TODO: should this be USS?
	subset_method_t *fp    = XALLOC(1, subset_method_t);
	fp->preprocess = (subset_preprocess_fp_t) uss_basecase_preprocess;
	fp->sample     = (subset_sample_fp_t)     uss_basecase_sample;
	fp->free       = (subset_free_fp_t)       uss_basecase_free;

	/*
	sample_t sam = {
		.internal = fp,
		.n        = n,
		.weights  = p_hat
	};
	*/

	/* Initialize sorted subset sampling */
	u->s = internal_sss_preprocess(fp, p_hat, n);

	/* Clean up temporary memory */
	XFREE(Bcnt);
	XFREE(prefix_sum);

	return u;
}

int uss_sample(uss_t *uss, int **output) {
	int idx, ridx,
		p_hat_samples = 0,
	   *p_hat_sampled = NULL,
		sampled = 0,
		alloced = uss->mu,
	   *samples = XALLOC(uss->mu, int),
		i;
	double p_hat;

	*output = samples;

	p_hat_samples = sss_sample(uss->s, &p_hat_sampled);

	for (i = 0; i < p_hat_samples; i++) {
		idx   = p_hat_sampled[i];
		ridx  = uss->index[idx];
		p_hat = uss->p_hat[idx];

		if ( xuni_rand() < uss->w[ridx]/p_hat ) {
			if ( alloced == sampled ) {
				uss->mu = alloced = alloced*2;
				uss->samples = samples = XREALLOC(samples, alloced, int);
			}

			samples[sampled] = ridx;
			sampled++;
		}
	}

	return sampled;
}

void uss_free(uss_t *u) {
	if ( NULL != u ) {
		if ( NULL != u->s ) {
			sss_free(u->s);
		}

		if ( NULL != u->index ) {
			XFREE(u->index);
		}

		if ( NULL != u->p_hat ) {
			XFREE(u->p_hat);
		}

		if ( NULL != u->samples ) {
			XFREE(u->samples);
		}

		XFREE(u);
	}
}
