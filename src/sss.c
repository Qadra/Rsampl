#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "defs.h"

#include "sample.h"
#include "sss.h"
#include "xutil.h"
#include "uss.h"
#include "uss_base.h"

// This is really a wrapper for the internal SSS.
sss_t *sss_preprocess(double *weights, int n) {
	subset_method_t *fp   = XALLOC(1, subset_method_t);
	fp->preprocess = (subset_preprocess_fp_t) uss_preprocess;
	fp->sample     = (subset_sample_fp_t)     uss_sample;
	fp->free       = (subset_free_fp_t)       uss_free;

	return internal_sss_preprocess(fp, weights, n);
}

static uint64_t pow_lookup[64] ={
	1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
	32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608,
	16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824,
	2147483648, 4294967296, 8589934592, 17179869184, 34359738368, 68719476736,
	137438953472, 274877906944, 549755813888, 1099511627776, 2199023255552,
	4398046511104, 8796093022208, 17592186044416, 35184372088832, 70368744177664,
	140737488355328, 281474976710656, 562949953421312, 1125899906842624,
	2251799813685248, 4503599627370496, 9007199254740992, 18014398509481984,
	36028797018963968, 72057594037927936, 144115188075855872, 288230376151711744,
	576460752303423488, 1152921504606846976, 2305843009213693952,
	4611686018427387904,
};

sss_t *internal_sss_preprocess(subset_method_t *recursive_method, double *weights, int n) {
	int i, j, idx, bucket_size,
		po2 = 1,
		L   = floor(log2(n) + 1);

	double mu = 0.,
		  *w  = weights,
	      *Q  = XALLOC(L, double),
	      *P  = XALLOC(L, double);

	sss_t *s     = XALLOC(1, sss_t);
	       s->L  = L;
	       s->n  = n;
	       s->Q  = Q;
	       s->P  = P;
	       s->w  = w;
	       s->fp = recursive_method;

	// DEBUG
	// printf("   idx " "  Q[i] " "  w[i] " " start " "   end " "  size " "    mu " "\n");

	for (j = 1, i = 0; j - 1. < n; j = i++, j = pow_lookup[i]) {
		// j is the 1-indexed index of how far we are into our weights.
		// Next power of two
		po2 = pow_lookup[i+1];
		// idx is now the 0-indexed version of j.
		idx = j - 1;

		// Maximal size of the current bucket
		bucket_size = po2 - j;

		// Check if our bucket is limited by the length of the array.
		bucket_size = likely(bucket_size + (j - 1) < n) ? 
			bucket_size : n - (j - 1);

		// Here w[idx] corresponds to p_2_k.
		Q[i] = 1. - pow(1. - w[idx], bucket_size);
		P[i] = w[idx];

		mu += w[idx] * bucket_size;
	}

	/*
	sample_t sam = {
		.n       = L,
		.weights = Q,
	};
	*/

	// TODO change to use of s->fp
	s->u       = s->fp->preprocess(Q, L);
	s->mu      = 2 * mu;
	s->samples = XALLOC(s->mu, int);

	return s;
}

int sss_sample(sss_t *sss, int **sampled) {

	int i, j, v, idx, next_power, count, po2, inner,
		n          = sss->n,
		sample_ptr = 0,
		alloced    = sss->mu,
	   *A;

	int n_samples = sss->fp->sample(sss->u, &A);

	//sampled    = sss->samples;
	int *output = XALLOC(sss->mu, int);

	bool first;
	
	double a, b, p2k, rand, num, denom;

	// DEBUG
	//printf("     i " "Q[idx] " "   p2k " "     j " "   idx\n");
	
	for (i = 0; i < n_samples; i++) {
		idx = A[i];

		// DEBUG
		//printf("Sampled: %d\n", idx);
		
		// Calculate power of two
		po2 = pow_lookup[idx]; //pow(2, idx)

		next_power = pow_lookup[idx+1] - 1; //pow(2, idx+1) - 1;

		// Limit next_power to the maximum amount of elements
		if ( next_power >= n ) {
			next_power = n;
		}

		p2k = sss->P[idx];

		assert(idx < 32); // Because we use int
		
		// The 2^k index into the array.
		j = po2 - 1; //pow(2, idx) - 1;

		// Amount of samples choosen
		count = 0;

		// Is it the first sample choosen in the bucket
		first = true;

		// DEBUG
		//printf("%6d %.4f %.4f %6d %6d\n", i, s->Q[idx], p2k, j, idx);
		
		denom = log(1. - p2k);
		while (1) {
			v = 0;

			if (p2k < 1.) {
				rand = xuni_rand();

				if ( first ) {
					first = false;
					rand = 1. - rand * sss->Q[idx];
				}

				num   = log(rand);

				if ( unlikely(isnan(num) || isinf(num)) ) {
					num = 0.;
				}

				if ( unlikely(denom == 0) ) {
					v = 0;
				} else {
					v = floor( num/denom );
				}

				j += v;
			}

			if ( j >= next_power ) {
				break;
			}

			if ( sample_ptr == alloced ) {
				sss->mu      = alloced = 2*alloced;
				output = XREALLOC(output, alloced, int);
			}

			assert(j >= 0);

			output[sample_ptr] = j;

			sample_ptr++;
			count++;
			j++;
		}

		assert(count > 0);
	}

	i = 0, j = 0;
	while (j < sample_ptr) {
		assert(output[j] >= 0);
		assert(output[j] < sss->n);

		inner = floor(log2(output[j]+1));

		assert( inner < sss->L );
		assert( inner >= 0 );

		a = sss->w[output[j]];
		b = sss->P[inner];

		assert( b > 0 );

		if (xuni_rand() <= a/b) {
			output[i] = output[j];
			i++;
		}
		j++;
	}

	//sample->samples = i;
	
	*sampled = output;

	return i;
}

void sss_free(sss_t *s) {
	if ( NULL != s ) {
		if ( NULL != s->samples ) {
			XFREE(s->samples);
		}

		if ( NULL != s->fp ) {
			s->fp->free(s->u);
			XFREE(s->fp);
		}

		if ( NULL != s->Q ) {
			XFREE(s->Q);
		}

		if ( NULL != s->P ) {
			XFREE(s->P);
		}

		XFREE(s);
	}
}
