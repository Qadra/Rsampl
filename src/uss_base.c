#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"

#include "sample.h"
#include "uss_base.h"
#include "xutil.h"

// Initialize O(n^2) USS solution, using O(n^2) space.
uss_basecase_t *uss_basecase_preprocess(double *weights, int n) {
	int i, j;
	int size;
	double prior, pr;
	double *warray;

	double mu         = 0.,
		   *w         = weights;
	uss_basecase_t *u = XALLOC(1, uss_basecase_t);
	u->X              = XALLOC(n, walkers_t*);
	u->n              = n;
	walkers_t **X     = u->X;

	// Calculate the expected amount of elements we sample
	for (i = 0; i < n; i++) {
		mu += w[i];
	}

	u->mu      = ceil(mu);
	u->sampled = XALLOC(u->mu, int);

	warray = XALLOC(1 + n, double);
	for (i = 0; i < n; i++) {
		prior = 1.;
		size  = n - i;

		for (j = i; j < n; j++) {
			pr          = prior * w[j];
			warray[j-i] = pr;
			prior      *= 1. - w[j];

			// DEBUG
			// printf("pr(%d, %d) = %f\n", i, j, pr);
		}

		// DEBUG
		// printf("pr(%d, inf) = %f\n", i, prior);

		// Prior is now probability of not sampling anything.
		warray[size] = prior;
		X[i] = walkers_preprocess(size + 1, warray);

	}
	XFREE(warray);

	return u;
}

int uss_basecase_sample(uss_basecase_t *uss_base, int **output) {
	int drawn,
		i          = 0,
		n          = uss_base->n,
		sample_ptr = 0,
		mu         = uss_base->mu,
		alloced    = mu,
	   *sampled    = uss_base->sampled;
	walkers_t *w;

	*output = sampled;

	// We expect to sample mu elements.
	while (i < n) {
		w     = uss_base->X[i];
		drawn = walkers_draw(w);
		i += drawn;

		if (i < n) {
			if (sample_ptr == alloced) {
				uss_base->mu                = alloced = 2*alloced;
				uss_base->sampled = sampled = XREALLOC(sampled, alloced, int);
			}

			sampled[sample_ptr] = i;
			sample_ptr++;
			i++;
		}
	}

	//sample->samples = sample_ptr;
	return sample_ptr;
}

void uss_basecase_free(uss_basecase_t *u) {
	int i;

	for (i = 0; i < u->n; i++) {
		walkers_free(u->X[i]);
	}

	if ( NULL != u->X ) {
		XFREE(u->X);
	}

	if ( NULL != u->sampled ) {
		XFREE(u->sampled);
	}

	if ( NULL != u ) {
		XFREE(u);
	}
}
