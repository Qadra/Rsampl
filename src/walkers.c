#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "xutil.h"
#include "walkers.h"

/**
 * Inspired by: http://explodecomputer.com/page.php?id=154
 *              https://github.com/wch/r-source/blob/ed66b715221d2720f5b334470335635bada520b1/src/main/random.c#L346
 **/
walkers_t *walkers_preprocess(const int n, void *restrict v) {
	double *q, sum = 0;
	int i, j, k;
	double *m = v;
	double *x = v;
	int *HL, *H, *L, *a;
	int intSize = n * sizeof(int), 
		doubleSize = n * sizeof(double), 
		structSize = sizeof(walkers_t);

	walkers_t *walkers = xmalloc(structSize);
	walkers->alias = a = xmalloc(intSize);
	walkers->probs = q = xmalloc(doubleSize);
	walkers->n = n;

	HL = xmalloc(intSize);
	H = HL - 1; 
	L = HL + n;

	for (i = 0; i < n; i++) {
		sum += m[i];
		x[i] = (double)m[i];
	}
	
	/* Normalize weights */
	for (i = 0; i < n; i++) {
		x[i] /= sum;
	}

	/** 
	 * Create the alias tables.
	 * The idea is that for HL[0] ... L-1 label the entries with q < 1
	 * and L ... H[n-1] label those >= 1.
	 * By rounding error we could have q[i] < 1. or > 1. for all entries.
	 **/
	for (i = 0; i < n; i++) {
		q[i] = x[i] * n;
		if (q[i] < 1.) {
			*++H = i;
		} else {
			*--L = i;
		}
	}

	if (H >= HL && L < HL + n) { /* So some q[i] are >= 1 and some < 1 */
		for (k = 0; k < n-1; k++) {
			i = HL[k];
			j = *L;
			a[i] = j;
			q[j] += q[i] - 1;

			if (q[j] < 1.) { 
				L++;
			}

			if (L >= HL + n) {
				break; /* now all are >= 1 */
			}
		}
	}

	for (i = 0; i < n; i++) {
		q[i] += i;
	}

	free(HL);

	return walkers;
}

int walkers_draw(walkers_t *restrict walkers) {
	double rU = (double) xuni_rand() * walkers->n;
	int i = (int) rU;
	return (rU < walkers->probs[i]) ? i : walkers->alias[i];
}

void walkers_free(walkers_t *restrict walkers) {
	if (NULL != walkers) {
		if (NULL != walkers->alias) {
			free(walkers->alias);
		}

		if (NULL != walkers->probs) {
			free(walkers->probs);
		}

		free(walkers);
	}
}
