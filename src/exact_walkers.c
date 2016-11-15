#define _DEFAULT_SOURCE

#include "defs.h"

#include "exact_walkers.h"
#include "xutil.h"

#include <immintrin.h>

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

exact_walkers_t *exact_walkers_preprocess(double *weights, int n, int k) {
	(void)k; // Unused k

	exact_walkers_t *this = XALLOC(1, exact_walkers_t);
	double *v = weights;
	walkers_t *w = walkers_preprocess(n, v);

	int words = ceil((double)n/64.0);

	this->words = words;
	this->bv = XALLOC(words, uint64_t);
	XMEMZERO(this->bv, words);
	this->sampled = XALLOC(n, int);

	this->w = w;
	this->n = n;

	return this;
}

static inline int __ew_get_and_set_bit(uint64_t *bv, int n) {
	int word =n >> 6,   // Division by 64
		idx = 0x3f & n; // 6 LSB determine index into word

	int pattern = 1 << idx;
	int res = bv[word] & pattern;
	bv[word] |= pattern;

	return res;
}

static inline int __ew_walkers_draw(walkers_t *restrict walkers) {
	double rU = (double) xuni_rand() * walkers->n;
	int i = (int) rU;
	return (rU < walkers->probs[i]) ? i : walkers->alias[i];
}

int exact_walkers_sample(exact_walkers_t *this, int k, int *sampled) {
	walkers_t *w = this->w;
	uint64_t *bv = this->bv;
	int q = 0; // The amount we have sampled.

	while (q < k) {
		int d = walkers_draw(w); // Random sample

		if (__ew_get_and_set_bit(bv, d) == 0) {
			sampled[q++] = d;
		} else {
			continue;
		}
	}

	//memset(bv, 0, sizeof(uint64_t) * this->words);

	//*
	for (int i = 0; i < k; i++) {
		int word = sampled[i] >> 6; // Division by 64
		bv[word] = 0;
	}
	// */

	return k;
}

void exact_walkers_free(exact_walkers_t *this) {
	if (this != NULL) {
		walkers_free(this->w);
		XFREE(this->bv);
		XFREE(this);
	}
}
