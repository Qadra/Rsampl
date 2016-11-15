#ifndef EXACTWALKER
#define EXACTWALKER

#include "walkers.h"
#include "sample.h"

#include <inttypes.h>

// TODO: Remove
#include <stdio.h>
#include <unistd.h>
// END TODO

typedef struct {
	walkers_t *w;
	int n, words;
	uint64_t *bv;
	int *sampled;
} exact_walkers_t;

exact_walkers_t *exact_walkers_preprocess(double *weights, int n, int k);
int exact_walkers_sample(exact_walkers_t *this, int k, int *sampled);
void exact_walkers_free(exact_walkers_t *this);

#endif
