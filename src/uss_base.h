#ifndef USS_BASECASE
#define USS_BASECASE

#include "walkers.h"
#include "sample.h"

typedef struct {
	int n;
	int mu;
	int *sampled;
	walkers_t **X;
} uss_basecase_t;

/**
 * Initialize the O(n^2) basecase solution for unsorted subset sampling
 */
uss_basecase_t *uss_basecase_preprocess(double *weights, int n);

/**
 * Sample a subset, you should free the returned array yourself
 */
int uss_basecase_sample(uss_basecase_t *uss_base, int **output);

/**
 * Free all allocated memory
 */
void uss_basecase_free(uss_basecase_t *u);

#endif
