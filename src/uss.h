#ifndef USS
#define USS

#include "sample.h"
#include "sss.h"
#include "uss_base.h"

typedef struct {
	// Number of weights.
	int n;

	// The sum of all weights, is also the expected amount of samples.
	int mu;

	// Preallocated space for samples.
	int *samples;

	// Used for reverse lookups of elements sampled with SSS.
	int *index;

	// Our SSS sample structure.
	double *p_hat;

	// The original weights.
	double *w;

	// Used for sorted sampling.
	sss_t *s;
} uss_t;

uss_t *uss_preprocess(double *weights, int n);
int uss_sample(uss_t *uss, int **output);
void   uss_free(uss_t *u);

#endif
