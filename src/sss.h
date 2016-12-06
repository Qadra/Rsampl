#ifndef H_SSS
#define H_SSS

#include "sample.h"

typedef struct {
	int     n;       // Amount of elements.
	int     L;       // ceil(log2(n))
	int     mu;      // Expected amount of elements we sample.
	double *Q;       // Pr[X_k = 1] for k = 0, ... , L, probability of sampling
	                 // from bucket X_k.
	double   *P;       // p_2^k, from above.
	double   *w;       // The original weights
	void     *u;       // Stores the recursive structure.
	int      *samples; // Used for storing sampled elements.
	subset_method_t *fp;
} sss_t;

/**
 * Initialize sorted subset sampling.
 *
 * The weights are used in the sampling procedure, so these must not be freed
 * until you are done.
 */
sss_t *sss_preprocess(double *weights, int n);

/**
 * Sample from our distribution.
 */
int sss_sample(sss_t *sss, int **sampled);

/**
 * Free allocated memory.
 *
 * The invoker is responsible for freeing 'w' after the last call to
 * sss_sample.
 */
void sss_free(sss_t *s);

/*** INTERNAL FUNCTIONS ***/
sss_t *internal_sss_preprocess(subset_method_t *recursive_method,
		double *weights, int n);

#endif
