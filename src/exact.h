#ifndef H_exact
#define H_exact

#include "sample.h"

#ifdef H_exact_list
typedef struct {
	unsigned int i;
	double       w;
	bool         e;
} tuple_t;
#endif

typedef struct {
	unsigned int n;
	double      *w;
	double      *p;
#ifdef H_exact_list
	tuple_t     *scales;
	int         round;
#else
	int         round;
	double      *scales;
	int         *scales_round;	
#endif
	unsigned int height;
	unsigned int leaves;
} exact_t;

exact_t *exact_preprocess(double *weights, int n, int k);
int exact_sample(exact_t *this, int k, int *sampled);
void     exact_free(exact_t *bin);

#endif
