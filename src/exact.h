#ifndef H_exact
#define H_exact

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

typedef struct {
	unsigned int n;
	double      *w;
	double      *p;
	int         round;
	double      *scales;
	int         *scales_round;	
	unsigned int height;
	unsigned int leaves;
} exact_t;

exact_t *exact_preprocess(int k, int n, double *w);
void *exact_sample(exact_t *e, int *sampled, int k, int n);
void  exact_free(exact_t *bin);

#endif
