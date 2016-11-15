#ifndef H_WALKER
#define H_WALKER

/**
 * This algorithm have redundancy of O(2n*log(n) + n), expected query time of 
 * O(1) and O(n) preproccessing time.
 */

#include <stdint.h>

typedef struct {
	int    n;               // Amount of classes and probabilities
	int    *alias;          // The classes to sample
	double *probs;          // The probabilities to sample
} walkers_t;

walkers_t *walkers_preprocess(
	const int        n,     // number of classes
	void *restrict x      // relative weights of each class
);

int walkers_draw(walkers_t *restrict walkers);

void walkers_free(walkers_t *restrict walkers);

#endif
