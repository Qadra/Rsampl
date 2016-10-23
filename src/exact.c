#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <Rmath.h>

#include <R_ext/RS.h>		/* for Calloc() */

#include "exact.h"

exact_t *exact_preprocess(int k, int n, double *w) {
	int i, j, idx, height, leaves;
	int nodes          = n-1;
	double *p;

	exact_t *e         = malloc( sizeof(exact_t) );
	e->n               = n;
	e->w               = w;
	e->p      = p      = Calloc(nodes, double);
	e->height = height = floor(log2(n)+1);
	e->leaves = leaves = nodes-((1 << (height-1))-1);

	e->round           = 0;
	e->scales          = Calloc(2*n, double);
	e->scales_round    = Calloc(2*n, int);

	Memzero(e->scales_round, 2 * n);
	Memzero(e->scales, 2 * n);

	idx = n-1;
	for (j = leaves, i = leaves*2; j > 0; i-=2, j--, idx--) {
		p[idx-1] = w[i-2] + w[i-1];
	}

	for (i = n-1; i > 2*leaves; i -= 2, idx--) {
		p[idx-1] = w[i] + w[i-1];
	}

	if ( leaves & 1 ) {
		p[idx-1] = w[i] + p[2*idx-1];
		idx--;
	}

	while (idx > 0) {
		p[idx-1] = p[2*idx-1] + p[2*idx];
		idx--;
	}

	return e;
}

void *exact_sample(exact_t *e, int *sampled, int k, int n) {
	int i, idx, level, index; 
	double rand, scale, left_child;
	int	key               = 0,
		height            = e->height,
		leaves            = e->leaves,
		nodes             = (1 << height),
		samples           = 0;
	double *p             = e->p,
	       *w             = e->w;
	int     round         = e->round;
	double *scales		  = e->scales;
	int    *scales_round  = e->scales_round;

	e->round             += 1;

	for (i = 0; i < k; i++) {
		idx   = 1;
		index = 0;
		scale = 0.;
		key   = idx-1;

		if (scales_round[key] == round) {
			scale = scales[key];
		}

  		rand = unif_rand() * (p[key] - scale);

		// Search binary tree
		while (idx < n) {
			key = idx-1;

			// Lookup left child node
			key = (idx << 1);

			// Find most significant bit aka. current level of the tree
			level  = sizeof(unsigned int)*8 - __builtin_clz(key) - 1;

			if ( unlikely(key >= n) ) {
				// Generally we always add the number of leaves outside the fully 
				// balanced binary tree (FBBT). If we have reached a left leaf 
				// outside the FBBT, we instead calculate the offset of the 
				// parent, which will be the offset we have to add to the index.
				if ( likely( nodes != n ) && level == height ) {
					index += (key/2) & ((1 << (level-1))-1);
				} else {
					index += leaves;
				}

				left_child = w[index];
			} else {
				left_child = p[key-1];
			}

			key -= 1;

			if (scales_round[key] == round) {
				scale = scales[key];
			} else {
				scale = 0.;
			}

			idx = (idx << 1);
			if (rand > (left_child-scale)) {
				rand -= (left_child-scale);
				idx++;

				if ( unlikely(idx-1 < n && idx == n) ) {
					// Find most significant bit aka. current level of the tree
					level  = sizeof(unsigned int)*8 - __builtin_clz(idx) - 1;

					if ( likely( nodes != n ) && level == height ) {
						index += (idx/2) & ((1 << (level-1))-1);
					} else {
						index += leaves;
					}
				}

				level = height - level - 1;
				level = (level > 0) ? level : 0;
				
				index += (1 << level);
			}
		}
	
		sampled[samples++] = index;
		key                = idx-1;

		// Insert scale of leaf node
		scales[key] = scale = w[index];
		scales_round[key] = round;
		
		// Update scale values according to the path visited in the tree
		while (idx > 1) {
			idx = idx/2;
			key = idx-1;

			if (scales_round[key] == round) {
				scales[key] += scale;
			} else {
				scales[key]       = scale;
				scales_round[key] = round;
			}
		}

	}
}

void exact_free(exact_t *e){
	if ( NULL != e ) {
		if (NULL != e->scales_round) Free(e->scales_round);
		if (NULL != e->scales) Free(e->scales);
		if (NULL != e->p) Free(e->p);

		Free(e);
	}
}
