#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "xutil.h"
#include "sample.h"
#include "exact.h"

exact_t *exact_preprocess(double *weights, int n, int k) {

	// Unused but required parameter for function prototype matching.
	(void)k;

	int i, j, idx, height, leaves;
	int nodes          = n-1;
	double *p, *w      = weights;

	exact_t *e         = XALLOC(1, exact_t);
	e->n               = n;
	e->w               = w;
	e->p      = p      = XALLOC(nodes, double);
	e->height = height = floor(log2(n)+1);
	e->leaves = leaves = nodes-((1 << (height-1))-1);

#ifdef H_exact_list
	e->scales          = XALLOC(k*(log2(n)+2), tuple_t);
#else
	e->round           = 0;
	e->scales          = XALLOC(2*n, double);
	e->scales_round    = XALLOC(2*n, int);

	memset(e->scales_round, 0, 2*n*sizeof(int));
	memset(e->scales, 0, 2*n*sizeof(double));
#endif

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

//	IACA_END

	return e;
}

int exact_sample(exact_t *this, int k, int *sampled) {
	int i, idx, level, index; 
	double rand, scale, left_child;
	exact_t *e            = this;
	int	key               = 0,
		n                 = this->n,
		height            = e->height,
		leaves            = e->leaves,
		nodes             = (1 << height),
		samples           = 0;
	double *p             = e->p,
	       *w             = e->w;

#ifdef H_exact_list
	tuple_t *t;
	tuple_t *list         = e->scales;
	e->round              = 0;
#else
	int     round         = e->round;
	double *scales		  = e->scales;
	int    *scales_round  = e->scales_round;

	e->round             += 1;
#endif

	for (i = 0; i < k; i++) {
		idx   = 1;
		index = 0;
		scale = 0.;
		key   = idx-1;

#ifdef H_exact_list
		scale = 0.;
#else
		if (scales_round[key] == round) {
			scale = scales[key];
		}
#endif
  		rand = XRANDFUN()*(p[key]-scale);

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

#ifdef H_exact_list
			scale = 0.;
#else
			if (scales_round[key] == round) {
				scale = scales[key];
			} else {
				scale = 0.;
			}
#endif

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
#ifdef H_exact_list
		t = list + e->round;
		t->i = index;
		t->w = scale = w[index];
		t->e = true;
		e->round++;
		w[index] = 0.;	
#else
		scales[key] = scale = w[index];
		scales_round[key] = round;
#endif
		
		// Update scale values according to the path visited in the tree
		while (idx > 1) {
			idx = idx/2;
			key = idx-1;
				
#ifdef H_exact_list
			t = list + e->round;
			t->i = key;
			t->w = scale;
			t->e = false;
			p[key] -= scale;
			e->round++;
#else
			if (scales_round[key] == round) {
				scales[key] += scale;
			} else {
				scales[key]       = scale;
				scales_round[key] = round;
			}
#endif
		}

	}

#ifdef H_exact_list
	for (i = 0; i < e->round; i++) {
		t = list + i;
		if (t->e) {
			w[t->i] += t->w;
		} else {
			p[t->i] += t->w;
		}
	}
#endif

	return samples;
}

void exact_free(exact_t *e){
	if ( NULL != e ) {
#ifndef H_exact_list 
		if (NULL != e->scales_round) {
			XFREE(e->scales_round);
		}
#endif
		if (NULL != e->scales) {
			XFREE(e->scales);
		}

		if (NULL != e->p) {
			XFREE(e->p);
		}
		XFREE(e);
	}
}
