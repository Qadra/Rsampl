#define _DEFAULT_SOURCE

#include <stdlib.h>  /* exit(), calloc(), malloc(), EXIT_FAILURE */
#include <stdio.h>   /* fprintf() */
#include <string.h>  /* memset() */
#include <cpuid.h>   /* __cpuid(), bit_POPCNT */
#include <stdbool.h> /* bool */
#include <math.h>

#include "defs.h"

#include "xutil.h"

unsigned int I1 = 1234;
unsigned int I2 = 5678;

extern inline double xuni_rand(void);
extern inline int    xceil_log2(unsigned long long x);

void shuffle_array(int n, double *w) {
	int i;

	for (i = n - 1; i > 0; i--) {
		int idx = (int)floor((double)xrand() * (i + 1));
		double t = w[i];
		w[i] = w[idx];
		w[idx] = t;
	}
}

void normalize_array(int n, double *w) {
	int i;
	double sum = 0.0;

	for (i = 0; i < n; i++) {
		sum += w[i];
	}

	for (i = 0; i < n; i++) {
		w[i] /= sum;
	}
}

double xrand(void) {
	return (double)random() / RAND_MAX;
}

void xerror(char *msg, int line, char *file) {
	fprintf(stderr, "Error %d@%s: %s\n", line, file, msg);
	exit(EXIT_FAILURE);
}

void xmemset(void *p, size_t size) {
	if (p != NULL) {
		memset(p, '\0', size);
	}
}

void *xcalloc(size_t nmemb, size_t size) {
	void *p;

	if (size == 0 || nmemb == 0) {
		return NULL;
	}

	p = calloc(nmemb, size);

	if (p == NULL) {
		xerror("Unable to c-allocate memory", __LINE__, __FILE__);
	}

	return p;
}

void *xmalloc(size_t size) {
	void *p;

	if (size == 0) {
		return NULL;
	}

	p = malloc(size);

	if (p == NULL) {
		xerror("Unable to m-allocate memory", __LINE__, __FILE__);
	}

	return p;
}

void *xrealloc(void *ptr, size_t size) {
	void *p;

	if (size == 0) {
		return NULL;
	}

	p = realloc(ptr, size);

	if (p == NULL) {
		xerror("Unable to m-allocate memory", __LINE__, __FILE__);
	}

	return p;
}
