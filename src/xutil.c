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
