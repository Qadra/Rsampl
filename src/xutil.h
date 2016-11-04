#ifndef H_XUTIL
#define H_XUTIL

#ifndef BYTE
#define BYTE 8
#endif

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

#include <stdlib.h>  /* random */
#include <math.h>    /* floor */

#include <immintrin.h>
#include <stdint.h>

inline double xfrac (double x) {
	return x - floor(x);
}

extern unsigned int I1, I2;

inline double xuni_rand(void) {
	I1 = 36969*(I1 & 0177777) + (I1>>16);
	I2 = 18000*(I2 & 0177777) + (I2>>16);
	return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
}

inline int xceil_log2(unsigned long long x) {
	static const unsigned long long t[6] = {
		0xFFFFFFFF00000000ull,
		0x00000000FFFF0000ull,
		0x000000000000FF00ull,
		0x00000000000000F0ull,
		0x000000000000000Cull,
		0x0000000000000002ull
	};

	int y = (((x & (x - 1)) == 0) ? 0 : 1);
	int j = 32;
	int i;

	for (i = 0; i < 6; i++) {
		int k = (((x & t[i]) == 0) ? 0 : j);
		y += k;
		x >>= k;
		j >>= 1;
	}

	return y;
}

double xrand(void);

void xerror(char *msg, int line, char *file);

void xmemset(void *p, size_t size); 

void *xcalloc(size_t nmemb, size_t size);

void *xmalloc(size_t size); 

void *xrealloc(void *ptr, size_t size);

void shuffle_array(int n, double *w);

void sort_array(int n, double *w);

void normalize_array(int n, double *w);
#endif
