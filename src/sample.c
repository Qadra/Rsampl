#include <stdbool.h>

#include "xutil.h"
#include "sample.h"

						  /* ALTERED IMPLEMENTATION */

void wrs_preprocess(method_t *method, double *weights, int n, int k) {
	method->data = method->preprocess(weights, n, k);
}

int wrs_sample(method_t *method, int k, int *sampled) {
	return method->sample(method->data, k, sampled);
}

void wrs_free(method_t *method) {
	method->free(method->data);
}
