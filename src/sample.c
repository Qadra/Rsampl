#include <stdbool.h>

#include "xutil.h"
#include "sample.h"

sample_t * create_sample(subset_t *fptr, const int n, 
		class_t *restrict classes) {
	sample_t *sample = xmalloc(sizeof(sample_t));

	sample->n          = n;
	sample->classes    = classes;
	sample->fptr       = fptr;
	sample->internal   = NULL;
	sample->redundancy = NULL;
	sample->k          = 0;
	sample->weights    = NULL;
	sample->samples    = 0;

	return sample;
}

void normalize_weights(sample_t *restrict sample, double *weights, int k) {
	int i, n = sample->n;
	double scale, sum = 0.;
	if (k != 0) {
		for (i = 0; i < n; i++) {
			sum += weights[i];
		}

		scale = sum/k;

		for (i = 0; i < n; i++) {
			weights[i] /= scale;
		}
	}
}

void preprocess_sample(sample_t *restrict sample, double *weights, int k) {
	sample->weights    = weights;
	sample->k          = k;
	sample->redundancy = sample->fptr->preprocess(sample);
}

void *draw_sample(sample_t *restrict sample) { 
	int i;
	int *indices     = sample->fptr->draw(sample);
	int samples      = sample->samples;
	class_t *classes = xmalloc(sizeof(class_t)*samples);

	for (i = 0; i < samples; i++) {
		classes[i] = sample->classes[indices[i]];
	}

	return classes;
}

void free_sample(sample_t *restrict sample) {
	if (NULL != sample) {
		if ( NULL != sample->fptr ) {
			sample->fptr->clean(sample->redundancy);	
		}
		free(sample);
	}
}

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
