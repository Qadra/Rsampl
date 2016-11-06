#ifndef H_SAMPLE
#define H_SAMPLE

#include <stdint.h>
#include <stdbool.h>

#ifndef CLASS_T
	#define CLASS_T int
#endif

typedef CLASS_T class_t;

typedef void *(*preprocess_fp_t)(double* weight, int n, int k);
typedef int   (*sample_fp_t)(uintptr_t *data, int k, int *sampled);
typedef void  (*free_fp_t)(uintptr_t *data);

typedef void *(*subset_preprocess_fp_t)(double* weight, int n);
typedef int   (*subset_sample_fp_t)(uintptr_t *data, int **sampled);
typedef void  (*subset_free_fp_t)(uintptr_t *data);

typedef struct {
	preprocess_fp_t preprocess;
	sample_fp_t     sample;
	free_fp_t       free;

	// Used to store the preprocessed structure.
	uintptr_t *data;
} method_t;

typedef struct {
	subset_preprocess_fp_t preprocess;
	subset_sample_fp_t     sample;
	subset_free_fp_t       free;

	// Used to store the preprocessed structure.
	uintptr_t *data;
} subset_method_t;

void wrs_preprocess(method_t *method, double *weights, int n, int k);
int wrs_sample(method_t *method, int k, int *sampled);
void wrs_free(method_t *method);

#endif
