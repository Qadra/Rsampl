#ifndef H_SAMPLE
#define H_SAMPLE

#include <stdint.h>
#include <stdbool.h>

#ifndef CLASS_T
	#define CLASS_T int
#endif

typedef CLASS_T class_t;

typedef void *(*preprocess)(void *);
typedef int  *(*draw)(void *);
typedef void  (*clean)(void *);

typedef void *(*preprocess_fp_t)(double* weight, int n, int k);
typedef int (*sample_fp_t)(uintptr_t *data, int k, int *sampled);
typedef void (*free_fp_t)();

typedef struct {
	draw draw;
	clean clean;
	preprocess preprocess;
	char *name;
} subset_t;

typedef struct {
	preprocess_fp_t preprocess;
	sample_fp_t sample;
	free_fp_t free;

	// Used to store the preprocessed structure.
	uintptr_t *data;
} method_t;

typedef struct {
	int       n;          // Amount of classes and probabilities
	int       k;          // Amount of samples expected
	int       samples;    // Amount of samples sampled
	class_t  *classes;    // The classes (types) to sample
	double   *weights;    // The weights to sample from
	void     *redundancy; // Preproccessed datastructure
	subset_t *fptr;       // Function pointers
	subset_t *internal;   // Function pointers
} sample_t;

sample_t * create_sample(subset_t *fptr, const int n, class_t *restrict classes);
void  preprocess_sample(sample_t *restrict sample, double *weights, int k);
void *draw_sample(sample_t *restrict sample); 
void  free_sample(sample_t *restrict sample);

void wrs_preprocess(method_t *method, double *weights, int n, int k);
int wrs_sample(method_t *method, int k, int *sampled);
void wrs_free(method_t *method);

#endif
