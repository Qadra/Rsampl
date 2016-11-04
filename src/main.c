#include "sample.h"
#include "exact.h"
#include "rstree.h"

#include <R.h>
#include <R_ext/RS.h> /* For Calloc etc. */
#include <Rinternals.h>

#include <stdbool.h>
#include <string.h>

void r_wrs_finalizer(SEXP ds_);
void wrs_set_method_t(const char *name, method_t *ptr);

// WRS without replacement
SEXP r_wrs_sample(SEXP W_, SEXP k_, SEXP methodstr_, SEXP ds_) {
	SEXP res;
	int k = asInteger(k_);

	bool has_struct = true;
	if (TYPEOF(ds_) == NILSXP || TYPEOF(ds_) != EXTPTRSXP) {
		has_struct = false;
	}

	method_t  method, // Used for allocating on stack.
			 *method_ptr;

	if(!has_struct) {
		double *W = REAL(W_); // The weights
		int n = XLENGTH(W_);  // Could be bigger than an int type.. <- TODO

		if (k > n) {
			error("sample size 'k' is larger than input set");
		}

		method_ptr = &method;
		wrs_set_method_t(CHAR(asChar(methodstr_)), method_ptr);

		// Preprocess the input into a tree datastructure
		wrs_preprocess(method_ptr, W, n, k);

	} else {
		method_ptr = R_ExternalPtrAddr(ds_);
	}

	// Create return output
	res = PROTECT(allocVector(INTSXP, k));

	// Do the sampling
	GetRNGstate();
	wrs_sample(method_ptr, k, INTEGER(res));
	PutRNGstate();

	// Cleanup allocated memory
	if (!has_struct) {
		wrs_free(method_ptr);
	}

	UNPROTECT(1);
	return res;
}

// Preprocessing for WRS without replacement
SEXP r_wrs_preprocess(SEXP W_, SEXP k_, SEXP methodstr_) {
	SEXP res;
	double *W = REAL(W_); // The weights
	int n = XLENGTH(W_),  // Could be bigger than an int type.. <- TODO
		k = asInteger(k_);

	if (k > n) {
		error("sample size 'k' is larger than input set");
	}

	method_t *method = Calloc(1, method_t);
	wrs_set_method_t(CHAR(asChar(methodstr_)), method);

	wrs_preprocess(method, W, n, k);

	SEXP structure = PROTECT(R_MakeExternalPtr(method, R_NilValue, R_NilValue));

	R_RegisterCFinalizerEx(structure, r_wrs_finalizer, TRUE);

	UNPROTECT(1);
	return structure;
}

void r_wrs_finalizer(SEXP ds_) {
	method_t *method = R_ExternalPtrAddr(ds_);

	R_ClearExternalPtr(ds_);

	wrs_free(method);
	Free(method);
}

void wrs_set_method_t(const char *name, method_t *method) {
	if (strncmp(name, "binary", 100) == 0) {
		method->preprocess = (preprocess_fp_t)exact_preprocess;
		method->sample     = (sample_fp_t)    exact_sample;
		method->free       = (free_fp_t)      exact_free;
	} else if (strncmp(name, "rstree", 100) == 0) {
		method->preprocess = (preprocess_fp_t)rstree_preprocess;
		method->sample     = (sample_fp_t)    rstree_sample;
		method->free       = (free_fp_t)      rstree_free;
	} else {
		error("Unknown method");
	}
}
