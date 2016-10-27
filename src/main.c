#include "exact.h"
#include "rstree.h"

#include <R.h>
#include <Rinternals.h>

// WRS without replacement
SEXP sample_impl(SEXP A_, SEXP k_) {
	SEXP res;
	double *A = REAL(A_); // The weights
	int n = XLENGTH(A_),  // Could in theory be bigger than an int type..
		k = asInteger(k_);

	if (k > n) {
		error("sample size 'k' is larger than input set");
	}

	// Preprocess the input into a tree datastructure
	//exact_t *E = exact_preprocess(k, n, A);
	rstree_t *E = rstree_preprocess(k, n, A);

	// Create return output
	PROTECT(res = allocVector(INTSXP, k));

	GetRNGstate();
	
	// Do the sampling
	//exact_sample(E, INTEGER(res), k, n);
	rstree_sample(E, INTEGER(res), k, n);
	PutRNGstate();

	// Cleanup allocated memory
	//exact_free(E);
	rstree_free(E);

	UNPROTECT(1);
	return res;
}
