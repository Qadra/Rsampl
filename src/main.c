#include "exact.h"

#include <R.h>
#include <Rinternals.h>

// WRS without replacement
SEXP sample_(SEXP A_, SEXP n_, SEXP k_) {
	SEXP res;
	double *A = REAL(A_);
	int n = asInteger(n_),
		k = asInteger(k_);

	exact_t *E = exact_preprocess(k, n, A);

	PROTECT(res = allocVector(INTSXP, k));

	exact_sample(E, INTEGER(res), k, n);

	exact_free(E);

	UNPROTECT(1);
	return res;
}

