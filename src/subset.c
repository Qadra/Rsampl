#include "sample.h"

#include "uss.h"
#include "uss_base.h"
#include "sss.h"

#include <R.h>
#include <R_ext/RS.h> /* For Calloc etc. */
#include <Rinternals.h>
#include <Rdefines.h>

#include <stdbool.h>
#include <string.h>

void r_subset_finalizer(SEXP r_structure);

SEXP r_subset_sample(SEXP W_, SEXP r_structure, SEXP draws_) {
	(void) draws_;
	int n = XLENGTH(W_),
		*output;

	sss_t *structure;

			  /*** Check if we have a preprocessed structure ***/

	bool has_struct = true;
	if (TYPEOF(r_structure) != EXTPTRSXP) {
		has_struct = false;
	}

	if (has_struct) {
		structure = R_ExternalPtrAddr(r_structure);
	} else {
		structure = sss_preprocess(REAL(W_), n);
	}

			 /*** Handle if several draws have been requested ***/

	int draws = 0;
	if (!isNull(draws_)) {
		draws = asInteger(draws_);
	}


					   /*** Handle the actual draws ***/
	SEXP res;

	int n_indices = sss_sample(structure, &output);
	res = PROTECT(allocVector(INTSXP, n_indices));
	Memcpy(INTEGER(res), output, n_indices);

	free(output);

	if (!has_struct) {
		sss_free(structure);
	}

	UNPROTECT(1);
	return res;
}

SEXP r_subset_preprocess(SEXP W_) {
	int n = XLENGTH(W_);

	sss_t *structure = sss_preprocess(REAL(W_), n);

	SEXP r_structure = PROTECT(R_MakeExternalPtr(structure, R_NilValue, R_NilValue));

	R_RegisterCFinalizerEx(r_structure, r_subset_finalizer, TRUE);

	UNPROTECT(1);
	return r_structure;
}

void r_subset_finalizer(SEXP r_structure) {
	sss_t *structure = R_ExternalPtrAddr(r_structure);

	R_ClearExternalPtr(r_structure);

	sss_free(structure);
}
