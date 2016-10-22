#include <R.h>
#include <Rinternals.h>

SEXP wrs_(SEXP A_, SEXP n_) {
	double *A = REAL(A_);
	int n =  asInteger(n_);

	return ScalarReal(A[n]);
}
