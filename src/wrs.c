#include "sample.h"
#include "exact.h"
#include "rstree.h"
#include "exact_walkers.h"

#include <R.h>
#include <R_ext/RS.h> /* For Calloc etc. */
#include <Rinternals.h>
#include <Rdefines.h>

#include <stdbool.h>
#include <string.h>

// For basename()
#include <libgen.h>
#include <stdlib.h>


void r_wrs_finalizer(SEXP ds_);
void wrs_set_method_t(const char *name, method_t *ptr);

// WRS without replacement
SEXP r_wrs_sample(SEXP W_, SEXP k_, SEXP methodstr_, SEXP ds_, SEXP draws_) {
	SEXP res;
	int k = asInteger(k_);

	bool has_struct = true;
	if (TYPEOF(ds_) == NILSXP || TYPEOF(ds_) != EXTPTRSXP) {
		has_struct = false;
	}

	method_t  method, // Used for allocating on stack.
			 *method_ptr;

	if(!has_struct) {
		// We assume that it has been check that W_ is a parameter.

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

	int draws = 0;

	//Rprintf("Draws is type: %s\n", type2char(TYPEOF(draws_)));
	if (!isNull(draws_) && isInteger(draws_)) {
		draws = asInteger(draws_);
		//Rprintf("Draws is size: %d\n", draws);
	}

	// Do the sampling
	if (draws == 0) { // Number of draws is unspecified
		res = PROTECT(allocVector(INTSXP, k));
		GetRNGstate();
		wrs_sample(method_ptr, k, INTEGER(res));
		PutRNGstate();
	} else { // Specific number of draws
		res = PROTECT(allocMatrix(INTSXP, k, draws));

		/* Currently I don't see this making sense.
		SEXP dimnames;

		// Specify dimension names
		PROTECT(dimnames = allocVector(VECSXP, 2));
		SET_VECTOR_ELT(dimnames, 0, GET_NAMES(k_));
		SET_VECTOR_ELT(dimnames, 1, GET_NAMES(draws_));

		setAttrib(res, R_DimNamesSymbol, dimnames);
		*/

		int *arr = INTEGER(res);

		for (int i = 0; i < draws; i++) {
			GetRNGstate();
			int *l_arr = &arr[i * k];
			wrs_sample(method_ptr, k, l_arr);
			PutRNGstate();
		}
	}

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
	} else if (strncmp(name, "srmethod", 100) == 0) {
		method->preprocess = (preprocess_fp_t)exact_walkers_preprocess;
		method->sample     = (sample_fp_t)    exact_walkers_sample;
		method->free       = (free_fp_t)      exact_walkers_free;
	} else {
		error("Parameter 'method' must be either 'binary', 'rstree' or 'srmethod'");
	}
}

int print_node(node_t *nodes, FILE *dot, int idx, int rank, int name, double *w) {
	node_t *node = &nodes[idx];
	if (node->type == NODE_TREE) {

		fprintf(dot, "%d [rank=%d, label=\""
				"Internal\n"
				"sum: %f\\l"
				"Inline index = %d\\l"
				"Left index = %d\\l"
				"Right index = %d\\l"
				"Parent index = %d\\l"
				"\"];\n",
				name, rank, node->internal.child_sum,
				node->index,
				node->internal.left_index,
				node->internal.right_index,
				node->parent_index
				);

		//fprintf(dot, "rankdir=LR;\n");

		int lname = print_node(nodes, dot, node->internal.left_index, rank + 1, name + 1, w);

		fprintf(dot, "%d -> %d;", name, name + 1);

		int rname = print_node(nodes, dot, node->internal.right_index, rank + 1, lname + 1, w);
		fprintf(dot, "%d -> %d;", name, lname + 1);

		return rname;
	} else {
		//arr_t *ar = &node->array;

		fprintf(dot,
				"%d [rank=%d, shape=box, label=\""
				"Array {\\N}\\n"
				"sum = %e\\l"
				"size = %d\\l"
				"p_hat = %e\\l"
				"Bucket index = %d\\l"
				"Sampled size = %d\\l"
				"Inline index = %d\\l"
				"Parent index = %d\\l"
				"Offset=%d\\l\"];\n",
				name, rank, node->array.sum, node->array.size,
				node->array.p_hat, node->array.bucket_index,
				node->array.sampled_size, node->index, 
				node->parent_index, node->array.offset);

		return name;
	}

	return 0;
}

SEXP print_rstree(SEXP filename_, SEXP ds_) {
	method_t *method_ptr = R_ExternalPtrAddr(ds_);
	rstree_t *ds = (rstree_t *)method_ptr->data;
	
	char *filename = strdup(CHAR(asChar(filename_)));

	FILE *dot = fopen(filename, "w");
	char *name = basename(filename);
	fprintf(dot, "digraph \"%s\" {", name);

	print_node(ds->tree, dot, 0, 0, 0, ds->w);

	fprintf(dot, "}");
	fclose(dot);

	free(filename);

	return R_NilValue;
}
