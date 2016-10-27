#ifndef RSTREE
#define RSTREE

#include <stdint.h>

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

enum NODE_TYPE {
	NODE_TREE,
	NODE_ARRAY
};

typedef union {
	uint64_t i;
	double   f;
} T;

struct node;

typedef struct internal {
	double child_sum;

	struct node *right;
	struct node *left;
	int right_index,
		left_index;
} internal_t;

typedef struct {
	double sum, sampled_sum, p_hat;
	int size, sampled_size, offset, bucket_index;
} arr_t;

typedef struct node {
	enum NODE_TYPE type;
	struct node *parent;
	int index,
		parent_index;

	union {
		internal_t internal;
		arr_t array;
	};
} node_t;

typedef struct {
	int size, allocated, offset;
	node_t **nodes;
} queue_t;

typedef struct {
	// The bucket sorted weights.
	double *w;

	// Stores the reverse lookup, to the correct index in the original array.
	int *reverse_map;

	// Used for faster undo operations.
	node_t **buckets;

	// Queue used or sampling
	queue_t queue;

	node_t * tree;
} rstree_t;

rstree_t *rstree_preprocess(int k, int n, double *w);
int rstree_sample(rstree_t *st, int *samples, int k, int n);
void rstree_free(rstree_t *st);

#endif
