#include "rstree.h"

#include "defs.h"
#include "xutil.h"

#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include <inttypes.h>

static int64_t mask = ~0xFFFFFFFFFFFFF,
		mask2 = 0x7FF;

static inline void __queue_insert(queue_t *queue, node_t *node) {
	if (queue->allocated < queue->size + 1) {
		queue->nodes = XREALLOC(queue->nodes, ((queue->allocated + 1) * 2), node_t*);
		queue->allocated = (queue->allocated + 1) * 2;
	}

	queue->nodes[queue->size] = node;
	queue->size++;
}

static inline void __queue_insert_fast(queue_t *queue, node_t *node) {
	assert(queue->allocated > queue->size);
	queue->nodes[queue->size] = node;
	queue->size++;
}

static inline int __bucket_index(double pi, int64_t L) {
	int64_t exp;
	T rep = {
		.f = pi
	};
	
	//int preexp;
	//frexp(pi, &preexp);

	exp = (rep.i & mask) >> 52;
	//int sign = (exp & invmask2) >> 11;

	exp &= mask2;
	exp -= 1022;

	if (exp == 1) {
		exp = 0;
	}

	exp = -exp;

	if (exp > L) {
		exp = L;
	}

	return exp;
}

// TODO: L should be adjusted so we rarely, if ever, need to visit the last
// bucket.
rstree_t *rstree_preprocess(double *weights, int n, int k) {
	// Parameter k must be there because of function pointers in sample.h
	(void)k;

	int L = 63,
		*bucket_sizes  = XALLOC(L + 1, int),
		*bucket_prefix = XALLOC(L + 1, int),
		*reverse_map   = XALLOC(n,     int),
		i;

	if (L > 63) {
		L = 63;
	};

	double *w = weights,
		   *w_shuffled = XALLOC(n,     double),
		   *bucket_max = XALLOC(L + 1, double);

	node_t *inline_tree;

	rstree_t *this = XALLOC(1, rstree_t);
	memset(this, 0, sizeof(rstree_t));

	this->w = w_shuffled;
	this->reverse_map = reverse_map;
	this->queue.allocated = xceil_log2(L);
	this->queue.nodes = XALLOC(this->queue.allocated, node_t*);

	// Clamp, so total sum equals 1.
	normalize_array(n, w);

	memset(bucket_sizes, 0, sizeof(int) * (L + 1));
	memset(bucket_max, 0, sizeof(double) * (L + 1));

	/**
	 * Calculate the size of our buckets + the prefix sum of the bucket sizes
	 */
	for (i = 0; i < n; i++) {
		int idx = __bucket_index(w[i], L);
		bucket_sizes[idx]++;

		if (bucket_max[idx] < w[i]) {
			bucket_max[idx] = w[i];
		}
	}

	int bsum = 0;
	for (i = 0; i <= L; i++) {
		bucket_prefix[i] = bsum;
		bsum += bucket_sizes[i];
	}

	/**
	 * Add the elements (weights) to their respective buckets in a single
	 * sequential array.
	 */

	// Use this for storing the amount of elements in each bucket.
	memset(bucket_sizes, 0, sizeof(int) * (L + 1));

	for (i = 0; i < n; i++) {
		int idx = __bucket_index(w[i], L);
		int offset = bucket_sizes[idx] + bucket_prefix[idx];

		w_shuffled[offset] = w[i];
		assert(offset >= 0 && offset < n);
		reverse_map[offset] = i;
		bucket_sizes[idx]++;
	}

	/**
	 * Add the buckets to a queue for construction.
	 */
	queue_t *qA = XALLOC(1, queue_t);
	XMEMZERO(qA, 1);

	queue_t *qB = XALLOC(1, queue_t);
	XMEMZERO(qB, 1);

	int buckets = 0;
	for (i = L; i >= 0; i--) {
		if (bucket_sizes[i] == 0) {
			continue;
		}
		buckets++;
	}
	this->buckets = XALLOC(buckets, node_t*);

	int tree_size = 2 * buckets - 1,
		inline_index = tree_size - 1;
	inline_tree = XALLOC(tree_size, node_t);

	/**
	 * Create the last elements of our tree
	 */
	for (i = L; i >= 0; i--) {
		if (bucket_sizes[i] == 0) {
			continue;
		}

		node_t *node = XALLOC(1, node_t);
		node->type = NODE_ARRAY;
		node->array.size = bucket_sizes[i];
		node->array.p_hat = bucket_max[i];

		node->array.sum  = bucket_sizes[i] * node->array.p_hat;
		node->array.sampled_size = 0;
		node->array.sampled_sum  = 0.0;
		node->array.offset = bucket_prefix[i];
		node->array.bucket_index = --buckets;
		node->index = INT_MAX;
		node->parent_index = INT_MAX;
		this->buckets[buckets] = node;

		__queue_insert(qA, node);
	}

	/**
	 * Build the rest of the tree, bottom-up.
	 */
	while (1) {
		node_t *a = qA->nodes[qA->offset++];
		// Either we have only one node left
		if (qA->size == qA->offset) {
			if (a->type == NODE_ARRAY && inline_index == 0) {
				// If this happens, we only have a single array in our tree.
				// This means we don't hit the else case here and thus our
				// node is not copied into the array. For this reason we need
				// to add it to the tree here.
				
				a->index = inline_index;
				this->buckets[a->array.bucket_index] = &inline_tree[inline_index];
				memcpy(&inline_tree[inline_index], a, sizeof(node_t));
				XFREE(a);
				a = &inline_tree[inline_index];
			}

			__queue_insert(qB, a);

		} else { // Or more than one
			if (a->index == INT_MAX) {
				int index = inline_index--;
				a->index = index;
				memcpy(&inline_tree[index], a, sizeof(node_t));

				XFREE(a);
				a = &inline_tree[index];

				if (a->type == NODE_ARRAY) {
					this->buckets[a->array.bucket_index] = &inline_tree[index];
				}
			}

			node_t *b = qA->nodes[qA->offset++];
			if (b->index == INT_MAX) {
				int index = inline_index--;
				b->index = index;
				memcpy(&inline_tree[index], b, sizeof(node_t));

				XFREE(b);
				b = &inline_tree[index];

				if (b->type == NODE_ARRAY) {
					this->buckets[b->array.bucket_index] = &inline_tree[index];
				}
			}

			node_t *parent = &inline_tree[inline_index];
			parent->index = inline_index--;

			parent->type = NODE_TREE;
			parent->internal.right = a;
			parent->internal.left = b;

			parent->internal.right_index = a->index;
			parent->internal.left_index = b->index;

			parent->internal.child_sum = a->internal.child_sum + b->internal.child_sum;
			parent->parent = NULL;
			parent->parent_index = INT_MAX;


			a->parent = parent;
			a->parent_index = parent->index;

			b->parent = parent;
			b->parent_index = parent->index;

			__queue_insert(qB, parent);
		}

		if (qA->size == qA->offset) {
			if (qB->size == 1) {
				break;
			}

			qA->offset = 0;
			qA->size = 0;
			queue_t *tmp = qA;
			qA = qB;
			qB = tmp;
		}
	}

	//assert(qB->size == 1);
	//memcpy(inline_tree, qB->nodes[0], sizeof(node_t));
	this->tree = inline_tree;

	XFREE(qA->nodes);
	XFREE(qA);
	XFREE(qB->nodes);
	XFREE(qB);
	XFREE(bucket_sizes);
	XFREE(bucket_prefix);
	XFREE(bucket_max);

	return this;
}

int rstree_sample(rstree_t *this, int k, int *sampled) {
	rstree_t *st = this;
	node_t *nodes = st->tree;
	double *w = st->w;
	int sampled_buckets = 0;
	
	queue_t *qA = &st->queue;

	//uint64_t rejects = 0;

	// We sample k numbers.
	for (int i = 0; i < k; i++) {
		node_t *n = &nodes[0];
		double rand = XRANDFUN() * n->internal.child_sum;
		double sample_weight = 0.0;

		do {
			if (n->type == NODE_TREE) {
				__queue_insert_fast(qA, n);

				node_t *left  = &nodes[n->internal.left_index],
					   *right = &nodes[n->internal.right_index];

				if (unlikely(rand > left->internal.child_sum)) {
					rand -= left->internal.child_sum;
					n = right;
				} else {
					n = left;
				}
			} else {
				arr_t *arr = &n->array;

				// Index into current "array"
				int randvar = (int)(rand/n->array.p_hat);
				int idx = arr->offset + arr->sampled_size + randvar;

				assert(idx >= 0);

				if (XRANDFUN() > w[idx]/n->array.p_hat) {
					// Cancel this sample.
					i--;
					//rejects += 1;
					goto REDO;
				}

				sampled[i] = st->reverse_map[idx];
				sample_weight = n->array.p_hat;

				// We now have one less element to sample here.
				arr->sampled_sum += sample_weight;
				arr->sum  -= sample_weight;
				arr->size -= 1;

				int sbuckets = sampled_buckets;

				/**
				 * Swap the sampled element to the front of the list
				 */

				// Amount of elements we have already sampled from this bucket.
				int sampled_size = arr->sampled_size;
				int offset = arr->offset; // Offset into weights for index 0.

				// Swap the weights
				double tmp = w[idx];
				w[idx] = w[offset + sampled_size];
				w[offset + sampled_size] = tmp;

				// Swap the reverse map
				int map = st->reverse_map[idx];
				st->reverse_map[idx] = st->reverse_map[offset + sampled_size];
				st->reverse_map[offset + sampled_size] = map;

				arr->sampled_size++; // All done

				assert(n->array.size >= 0);

				/**
				 * Swap around the sampled buckets, so we can reset the tree
				 */

				// If we haven't sampled from this bucket yet.
				if (sbuckets - 1 < arr->bucket_index) {

					// If we actually need to move it, or we can just->increment
					// the size of our samples.
					if (sbuckets < n->array.bucket_index) {

						// n2 is the element after our sampled list->
						node_t *n2 = st->buckets[sbuckets];

						// Swap bucket indexes.
						n2->array.bucket_index = n->array.bucket_index;
						st->buckets[n2->array.bucket_index] = n2;
						n->array.bucket_index = sbuckets;
						st->buckets[sbuckets] = n;
					}

					sampled_buckets++;
				}

				break;
			}
		} while(1);

		for (int j = 0; j < qA->size; j++) {
			qA->nodes[j]->internal.child_sum -= sample_weight;
		}
REDO:
		qA->size = 0;
	}

	for (int j = 0; j < sampled_buckets; j++) {
		arr_t *arr = &st->buckets[j]->array;

		arr->size += arr->sampled_size;
		arr->sampled_size = 0;
		double add_weight = arr->sampled_sum;

		arr->sum += add_weight;
		arr->sampled_sum = 0.0;

		int pindex = st->buckets[j]->parent_index;

		node_t *parent;
		while (pindex != INT_MAX) {
			parent = &nodes[pindex];
			parent->internal.child_sum += add_weight;

			pindex = parent->parent_index;
		}
	}

	//printf("Rejects: %"PRIu64"\n", rejects);

	return 0;
}

void rstree_free(rstree_t *st) {
	if (st != NULL) {
		XFREE(st->reverse_map);
		XFREE(st->w);
		XFREE(st->buckets);
		XFREE(st->queue.nodes);
		XFREE(st->tree);

		XFREE(st);
	}
}
