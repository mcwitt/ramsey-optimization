/* Code adapted from Numerical Recipes in C (2nd & 3rd Eds.) */

#define QSELECT_MAXLEN 1024

typedef int elem_t;

double qselect(unsigned long k, unsigned long n, elem_t arr[]);
/*
 * Given k in [0..n-1] returns an array value from arr[0..n-1] such that k
 * array values are less than or equal to the one returned. The input array
 * will be rearranged to have this value in location arr[k], with all smaller
 * elements movd to arr[0..k-1] (in arbitrary order) and all larger elements in
 * arr[k+1..n-1] (also in arbitrary order).
 */

unsigned long qselect_index(unsigned long k, unsigned long n, elem_t arr[]);
/* Returns the index of the kth smallest value in the array arr[0..n-1] */

