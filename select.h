/* Code adapted from Numerical Recipes in C (2nd Ed.) */

double select(unsigned long k, unsigned long n, double arr[]);
/*
 * Returns the kth smallest value in the array arr[1..n]. The input array will
 * be rearranged to have this value in location arr[k], with all smaller
 * elements moved to arr[1..k-1] (in arbitrary order) and all larger elements
 * in arr[k+1..n] (also in arbitrary order).
 */
