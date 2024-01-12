#ifndef _UTILITY_H
#define _UTILITY_H

// Prefix sum (Thread parallel)
void scan(IT *in, IT *out, IT N);

/* return the coordinate of the first occurrence of val in A */
IT lower_bound(IT *A, IT L, IT R, IT val); 

#endif