#include "../include/malloc2D.h"

#include <stdlib.h>

#include "../include/CSR.h"

/* Single contiguous memory allocation of a 2D array */
NT **malloc2D_NT(IT imax, IT jmax) {
    NT **x = (NT **)malloc(imax * sizeof(NT *) + imax * jmax * sizeof(NT));

    x[0] = (NT *)x + imax;

    for (IT i = 1; i < imax; i++) {
        x[i] = x[i - 1] + jmax;
    }
    return x;
}

IT **malloc2D_IT(IT imax, IT jmax) {
    IT **x = (IT **)malloc(imax * sizeof(IT *) + imax * jmax * sizeof(IT));

    x[0] = (IT *)x + imax;

    for (IT i = 1; i < imax; i++) {
        x[i] = x[i - 1] + jmax;
    }
    return x;
}