#ifndef _UTILITY_H
#define _UTILITY_H
#include "CSR.h"
#include </opt/homebrew/opt/libomp/include/omp.h>

// Prefix sum (Thread parallel)
void scan(IT *in, IT *out, IT N) {
    if (N < (1 << 17)) {
        out[0] = 0;
        for (IT i = 0; i < N - 1; ++i) {
            out[i + 1] = out[i] + in[i];
        }
    } else {
        int tnum = 64;

        IT each_n = N / tnum;
        IT *partial_sum = (IT *)malloc(sizeof(IT) * (tnum));
#pragma omp parallel num_threads(tnum)
        {
            int tid = omp_get_thread_num();
            IT start = each_n * tid;
            IT end = (tid < tnum - 1) ? start + each_n : N;
            out[start] = 0;
            for (IT i = start; i < end - 1; ++i) {
                out[i + 1] = out[i] + in[i];
            }
            partial_sum[tid] = out[end - 1] + in[end - 1];
#pragma omp barrier

            IT offset = 0;
            for (int i = 0; i < tid; ++i) {
                offset += partial_sum[i];
            }
            for (IT i = start; i < end; ++i) {
                out[i] += offset;
            }
        }
        free(partial_sum);
    }
}

/* return the coordinate of the first occurrence of val in A */
IT lower_bound(IT *A, IT L, IT R, IT val) {
    IT l = L, r = R;
    while (l < r) {
        IT m = (l + r) / 2;
        if (val <= A[m]) {
            r = m - 1;
        } else {
            l = m + 1;
        }
    }
    return val <= A[l] ? l : l + 1;
}

#endif