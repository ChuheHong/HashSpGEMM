#ifndef _CSR_H_
#define _CSR_H_

#include <stdlib.h>
#include <sys/time.h>
#define IT int
#define NT float

typedef struct {
    IT num_rows;
    IT num_cols;
    IT nnz;
    IT *row_ptr;
    IT *col_idx;
    NT *values;
} CSRMatrix;

void PrepareSparseMatrix(NT *matrix, IT num_rows, IT num_cols, NT ratio);
void to_CSR(NT *matrix, CSRMatrix *csr_matrix, IT m, IT n, IT nnz);

#endif