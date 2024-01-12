#include "../include/CSR.h"

void PrepareSparseMatrix(NT *matrix, IT num_rows, IT num_cols, NT ratio) {
    IT nozero = ratio * num_rows * num_cols;
    srand((unsigned)time(NULL));
    for (IT h = 0; h < nozero; h++) {
        IT i = rand() % num_rows;
        IT j = rand() % num_cols;
        if (matrix[i * num_cols + j] != 0) {
            h--;
            continue;
        }
        NT t = (NT)(rand() % 10);
        matrix[i * num_cols + j] = (NT)t;
    }
}

void to_CSR(NT *matrix, CSRMatrix *csr_matrix, IT m, IT n, IT nnz) {
    csr_matrix->num_rows = m;
    csr_matrix->num_cols = n;
    csr_matrix->nnz = nnz;
    csr_matrix->row_ptr = (IT *)malloc((m + 1) * sizeof(IT));
    csr_matrix->col_idx = (IT *)malloc(nnz * sizeof(IT));
    csr_matrix->values = (NT *)malloc(nnz * sizeof(NT));
    IT size = 0;
    for (IT row = 0; row < m; row++) {
        csr_matrix->row_ptr[row] = size;
        for (IT col = 0; col < n; col++) {
            NT num = matrix[row * n + col];
            // if the number is non-zero, store it in the CSR matrix
            if (num != 0) {
                csr_matrix->values[size] = num;
                csr_matrix->col_idx[size] = col;
                size++;
            }
        }
    }
    csr_matrix->row_ptr[m] = size;
}