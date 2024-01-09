#include "/opt/homebrew/opt/libomp/include/omp.h"
#include <stdio.h>
#include <stdlib.h>

#include "CSR.h"
#include "hash_table.h"
#include "utility.h"

void HashSpGEMM() {
    /* prepare the input matrix a and b */
    IT m, k, n, rows;
    m = k = n = rows = 8192;
    NT *ma = (NT *)malloc(n * m * sizeof(NT));
    NT *mb = (NT *)malloc(m * n * sizeof(NT));
    PrepareSparseMatrix(ma, n, m, 0.01);
    PrepareSparseMatrix(mb, n, m, 0.01);
    CSRMatrix a, b;
    to_CSR(ma, &a, m, n, m * n * 0.01);
    to_CSR(mb, &b, m, n, m * n * 0.01);
    free(ma);
    free(mb);

    /* prepare the output matrix c*/
    CSRMatrix c;
    c.num_rows = a.num_rows;
    c.num_cols = b.num_cols;

    /* prepare the HashTable accumulator */
    HashTable hash_table;
    hash_table.total_intprod = 0;
    hash_table.thread_num = omp_get_max_threads();
    printf("Num of thread: %d\n", hash_table.thread_num);
    hash_table.row_nnz = (IT *)malloc(rows * sizeof(IT));
    hash_table.rows_offset =
        (IT *)malloc((hash_table.thread_num + 1) * sizeof(IT));
    hash_table.table_id = (IT *)malloc(rows * sizeof(IT));
    hash_table.local_hash_table_id =
        (IT **)malloc(hash_table.thread_num * sizeof(IT *));
    hash_table.local_hash_table_val =
        (NT **)malloc(hash_table.thread_num * sizeof(NT *));

    double parallel_time_begin = omp_get_wtime();
    set_intprod_num(&hash_table, a.row_ptr, a.col_idx, b.row_ptr, c.num_rows);
    set_rows_offset(&hash_table, c.num_rows);
    set_table_id(&hash_table, c.num_rows, c.num_cols, MIN_HT_S);
    create_local_hash_table(&hash_table, c.num_cols);
    c.row_ptr = (IT *)malloc((c.num_rows + 1) * sizeof(IT));

    // symbolic phase
    hash_symbolic_phase(&hash_table, a.row_ptr, a.col_idx, b.row_ptr,
                        b.col_idx);
    // the nnz of each row is accurete
    scan(hash_table.row_nnz, c.row_ptr, c.num_rows + 1);
    c.nnz = c.row_ptr[c.num_rows];
    c.col_idx = (IT *)malloc(c.nnz * sizeof(IT));
    c.values = (NT *)malloc(c.nnz * sizeof(NT));
    // numeric phase
    hash_numeric_phase(&hash_table, a.row_ptr, a.col_idx, a.values, b.row_ptr,
                       b.col_idx, b.values, c.row_ptr, c.col_idx, c.values);
    double parallel_time_end = omp_get_wtime();

    // show the results
    printf("parallel time: %lfs\n", parallel_time_end - parallel_time_begin);
    double up = (double)(0.01 * 0.01 * m * k * n);
    double down = (parallel_time_end - parallel_time_begin) * 1e9;
    printf("GFlops: %lf\n", up / down);
}

void test() {
    int n, k, m;
    n = k = m = 256;
    NT *matrix = (NT *)malloc(n * m * sizeof(NT));
    PrepareSparseMatrix(matrix, n, m, 0.01);
    // for (int i = 0; i < n * m; ++i)
    //     if (matrix[i] > 0) printf("%f\n", matrix[i]);
    CSRMatrix csr_matrix;
    to_CSR(matrix, &csr_matrix, m, n, m * n * 0.01);
    printf("%d %d %d\n", csr_matrix.num_rows, csr_matrix.num_cols,
           csr_matrix.nnz);
    printf("%d %d %f\n", csr_matrix.row_ptr[2], csr_matrix.col_idx[2],
           csr_matrix.values[2]);
}

int main(int argc, char *argv[]) {
    HashSpGEMM();
    return 0;
}
