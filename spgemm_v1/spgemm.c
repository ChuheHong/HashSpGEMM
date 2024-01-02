// We implement matrix multiplication using the Gustavson method.
// And here we use the OpenMP library which is more convenient to our accomplishment.
// If you want to run our program successfully, you should use the next shown command:
// module load GCC/9.3.0
// gcc -fopenmp spgemm.c -o spgemm
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "matrix.h"

// C = A * B
void sparse_matrix_mul(struct CSRMatrix matrix_a, struct CSRMatrix matrix_b,
                       struct COOMatrix *matrix_c, int threads);
// transfer the COO format to the CSR format
void coo_to_csr(struct COOMatrix *coo_matrix, struct CSRMatrix *csr_matrix);
// read the matrix of COO format
void coo_read_matrix(struct COOMatrix *matrix, FILE *f);

int main(int argc, char const *argv[])
{
    struct COOMatrix coo_matrix_a;
    struct COOMatrix coo_matrix_b;
    coo_read_matrix(&coo_matrix_a, fopen("datasets/494_bus.mtx", "r"));
    coo_read_matrix(&coo_matrix_b, fopen("datasets/494_bus.mtx", "r"));
    struct CSRMatrix csr_matrix_a;
    struct CSRMatrix csr_matrix_b;
    coo_to_csr(&coo_matrix_a, &csr_matrix_a);
    coo_to_csr(&coo_matrix_b, &csr_matrix_b);
    struct COOMatrix coo_matrix_c;
    sparse_matrix_mul(csr_matrix_a, csr_matrix_b, &coo_matrix_c, 10);
    
    return 0;
}

void sparse_matrix_mul(struct CSRMatrix matrix_a, struct CSRMatrix matrix_b,
                       struct COOMatrix *matrix_c, int threads)
{
    matrix_c->rows = matrix_a.rows;
    matrix_c->cols = matrix_b.cols;
    matrix_c->value = (float *)malloc(matrix_c->rows * matrix_c->cols * sizeof(float));
    matrix_c->nnz = 0;
#pragma omp parallel for num_threads(threads)
    // traverse the row of A parallelly
    for (int m = 0; m < matrix_a.rows; m++)
    {
        int row_start = matrix_a.row_ptr[m];
        int row_end = matrix_a.row_ptr[m + 1];
        // traverse the mth row of A
        for (int idx1 = row_start; idx1 < row_end; idx1++)
        {
            // the kth column of A
            int k = matrix_a.col_index[idx1];
            // traverse the kth row of B
            for (int idx2 = matrix_b.row_ptr[k]; idx2 < matrix_b.row_ptr[k + 1]; idx2++)
            {
                // the row of A plus the column of B
                matrix_c->value[m * matrix_c->cols + matrix_b.col_index[idx2]] += matrix_a.value[idx1] * matrix_b.value[idx2];
            }
        }
    }
    // optimize
    // for (int i = 0; i < matrix_c.rows * matrix_c.cols; i++)
    //     if (matrix_c.value[i] != 0)
    //         matrix_c.nnz++;
}

void coo_read_matrix(struct COOMatrix *matrix, FILE *f)
{
#define MAXLEN 511
    for (char s[MAXLEN]; fgets(s, MAXLEN, f) && !feof(f);)
        if (s[0] != '%')
        {
            int rows, cols, nnz;
            sscanf(s, "%d%d%d", &rows, &cols, &nnz);
            // create the COO matrix
            matrix->rows = rows;
            matrix->cols = cols;
            matrix->value = (float *)malloc(rows * cols * sizeof(float));
            matrix->nnz = nnz;
            for (int i = 0; i < rows * cols; i++)
                matrix->value[i] = 0;
            for (int line = 0; line < nnz; ++line)
            {
                int r, c;
                float v;
                fscanf(f, "%d%d%f", &r, &c, &v);
                if (0 <= r && r < rows && 0 <= c && c < cols)
                {
                    matrix->value[(r - 1) * cols + (c - 1)] = v;
                }
            }
            break;
        }
#undef MAXLEN
}

void coo_to_csr(struct COOMatrix *coo_matrix, struct CSRMatrix *csr_matrix)
{
    csr_matrix->rows = coo_matrix->rows;
    csr_matrix->cols = coo_matrix->cols;
    csr_matrix->nnz = coo_matrix->nnz;
    csr_matrix->row_ptr = (int *)malloc((csr_matrix->rows + 1) * sizeof(int));
    csr_matrix->col_index = (int *)malloc(csr_matrix->nnz * sizeof(int));
    csr_matrix->value = (float *)malloc(csr_matrix->nnz * sizeof(float));
    int size = 0;
    for (int row = 0; row < coo_matrix->rows; row++)
    {
        csr_matrix->row_ptr[row] = size;
        for (int col = 0; col < coo_matrix->cols; col++)
        {
            float num = coo_matrix->value[row * coo_matrix->cols + col];
            // if the number is non-zero, store it in the CSR matrix
            if (num != 0)
            {
                csr_matrix->value[size] = num;
                csr_matrix->col_index[size] = col;
                size++;
            }
        }
    }
    csr_matrix->row_ptr[csr_matrix->rows] = size;
}