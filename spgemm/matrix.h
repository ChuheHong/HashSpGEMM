#ifndef MATRIX_H
#define MATRIX_H

struct COOMatrix
{
    float *value;
    int rows;
    int cols;
    int nnz;
};

struct CSRMatrix
{
    float *value;
    int rows;
    int cols;
    int nnz;
    int *row_ptr;
    int *col_index;
};



#endif