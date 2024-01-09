#ifndef _CSC_H_
#define _CSC_H_

#define IT int
#define NT float

typedef struct {
    IT num_rows;
    IT num_cols;
    IT nnz;
    IT *row_idx;
    IT *col_ptr;
    NT *values;
} CSCMatrix;

#endif