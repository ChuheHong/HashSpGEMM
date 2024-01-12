#include "../include/hash_table.h"

#include </opt/homebrew/opt/libomp/include/omp.h>
#include <stdlib.h>

#include "../include/CSR.h"
#include "../include/malloc2D.h"
#include "../include/utility.h"

/* Compute the total intermediate product and the "max" nnz of each row */
void set_intprod_num(HashTable *hash_table, const IT *arpt, const IT *acol,
                     const IT *brpt, const IT rows) {
#pragma omp parallel
    {
        IT each_int_prod = 0;
#pragma omp for
        for (IT i = 0; i < rows; ++i) {
            IT nnz_per_row = 0;
            // j represents to the j-th entry of row i
            for (IT j = arpt[i]; j < arpt[i + 1]; ++j) {
                // acol[j] represents to the column of j-th entry and the
                // corresponding (acol[j]) row of matrix b, so "brpt[acol[j] +
                // 1] - brpt[acol[j]]" can compute the corresponding nums of
                // entries of b
                nnz_per_row += brpt[acol[j] + 1] - brpt[acol[j]];
            }
            hash_table->row_nnz[i] = nnz_per_row;
            each_int_prod += nnz_per_row;
        }
#pragma omp atomic
        hash_table->total_intprod += each_int_prod;
    }
}

/* Get total number of floating operations and average
 * then, use it for assigning rows to thread as the amount of work is equally
 * distributed
 */
void set_rows_offset(HashTable *hash_table, const IT rows) {
    IT *ps_row_nnz = (IT *)malloc((rows + 1) * sizeof(IT));

    /* Prefix sum of #intermediate products */
    scan(hash_table->row_nnz, ps_row_nnz, rows + 1);

    IT average_intprod =
        (hash_table->total_intprod + hash_table->thread_num - 1) /
        hash_table->thread_num;
    // long long int average_intprod = total_intprod / thread_num;

    /* Search end point of each range */
    hash_table->rows_offset[0] = 0;
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        // assign the begin and end of the rows for the rows_offset
        IT end_itr =
            lower_bound(ps_row_nnz, 0, rows, average_intprod * (tid + 1));
        hash_table->rows_offset[tid + 1] = end_itr;
    }
    hash_table->rows_offset[hash_table->thread_num] = rows;
    free(ps_row_nnz);
}

/*
 * Precompute how many entries each row requires for the hash table
 * the size is 2^table_id
 */
void set_table_id(HashTable *hash_table, const IT rows, const IT cols,
                  const IT min) {
    IT i;
#pragma omp parallel for
    for (i = 0; i < rows; ++i) {
        IT j;
        IT nnz_per_row = hash_table->row_nnz[i];
        if (nnz_per_row > cols) nnz_per_row = cols;
        if (nnz_per_row == 0) {
            hash_table->table_id[i] = 0;
        } else {
            j = 0;
            while (nnz_per_row > (min << j)) {
                j++;
            }
            hash_table->table_id[i] = j + 1;
        }
    }
}

void create_local_hash_table(HashTable *hash_table, const IT cols) {
#pragma omp parallel
    {
        IT ht_size = 0;
        int tid = omp_get_thread_num();
        /* Get max size of hash table */
        for (IT j = hash_table->rows_offset[tid];
             j < hash_table->rows_offset[tid + 1]; ++j) {
            if (ht_size < hash_table->row_nnz[j])
                ht_size = hash_table->row_nnz[j];
        }
        /* the size of hash table is aligned as 2^n */
        if (ht_size > 0) {
            if (ht_size > cols) ht_size = cols;
            int k = MIN_HT_S;
            while (k < ht_size) {
                k <<= 1;
            }
            ht_size = k;
        }

        hash_table->local_hash_table_id[tid] =
            (IT *)malloc(ht_size * sizeof(IT));
        hash_table->local_hash_table_val[tid] =
            (NT *)malloc(ht_size * sizeof(NT));
    }
}

void hash_symbolic_phase(HashTable *hash_table, IT *arpt, IT *acol, IT *brpt,
                         IT *bcol) {
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        IT start_row = hash_table->rows_offset[tid];
        IT end_row = hash_table->rows_offset[tid + 1];
        IT *check = hash_table->local_hash_table_id[tid];
        for (IT i = start_row; i < end_row; ++i) {
            IT nnz = 0;
            // the table_id represents the num of << (multiply by 2) of the
            // MIN_HT_S
            IT table_id = hash_table->table_id[i];

            if (table_id > 0) {
                // determine hash table size for i-th row
                IT ht_size = MIN_HT_S << (table_id - 1);
                // initialize hash table
                for (IT j = 0; j < ht_size; ++j) {
                    check[j] = -1;
                }

                for (IT j = arpt[i]; j < arpt[i + 1]; ++j) {
                    IT t_acol = acol[j];
                    for (IT k = brpt[t_acol]; k < brpt[t_acol + 1]; ++k) {
                        IT key = bcol[k];
                        IT hash = (key * HASH_SCAL) & (ht_size - 1);
                        // Loop for hash prohash_tableg
                        while (1) {
                            // if the key is already inserted, it's ok
                            if (check[hash] == key) {
                                break;
                            }
                            // if the key has not been inserted yet, then it's
                            // added.
                            else if (check[hash] == -1) {
                                check[hash] = key;
                                nnz++;
                                break;
                            }
                            // linear prohash_tableg: check next entry
                            else {
                                // hash = (hash + 1) % ht_size
                                hash = (hash + 1) & (ht_size - 1);
                            }
                        }
                    }
                }
            }
            hash_table->row_nnz[i] = nnz;
        }
    }
}

void hash_numeric_phase(HashTable *hash_table, IT *arpt, IT *acol, NT *aval,
                        IT *brpt, IT *bcol, NT *bval, IT *crpt, IT *ccol,
                        NT *cval) {
#pragma omp parallel
    {
        IT tid = omp_get_thread_num();
        IT start_row = hash_table->rows_offset[tid];
        IT end_row = hash_table->rows_offset[tid + 1];

        IT *ht_check = hash_table->local_hash_table_id[tid];
        NT *ht_value = hash_table->local_hash_table_val[tid];

        for (IT i = start_row; i < end_row; ++i) {
            IT table_id = hash_table->table_id[i];
            if (table_id > 0) {
                IT offset = crpt[i];
                IT ht_size = MIN_HT_N << (table_id - 1);
                for (IT j = 0; j < ht_size; ++j) {
                    ht_check[j] = -1;
                }
                for (IT j = arpt[i]; j < arpt[i + 1]; ++j) {
                    IT t_acol = acol[j];
                    NT t_aval = aval[j];
                    for (IT k = brpt[t_acol]; k < brpt[t_acol + 1]; ++k) {
                        NT t_val = t_aval * bval[k];
                        IT key = bcol[k];
                        IT hash = (key * HASH_SCAL) & (ht_size - 1);
                        // Loop for hash probing
                        while (1) {
                            // key is already inserted
                            if (ht_check[hash] == key) {
                                ht_value[hash] += t_val;
                                break;
                            }
                            // insert new entry
                            else if (ht_check[hash] == -1) {
                                ht_check[hash] = key;
                                ht_value[hash] = t_val;
                                break;
                            } else {
                                // (hash + 1) % ht_size
                                hash = (hash + 1) & (ht_size - 1);
                            }
                        }
                    }
                }
                // insert the value to the result matrix
                IT index = 0;
                for (IT j = 0; j < ht_size; ++j) {
                    if (ht_check[j] != -1) {
                        ccol[index + offset] = ht_check[j];
                        cval[index + offset] = ht_value[j];
                        index++;
                    }
                }
            }
        }
    }
}