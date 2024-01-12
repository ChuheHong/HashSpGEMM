#ifndef _HASH_TABLE_H_
#define _HASH_TABLE_H_

#define IT int
#define NT float
#define HASH_SCAL 107  // Set disjoint number to hash table size (=2^n)
#define MIN_HT_S 8     // minimum hash table size per row in symbolic phase
#define MIN_HT_N 8     // minimum hash table size per row in numeric phase

typedef struct {
    long long int total_intprod;  // total intermediate product
    IT thread_num;

    IT *row_nnz;  // the number of flop or non-zero elements of output matrix,
                  // size of rows
    IT *rows_offset;  // allocate the rows for the thread, size of (thread + 1)
    IT *table_id;  // MIN_HT_S << table_id, to precompute how many entries each
                   // row requires for the hash table, size of rows
    IT **local_hash_table_id;   // size of thread
    NT **local_hash_table_val;  // size of thread
} HashTable;

/* Compute the total intermediate product and the "max" nnz of each row */
void set_intprod_num(HashTable *hash_table, const IT *arpt, const IT *acol,
                     const IT *brpt, const IT rows);

/* Get total number of floating operations and average
 * then, use it for assigning rows to thread as the amount of work is equally
 * distributed
 */
void set_rows_offset(HashTable *hash_table, const IT rows);
/*
 * Precompute how many entries each row requires for the hash table
 * the size is 2^table_id
 */
void set_table_id(HashTable *hash_table, const IT rows, const IT cols,
                  const IT min);

void create_local_hash_table(HashTable *hash_table, const IT cols);

void hash_symbolic_phase(HashTable *hash_table, IT *arpt, IT *acol, IT *brpt,
                         IT *bcol);

void hash_numeric_phase(HashTable *hash_table, IT *arpt, IT *acol, NT *aval,
                        IT *brpt, IT *bcol, NT *bval, IT *crpt, IT *ccol,
                        NT *cval);
                        
#endif