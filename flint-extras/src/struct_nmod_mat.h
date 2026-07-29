#include "pml.h"
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>

/* Toeplitz matrices stored as
 *     a3 a2 a1 a0
 *     a4 a3 a2 a1
 *     a5 a4 a3 a2
 */
typedef struct
{
    long nrows,ncols;
    ulong n;
    nn_ptr data;
} struct_nmod_toeplitz_struct;

typedef struct_nmod_toeplitz_struct struct_nmod_toeplitz_t[1];

/* Initialization */
void struct_nmod_toeplitz_init(struct_nmod_toeplitz_t mat,
                                long nrows, 
                                long ncols,
                                ulong N);

void struct_nmod_toeplitz_set(struct_nmod_toeplitz_t mat,
                               nn_ptr data);

void struct_nmod_toeplitz_init_set(struct_nmod_toeplitz_t mat,
                                    long nrows, 
                                    long ncols,
                                    ulong N,
                                    nn_ptr data);
    
/* Getters */
slong struct_nmod_toeplitz_nrows(struct_nmod_toeplitz_t mat);
slong struct_nmod_toeplitz_ncols(struct_nmod_toeplitz_t mat);

/* Computes the dense representation of a Toeplitz matrix */
slong struct_nmod_toeplitz_dense(struct_nmod_toeplitz_t mat, 
                                  nmod_mat_t dense_mat);

/* Computes the matrix-vector product with a Toeplitz matrix */
slong struct_nmod_toeplitz_pdt_vec(struct_nmod_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res);
